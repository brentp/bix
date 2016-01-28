// Tabix queries for go
package bix

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"unsafe"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/tabix"
	"github.com/brentp/cgotbx"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/irelate/parsers"
	"github.com/brentp/vcfgo"
)

const (
	VCF = iota
	BED
	GENERIC
)

type location struct {
	chrom string
	start int
	end   int
}

func (s location) RefName() string {
	return s.chrom
}
func (s location) Start() int {
	return s.start
}
func (s location) End() int {
	return s.end
}

// Bix provides read access to tabix files.
type Bix struct {
	*tabix.Index
	bgzf *bgzf.Reader
	path string

	VReader *vcfgo.Reader

	file *os.File
	buf  *bufio.Reader
	cgo  *cgotbx.Tbx
}

// create a new bix that does as little as possible from the old bix
func newShort(old *Bix) (*Bix, error) {
	tbx := &Bix{
		Index:   old.Index,
		path:    old.path,
		VReader: old.VReader,
	}
	var err error
	tbx.file, err = os.Open(tbx.path)
	if err != nil {
		return nil, err
	}
	tbx.bgzf, err = bgzf.NewReader(tbx.file, 1)
	if err != nil {
		return nil, err
	}
	return tbx, nil
}

// New returns a &Bix
func New(path string, workers ...int) (*Bix, error) {
	f, err := os.Open(path + ".tbi")
	if err != nil {
		return nil, err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer gz.Close()

	idx, err := tabix.ReadFrom(gz)
	if err != nil {
		return nil, err
	}
	n := 1
	if len(workers) > 0 {
		n = workers[0]
	}

	b, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	bgz, err := bgzf.NewReader(b, n)
	if err != nil {
		return nil, err
	}

	/*
		c := cache.NewFIFO(4)
		bgz.SetCache(c)
	*/

	buf := bufio.NewReaderSize(bgz, 16384)
	h := make([]string, 0)
	tbx := &Bix{bgzf: bgz, path: path, file: b}

	l, err := buf.ReadString('\n')
	if err != nil {
		return tbx, err
	}

	for i := 0; i < int(idx.Skip) || rune(l[0]) == idx.MetaChar; i++ {
		h = append(h, l)
		l, err = buf.ReadString('\n')
		if err != nil {
			return tbx, err
		}
	}
	header := strings.Join(h, "")

	if len(h) > 0 && strings.HasSuffix(tbx.path, ".vcf.gz") {
		var err error
		h := strings.NewReader(header)

		tbx.VReader, err = vcfgo.NewReader(h, true)
		if err != nil {
			return nil, err
		}
	}
	tbx.cgo, err = cgotbx.New(path)
	if err != nil {
		return nil, err
	}

	tbx.buf = buf
	tbx.Index = idx
	return tbx, nil
}

func (b *Bix) Close() error {
	if b.bgzf != nil {
		b.bgzf.Close()
	}
	if b.file != nil {
		b.file.Close()
	}
	return nil
}

func (tbx *Bix) toPosition(toks [][]byte) interfaces.Relatable {
	isVCF := tbx.VReader != nil

	if isVCF {
		v := tbx.VReader.Parse(toks)
		return interfaces.AsRelatable(v)
	} else {
		g, _ := newgeneric(toks, int(tbx.Index.NameColumn-1), int(tbx.Index.BeginColumn-1),
			int(tbx.Index.EndColumn-1), tbx.Index.ZeroBased)
		return g
	}
}

/*
func (tbx *Bix) GetHeaderType(field string) string {
	if tbx.VReader == nil {
		return ""
	}
	hdr := tbx.VReader.Header

	info, ok := hdr.Infos[field]
	if !ok {
		return ""
	}
	return info.Type
}*/

func unsafeString(b []byte) string {
	return *(*string)(unsafe.Pointer(&b))
}

// return an interval using the info from the tabix index
func newgeneric(fields [][]byte, chromCol int, startCol int, endCol int, zeroBased bool) (*parsers.Interval, error) {
	s, err := strconv.Atoi(unsafeString(fields[startCol]))
	if err != nil {
		return nil, err
	}
	if !zeroBased {
		s -= 1
	}
	e, err := strconv.Atoi(unsafeString(fields[endCol]))
	if err != nil {
		return nil, err
	}
	return parsers.NewInterval(string(fields[chromCol]), uint32(s), uint32(e), fields, uint32(0), nil), nil
}

func (tbx *Bix) ChunkedReader(r tabix.Record) (io.Reader, error) {
	chunks, err := tbx.Chunks(r)
	if err == index.ErrNoReference {
		l := location{r.RefName(), r.Start(), r.End()}
		if strings.HasPrefix(l.chrom, "chr") {
			l.chrom = l.chrom[3:]
			chunks, err = tbx.Chunks(l)
		} else {
			l.chrom = "chr" + l.chrom
			chunks, err = tbx.Chunks(l)
		}
	}
	if err == index.ErrInvalid {
		return index.NewChunkReader(tbx.bgzf, []bgzf.Chunk{})
	} else if err == index.ErrNoReference {
		log.Printf("chromosome %s not found in %s\n", r.RefName(), tbx.path)
		return index.NewChunkReader(tbx.bgzf, []bgzf.Chunk{})
	} else if err != nil {
		return nil, err
	}
	/*
		if len(chunks) == 1 && chunks[0].Begin.File == chunks[0].End.File {
			//fmt.Fprintf(os.Stderr, "\n%+v\n%+v\n", r, chunks)
			tbx.bgzf.Seek(chunks[0].Begin)
			r := io.LimitReader(tbx.bgzf, int64(chunks[0].End.Block))
			return r, nil
		}*/

	cr, err := index.NewChunkReader(tbx.bgzf, chunks)
	if err != nil {
		return nil, err
	}
	return cr, nil
}

// bixerator meets interfaces.RelatableIterator
type bixerator struct {
	rdr io.Reader
	buf *bufio.Reader
	tbx *Bix

	region interfaces.IPosition
}

func makeFields(line []byte) [][]byte {
	fields := make([][]byte, 9)
	copy(fields[:8], bytes.SplitN(line, []byte{'\t'}, 8))
	s := 0
	for i, f := range fields {
		if i == 7 {
			break
		}
		s += len(f) + 1
	}
	e := bytes.IndexByte(line[s:], '\t')
	if e == -1 {
		e = len(line)
	} else {
		e += s
	}

	fields[7] = line[s:e]
	if len(line) > e+1 {
		fields[8] = line[e+1:]
	} else {
		fields = fields[:8]
	}

	return fields
}

func (b bixerator) Next() (interfaces.Relatable, error) {

	for {
		line, err := b.buf.ReadBytes('\n')

		if err == io.EOF && len(line) == 0 {
			return nil, io.EOF
		} else if err != nil {
			return nil, err
		}
		if len(line) == 0 {
			return nil, io.EOF
		}
		if line[len(line)-1] == '\n' {
			line = line[:len(line)-1]
		}
		//var in bool = true
		var toks [][]byte
		/*
			if b.region != nil {
				var err error

				in, err, toks = b.inBounds(line)
				if err != nil {
					return nil, err
				}
			} else {
		*/
		if b.tbx.VReader != nil {
			toks = makeFields(line)
		} else {
			toks = bytes.Split(line, []byte{'\t'})
		}
		//}

		//if in {
		return b.tbx.toPosition(toks), nil
		//}
	}
	return nil, io.EOF
}

func (b bixerator) Close() error {
	if b.rdr != nil {
		if c, ok := b.rdr.(io.ReadCloser); ok {
			c.Close()
		}
	}
	return b.tbx.Close()
}

var _ interfaces.RelatableIterator = bixerator{}

func (tbx *Bix) Query(region interfaces.IPosition) (interfaces.RelatableIterator, error) {
	tbx2 := &Bix{VReader: tbx.VReader, Index: tbx.Index}
	/*
		tbx2, err := newShort(tbx)
		if err != nil {
			return nil, err
		}
		if region == nil {
			var l string
			var err error
			buf := bufio.NewReader(tbx2.bgzf)
			l, err = buf.ReadString('\n')
			for i := 0; i < int(tbx2.Index.Skip) || rune(l[0]) == tbx2.Index.MetaChar; i++ {
				l, err = buf.ReadString('\n')
				if err != nil {
					return nil, err
				}
			}
			if tbx2.Index.Skip == 0 && rune(l[0]) != tbx2.Index.MetaChar {
				buf = bufio.NewReader(io.MultiReader(strings.NewReader(l), buf))
			}
			return bixerator{nil, buf, tbx2, region}, nil
		}*/
	//cr, err := tbx2.ChunkedReader(location{chrom: region.Chrom(), start: int(region.Start()), end: int(region.End())})
	cr, err := tbx.cgo.Get(region.Chrom(), int(region.Start()), int(region.End()))

	// return an empty iterator
	if err != nil {
		if cr != nil {
			tbx2.Close()
			if c, ok := cr.(io.ReadCloser); ok {
				c.Close()
			}
		}
		return nil, err
	}
	if cr == nil {
		return nil, nil
	}
	return bixerator{cr, bufio.NewReaderSize(cr, 8192), tbx2, region}, nil
}

func (tbx *Bix) AddInfoToHeader(id, number, vtype, desc string) {
	if tbx.VReader == nil {
		return
	}
	tbx.VReader.AddInfoToHeader(id, number, vtype, desc)
}

func (tbx *Bix) GetHeaderType(field string) string {
	if tbx.VReader == nil {
		return ""
	}
	return tbx.VReader.GetHeaderType(field)
}

func (b *bixerator) inBounds(line []byte) (bool, error, [][]byte) {

	var readErr error
	line = bytes.TrimRight(line, "\r\n")
	var toks [][]byte
	if b.tbx.VReader != nil {
		toks = makeFields(line)
	} else {
		toks = bytes.Split(line, []byte{'\t'})
	}

	s, err := strconv.Atoi(unsafeString(toks[b.tbx.BeginColumn-1]))
	if err != nil {
		return false, err, toks
	}

	pos := s
	if !b.tbx.ZeroBased {
		pos -= 1
	}
	if pos >= int(b.region.End()) {
		return false, io.EOF, toks
	}

	if b.tbx.EndColumn != 0 {
		e, err := strconv.Atoi(unsafeString(toks[b.tbx.EndColumn-1]))
		if err != nil {
			return false, err, toks
		}
		if e < int(b.region.Start()) {
			return false, readErr, toks
		}
		return true, readErr, toks
	} else if b.tbx.VReader != nil {
		start := int(b.region.Start())
		alt := strings.Split(string(toks[4]), ",")
		lref := len(toks[3])
		if start >= pos+lref {
			for _, a := range alt {
				if a[0] != '<' || a == "<CN0>" {
					e := pos + lref
					if e > start {
						return true, readErr, toks
					}
				} else if strings.HasPrefix(a, "<DEL") || strings.HasPrefix(a, "<DUP") || strings.HasPrefix(a, "<INV") || strings.HasPrefix(a, "<CN") {
					info := string(toks[7])
					if idx := strings.Index(info, ";END="); idx != -1 {
						v := info[idx+5 : idx+5+strings.Index(info[idx+5:], ";")]
						e, err := strconv.Atoi(v)
						if err != nil {
							return false, err, toks
						}
						if e > start {
							return true, readErr, toks
						}
					} else {
						log.Println("no end:", b.tbx.path, string(toks[0]), pos, string(toks[3]), a)
					}
				}
			}
		} else {
			return true, readErr, toks
		}
		return false, readErr, toks
	}
	return false, readErr, toks

}

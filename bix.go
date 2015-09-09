// Tabix queries for go
package bix

import (
	"bufio"
	"compress/gzip"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/tabix"
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
	bgzf   *bgzf.Reader
	path   string
	Header string

	vReader *vcfgo.Reader
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

	idx, err := tabix.ReadTabix(gz)
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
	bgz.Blocked = false
	if err != nil {
		return nil, err
	}

	/*
		c := cache.NewFIFO(4)
		bgz.SetCache(c)
	*/

	buf := bufio.NewReaderSize(bgz, 16384/2)
	h := make([]string, 0)
	tbx := &Bix{bgzf: bgz, path: path}

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
	tbx.Header = strings.Join(h, "")

	if len(h) > 0 && strings.HasSuffix(tbx.path, ".vcf.gz") {
		rdr := strings.NewReader(tbx.Header)
		var err error
		tbx.vReader, err = vcfgo.NewReader(rdr, true)
		if err != nil {
			return nil, err
		}
	}

	tbx.Index = idx
	return tbx, nil
}

// bixReader is return from Bix.Query() it meets the io.Reader interface
type bixReader struct {
	startCol int
	endCol   int
	maxCol   int
	add      int

	start int
	end   int

	rdr   *bufio.Reader
	isVCF bool

	buf []byte
}

func newReader(tbx *Bix, cr *index.ChunkReader, start, end int) (io.Reader, error) {
	r := &bixReader{}
	r.startCol = int(tbx.BeginColumn - 1)
	r.endCol = int(tbx.EndColumn - 1)
	r.maxCol = r.endCol
	if r.startCol > r.endCol {
		r.maxCol = r.startCol
	}
	r.start = start
	r.end = end
	r.maxCol += 1
	if !tbx.ZeroBased {
		r.add = -1
	}
	r.rdr = bufio.NewReader(cr)
	r.buf = make([]byte, 0)
	r.isVCF = strings.HasSuffix(tbx.path, ".vcf.gz")
	if r.isVCF {
		r.maxCol = 9
	}

	return r, nil
}

func (r *bixReader) inBounds(line string) (bool, error, []string) {

	var readErr error
	toks := strings.SplitN(line, "\t", r.maxCol+1)

	s, err := strconv.Atoi(toks[r.startCol])
	if err != nil {
		return false, err, toks
	}

	pos := s + r.add
	if pos >= r.end {
		return false, io.EOF, toks
	}

	if r.endCol != -1 {
		e, err := strconv.Atoi(toks[r.endCol])
		if err != nil {
			return false, err, toks
		}
		if e < r.start {
			return false, readErr, toks
		}
	} else if r.isVCF {
		alt := strings.Split(toks[4], ",")
		lref := len(toks[3])
		//log.Println(maxEnd, r.start, lref, toks[3], toks[4])
		if r.start >= pos+lref {
			for _, a := range alt {
				if a[0] != '<' {
					e := pos + lref
					if e > r.start {
						//log.Println("true")
						return true, readErr, toks
					}
				} else if strings.HasPrefix(a, "<DEL") || strings.HasPrefix(a, "<DUP") || strings.HasPrefix(a, "<INV") || strings.HasPrefix(a, "<CN") {
					info := toks[7]
					if idx := strings.Index(info, ";END="); idx != -1 {
						v := info[idx+5 : idx+5+strings.Index(info[idx+5:], ";")]
						e, err := strconv.Atoi(v)
						if err != nil {
							return false, err, toks
						}
						if e > r.start {
							return true, readErr, toks
						}
					} else {
						log.Println("no end")
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

func (tbx *Bix) Get(q interfaces.IPosition) []interfaces.IPosition {
	overlaps := make([]interfaces.IPosition, 0)
	chunkReader, err := tbx.chunkedReader(location{q.Chrom(), int(q.Start()), int(q.End())})
	if err != nil {
		return overlaps
	}
	br := chunkReader.(*bixReader)
	isVCF := tbx.vReader != nil
	isBED := tbx.Index.NameColumn == 1 && tbx.Index.BeginColumn == 2 && tbx.Index.EndColumn == 3
	for {
		line, err := br.rdr.ReadString('\n')

		if err != nil && err != io.EOF {
			log.Println(err)
			break
		}
		if len(line) == 0 || (err == io.EOF && line[len(line)-1] != '\n') {
			break
		}
		if err != nil {
			log.Println(err)
		}
		in, rerr, toks := br.inBounds(line)
		if !in {
			if rerr == io.EOF {
				break
			}
			continue
		}

		if isVCF {
			v := tbx.vReader.Parse(toks)
			if v != nil {
				overlaps = append(overlaps, v)
			}
		} else if isBED {
			v, e := parsers.IntervalFromBedLine(line)
			if e != nil {
				overlaps = append(overlaps, v)
			}
		} else {

			g, e := newgeneric(toks, int(tbx.Index.NameColumn-1), br.startCol, br.endCol, tbx.Index.ZeroBased)
			if e != nil {
				overlaps = append(overlaps, g)
			}

		}
		if err == io.EOF || rerr == io.EOF {
			break
		}
	}
	return overlaps

}

// return an interval using the info from the tabix index
func newgeneric(fields []string, chromCol int, startCol int, endCol int, zeroBased bool) (*parsers.Interval, error) {
	s, err := strconv.Atoi(fields[startCol])
	if err != nil {
		return nil, err
	}
	if !zeroBased {
		s -= 1
	}
	e, err := strconv.Atoi(fields[endCol])
	if err != nil {
		return nil, err
	}

	return parsers.NewInterval(fields[chromCol], uint32(s), uint32(e), fields, uint32(0), nil), nil
}

// Read from a BixReader (will contain only overlapping intervals).
func (r *bixReader) Read(p []byte) (int, error) {
	var readErr error
	var line string

	for len(r.buf) < len(p) {
		line, readErr = r.rdr.ReadString('\n')
		if readErr == io.EOF {
			if len(line) == 0 {
				break
			}
		} else if readErr != nil {
			log.Println(readErr)
			break
		}
		var in bool
		var err error
		if in, err, _ = r.inBounds(line); in {
			r.buf = append(r.buf, line...)
		}
		if err == io.EOF {
			readErr = err
			break
		} else if err != nil {
			log.Println(err)
		}
	}
	if len(r.buf) >= len(p) {
		copy(p, r.buf[:len(p)])
		r.buf = r.buf[len(p):]
		if len(r.buf) > 0 && readErr == io.EOF {
			readErr = nil
		}
		return len(p), readErr
	}
	copy(p, r.buf)
	l := len(r.buf)
	r.buf = r.buf[:0]
	return l, readErr
}

func (tbx *Bix) chunkedReader(l location) (io.Reader, error) {
	chunks, err := tbx.Chunks(l)
	if err != nil {
		return nil, err
	}

	chunkReader, err := index.NewChunkReader(tbx.bgzf, chunks)
	if err != nil {
		return nil, err
	}
	rdr, err := newReader(tbx, chunkReader, l.start, l.end)
	if err != nil {
		return nil, err
	}
	return rdr, nil
}

// Query a Bix struct with genomic coordinates. Returns an io.Reader.
func (tbx *Bix) Query(chrom string, start int, end int, printHeader bool) (io.Reader, error) {

	rdr, err := tbx.chunkedReader(location{chrom: chrom, start: start, end: end})
	if err != nil {
		return nil, err
	}
	if printHeader {
		rdr.(*bixReader).buf = []byte(tbx.Header)
	}
	return rdr, err
}

// Close closes the files associate with a Bix struct.
func (tbx *Bix) Close() {
	tbx.bgzf.Close()
}

func (tbx *Bix) AddInfoToHeader(id, number, vtype, desc string) {
	if tbx.vReader == nil {
		return
	}
	tbx.vReader.AddInfoToHeader(id, number, vtype, desc)
}

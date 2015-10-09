// Tabix queries for go
package bix

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"strconv"
	"strings"
	"unsafe"

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

	vReader    *vcfgo.Reader
	lastChunks []bgzf.Chunk
	cache      []interfaces.IPosition

	file *os.File
}

func chunksEqual(a, b []bgzf.Chunk) bool {
	if len(a) != len(b) || len(a) == 0 || len(b) == 0 {
		return false
	}
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			return false
		}
	}
	return true
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
	bgz.Blocked = false
	if err != nil {
		return nil, err
	}

	/*
		c := cache.NewFIFO(4)
		bgz.SetCache(c)
	*/

	buf := bufio.NewReaderSize(bgz, 16384)
	h := make([]string, 0)
	tbx := &Bix{bgzf: bgz, path: path, cache: make([]interfaces.IPosition, 0, 4000), file: b}

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

func (b *Bix) Close() error {
	//tbx.bgzf.Close()
	return b.file.Close()
}

func (tbx *Bix) toPosition(toks [][]byte) interfaces.Relatable {
	isVCF := tbx.vReader != nil

	if isVCF {
		v := tbx.vReader.Parse(toks)
		return interfaces.AsRelatable(v)
	} else {

		g, _ := newgeneric(toks, int(tbx.Index.NameColumn-1), int(tbx.Index.BeginColumn-1),
			int(tbx.Index.EndColumn-1), tbx.Index.ZeroBased)
		return g
	}
}

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

type frdr struct {
	io.Reader
	h   io.Reader
	c   *index.ChunkReader
	bix *Bix
}

func (x frdr) Close() error {
	if x.bix != nil {
		b := x.bix
		x.bix = nil

		return b.Close()

	}
	return nil

}

func (tbx *Bix) ChunkedReader(r tabix.Record) (io.ReadCloser, error) {
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
	if err != nil {
		return nil, err
	}
	chunkReader, err := index.NewChunkReader(tbx.bgzf, chunks)
	if err != nil {
		return nil, err
	}

	if tbx.Header != "" {
		h := strings.NewReader(tbx.Header)
		b := io.MultiReader(h, chunkReader)
		return frdr{b, h, chunkReader, tbx}, nil
	} else {
		return chunkReader, nil
	}
}

// bixerator meets interfaces.RelatableIterator
type bixerator struct {
	rdr io.ReadCloser
	buf *bufio.Reader
	tbx *Bix
}

func (b bixerator) Next() (interfaces.Relatable, error) {
	line, err := b.buf.ReadBytes('\n')
	if err == io.EOF && len(line) == 0 {
		return nil, io.EOF
	} else if err != nil {
		return nil, err
	}
	toks := bytes.Split(line, []byte{'\t'})
	return b.tbx.toPosition(toks), nil
}

func (b bixerator) Close() error {
	return b.rdr.Close()
}

var _ interfaces.RelatableIterator = bixerator{}

func (tbx *Bix) Query(region interfaces.IPosition) (interfaces.RelatableIterator, error) {
	cr, err := tbx.ChunkedReader(location{chrom: region.Chrom(), start: int(region.Start()), end: int(region.End())})
	if err != nil {
		return nil, err
	}
	return bixerator{cr, bufio.NewReader(cr), tbx}, nil
}

func (tbx *Bix) AddInfoToHeader(id, number, vtype, desc string) {
	if tbx.vReader == nil {
		return
	}
	tbx.vReader.AddInfoToHeader(id, number, vtype, desc)
}

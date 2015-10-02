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

// bixReader is return from Bix.Query() it meets the io.Reader interface
type bixReader struct {
	startCol int
	endCol   int
	maxCol   int
	add      int

	start int
	end   int

	cr    *index.ChunkReader
	rdr   *bufio.Reader
	isVCF bool

	buf []byte
}

func newReader(tbx *Bix, cr *index.ChunkReader, start, end int) (io.ReadCloser, error) {
	r := &bixReader{}
	r.startCol = int(tbx.BeginColumn - 1)
	r.endCol = int(tbx.EndColumn - 1)
	r.maxCol = r.endCol
	if r.startCol > r.endCol {
		r.maxCol = r.startCol
	}
	r.start = start
	r.end = end
	r.maxCol += 2
	if !tbx.ZeroBased {
		r.add = -1
	}
	r.rdr = bufio.NewReader(cr)
	r.cr = cr
	r.buf = make([]byte, 0)
	r.isVCF = strings.HasSuffix(tbx.path, ".vcf.gz")
	if r.isVCF {
		r.maxCol = 9
	}

	return r, nil
}

func (r *bixReader) Close() error {
	return r.cr.Close()
}

func (r *bixReader) inBounds(line []byte) (bool, error, [][]byte) {

	var readErr error
	line = bytes.TrimRight(line, "\r\n")
	toks := bytes.Split(line, []byte{'\t'}) //, r.maxCol+1)

	s, err := strconv.Atoi(unsafeString(toks[r.startCol]))
	if err != nil {
		return false, err, toks
	}

	pos := s + r.add
	if pos >= r.end {
		return false, io.EOF, toks
	}

	if r.endCol != -1 {
		e, err := strconv.Atoi(unsafeString(toks[r.endCol]))
		if err != nil {
			return false, err, toks
		}
		if e < r.start {
			return false, readErr, toks
		}
		return true, readErr, toks
	} else if r.isVCF {
		alt := strings.Split(string(toks[4]), ",")
		lref := len(toks[3])
		if r.start >= pos+lref {
			for _, a := range alt {
				if a[0] != '<' {
					e := pos + lref
					if e > r.start {
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

func (tbx *Bix) fillCache(br *bixReader) {
	// this is a naive cache. It just parses the entire chunk.
	// Even if we only read up to the requested location, then
	// a Future Bix.Get() call would have a new reader and it would
	// have to bgzf.Read() again through the block, past what we already
	// knew. We could do one big Read(), then just parse each time, but
	// this version actually works quite well.
	if len(tbx.cache) != 0 {
		return
	}
	var err error
	for err == nil {
		line, err := br.rdr.ReadBytes('\n')
		if err != nil && err != io.EOF {
			log.Println(err)
			return
		}
		if len(line) == 0 || (err == io.EOF && line[len(line)-1] != '\n') {
			return
		}
		toks := bytes.SplitN(line, []byte{'\t'}, br.maxCol)
		ip := tbx.toPosition(toks)
		tbx.cache = append(tbx.cache, ip)
	}

}

func (tbx *Bix) Get(q interfaces.IPosition) []interfaces.IPosition {
	overlaps := make([]interfaces.IPosition, 0)
	chunkReader, err := tbx.ChunkedReader(location{q.Chrom(), int(q.Start()), int(q.End())}, false)

	if err != nil {
		return overlaps
	}
	br := chunkReader.(*bixReader)

	var k int
	sizeCutoff := uint32(65536 / 2)

	if len(tbx.cache) == 0 && q.End()-q.Start() < sizeCutoff {
		tbx.fillCache(br)
	}
	for {
		if q.End()-q.Start() < sizeCutoff {
			if k >= len(tbx.cache) {
				break
			}
			v := tbx.cache[k]
			if v != nil {
				if interfaces.OverlapsPosition(v, q) {
					overlaps = append(overlaps, v)
				} else if v.Start() >= q.End() {
					break
				}
			}
			k += 1
		} else {
			line, err := br.rdr.ReadBytes('\n')

			if err != nil && err != io.EOF {
				log.Println(err)
				break
			}
			if len(line) == 0 || (err == io.EOF && line[len(line)-1] != '\n') {
				break
			}
			k += 1
			if err != nil {
				log.Println(err)
			}
			in, rerr, toks := br.inBounds(line)

			if !in {
				if rerr == io.EOF {
					err = io.EOF
					break
				}
				continue
			}
			ip := tbx.toPosition(toks)
			if ip != nil {
				overlaps = append(overlaps, ip)
			}
			if err == io.EOF {
				break
			}
		}
	}
	return overlaps
}

func (tbx *Bix) toPosition(toks [][]byte) interfaces.IPosition {
	isVCF := tbx.vReader != nil

	if isVCF {
		v := tbx.vReader.Parse(toks)
		return v
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

// Read from a BixReader (will contain only overlapping intervals).
func (r *bixReader) Read(p []byte) (int, error) {
	var readErr error
	var line []byte

	for len(r.buf) < len(p) {
		line, readErr = r.rdr.ReadBytes('\n')
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

func (tbx *Bix) ChunkedReader(r tabix.Record, justChunk bool) (io.ReadCloser, error) {
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

	if !chunksEqual(chunks, tbx.lastChunks) {
		(*tbx).lastChunks = chunks
		tbx.cache = tbx.cache[:0]
	}

	chunkReader, err := index.NewChunkReader(tbx.bgzf, chunks)
	if err != nil {
		return nil, err
	}
	if justChunk {

		if tbx.Header != "" {
			h := strings.NewReader(tbx.Header)
			b := bufio.NewReader(io.MultiReader(h, chunkReader))
			return frdr{b, h, chunkReader, tbx}, nil
		} else {
			return chunkReader, nil
		}
	}
	rdr, err := newReader(tbx, chunkReader, r.Start(), r.End())
	if err != nil {
		return nil, err
	}
	return rdr, nil
}

// Query a Bix struct with genomic coordinates. Returns an io.Reader.
func (tbx *Bix) Query(chrom string, start int, end int, printHeader bool) (io.ReadCloser, error) {

	rdr, err := (*tbx).ChunkedReader(location{chrom: chrom, start: start, end: end}, false)
	if err != nil {
		return nil, err
	}
	if printHeader {
		rdr.(*bixReader).buf = []byte(tbx.Header)
	}
	return rdr, err
}

func (tbx *Bix) AddInfoToHeader(id, number, vtype, desc string) {
	if tbx.vReader == nil {
		return
	}
	tbx.vReader.AddInfoToHeader(id, number, vtype, desc)
}

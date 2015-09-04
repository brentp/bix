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
	if err != nil {
		return nil, err
	}

	tbx := &Bix{bgzf: bgz, path: path}
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
			return 0, readErr
		}

		toks := strings.SplitN(line, "\t", r.maxCol+1)

		s, err := strconv.Atoi(toks[r.startCol])
		if err != nil {
			return 0, err
		}
		pos := s + r.add
		if pos >= r.end {
			readErr = io.EOF
			break
		}

		if r.endCol != -1 {
			e, err := strconv.Atoi(toks[r.endCol])
			if err != nil {
				return 0, err
			}
			if e < r.start {
				continue
			}
		} else if r.isVCF {
			alt := strings.Split(toks[4], ",")
			lref := len(toks[5])
			maxEnd := pos + 2
			anyLess := true
			if maxEnd <= r.start {
				anyLess = false
				for _, a := range alt {
					if a[0] != '<' {
						e := pos + lref
						if e >= r.start {
							anyLess = true
							break
						}
					} else if strings.HasPrefix(a, "<DEL") || strings.HasPrefix(a, "<DUP") || strings.HasPrefix(a, "<INV") || strings.HasPrefix(a, "<CN") {
						info := toks[7]
						if idx := strings.Index(info, ";END="); idx != -1 {
							v := info[idx+5 : idx+5+strings.Index(info[idx+5:], ";")]
							e, err := strconv.Atoi(v)
							if err != nil {
								return 0, err
							}
							if e > maxEnd {
								maxEnd = e
							}
						} else {
							log.Println("no end")
						}
					} else {
						if pos+lref > maxEnd {
							maxEnd = pos + lref
						}
					}
					if maxEnd > r.start {
						anyLess = true
						break
					}
				}
			}
			if !anyLess {
				continue
			}
		}

		r.buf = append(r.buf, line...)
	}
	if len(r.buf) >= len(p) {
		copy(p, r.buf[:len(p)])
		r.buf = r.buf[len(p):]
		return len(p), nil
	}
	copy(p, r.buf)
	return len(r.buf), readErr
}

// Query a Bix struct with genomic coordinates. Returns an io.Reader.
func (tbx *Bix) Query(chrom string, start int, end int) (io.Reader, error) {
	chunks, err := tbx.Chunks(location{chrom: chrom, start: start, end: end})
	if err != nil {
		return nil, err
	}

	chunkReader, err := index.NewChunkReader(tbx.bgzf, chunks)
	if err != nil {
		return nil, err
	}
	return newReader(tbx, chunkReader, start, end)
}

// Close closes the files associate with a Bix struct.
func (tbx *Bix) Close() {
	tbx.bgzf.Close()
}

// Tabix queries for go
package bix

import (
	"bytes"
	"encoding/binary"
	"io"
	"strings"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/csi"
	"github.com/biogo/hts/tabix"
)

// Index unifies CSI and tabix.
type Index interface {
	Chunks(string, int, int) ([]bgzf.Chunk, error)
	NameColumn() int
	BeginColumn() int
	EndColumn() int
	ZeroBased() bool
	MetaChar() rune
	Skip() int
}

type tIndex struct{ *tabix.Index }

func (t tIndex) NameColumn() int {
	return int(t.Index.NameColumn)
}
func (t tIndex) BeginColumn() int {
	return int(t.Index.BeginColumn)
}

func (t tIndex) EndColumn() int {
	return int(t.Index.EndColumn)
}

func (t tIndex) ZeroBased() bool {
	return t.Index.ZeroBased
}

func (t tIndex) MetaChar() rune {
	return t.Index.MetaChar
}

func (t tIndex) Skip() int {
	return int(t.Index.Skip)
}

type cIndex struct {
	*csi.Index
	chroms      []string
	nameColumn  int
	beginColumn int
	endColumn   int
	metaChar    rune
	zeroBased   bool
	skip        int
}

func (c cIndex) NameColumn() int {
	return c.nameColumn
}
func (c cIndex) BeginColumn() int {
	return c.beginColumn
}

func (c cIndex) EndColumn() int {
	return c.endColumn
}

func (c cIndex) ZeroBased() bool {
	return c.zeroBased
}

func (c cIndex) MetaChar() rune {
	return c.metaChar
}

func (c cIndex) Skip() int {
	return c.skip
}

// StripChr removes the "chr" prefix if it is present
func stripChr(c string) string {
	if strings.HasPrefix(c, "chr") {
		return c[3:]
	}
	return c
}

func (c cIndex) Chunks(chrom string, start int, end int) ([]bgzf.Chunk, error) {
	idx := -1
	chrom = stripChr(chrom)
	for i, ichrom := range c.chroms {
		if stripChr(ichrom) == chrom {
			idx = i
			break
		}

	}
	return c.Index.Chunks(idx, start, end), nil
}

func NewCSI(r io.Reader) (cIndex, error) {
	c, err := csi.ReadFrom(r)
	ci := cIndex{Index: c}
	if err != nil {
		return ci, err
	}
	aux := c.Auxilliary

	ci.nameColumn = int(binary.LittleEndian.Uint32(aux[4:8]))
	ci.beginColumn = int(binary.LittleEndian.Uint32(aux[8:12]))
	ci.endColumn = int(binary.LittleEndian.Uint32(aux[12:16]))
	ci.metaChar = rune(binary.LittleEndian.Uint32(aux[16:20]))
	ci.skip = int(binary.LittleEndian.Uint32(aux[20:24]))

	L := int(binary.LittleEndian.Uint32(aux[24:28]))
	// ends in trailing '0' so we take all but the last one byte.
	chroms := bytes.Split(aux[28:28+L-1], []byte{0})
	ci.chroms = make([]string, len(chroms))
	for i, chrom := range chroms {
		ci.chroms[i] = string(chrom)
	}
	return ci, nil
}

package bix

import (
	"compress/gzip"
	"os"
	"testing"

	. "gopkg.in/check.v1"
)

func Test(t *testing.T) { TestingT(t) }

type BixSuite struct{}

var _ = Suite(&BixSuite{})

func (s *BixSuite) TestReadCSI(c *C) {
	f, err := os.Open("tests/csitest.bed.gz.csi")
	if err != nil {
		c.Fatal(err)
	}
	g, err := gzip.NewReader(f)
	if err != nil {
		c.Fatal(err)
	}

	cs, err := NewCSI(g)
	if err != nil {
		c.Fatal(err)
	}
	c.Check(cs.NameColumn(), Equals, 1)
	c.Check(cs.BeginColumn(), Equals, 2)
	c.Check(cs.EndColumn(), Equals, 3)
	c.Check(cs.Skip(), Equals, 0)
	c.Check(cs.NumRefs(), Equals, len(cs.chroms))
	c.Check(cs.MetaChar(), Equals, '#')
	c.Check(cs.chroms, DeepEquals, []string{"1", "2", "3", "4", "5", "6", "9"})

}

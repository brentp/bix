package bix

import (
	"io"
	"testing"

	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/irelate/parsers"
)

func TestHeaderOnlyVCF(t *testing.T) {
	for _, kind := range []string{"tbi", "csi"} {
		t.Run(kind, func(t *testing.T) {
			b, err := New("tests/empty-" + kind + "/annotation.vcf.gz")
			if err != nil {
				t.Fatal(err)
			}
			defer b.Close()
			if b.VReader == nil || b.GetHeaderType("SCORE") != "Float" {
				t.Fatal("VCF header was not read")
			}
			region := parsers.NewInterval("chr1", 999, 1000, nil, 0, nil)
			queries := []struct {
				name string
				run  func() (interfaces.RelatableIterator, error)
			}{
				{"region", func() (interfaces.RelatableIterator, error) { return b.Query(region) }},
				{"all", func() (interfaces.RelatableIterator, error) { return b.Query(nil) }},
				{"fast", func() (interfaces.RelatableIterator, error) { return b.FastQuery(region) }},
			}
			for _, query := range queries {
				t.Run(query.name, func(t *testing.T) {
					iter, err := query.run()
					if err != nil {
						t.Fatal(err)
					}
					defer iter.Close()
					if _, err := iter.Next(); err != io.EOF {
						t.Fatalf("Next() error = %v, want EOF", err)
					}
				})
			}
		})
	}
}

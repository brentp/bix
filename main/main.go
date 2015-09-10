package main

import (
	"fmt"
	"log"
	"os"
	"strconv"

	"github.com/brentp/bix"
	"github.com/brentp/vcfgo"
)

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

type loc struct {
	chrom string
	start int
	end   int
}

func (s loc) Chrom() string {
	return s.chrom
}
func (s loc) Start() uint32 {
	return uint32(s.start)
}
func (s loc) End() uint32 {
	return uint32(s.end)
}

func main() {

	f := os.Args[1]
	tbx, err := bix.New(f)
	check(err)

	chrom := os.Args[2]

	s, err := strconv.Atoi(os.Args[3])
	check(err)

	e, err := strconv.Atoi(os.Args[4])
	check(err)

	vals := tbx.Get(loc{chrom, s, e})
	for i := 0; i < len(vals); i++ {
		fmt.Println(vals[i].(*vcfgo.Variant).String()[:40])
	}
	tbx.Close()

}

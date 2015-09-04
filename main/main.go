package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"

	"github.com/brentp/bix"
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

func main() {

	f := os.Args[1]
	tbx, err := bix.New(f)
	check(err)

	chrom := os.Args[2]

	s, err := strconv.Atoi(os.Args[3])
	check(err)

	e, err := strconv.Atoi(os.Args[4])
	check(err)

	rdr, err := tbx.Query(chrom, s, e, true)
	check(err)
	bufr := bufio.NewReader(rdr)
	for {
		v, err := bufr.ReadString('\n')
		if err == io.EOF {
			break
		}
		check(err)
		fmt.Println(v[:min(len(v), 20)])
	}
	tbx.Close()

}

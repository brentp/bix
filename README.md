bix
===

Tabix files in go.

[![GoDoc] (https://godoc.org/github.com/brentp/bix?status.png)](https://godoc.org/github.com/brentp/bix)


`Bix`is a pure go tabix reader. Tabix has a minimum resolution of 16K bases. So, for dense annotations
like ExAC, repeated queries will have to re-parse the same intervals. To mitigate this, `Bix` caches
all annotations for the blocks from each query and re-uses them for later querires to the same block.
In practice, this results in 10X-100X speedup for any type of data.




```go
tbx, err := bix.New(f)

// Query returns an io.Reader
ph := true // print header ?
rdr, err := tbx.Query(chrom, start, end, ph)
buf := bufio.NewReader(rdr)
for {
	line, err := bufr.ReadString('\n')
	if err == io.EOF {
		break
	}
	fmt.Println(line)
}
// or

intervals := tbx.Get(chrom, start, end)


```

where intervals will be either a vcfgo.Variant or an Interval object with a Chrom(), Start(), and End() method.

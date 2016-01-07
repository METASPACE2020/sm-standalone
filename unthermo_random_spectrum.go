package main

import (
	"bitbucket.org/lomereiter/ms/unthermo"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"time"
)

var fileName string

func init() {
	flag.StringVar(&fileName, "raw", "small.RAW", "name of the subject RAW file")
	flag.Parse()
}

func main() {
	seed := time.Now().UnixNano()
	r := rand.New(rand.NewSource(seed))

	file, err := unthermo.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	i := r.Intn(file.NScans())
	spectrum := file.Scan(i + 1).Spectrum()
	for _, peak := range spectrum {
		fmt.Println(peak.Mz, peak.I)
	}
}

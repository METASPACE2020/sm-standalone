package main

import (
	"bitbucket.org/lomereiter/ms/unthermo"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"time"
)

func main() {

	var fileName string
	var scannumber int

	flag.IntVar(&scannumber, "sn", -1, "the scan number")
	flag.StringVar(&fileName, "raw", "small.RAW", "name of the subject RAW file")
	flag.Parse()

	file, err := unthermo.Open(fileName)
	if scannumber == -1 {
		seed := time.Now().UnixNano()
		r := rand.New(rand.NewSource(seed))
		scannumber = r.Intn(file.NScans()) + 1
	}

	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	spectrum := file.Scan(scannumber).Spectrum()
	for _, peak := range spectrum {
		fmt.Println(peak.Mz, peak.I)
	}
}

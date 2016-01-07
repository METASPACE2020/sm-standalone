package main

import (
	"bitbucket.org/lomereiter/ms/unthermo"
	"flag"
	"fmt"
	"log"
)

var fileName string

func init() {
	flag.StringVar(&fileName, "raw", "small.RAW", "name of the subject RAW file")
	flag.Parse()
}

func main() {
	file, err := unthermo.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	spectrum := file.ComputeMeanSpectrum()
	for _, peak := range spectrum {
		fmt.Println(peak.Mz, peak.I)
	}
}

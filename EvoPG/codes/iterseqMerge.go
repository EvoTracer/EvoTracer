package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Printf("Lack of input file!\n")
		return
	}

	seqFileName := os.Args[1]
	seqFile, seqopenErr := os.Open(seqFileName)
	if seqopenErr != nil {
		fmt.Printf("Error happened in opening %s!\n", seqFileName)
		return
	}
	defer seqFile.Close()

	seqReader := bufio.NewReader(seqFile)
	for {
		line, readErr := seqReader.ReadString('\n')
		line = strings.Trim(line, "\r\n")
		line = strings.Trim(line, "\n")
		if strings.Contains(line, ">") {
			fmt.Print(line, "\n")
		} else if len(line) != 0 {
			fmt.Print(line, "\n")
		}
		if readErr == io.EOF {
			fmt.Print(line, "\n")
			break
		}
	}
}

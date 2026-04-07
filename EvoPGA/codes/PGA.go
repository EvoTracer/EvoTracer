// PGA.go
/* To annotate PGAG or RAST gene annotation table with Salmonella PG UID for each protein gene. Here, we use 26-genome Salmonella PG to make the annotation.
   Procedure: (1) Making the alignment between each proteome of the 26 Salmonella strains and the proteome to be annotated;
              (2) Finding the mutual best alignments for each strain and the one to be annotated;
              (3) Using "26_PG.txt" to annotate the genes of targeted strains with PG UIDs; Some genes could not be annotated because of no PG cluster was
                  mapped or pseudogenes / non-coding RNAs.
   USAGE:     Copy the compiled program, "26_PG.txt" and "S.genus.AG_PGAG.gene.tab.txt" or "Salmonella.genus.ancient.orthologous.genome_RAST.tab.txt" into
              ~/BactCG1.0/results/out_mutbest_filt/, running the command below.

   Command:   ./PGA  <26_PG.txt>    <PGAG|RAST_tab.txt>   <PGAG|RAST>   <Strain_to_be_annotated>
*/

package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

func main() {
	if len(os.Args) < 5 {
		fmt.Printf("Errors in inputing files!\n")
		return
	}
	pgmapmatFileName := os.Args[1]
	pgmapFile, openErr1 := os.Open(pgmapmatFileName)
	if openErr1 != nil {
		fmt.Printf("Error happened in opening %s!\n", pgmapmatFileName)
		return
	}
	defer pgmapFile.Close()

	pgmapReader := bufio.NewReader(pgmapFile)
	pg := make(map[string]string)
	for {
		pgmap, readErr1 := pgmapReader.ReadString('\n')
		pgmap = seqTrim(pgmap)
		if len(pgmap) > 0 {
			pgmapArray := strings.Split(pgmap, "\t")
			//fmt.Printf("%d\n",len(pgmapArray))
			for i := 2; i < len(pgmapArray)-1; i++ {
				if !(pgmapArray[i] == "-") {
					pg[pgmapArray[i]] = pgmapArray[0]
					//fmt.Printf("%s\t%s\n",pgmapArray[i],pg[pgmapArray[i]])
				}
			}
		}
		if readErr1 == io.EOF {
			break
		}
	}

	var flag int
	if strings.Contains(os.Args[3], "PGAG") {
		flag = 2
	} else if strings.Contains(os.Args[3], "RAST") {
		flag = 1
	}

	strain := os.Args[4]

	annMap := make(map[string]string)

	dir := "/home/yaozikun/Public_Dir/Bact_PGA/BactPGA/test/seq"
	entries := os.ReadDir(dir)

	count := 0
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()
		if filepath.Ext(name) != ".faa" {
			continue
		}
		base := name[:len(name)-4]
		if _, err := strconv.Atoi(base); err == nil {
			count++
		}
	}

	for i := 1; i <= count; i++ {
		mutbestFileName := strain + ".vs." + strconv.Itoa(i) + ".mutbest.filt.txt"
		mutbestFile, openErr2 := os.Open(mutbestFileName)
		if openErr2 != nil {
			fmt.Printf("Error happened in opening %s!\n", mutbestFileName)
			return
		}
		defer mutbestFile.Close()

		mutbestReader := bufio.NewReader(mutbestFile)

		for {
			mutbest, readErr2 := mutbestReader.ReadString('\n')
			mutbest = seqTrim(mutbest)
			if len(mutbest) > 0 {
				mutbestArray := strings.Split(mutbest, "\t")
				if !(len(annMap[mutbestArray[0]]) > 0) {
					if len(pg[mutbestArray[1]]) > 0 {
						annMap[mutbestArray[0]] = pg[mutbestArray[1]]
						//fmt.Printf("%s\t%s\n",mutbestArray[0],annMap[mutbestArray[0]])
					}
				}

			}
			if readErr2 == io.EOF {
				break
			}
		}
	}

	agFileName := os.Args[2]
	agFile, openErr3 := os.Open(agFileName)
	if openErr3 != nil {
		fmt.Printf("Error happened in opening %s!\n", agFileName)
		return
	}
	defer agFile.Close()

	agReader := bufio.NewReader(agFile)
	for {
		ag, readErr3 := agReader.ReadString('\n')
		ag = seqTrim(ag)
		if len(ag) > 0 {
			agArray := strings.Split(ag, "\t")
			for i := 0; i <= flag; i++ {
				fmt.Printf("%s\t", agArray[i])
			}
			if len(annMap[agArray[flag]]) > 0 {
				fmt.Printf("%s", annMap[agArray[flag]])
			}
			for i := flag + 1; i < len(agArray); i++ {
				fmt.Printf("\t%s", agArray[i])
			}
			fmt.Printf("\n")
		}
		if readErr3 == io.EOF {
			break
		}
	}
}

func seqTrim(s string) string {
	s = strings.Trim(s, "\r\n")
	s = strings.Trim(s, "\n")
	//   s = strings.Replace(s," ","",-1)
	return s
}

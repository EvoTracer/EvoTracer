//gbParse.go
//The program is used to retrieve all the annotated protein encoding genes (GeneID, GI_acc, Protein_ID), genome coordinate, strand, and protein sequences.

package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"regexp"
	"strings"
)

func main() {
	if len(os.Args) <= 1 {
		fmt.Printf("Please input the gb file that is going to be parsed!\n")
		return
	}
	gbFileName := os.Args[1]
	gbFile, openErr := os.Open(gbFileName)
	if openErr != nil {
		fmt.Printf("Error happened in opening %s!\n", gbFileName)
		return
	}
	defer gbFile.Close()

	inputReader := bufio.NewReader(gbFile)
	var seq, molType, geneAcc, gene, locusTag, product, geneID string
	var flag, flag0, flag2 int
	//	var gi,gene,prot,geneStart,geneEnd,strand string
	for {
		line, readErr := inputReader.ReadString('\n')
		pat01, _ := regexp.MatchString("  rRNA  ", line)
		pat02, _ := regexp.MatchString("  tRNA  ", line)
		pat03, _ := regexp.MatchString("  misc_feature  ", line)
		pat04, _ := regexp.MatchString("  ncRNA  ", line)
		patA, _ := regexp.MatchString("/gene=", line)
		patB, _ := regexp.MatchString("/locus_tag=", line)
		patC, _ := regexp.MatchString("/product=", line)
		patD, _ := regexp.MatchString("/EC_number=", line)
		patF, _ := regexp.MatchString("                     /pseudo", line)

		pat1, _ := regexp.MatchString("  CDS  ", line)
		patE, _ := regexp.MatchString("/translation=", line)
		//pat3, _ := regexp.MatchString("/db_xref=\"GeneID:",line)
		//pat4, _ := regexp.MatchString("/db_xref=\"GI:",line)
		//pat5, _ := regexp.MatchString("/protein_id=\"",line)
		if pat01 {
			molType = "rRNA"
			line = strings.Replace(line, "rRNA", "", 1)
			geneAcc = seqTrim(line)
			flag0 = 1
			gene = ""
			geneID = ""
			locusTag = ""
			product = ""
		} else if pat02 {
			molType = "tRNA"
			line = strings.Replace(line, "tRNA", "", 1)
			geneAcc = seqTrim(line)
			flag0 = 1
			gene = ""
			geneID = ""
			locusTag = ""
			product = ""
		} else if pat03 {
			molType = "misc_feature"
			line = strings.Replace(line, "misc_feature", "", 1)
			geneAcc = seqTrim(line)
			flag0 = 1
			gene = ""
			geneID = ""
			locusTag = ""
			product = ""
		} else if pat04 {
			molType = "ncRNA"
			line = strings.Replace(line, "ncRNA", "", 1)
			geneAcc = seqTrim(line)
			flag0 = 1
			gene = ""
			geneID = ""
			locusTag = ""
			product = ""
		} else if flag0 == 1 {
			if patA {
				gene = strings.Replace(line, "/gene=\"", "", -1)
				gene = strings.Replace(gene, "\"", "", -1)
				gene = seqTrim(gene)
			}
			if patB {
				locusTag = strings.Replace(line, "/locus_tag=\"", "", -1)
				locusTag = strings.Replace(locusTag, "\"", "", -1)
				locusTag = seqTrim(locusTag)
			}
			if patC {
				product = strings.Replace(line, "                     /product=\"", "", -1)
				product = strings.Replace(product, "\"", "", -1)
				product = strings.Replace(product, "\n", "", -1)
				flag0 = 0
				fmt.Printf("%s\t%s\t%s\t%s\t%s\t%s\t\n", geneAcc, molType, locusTag, geneID, gene, product)
				gene = ""
				geneID = ""
				locusTag = ""
				product = ""
			}
		}
		if pat1 {
			molType = "CDS"
			line = strings.Replace(line, "CDS", "", 1)
			line = seqTrim(line)
			geneAcc = seqTrim(line)
			flag0 = 0
			flag = 1
			gene = ""
			geneID = ""
			locusTag = ""
			product = ""
		} else if flag == 1 {
			if patA {
				gene = strings.Replace(line, "/gene=\"", "", -1)
				gene = strings.Replace(gene, "\"", "", -1)
				gene = seqTrim(gene)
			}
			if patB {
				locusTag = strings.Replace(line, "/locus_tag=\"", "", -1)
				locusTag = strings.Replace(locusTag, "\"", "", -1)
				locusTag = seqTrim(locusTag)
			}
			if patC {
				product = strings.Replace(line, "                     /product=\"", "", -1)
				product = strings.Replace(product, "\"", "", -1)
				product = strings.Replace(product, "\n", "", -1)
			}
			if patD {
				geneID = strings.Replace(line, "/EC_number=\"", "", -1)
				geneID = strings.Replace(geneID, "\"", "", -1)
				geneID = seqTrim(geneID)
			}
			if patE {
				line = strings.Replace(line, "/translation=\"", "", -1)
				line = seqTrim(line)
				flag = 0
				if strings.Contains(line, "\"") {
					seq = strings.Replace(line, "\"", "", 1)
					fmt.Printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", geneAcc, molType, locusTag, geneID, gene, product, seq)
					gene = ""
					geneID = ""
					locusTag = ""
					product = ""
					seq = ""
				} else {
					flag2 = 1
					seq = line
				}
			}
			if patF {
				fmt.Printf("%s\t%s\t%s\t%s\t%s\tpseudo\n", geneAcc, molType, locusTag, geneID, gene)
				gene = ""
				geneID = ""
				locusTag = ""
				product = ""
				seq = ""
				flag = 0
			}
		} else if flag2 == 1 {
			line = seqTrim(line)
			if strings.Contains(line, "\"") {
				line = strings.Replace(line, "\"", "", -1)
				seq += line
				fmt.Printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", geneAcc, molType, locusTag, geneID, gene, product, seq)
				flag2 = 0
				gene = ""
				geneID = ""
				locusTag = ""
				product = ""
				seq = ""
			} else {
				seq += line
			}
		}
		if readErr == io.EOF {
			break
		}
	}
}

func seqTrim(s string) string {
	s = strings.Trim(s, "\r\n")
	s = strings.Trim(s, "\n")
	s = strings.Replace(s, " ", "", -1)
	return s
}

/* func substrPrint (s string, l int) string {
   var s0 string
   n := (len(s) - len(s)%l)/l
   for i:=0;i<n;i++ {
   	  s0 += s[i*l:(i+1)*l] + "\n"
   }

   if len(s)%l != 0 {
   	  s0 += s[l*n:len(s)] + "\n"
   }
   return s0
} */

package cmd

import (
"bufio"
"fmt"
"os"
"path/filepath"
"strings"
"github.com/spf13/cobra"
)

var (
inputCdhitRun   string
inputCgFileRun  string
baseOutDirRun   string
)

var getFaCmd = &cobra.Command{
Use:   "getfa",
Short: "Extract sequences from CD-HIT fasta to independent .fa files per family",
Run: func(cmd *cobra.Command, args []string) {
err := GetFaFileAll(inputCdhitRun, inputCgFileRun, baseOutDirRun)
if err != nil {
fmt.Println("Error:", err)
}
},
}

func init() {
rootCmd.AddCommand(getFaCmd)
getFaCmd.Flags().StringVarP(&inputCdhitRun, "cdhit", "c", "output/CG_results/1.cd-hit_output/1.cd-hit_fatsa", "Directory with cdhit FASTA")
getFaCmd.Flags().StringVarP(&inputCgFileRun, "cgfile", "g", "output/CG_results/CG_ALL.txt", "CG core gene table file")
getFaCmd.Flags().StringVarP(&baseOutDirRun, "outdir", "o", "output/CG_results/all-strain-together", "Base output directory")
}

func GetFaFileAll(inputCdhitDir, cgFile, baseOutDir string) error {
fmt.Println("Running Post-process Step 2: get_fa_file-all...")
outDir := filepath.Join(baseOutDir, "2.result")
os.MkdirAll(outDir, 0755)

seqDict := make(map[string]string)
files, err := os.ReadDir(inputCdhitDir)
if err != nil {
return fmt.Errorf("failed to read cdhit dir: %v", err)
}

for _, f := range files {
if !f.IsDir() && strings.HasSuffix(f.Name(), ".fasta") {
indexName := strings.TrimSuffix(f.Name(), ".fasta")
path := filepath.Join(inputCdhitDir, f.Name())
file, err := os.Open(path)
if err != nil {
continue
}
scanner := bufio.NewScanner(file)
var currentSeqName string
for scanner.Scan() {
line := scanner.Text()
if strings.HasPrefix(line, ">") {
currentSeqName = ">" + indexName + "---" + strings.TrimPrefix(line, ">")
seqDict[currentSeqName] = ""
} else if line != "" {
seqDict[currentSeqName] += strings.TrimSpace(line)
}
}
file.Close()
}
}

cgLineOutput := make([][]string, 0)
cgf, err := os.Open(cgFile)
if err != nil {
return fmt.Errorf("failed to open CG file: %v", err)
}
defer cgf.Close()

scanner := bufio.NewScanner(cgf)
if !scanner.Scan() {
return nil
}
firstLine := scanner.Text()
headers := strings.Split(firstLine, "\t")
var indexList []string
for _, h := range headers {
indexList = append(indexList, strings.TrimSpace(h)+"---")
}

for scanner.Scan() {
line := scanner.Text()
if line == "" {
continue
}
parts := strings.Split(line, "\t")
var outLine []string
for i := 0; i < len(parts) && i < len(indexList); i++ {
outLine = append(outLine, indexList[i]+strings.TrimSpace(parts[i]))
}
cgLineOutput = append(cgLineOutput, outLine)
}

for _, row := range cgLineOutput {
if len(row) == 0 {
continue
}
outFileName := strings.TrimPrefix(row[0], indexList[0])
outPath := filepath.Join(outDir, outFileName+".fa")
outf, err := os.Create(outPath)
if err != nil {
continue
}
writer := bufio.NewWriter(outf)
for _, el := range row {
key := ">" + el
if val, ok := seqDict[key]; ok {
writer.WriteString(key + "\n" + val + "\n")
}
}
writer.Flush()
outf.Close()
}
return nil
}

func GetSnpMegaAll(inputMegaDir, outputDir string) error {
fmt.Println("Running Post-process Step 4: get_SNP_mega-all...")
os.MkdirAll(outputDir, 0755)

files, err := os.ReadDir(inputMegaDir)
if err != nil {
return fmt.Errorf("failed to read mega dir: %v", err)
}

seqs := make(map[string]string)
var firstLen int
var indexList []string

for _, f := range files {
if !strings.HasSuffix(f.Name(), ".meg") {
continue
}

path := filepath.Join(inputMegaDir, f.Name())
file, err := os.Open(path)
if err != nil {
continue
}

scanner := bufio.NewScanner(file)
isFirst := (len(indexList) == 0)

for scanner.Scan() {
line := strings.TrimSpace(scanner.Text())
if line == "" || strings.HasPrefix(line, "#") {
continue
}

if strings.HasPrefix(line, ">") {
name := strings.TrimPrefix(line, ">")
if isFirst {
indexList = append(indexList, name)
}
}
}
file.Close()
break // just need the index list from first file
}

for _, id := range indexList {
seqs[id] = ""
}

for _, f := range files {
if !strings.HasSuffix(f.Name(), ".meg") {
continue
}

path := filepath.Join(inputMegaDir, f.Name())
file, err := os.Open(path)
if err != nil {
continue
}

scanner := bufio.NewScanner(file)
var currentSeqId string
for scanner.Scan() {
line := strings.TrimSpace(scanner.Text())
if line == "" || strings.HasPrefix(line, "#") {
continue
}
if strings.HasPrefix(line, ">") {
currentSeqId = strings.TrimPrefix(line, ">")
} else {
if _, ok := seqs[currentSeqId]; ok {
seqs[currentSeqId] += line
}
}
}
file.Close()
}

if len(indexList) > 0 {
firstLen = len(seqs[indexList[0]])
}

outPath := filepath.Join(outputDir, "all_core_gene.meg")
outf, err := os.Create(outPath)
if err != nil {
return err
}
defer outf.Close()

writer := bufio.NewWriter(outf)
writer.WriteString(fmt.Sprintf("#MEGA\n!Title SNPMega;\n!Format DataType=DNA indel=-;\n\n"))

for _, id := range indexList {
writer.WriteString(fmt.Sprintf("#%s\n", id))
writer.WriteString(seqs[id] + "\n")
}
writer.Flush()

fmt.Printf("Step 4 output to %s with length %d\n", outPath, firstLen)
return nil
}

#!/usr/bin/env python3
"""
gbParse.py

Extract annotated protein-coding genes (CDS), rRNAs, tRNAs, ncRNAs, and misc_feature records from a GenBank file, including:

- genomic coordinates (geneAcc)
- feature type (molType)
- locus_tag
- GeneID (taken from /EC_number qualifier)
- gene name
- product description
- protein sequence (for CDS only)

Pseudo-genes (marked with /pseudo) are also reported.

Usage:
    python gbParse.py <genbank.gb> [rootDir]
"""

import os
import re
import sys

def seq_trim(s: str) -> str:
    return s.replace(" ", "").replace("\n", "").replace("\r", "")

def parse_genbank(gb_path: str):
    mol_type = ""
    gene_acc = ""
    gene = ""
    locus_tag = ""
    product = ""
    gene_id = ""
    seq = ""
    flag = 0  
    flag0 = 0  
    flag2 = 0  

    with open(gb_path, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')

            if re.search(r'\s+CDS\s+', line):
                mol_type = "CDS"
                gene_acc = seq_trim(line.replace("CDS", ""))
                flag = 1
                flag0 = 0
                gene = ""
                locus_tag = ""
                product = ""
                gene_id = ""
                seq = ""
            elif re.search(r'\s+rRNA\s+', line):
                mol_type = "rRNA"
                gene_acc = seq_trim(line.replace("rRNA", ""))
                flag0 = 1
                flag = 0
                gene = ""
                locus_tag = ""
                product = ""
                gene_id = ""
            elif re.search(r'\s+tRNA\s+', line):
                mol_type = "tRNA"
                gene_acc = seq_trim(line.replace("tRNA", ""))
                flag0 = 1
                flag = 0
                gene = ""
                locus_tag = ""
                product = ""
                gene_id = ""
            elif re.search(r'\s+misc_feature\s+', line):
                mol_type = "misc_feature"
                gene_acc = seq_trim(line.replace("misc_feature", ""))
                flag0 = 1
                flag = 0
                gene = ""
                locus_tag = ""
                product = ""
                gene_id = ""
            elif re.search(r'\s+ncRNA\s+', line):
                mol_type = "ncRNA"
                gene_acc = seq_trim(line.replace("ncRNA", ""))
                flag0 = 1
                flag = 0
                gene = ""
                locus_tag = ""
                product = ""
                gene_id = ""

            if flag0 == 1:
                if line.startswith("                     /gene="):
                    gene = line.split('"')[1]
                elif line.startswith("                     /locus_tag="):
                    locus_tag = line.split('"')[1]
                elif line.startswith("                     /product="):
                    product = line.split('"')[1]
                    flag0 = 0
                    print(f"{gene_acc}\t{mol_type}\t{locus_tag}\t{gene_id}\t{gene}\t{product}")

            elif flag == 1:
                if line.startswith("                     /gene="):
                    gene = line.split('"')[1]
                elif line.startswith("                     /locus_tag="):
                    locus_tag = line.split('"')[1]
                elif line.startswith("                     /product="):
                    product = line.split('"')[1]
                elif line.startswith("                     /EC_number="):
                    gene_id = line.split('"')[1]
                elif line.startswith("                     /pseudo"):
                    print(f"{gene_acc}\t{mol_type}\t{locus_tag}\t{gene_id}\t{gene}\tpseudo")
                    flag = 0
                    gene = ""
                    locus_tag = ""
                    product = ""
                    gene_id = ""
                    seq = ""
                elif line.startswith("                     /translation="):
                    seq = line.split('"')[1]
                    if seq.endswith('"'):
                        seq = seq[:-1]
                        print(f"{gene_acc}\t{mol_type}\t{locus_tag}\t{gene_id}\t{gene}\t{product}\t{seq}")
                        flag = 0
                        seq = ""
                    else:
                        flag2 = 1
                elif flag2 == 1:
                    if line.endswith('"'):
                        seq += line[:-1]
                        print(f"{gene_acc}\t{mol_type}\t{locus_tag}\t{gene_id}\t{gene}\t{product}\t{seq}")
                        flag2 = 0
                        flag = 0
                        seq = ""
                    else:
                        seq += seq_trim(line)

def main():
    root = "/home/yaozikun/Public_Dir/Bact_PGA/BactPGA/input"
    if len(sys.argv) >= 3:
        root = sys.argv[2]

    gb_path = os.path.join(root, sys.argv[1])
    if not os.path.isfile(gb_path):
        print(f"Error: Cannot open {gb_path}")
        sys.exit(1)

    parse_genbank(gb_path)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python gbParse.py <genbank.gb> [rootDir]")
        sys.exit(1)
    main()
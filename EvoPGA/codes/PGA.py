#!/usr/bin/env python3
"""
PGA.py
======

Python reimplementation of PGA.go

Purpose:
    Annotate PGAG or RAST gene annotation table with Salmonella PG UID for each protein gene.
    Uses 26-genome Salmonella PG for annotation.

Procedure:
    1. Align each proteome of 26 Salmonella strains with the target strain.
    2. Find mutual best alignments.
    3. Use "26_PG.txt" to annotate genes with PG UIDs.

Usage:
    python PGA.py <26_PG.txt> <PGAG_or_RAST_tab.txt> PGAG|RAST <StrainName> [rootDir]

Example:
    python PGA.py 26_PG.txt S.genus.AG_PGAG.gene.tab.txt PGAG StrainName
"""

import os
import sys 

def seq_trim(s: str) -> str:
    return s.strip('\r\n')

def load_pg_map(pg_file: str) -> dict:
    pg = {}
    with open(pg_file, 'r') as f:
        for line in f:
            line = seq_trim(line)
            if not line:
                continue
            parts = line.split('\t')
            pg_uid = parts[0]
            for i in range(2, len(parts) - 1):
                gene = parts[i]
                if gene != '-':
                    pg[gene] = pg_uid
    return pg

def load_mutbest_files(strain: str, count: int, root: str) -> dict:
    ann_map = {}
    for i in range(1, count + 1):
        filename = f"{strain}.vs.{i}.best.filt.txt"
        filepath = os.path.join(root, filename)
        if not os.path.isfile(filepath):
            print(f"Error: Cannot open {filepath}")
            sys.exit(1)

        with open(filepath, 'r') as f:
            for line in f:
                line = seq_trim(line)
                if not line:
                    continue
                parts = line.split('\t')
                query = parts[0]
                subject = parts[1]
                if query not in ann_map:
                    ann_map[query] = subject
    return ann_map

def count_faa_files(seq_dir: str) -> int:
    count = 0
    for fname in os.listdir(seq_dir):
        if fname.endswith('.faa'):
            try:
                int(fname[:-4])
                count += 1
            except ValueError:
                continue
    return count

def annotate_and_print(ag_file: str, ann_map: dict, pg: dict, flag: int):
    """Read annotation file and print annotated output."""
    with open(ag_file, 'r') as f:
        for line in f:
            line = seq_trim(line)
            if not line:
                continue
            parts = line.split('\t')
            gene = parts[flag]
            pg_uid = pg.get(ann_map.get(gene), '')

            print('\t'.join(parts[:flag + 1]), end='')
            print(f"\t{pg_uid}", end='')
            if len(parts) > flag + 1:
                print(f"\t" + '\t'.join(parts[flag + 1:]))
            else:
                print()

def main():
    if len(sys.argv) < 5:
        print("Usage: python PGA.py <255_PG.txt> <PGAG_or_RAST_tab.txt> PGAG|RAST <StrainName> [rootDir]")
        sys.exit(1)

    root = "/home/yaozikun/Public_Dir/Bact_PGA/BactPGA/input"
    if len(sys.argv) >= 6:
        root = sys.argv[5]

    pg_file = os.path.join(root, sys.argv[1])
    ag_file = os.path.join(root, sys.argv[2])
    strain = sys.argv[4]
    mode = sys.argv[3].upper()

    if mode == "PGAG":
        flag = 2
    elif mode == "RAST":
        flag = 1
    else:
        print("Third argument must be PGAG or RAST")
        sys.exit(1)

    seq_dir = os.path.join(root, "seq")
    if not os.path.isdir(seq_dir):
        print(f"Error: seq directory not found at {seq_dir}")
        sys.exit(1)

    count = count_faa_files(seq_dir)

    pg = load_pg_map(pg_file)
    ann_map = load_mutbest_files(strain, count, root)
    annotate_and_print(ag_file, ann_map, pg, flag)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Batch replace text in files under a directory using a mapping list.

Mapping file format (one rule per line):
1) Tab-separated: old<TAB>new
2) Comma-separated: old,new

Notes:
- Empty lines and lines starting with # are ignored.
- Rules are applied in order.
- By default, only UTF-8 text files are processed.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import List, Sequence, Tuple


Rule = Tuple[str, str]

# Direct configuration: edit these values and run the script directly.
TARGET_DIR = "/home/wangyiwei/1DGR/MyResults/Final_Results"
MAPPING_FILE = "/home/wangyiwei/1DGR/OneDGR/replace_rules.example.txt"
ENCODING = "utf-8"
ERRORS = "strict"  # strict | ignore | replace
INCLUDE_EXT: List[str] | None = [".txt"]  # None means all files
EXCLUDE_DIRS: List[str] = [".git"]
BACKUP = False
DRY_RUN = False


def load_rules(mapping_file: Path, encoding: str, errors: str) -> List[Rule]:
    rules: List[Rule] = []
    with mapping_file.open("r", encoding=encoding, errors=errors) as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            if "\t" in line:
                old, new = line.split("\t", 1)
            elif "," in line:
                old, new = line.split(",", 1)
            else:
                parts = line.split(maxsplit=1)
                if len(parts) == 2:
                    old, new = parts
                else:
                    raise ValueError(
                        f"Invalid rule at line {line_no}: use TAB/comma/or space separator"
                    )

            old = old.strip()
            new = new.strip()

            if old == "":
                raise ValueError(f"Invalid rule at line {line_no}: old text is empty")

            rules.append((old, new))

    if not rules:
        raise ValueError("No valid replacement rules found in mapping file")
    return rules


def is_binary_file(path: Path) -> bool:
    try:
        with path.open("rb") as f:
            sample = f.read(2048)
        return b"\x00" in sample
    except OSError:
        return True


def should_process(path: Path, include_ext: Sequence[str] | None) -> bool:
    if include_ext is None:
        return True
    return path.suffix in include_ext


def apply_rules(content: str, rules: Sequence[Rule]) -> Tuple[str, int]:
    updated = content
    total_hits = 0
    for old, new in rules:
        hits = updated.count(old)
        if hits:
            updated = updated.replace(old, new)
            total_hits += hits
    return updated, total_hits


def normalize_extensions(exts: Sequence[str] | None) -> List[str] | None:
    if exts is None:
        return None
    normalized: List[str] = []
    for ext in exts:
        normalized.append(ext if ext.startswith(".") else f".{ext}")
    return normalized


def main() -> int:
    target_dir = Path(TARGET_DIR).resolve()
    mapping_file = Path(MAPPING_FILE).resolve()

    if not target_dir.is_dir():
        raise FileNotFoundError(f"Target directory not found: {target_dir}")
    if not mapping_file.is_file():
        raise FileNotFoundError(f"Mapping file not found: {mapping_file}")

    include_ext = normalize_extensions(INCLUDE_EXT)
    exclude_dirs = set(EXCLUDE_DIRS)
    rules = load_rules(mapping_file, ENCODING, ERRORS)

    scanned_files = 0
    changed_files = 0
    skipped_binary = 0
    skipped_decode = 0
    total_replacements = 0

    for path in target_dir.rglob("*"):
        if not path.is_file():
            continue
        if any(part in exclude_dirs for part in path.parts):
            continue
        if path == mapping_file:
            continue
        if not should_process(path, include_ext):
            continue

        scanned_files += 1

        if is_binary_file(path):
            skipped_binary += 1
            continue

        try:
            original = path.read_text(encoding=ENCODING, errors=ERRORS)
        except UnicodeDecodeError:
            skipped_decode += 1
            continue

        updated, hits = apply_rules(original, rules)
        if hits == 0:
            continue

        changed_files += 1
        total_replacements += hits

        if DRY_RUN:
            print(f"[DRY-RUN] {path} | replacements={hits}")
            continue

        if BACKUP:
            backup_path = path.with_suffix(path.suffix + ".bak")
            shutil.copy2(path, backup_path)

        path.write_text(updated, encoding=ENCODING, errors=ERRORS)
        print(f"[UPDATED] {path} | replacements={hits}")

    print("\nSummary")
    print(f"- scanned_files: {scanned_files}")
    print(f"- changed_files: {changed_files}")
    print(f"- total_replacements: {total_replacements}")
    print(f"- skipped_binary: {skipped_binary}")
    print(f"- skipped_decode_error: {skipped_decode}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

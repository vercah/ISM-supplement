#!/usr/bin/env python3

from pathlib import Path


FASTA_SUFFIXES = (
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fa",
    ".fasta",
    ".fna",
)


def canonical_sample_name(value):
    name = Path(value.strip()).name
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem

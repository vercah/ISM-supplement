#!/usr/bin/env python3

import argparse
import random


def parse_args():
    parser = argparse.ArgumentParser(
        description="Deterministically shuffle input lines using a seed."
    )
    parser.add_argument("--seed", type=int, required=True, help="Base integer seed.")
    parser.add_argument("input", help="Input text file.")
    parser.add_argument("output", help="Output text file.")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.input, "r", encoding="utf-8", newline="") as infile:
        lines = infile.readlines()

    rng = random.Random(args.seed)
    rng.shuffle(lines)

    with open(args.output, "w", encoding="utf-8", newline="") as outfile:
        outfile.writelines(lines)


if __name__ == "__main__":
    main()

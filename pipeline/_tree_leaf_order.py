#!/usr/bin/env python3

import argparse
import re


def extract_newick_leaf_order(in_tree_fn):
    with open(in_tree_fn) as f:
        nw = f.read()

    leaves = re.findall(r"(?<=[(,])\s*([^():,;\s]+)\s*:", nw)

    assert len(leaves) > 0, f"No leaves found in Newick file: {in_tree_fn}"

    for leaf in leaves:
        assert leaf != "", f"Empty leaf name found in Newick file: {in_tree_fn}"
        assert leaf != "merge_root", (
            f"Invalid leaf name 'merge_root' found in Newick file: {in_tree_fn}"
        )
        print(leaf)


def main():
    parser = argparse.ArgumentParser(
        description="Extract left-to-right leaf order from a Newick tree."
    )
    parser.add_argument("in_tree_fn", metavar="input.tree.nw")

    args = parser.parse_args()
    extract_newick_leaf_order(args.in_tree_fn)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import os
import re
import sys


HEADER = [
    "dataset",
    "method",
    "k",
    "N",
    "type",
    "n_bit_changes",
    "n_runs",
    "n_cylinder_runs",
]
RUN_KEYS = ["n_bit_changes", "n_runs", "n_cylinder_runs"]
RUN_PATH_RE = re.compile(
    r"^10_runs/(?P<dataset>.+)_(?P<method>[^_]+)_k(?P<k>\d+)_N(?P<N>\d+)_(?P<type>[^_/]+)\.runs$"
)


class InputError(ValueError):
    pass


def parse_run_path(path):
    normalized = os.path.normpath(path)
    match = RUN_PATH_RE.match(normalized)
    if not match:
        raise InputError(
            f"Malformed .runs filename {path!r}; expected "
            "10_runs/{dataset}_{method}_k{k}_N{N}_{type}.runs"
        )
    return match.groupdict()


def parse_run_file(path):
    values = {}
    with open(path) as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                raise InputError(f"{path}:{line_number}: empty lines are not allowed")

            fields = line.split("\t")
            if len(fields) != 2:
                raise InputError(
                    f"{path}:{line_number}: expected two tab-separated fields"
                )

            key, value = fields
            if key not in RUN_KEYS:
                raise InputError(f"{path}:{line_number}: unknown key {key!r}")
            if key in values:
                raise InputError(f"{path}:{line_number}: duplicate key {key!r}")
            if not value.isdigit():
                raise InputError(
                    f"{path}:{line_number}: value for {key!r} must be an integer"
                )
            values[key] = value

    missing = [key for key in RUN_KEYS if key not in values]
    if missing:
        raise InputError(f"{path}: missing keys: {', '.join(missing)}")

    return values


def write_tsv(paths):
    rows = []
    for path in sorted(paths):
        metadata = parse_run_path(path)
        values = parse_run_file(path)
        row = [metadata[column] for column in HEADER[:5]]
        row.extend(values[key] for key in RUN_KEYS)
        rows.append(row)

    print("\t".join(HEADER))
    for row in rows:
        print("\t".join(row))


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate final .runs files into a single TSV."
    )
    parser.add_argument("runs", nargs="+", help="10_runs/*.runs files")
    args = parser.parse_args()

    try:
        write_tsv(args.runs)
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except InputError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

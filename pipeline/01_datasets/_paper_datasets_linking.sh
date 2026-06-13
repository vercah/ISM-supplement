#! /usr/bin/env bash

set -e
set -o pipefail
set -u

./link_genomes.sh ngono
./link_genomes.sh ngono-spneumo
./link_genomes.sh diverse

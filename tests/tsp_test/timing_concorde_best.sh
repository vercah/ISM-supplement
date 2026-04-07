#! /usr/bin/env bash

set -e
set -o pipefail
set -u

d=$(mktemp -d)
p=$(pwd)

(cd $d; time $p/../../bin/concorde $p/instance_ngono_best/instance.tsp)


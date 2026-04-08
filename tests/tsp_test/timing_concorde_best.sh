#! /usr/bin/env bash

set -e
set -o pipefail
set -u

d=$(mktemp -d)
p=$(pwd)

#time ../../bin/concorde instance_ngono_best/instance.tsp
(cd $d; time $p/../../bin/concorde $p/instance_ngono_best/instance.tsp)


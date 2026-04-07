#! /usr/bin/env bash

set -e
set -o pipefail
set -u

time ../../bin/concorde instance_ngono_best/instance.tsp


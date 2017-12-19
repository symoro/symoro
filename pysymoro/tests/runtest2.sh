#!/bin/sh

DIR="$(dirname $0)"

PYTHONPATH="$DIR"/../..:"$PYTHONPATH" python2 "$DIR"/run_all_tests.py


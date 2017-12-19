#!/bin/sh

DIR="$(dirname $0)"

PYTHONPATH="$DIR"/../..:"$PYTHONPATH" python3 "$DIR"/run_all_tests.py


#!/bin/sh

DIR="$(dirname $0)"

PYTHONPATH="$DIR":"$PYTHONPATH" python3 "$DIR"/bin/symoro-bin.py


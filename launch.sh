#!/bin/sh

DIR="$(dirname $0)"

PYTHONPATH="$DIR":"$PYTHONPATH" python "$DIR"/bin/symoro-bin.py


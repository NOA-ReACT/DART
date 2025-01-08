#!/bin/bash

set -e

# Use compiler from $CC if set, otherwise use icx if in path, otherwise use gcc
if [ -z "$CC" ]; then
    if command -v icx >/dev/null 2>&1; then
        CC=icx
        echo "Using Intel compiler (override with CC environment variable)"
    else
        CC=gcc
        echo "Using GCC compiler (override with CC environment variable)"
    fi
else
    echo "Using $CC compiler"
fi

CMD="$CC -c coda_fortran.c -o coda_fortran.o -I /home/thgeorgiou/.local/coda/include/"
echo $CMD
$CMD
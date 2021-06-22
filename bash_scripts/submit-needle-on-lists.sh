#!/bin/bash

# prendi riga t-esima di $1 e chiama needle:
#./needle-on-list.sh lists/list_t.dat

if [ "$#" -ne 1 ]; then
    echo "Usage $0 <list of lists file>"
	exit 1
fi

nthline=$(cat "$1"|sed "${SGE_TASK_ID}q;d")
./needle-on-list.sh "$nthline"


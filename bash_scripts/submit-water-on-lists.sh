#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage $0 <list of lists file>"
	exit 1
fi

nthline=$(cat "$1"|sed "${SGE_TASK_ID}q;d")
./water-on-list.sh "$nthline"

#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage $0 <list to split>"
	exit 1
fi

listdir=lists
fnlist=list_
fnlistoflists=list_of_lists.dat
#filesperlist=20
filesperlist=5000
count=1
listnum=1
tot=0

mkdir -p ${listdir}
rm -f ${listdir}/${fnlist}*.dat

echo -en "Writing to ${listdir}/${fnlist}${listnum}.dat ... "
while read l; do
	if [ $count -gt $filesperlist ]; then
		echo "($((count-1)))"
		listnum=$((listnum+1))
		count=1
		echo -en "Writing to ${listdir}/${fnlist}${listnum}.dat ... "
	fi
	echo "$l" >>${listdir}/${fnlist}${listnum}.dat
	tot=$((tot+1))
	count=$((count+1))
#done <<< "$(cat "$1" |head -n 100)"
done <<< "$(cat "$1")"
echo "($((count-1)))"
echo "Done writing lists"
echo "Wrote $tot lines in $listnum files"

ls ${listdir}/${fnlist}*.dat >$fnlistoflists
echo "Wrote $fnlistoflists list file with $listnum lines"


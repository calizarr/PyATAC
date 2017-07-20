#!/bin/bash
IFS=,
shopt -s nullglob
while read sample protocol ab size type rep ind1name ind1seq
do    
    for file in $1/$sample*
    do
        ext=${file#*.}
        mv "$file" "$1/$sample-$protocol-$ab-$size-$type-$rep.$ext"
    done
done < $2

#!/usr/bin/env bash
for i in `seq 0 125`
do
    cd config$i
    max=`ls | sed 's/e_//' | sort -n | tail -1`
    for j in `seq 0 "${max}"`
    do
        dist.pl s/POSCAR e_${j}/POSCAR;
    done
    cd ..
done

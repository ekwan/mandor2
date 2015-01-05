#!/bin/bash
for i in *.xyz; do
    cp $i old_xyz
    awk -f translate.awk $i > temp.xyz
    sleep 1
    echo $i
    mv temp.xyz $i
done

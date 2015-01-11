#!/bin/bash

# takes a tinker xyz file ($1 is the name of the xyz file with extension)
# and prints out atom numbers and their surface tensions (kcal/A^2)

# usage: ./getSurfaceTension.sh Ala.xyz

awk -f convert.awk $1 > temp_geometry.txt
cat header.txt temp_geometry.txt > input.dat
./omnisol < input.dat > omnisol.txt
#cat omnisol.txt
awk '/kcal/,/CS Contribution/ {if (NF==4 && $1 + 0 > 0) {printf "surfaceTension  %-4d %12.8f\n", $1, $4/$3}}' omnisol.txt

rm temp_geometry.txt
rm input.dat
rm omnisol.txt

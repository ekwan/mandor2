#!/bin/bash

cd $1
./analyze $2.xyz -k $2 D > $2.txt1 2> $2.txt2
cat $2.txt2 $2.txt1 > $2.txt
rm $2.txt1 $2.txt2


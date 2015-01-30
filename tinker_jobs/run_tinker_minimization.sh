#!/bin/bash

cd $1
./minimize $2.xyz -k $2 0.2 > $2.out
mv $2.xyz_2 $2_minimized.xyz

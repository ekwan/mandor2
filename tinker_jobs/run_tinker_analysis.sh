#!/bin/bash

cd $1
analyze $2.xyz -k $2 D > $2.txt


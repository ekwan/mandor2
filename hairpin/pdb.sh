#!/bin/bash

for i in *.xyz; do
    xyzpdb $i ../../../tinker/tinker/params/amoebapro13.prm
done

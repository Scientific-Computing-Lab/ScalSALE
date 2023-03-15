#!/bin/bash
module load scr/1.2.0-openmpi-1.10.4-intel-2017
for i in `scr_index --prefix=$1 | awk 'NR != 1 {print $NF}'`; do
    echo Removing: $i;
    scr_index --prefix=$1 --remove=$i;
done

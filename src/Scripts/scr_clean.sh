#!/bin/bash
scr_prefix=`grep '^SCR_PREFIX' ../CR/scr_conf.conf  | awk '{print $1}' | awk -F= '{print $2}'`
./scr_clear_index.sh `eval "echo $scr_prefix"`

for d in `grep '^STORE' ../CR/scr_conf.conf  | awk '{print $1}' | awk -F= '{print $2}'`; do
    echo Removing: `eval "echo $d/${USER}"`;
    rm -rf `eval "echo $d/${USER}"`

done

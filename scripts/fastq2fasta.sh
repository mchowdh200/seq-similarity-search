#!/bin/env bash

input=$1
output=$input.fa.gz

zcat $input |
    paste -d '\t' - - - - |
    cut -f 1,2 |
    sed 's/@/>/'g |
    awk -F '\t' '{print $1"\n"$2}'|
    bgzip -c > $output

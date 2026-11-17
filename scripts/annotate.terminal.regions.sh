#!/usr/bin/env bash

# Usage: ./annotate_terminal_regions.sh 1 100kb.windows.bed intervals.bed 500000

chr_windows="$1"
bedfile="$2"
terminal="$3"

awk -v T="$terminal" '
    NR==FNR { size[$1]=$2; next }
    {
        chr=$1; start=$2; end=$3;
        fivep_end = T;
        threep_start = size[chr] - T;

        if (end <= fivep_end) {
            region = "5p";
        } else if (start >= threep_start) {
            region = "3p";
        } else {
            region = "IT";
        }

        print $0, region;
    }
' "$chr_windows" "$bedfile"

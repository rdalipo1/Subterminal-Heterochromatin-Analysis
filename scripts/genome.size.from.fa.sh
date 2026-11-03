#!/usr/bin/env bash

# Usage: ./genome.size.from.fa.sh input.fa output.txt


input="$1"
output="$2"

awk '/^>/ {
    if (seq) {
        print name "\t" length(seq)
    }
    # grab header name (remove ">")
    name = substr($1, 2)
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    if (seq) {
        print name "\t" length(seq)
    }
}' "$input" > "$output"

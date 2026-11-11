#!/usr/bin/env bash


FASTA="$1"
KMERS="$2"
OUTDIR="$3"

# Create output directory
mkdir -p "$OUTDIR"

# Process each kmer individually
while IFS= read -r KMER; do
    # Skip empty lines
    [[ -z "$KMER" ]] && continue
    
    echo "Processing k-mer: $KMER"
    
    BED_OUT="${OUTDIR}/${KMER}.bed"
    
    # Run fastaRegexFinder.py for this specific kmer
    python3 ./fastaRegexFinder.py \
        -f "$FASTA" \
        -r "$KMER" \
        | awk -v kmer="$KMER" 'BEGIN{OFS="\t"} {print $1, $2, $3, kmer}' \
        > "$BED_OUT"
    
    echo "  -> Written to $BED_OUT"
    
done < "$KMERS"

echo "All k-mers processed"

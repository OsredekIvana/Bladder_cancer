#!/bin/bash

## Pileup
for BAM in ./*_assembly.bam
do
    ID=$(basename "${BAM}" .bam)
    FASTA="${ID}.fasta"

    if [[ ! -f "${FASTA}" ]]; then
        echo " FASTA file not found for ${BAM}, expected ${FASTA}" >&2
        continue
    fi

    echo "Processing: ${ID}"
    modkit pileup \
        --combine-strands --cpg --ref "${FASTA}" \
        --ignore h --threads 8 --log-filepath "${ID}.log" \
        "${BAM}" "${ID}.pileup.bed"

    bgzip "${ID}.pileup.bed"
    tabix -p bed "${ID}.pileup.bed.gz"
done

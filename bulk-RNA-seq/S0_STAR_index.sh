#!/usr/bin/bash

# run STAR to make index for organism of interest
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir ref/star_index \
     --genomeFastaFiles ref/test_genome.fa \
     --sjdbGTFfile ref/test_genes.gtf \
     --sjdbOverhang 99     # readLength-1 for 100bp reads; adjust to your read length

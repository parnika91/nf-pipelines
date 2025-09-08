First run S0_STAR_index.sh with the fasta and GTF file of the organism(s) of interest.

Do this only once. Index can be reused for other experiments.

Next, run the command in S1_nf-run on the terminal. It maps fastq files from data/fastq/ with STAR, quantifies genes with featureCounts and analyses differentially expressed genes with EdgeR.

--------------------------------------

Notes on indexing with STAR:

1. Why index a genome?

"Like the index at the end of a book, an index of a large DNA sequence allows one to rapidly find shorter sequences embedded within it." -[How to map billions of short reads onto genomes](10.1038/nbt0509-455)

This is why we need a gene annotation file such as a .GTF and the actual DNA sequence from the .fa file.


2. It is important that chromosome names in the fasta and GTF files match exactly. For example,

.gtf file:

>**chr1**  test  gene    1   60   .   +   .   gene_id "GENE1";
>**chr1**  test  exon    1   60   .   +   .   gene_id "GENE1"; transcript_id "TX1";
>
>**chr2**  test  gene    10  80   .   -   .   gene_id "GENE2";
>**chr2**  test  exon    10  80   .   -   .   gene_id "GENE2"; transcript_id "TX2";

.fa file:

>**chr1**
ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
>**chr2**
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA

---------------------------------------


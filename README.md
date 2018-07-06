# 10xrnaseq
Some codes for processing single-cell RNA-seq data

10x Genomics supplies massive single-cell NGS data.
To process the data efficiently some codes were prepared.

These programs process following datasets.

1. Paired end having barcode in R1 and target sequences in R2 in fastq (or fastq.gz) format.

2. The order of sequences were common in R1 and R2.

3. In the first line of each sequence in R1/R2, plate index is shown after last colon, such as "@ST-E00127:770:HMG5JCCXY:1:1101:5315:1344 1:N:0:NAGTACTG".

4. STAR is supposed as alignment software and featureCount as counting software.

## Procedure

1. All reads will be renamed by correct barcodes. Frequency of indexes are counted and 
ount_and_rename --R1 [R1.fastq.gz] --R2 [R2.fastq.gz] -o [renamed_fastq.gz] --mismatches [num]
###OPTIONS
- --R1, -1 [filename] fastq file having index
- --R2, -2 [filename] fastq file having reads
- --barcode-size, -b [number] sequence of R1 is split into [cell barcode] + [UMI]. this option gives the length of cell index.
- -k [number] mismatch tolerance (default 1)
- -m [number] threshold of read counts as available (default 500)
- -o [filename] fastq.gz file of renamed sequences

2. Align renamed sequence on genome. (STAR is supposed but other aligners are assumed available.)

3. Split SAM/BAM file into single cells.
python split_sam.py
###OPTIONS
- -i [filename] BAM/SAM filename
- -o [directory] Split SAM file container
- -n [number]

4. Use featureCount

5.
normalize_fc

###OPTIONS

-i [filename...] output of featureCounts
-o prefix (<prefix>.cnt and <prefix>.tpm will be genrated)

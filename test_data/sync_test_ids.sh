#!/bin/bash

ID=$1
echo $ID

gunzip -dc Yale-${ID}_R2* | grep \@ | perl -pe 's/^\@//gi' | sort > R2ids.txt
gunzip -dc Yale-${ID}_R1* | grep \@ | perl -pe 's/^\@//gi' | sort > R1ids.txt
comm R?ids.txt  -12 > R12ids.txt

seqtk subseq Yale-${ID}_R1_testdata.fastq.gz R12ids.txt \
    | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | \
    gzip > Yale-${ID}_syncR1_testdata.fastq.gz

seqtk subseq Yale-${ID}_R2_testdata.fastq.gz R12ids.txt \
    | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | \
    gzip > Yale-${ID}_syncR2_testdata.fastq.gz

rm R1ids.txt R2ids.txt R12ids.txt

mkdir unsorted
mv Yale-${ID}_R1_testdata.fastq.gz unsorted/Yale-${ID}_R1_testdata.fastq.gz
mv Yale-${ID}_R2_testdata.fastq.gz unsorted/Yale-${ID}_R2_testdata.fastq.gz

mv Yale-${ID}_syncR1_testdata.fastq.gz Yale-${ID}_R1_testdata.fastq.gz
mv Yale-${ID}_syncR2_testdata.fastq.gz Yale-${ID}_R2_testdata.fastq.gz

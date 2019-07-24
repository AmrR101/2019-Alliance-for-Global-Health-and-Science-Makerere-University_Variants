#!/bin/bash

## assumes star version 2.7.0e

start=`date +%s`
echo $HOSTNAME

outpath="References"
[[ -d ${outpath} ]] || mkdir ${outpath}

cd ${outpath}
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
FASTA="../GRCh38.primary_assembly.genome.fa"

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
gunzip gencode.v29.primary_assembly.annotation.gtf.gz
GTF="../gencode.v29.primary_assembly.annotation.gtf"

mkdir star.overlap100.gencode.v29
cd star.overlap100.gencode.v29

call="STAR
     --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir . \
     --sjdbOverhang 100 \
     --sjdbGTFfile ${GTF} \
     --genomeFastaFiles ${FASTA}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
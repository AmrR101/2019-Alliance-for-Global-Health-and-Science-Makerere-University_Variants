#!/bin/bash

## assumes bwa is on the path

start=`date +%s`
echo $HOSTNAME

outpath="Reference"
[[ -d ${outpath} ]] || mkdir ${outpath}

cd ${outpath}
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/789/395/GCF_000789395.1_ASM78939v1/GCF_000789395.1_ASM78939v1_genomic.fna.gz
gunzip GCF_000789395.1_ASM78939v1_genomic.fna.gz
FASTA="GCF_000789395.1_ASM78939v1_genomic.fna"

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/789/395/GCF_000789395.1_ASM78939v1/GCF_000789395.1_ASM78939v1_genomic.gff.gz
gunzip GCF_000789395.1_ASM78939v1_genomic.gff.gz
GFF="GCF_000789395.1_ASM78939v1_genomic.gff"

call="bwa index ${FASTA}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime

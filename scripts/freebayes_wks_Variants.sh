#!/bin/bash

## assumes freebayes is available on the Path

start=`date +%s`
echo $HOSTNAME

outpath='03-Freebayes'
[[ -d ${outpath} ]] || mkdir ${outpath}

outfile=${outpath}/freebayes.vcf

mapfasta=./Reference/GCF_000789395.1_ASM78939v1_genomic.fna

BAMLIST="bamlist.txt"

call="freebayes \
  --bam-list ${BAMLIST} \
  --fasta-reference ${mapfasta} \
  --vcf  ${outfile} \
  --ploidy 1 \
  --use-best-n-alleles 6 \
  --max-complex-gap 75 \
  --min-mapping-quality 30 \
  --min-base-quality 20 \
  --min-supporting-allele-qsum 0 \
  --genotype-variant-threshold 0"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

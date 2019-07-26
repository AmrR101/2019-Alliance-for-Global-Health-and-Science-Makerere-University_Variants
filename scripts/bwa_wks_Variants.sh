#!/bin/bash
#

start=`date +%s`
echo $HOSTNAME


THREADS=4
MAPTHREADS=$(expr ${THREADS} - 2)
SORTTHREADS=$(expr ${THREADS} - ${MAPTHREADS})

inpath=01-HTS_Preproc

outpath='02-BWA'
[[ -d ${outpath} ]] || mkdir ${outpath}

mapfasta=./Reference/GCF_000789395.1_ASM78939v1_genomic.fna

for sample in `cat samples.txt`
do

    r1=${inpath}/${sample}/${sample}_R1.fastq.gz
    r2=${inpath}/${sample}/${sample}_R2.fastq.gz
    se=${inpath}/${sample}/${sample}_SE.fastq.gz

    [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

    echo "SAMPLE: ${sample}"

    output=${outpath}/${sample}/${sample}_bwa.bam

    call="bwa mem -t ${MAPTHREADS} \
      -R '@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tDS:Paired' \
      ${mapfasta} ${r1} ${r2} | \
      samtools sort -m 768M --threads ${SORTTHREADS} -o ${output}-pe -"
    echo $call
    eval $call

    call="bwa mem -t ${MAPTHREADS} \
      -R '@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tDS:Paired' \
      ${mapfasta} ${se} | \
      samtools sort -m 768M --threads ${SORTTHREADS} -o ${output}-se -"
    echo $call
    eval $call

    call="samtools merge -f -@ ${THREADS} ${output} ${output}-pe ${output}-se"
    echo $call
    eval $call

    call="samtools index -@ ${THREADS} ${output}"
    echo $call
    eval $call

    call="samtools idxstats ${output} > ${output}.idxstats"
    echo $call
    eval $call

    call="samtools flagstat -@ ${THREADS} ${output} > ${output}.flagstat"
    echo $call
    eval $call

    call="samtools stats -@ ${THREADS} ${output} > ${output}.stats"
    echo $call
    eval $call

done

end=`date +%s`

runtime=$((end-start))

echo $runtime

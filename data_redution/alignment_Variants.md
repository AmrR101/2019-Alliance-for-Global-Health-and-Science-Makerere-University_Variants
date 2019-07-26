# Alignment to Read Counts & Visualization in IGV

This document assumes [preproc htstream](./preproc_htstream_Variants.md) has been completed.

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can copy over from the flash drive

<img src="alignment_figures/wkflow_2.png" alt="workflow flowchart" width="600px"/>

## Alignment vs Assembly

**Given sequence data**,

_Assembly seeks to put together the puzzle without knowing what the picture is._  

_Mapping tries to put together the puzzle pieces directly onto an image of the picture._

In mapping the question is more, given a small chunk of sequence, where in the genome did this piece most likely come from.

The goal then is to find the match(es) with either the “best” edit distance (smallest), or all matches with edit distance less than max edit dist. Main issues are:

* Large search space
* Regions of similarity (aka repeats)
* Gaps (INDELS)
* Complexity (RNA, transcripts)

#### Considerations
* Placing reads in regions that do not exist in the reference genome (reads extend off the end) [ mitochondrial, plasmids, structural variants, etc.].
* Sequencing errors and variations: alignment between read and true source in genome may have more differences than alignment with some other copy of repeat.
* What if the closest fully sequenced genome is too divergent? (3% is a common alignment capability)
* Placing reads in repetitive regions: Some algorithms only return 1 mapping; If multiple: map quality = 0
* Algorithms that use paired-end information => might prefer correct distance over correct alignment.

In RNAseq data, you must also consider effect of splice junctions, reads may span an intron.

<img src="alignment_figures/alignment_figure1.png" alt="alignment_figure1" width="600px"/>

### Aligners
Many [alignment algorithms](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software
) to choose from.
* Aligners (global and local)
  * bwa-mem
  * Bowtie2 (global default, local possible)
* Spliced Aligners
  * STAR
  * HiSAT2 (formerly Tophat [Bowtie2])
  * GMAP - GSNAP
  * SOAPsplice
  * MapSplice
* Pseudo-Aligners (transcriptomes)
  * salmon
  * kalisto

#### Pseudo-aligners
* Quasi-mapping
* Probabilistic
* Map to transcripts, not genome
* Does transcript quantifications (or gene)
* Blazing FAST and can run on most laptops
* Experience suggests differences between “traditional” mappers are in the low abundance genes.

### Genome and Genome Annotation

Genome sequence fasta files and annotation (gff, gtf) files go together! These should be identified at the beginning of analysis.
* Genome fasta files should include all primary chromosomes, unplaced sequences
and un-localized sequences, as well as any organelles. Should bet contain any
contigs that represent patches, or alternative haplotypes.
* If you expect contamination, or the presence of additional sequence/genome,
add the sequence(s) to the genome fasta file.
* Annotation file should be GTF (preferred), and should be the most
comprehensive you can find.
  * Chromosome names in the GTF must match those in the fasta file (they don’t always do).

### Alignment concepts

* Multimappers:
  * Reads that align equally well to more than one reference location.
  * Generally, multimappers are discounted in variant detection, and are often discounted in counting
applications (like RNA-Seq ... would “cancel” out anyway).
  * Note: multimapper “rescue” in some algorithms (RSEM, Express?).
* Duplicates:
  * Reads or read pairs arising from the same original library fragment, either during library preparation (PCR duplicates).
  * Generally, duplicates can only be detected reliably with paired-end sequencing. If PE, they’re discounted in variant detection, and discounted in counting applications (like RNA-Seq).
* Clipping vs Splicing  
<img src="alignment_figures/alignment_figure3.png" alt="alignment_figure3" width="400px"/>  
* Inner length, insert size, fragment length  
<img src="alignment_figures/alignment_figure4.jpg" alt="alignment_figure4" width="400px"/>  
*From https://www.biostars.org/p/106291/

## Indexing a Reference sequence and annotation

**1\.** First lets make sure we are where we are supposed to be and create a
References directory.

    cd ~/variant_example
    mkdir Reference
    cd Reference

---
**2\.** To align our data we will need the genome (fasta) and annotation (gff) for Clostridium baratii str. Sullivan. There are many places to find them, but we are going to get them from the [NCBI](ftp://ftp.ncbi.nlm.nih.gov).

We need to first get the url for the annotation gff.

TODO: Add fig  

Then the url for the genome fasta.

TODO: Add fig

**3\.** We are going to use an aligner called ['BWA MEM'](https://github.com/lh3/bwa) to align the sequence reads, but first we need to index the genome for BWA. Lets pull down a shell script to index the genome.

  cd ~/variant_example
  wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-Alliance-for-Global-Health-and-Science-Makerere-University_Variants/master/scripts/bwa_index_wks_Variants.sh
  less bwa_index_wks_Variants.sh

 When you are done, type "q" to exit.

 ```bash
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
 ```

1. The script uses wget to download the fasta and GFF files from NCBI using the links you found earlier.
1. Uncompresses them using gunzip.
1. Creates the bwa index

  cd ~/variant_example
  bash bwa_index_wks_Variants.sh > scriptout/bwa_index.out 2> scriptout/bwa_index.err

This step will take a seconds.

**IF** for some reason it didn't finish, is corrupted, or you missed the session, you can copy over from the flash drive.

## Alignments

**1\.** We are now ready to try an alignment. Let's create an output directory for BWA:

    cd ~/variant_example/HTS_testing

and let's run bwa on the pair of hstream streamed test files we created earlier.

Then run the bwa command (Its on multiple lines for readability)

  bwa mem \
  -t 2 \
  -R '@RG\tID:sample1.streamed\tSM:sample1.streamed\tPL:ILLUMINA\tDS:Paired' \
  ../Reference/GCF_000789395.1_ASM78939v1_genomic.fna  \
  sample1.streamed_R1.fastq.gz sample1.streamed_R2.fastq.gz | \
  samtools sort -m 768M --threads 2 -o sample1.streamed_bwa.bam -


In the command, we are telling bwa to map reads to the genome using 2 threads, and providing a RG tag (see [filetype.md](./filetypes.md)), __piping__ to samtools sort (max mem 768M and threads 2), the name for the output file will be sample1.streamed_bwa.bam.

##  Now let's take a look at an alignment in IGV. TODO: GET IGV

**1\.** We first need to index the bam file, will use 'samtools' for this step, which is a program to manipulate SAM/BAM files. Take a look at the options for samtools and 'samtools index'.

    samtools
    samtools index

We need to index the BAM file:

    cd ~/variant_example/HTS_testing
    samtools index sample1.streamed_bwa.bam

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can copy over from the flash drive

---
**2\.** Now we are ready to use IGV. Go to the [IGV page at the Broad Institute](http://software.broadinstitute.org/software/igv/).

<img src="alignment_figures/index_igv1.png" alt="index_igv1" width="600px"/>

And then navigate to the download page, [IGV download](http://software.broadinstitute.org/software/igv/download)
<img src="alignment_figures/index_igv2.png" alt="index_igv2" width="600px"/>

Here you can download IGV for your respective platform (Window, Mac OSX, Linux), but we are going to use the web application they supply, [IGV web app](https://igv.org/app/).

<img src="alignment_figures/index_igv3.png" alt="index_igv3" width="600px"/>

---
**3\.** The first thing we want to do is load the Human genome. Click on "Genomes" in the menu and choose "Human (GRCh38/hg38)".

<img src="alignment_figures/index_igv4.png" alt="index_igv4" width="600px"/>

---
**4\.** Now let's load the alignment bam and index files. Click on "Tracks" and choose "Local File ...".

<img src="alignment_figures/index_igv5.png" alt="index_igv5" width="600px"/>

Navigate to where you transferred the bam and index file and select them both.

<img src="alignment_figures/index_igv6.png" alt="index_igv6" width="600px"/>

Now your alignment is loaded. Any loaded bam file aligned to a genome is called a "track".

<img src="alignment_figures/index_igv7.png" alt="index_igv7" width="600px"/>

---
**5\.** Lets take a look at the alignment associated with the gene __HBB__:

<img src="alignment_figures/index_igv8.png" alt="index_igv8" width="600px"/>

<img src="alignment_figures/index_igv9.png" alt="index_igv9" width="600px"/>

You can zoom in by clicking on the plus sign (top right) or zoom out by clicking the negative sign. You also may have to move around by clicking and dragging in the BAM track window.

You can also zoom in by clicking and dragging across the number line at the top. That section will highlight, and when you release the button, it will zoom into that section.

<img src="alignment_figures/index_igv10.png" alt="index_igv10" width="600px"/>

<img src="alignment_figures/index_igv11.png" alt="index_igv11" width="600px"/>

Reset the window by searching for HBB again. And zoom in 1 step.

<img src="alignment_figures/index_igv12.png" alt="index_igv12" width="600px"/>

---
**6\.** See that the reads should be aligning within the exons in the gene. This makes sense, since RNA-Seq reads are from exons. Play with the settings on the right hand side a bit.

<img src="alignment_figures/index_igv13.png" alt="index_igv13" width="600px"/>

<img src="alignment_figures/index_igv14.png" alt="index_igv14" width="600px"/>

<img src="alignment_figures/index_igv15.png" alt="index_igv15" width="600px"/>

<img src="alignment_figures/index_igv16.png" alt="index_igv16" width="600px"/>

<img src="alignment_figures/index_igv17.png" alt="index_igv17" width="600px"/>

---

## Running BWA on the experiment

**1\.** We can now run BWA across all samples on the real data using a shell script, [bwa_wks_Variants.sh](../scripts/bwa_wks_Variants.sh), that we should take a look at now.

    cd ~/variant_example  # We'll run this from the main directory
    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-Alliance-for-Global-Health-and-Science-Makerere-University_Variants/master/scripts/bwa_wks_Variants.sh
    less bwa_wks_Variants.sh

When you are done, type "q" to exit.


```bash
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

    output=${outpath}/${sample}/${sample}-_bwa.bam

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
```

Because in preprocessing we overlapped reads, we have both single-end and paired-end reads. We need to process these separately and then merge (using samtools merge) the results. In addition to mapping and merging, we will also use samtools to index, idxstats, flagstats and stats. Take a quick look at the help documentation for these apps.


  samtools merge
  samtools index
  samtools idxstats
  samtools flagstats
  samtools stats


After looking at the script, lets run it.

    bash bwa_wks_Variants.sh > scriptout/bwa.out 2> scriptout/bwa.err # moment of truth!

Takes about 15 minutes to run all samples.

## Quality Assurance - Mapping statistics as QA/QC.

**1\.** Once your jobs have finished successfully (check the error and out logs like we did in the previous exercise), use a script of ours, [star_stats.sh](../scripts/star_stats.sh) to collect the alignment stats. Don't worry about the script's contents at the moment. For now:

    cd ~/variant_example  # We'll run this from the main directory
    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-Alliance-for-Global-Health-and-Science-Makerere-University_Variants/master/scripts/bwa_stats_wks_Variants.R

    R CMD BATCH bwa_stats_wks_Variants.R
    cat bwa_stats.txt

Open in excel (or excel like application), you may have to move the header column 1 cell to the right, and lets review.

Are all samples behaving similarly? Discuss ...


Once we are satisfied the data look good, lets remove the pre-merged bam files that at in pe and se.


  cd ~/variant_example  # We'll run this from the main directory
  rm -rf 02-BWA/*/*-pe 02-BWA/*/*-se


## Scripts

shell script for indexing the genome

[bwa_index_wks_Variants.sh](../scripts/bwa_index_wks_Variants.sh)

shell script for mapping using bash loop and bwa.

[bwa_wks_Variants.sh](../scripts/bwa_wks_Variants.sh)

shell script to produce summary mapping table

[bwa_stats_wks_Variants.sh](../scripts/bwa_stats_wks_Variants.sh)

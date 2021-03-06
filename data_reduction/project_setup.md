# Data Setup

Let's set up a project directory for the workshop, and talk a bit about project philosophy..

**1\.** First, create a directory for the example project in your home directory:

    mkdir -p ~/variant_example

---

**2a\.** Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets copy the workshop data into this folder.

    cd ~/variant_example

Now copy over the data in 00-RawData from the flash drive.

    cd 00-RawData/

This directory now contains a folder for each sample and the fastq files for each sample are in the sample folders.

**2b\.** lets create a sample sheet for the project, store sample names in a file called samples.txt

    ls > ../samples.txt
    cat ../samples.txt

---
**3a\.** Now, take a look at the raw data directory.

    ls ~/variant_example/00-RawData


**3b\.** You will see a list of the contents of each directory.

    ls *

**3c\.** Lets get a better look at all the files in all of the directories.

    ls -lah */*

---

**5\.** Pick a directory and go into it. View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files):

    cd sample1/
    zless sample1_R1.fastq.gz

Make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:

    gunzip -c sample1_R1.fastq.gz | wc -l

Divide this number by 4 and you have the number of reads in this file. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    gunzip -c sample1_R1.fastq.gz  | head -2 | tail -1

Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block. Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):

    echo -n [sequence] | wc -c

This will give you the length of the read.

Also can do the bash one liner:

    echo -n $(gunzip -c sample1_R1.fastq.gz  | head -2 | tail -1) | wc -c

See if you can figure out how this command works.

---

**5\.** Now go back to your 'variant_example' directory and create two directories called 'scriptout' and '01-HTS_Preproc':

    cd ~/variant_example
    mkdir scriptout
    mkdir 01-HTS_Preproc

The results of all our scripts will output to .out and .err files into the scriptout folder. The results of our preprocessing steps will be put into the 01-HTS_Preproc directory. The next step after that will go into a "02-..." directory, etc. You can collect scripts that perform each step, and notes and metadata relevant for each step, in the directory for that step. This way anyone looking to replicate your analysis has limited places to search for the commands you used. In addition, you may want to change the permissions on your original 00-RawData directory to "read only", so that you can never corrupt your raw data. (We won't worry about this here, because we've linked in sample folders).

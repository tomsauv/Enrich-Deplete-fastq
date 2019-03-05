# Mapping-tid-bits

Sometimes, we need to ENRICH and/or DEPLETE paired-ends Illumina data set for a particular genomic compartments<br/>

Below are example command lines using BWA-mem + Samtools to do so<br/>

Raw Illumina data: Raw_R1.fastq and Raw_R2.fastq<br/>
File with scaffolds: query.fasta<br/>

BWA-mem v0.7.15-r1140 and Samtools v1.6 are in the path on a linux platform with 40 cores<br/>

**1) Index the file containing scaffolds**

*bwa index query.fasta*

**2) Run the actual read mapping and pipe results to samtools (adjust -t parameters to your number of cores)**

*bwa mem -t 40 query.fasta Raw_R1.fastq Raw_R2.fastq | samtools sort > query.sorted.bam*

**3) Index the sorted bam file**

*samtools index query.sorted.bam*

# ENRICH (keep matching reads)
**Enrich by pulling out reads mapping to the query (with -F flag and integer 4)**

*samtools view -b -F 4 query.sorted.bam > enriched.sorted.bam*

**Extract the reads to fastq file**

*samtools fastq enriched.sorted.bam -1 enriched_R1.fastq -2 enriched_R2.fastq -0 /dev/null -n*

# and/or DEPLETE (keep non-matching reads)
**Deplete by pulling out reads not mapping to query (with -f flag and integer 4)**

*samtools view -b -f 4 query.sorted.bam > depleted.sorted.bam*

**Extract reads to fastq file**

*samtools fastq depleted.sorted.bam -1 depleted_R1.fastq -2 depleted_R2.fastq -0 /dev/null -n*

# COUNT reads in R1 fastq files to check
**Raw (initial) R1 fastq file**<br/>

*echo "Raw (R1)" > counts.txt*<br/>
*echo $(cat R1.fastq|wc -l)/4|bc >> counts.txt*<br/>

**Enriched R1 fastq file**<br/>

*echo "Enriched (R1)" >> counts.txt*<br/>
*echo $(cat enriched_R1.fastq|wc -l)/4|bc >> counts.txt*<br/>

**Depleted R1 fastq file**<br/>

*echo "Depleted (R1)" >> counts.txt*<br/>
*echo $(cat depleted_R1.fastq|wc -l)/4|bc >> counts.txt*<br/>

**Display counts on screen**<br/>

*head counts.txt*<br/>

When summed, enriched and depleted reads counts should equal those of the raw R1 fastq file





# Mapping-tid-bits

Sometimes, we need to ENRICH and/or DEPLETE paired-ends Illumina data sets for particular scaffolds<br/>

Below are example command lines using BWA-mem + Samtools to do so<br/>

Here, BWA-mem v0.7.15-r1140 and Samtools v1.6 are in the path on a linux platform with 40 cores<br/>

**1) Index the file containing scaffolds**

*bwa index query.fasta*

**2) Run the actual read mapping and pipe results to samtools (adjust -t parameters to your number of cores)**

*bwa mem -t 40 query.fasta R1.fastq.gz R2.fastq.gz | samtools sort > query.sorted.bam*

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

# COUNT reads in fastq files to check
**Initial fastq file (R1)**<br/>

*echo "Raw (R1)" > counts.txt*<br/>
*echo $(zcat R1.fastq.gz|wc -l)/4|bc >> counts.txt*<br/>

**Initial fastq file (R2)**<br/>

*echo "Raw (R2)" >> counts.txt*<br/>
*echo $(zcat R2.fastq.gz|wc -l)/4|bc >> counts.txt*<br/>

**Enriched fastq file (contains R1+R2)**<br/>

*echo "Enriched (R1+R2)" >> counts.txt*<br/>
*echo $(cat enriched.fastq|wc -l)/4|bc >> counts.txt*<br/>

**Depleted fastq file (contains R1+R2)**<br/>

*echo "Depleted (R1+R2)" >> counts.txt*<br/>
*echo $(cat depleted.fastq|wc -l)/4|bc >> counts.txt*<br/>

**Display counts on screen**<br/>

*head counts.txt*<br/>

When summed, enriched and depleted reads counts should equal those of the initial fastq files (R1+R2) 





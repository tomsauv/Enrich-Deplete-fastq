# Mapping workflow to sort Paired-ends fastq<br/>

Sometimes, we need to ENRICH and/or DEPLETE paired-ends Illumina data set for a particular target<br/>

Below are example command lines using BWA-mem + Samtools to do so<br/>

Raw Illumina data: *raw_R1.fastq* and *raw_R2.fastq*<br/>
File with scaffolds: *query.fasta*<br/>

BWA-mem v0.7.15-r1140 and Samtools v1.6 are in the path on a linux platform with 40 cores<br/>

**1) Index the file containing scaffolds**<br/>

*bwa index query.fasta*<br/>

**2) Run the actual read mapping and pipe results to samtools (adjust -t parameters to your number of cores)**<br/>

*bwa mem -t 40 query.fasta raw_R1.fastq raw_R2.fastq | samtools sort > query.sorted.bam*<br/>

**3) Index the sorted bam file**<br/>

*samtools index query.sorted.bam*<br/>

# ENRICH target (skip unmapped reads)<br/>
**Enrich by keeping mapped reads (-f ) and skipping unmapped reads (-F 4)**<br/>

*samtools view -b -F 4 query.sorted.bam > enriched.sorted.bam*<br/>

(*Use flags -f 2 -F 2052* if you want to only keep proper pairs and thus exclude singletons, unmapped reads and secondary alignments)<br/>

**Extract the reads to fastq file**<br/>

*samtools fastq enriched.sorted.bam -1 enriched_R1.fastq -2 enriched_R2.fastq -0 /dev/null -n*<br/>

# and/or DEPLETE target (keep unmapped reads)<br/>
**Deplete by keeping unmapped reads (-f 4) and skipping mapped reads (-F 3)**<br/>

*samtools view -b -f 4 -F 3 query.sorted.bam > depleted.sorted.bam*<br/>

(*-F 3* is to be thorough, but it may be omitted)<br/>

**Extract reads to fastq file**<br/>

*samtools fastq depleted.sorted.bam -1 depleted_R1.fastq -2 depleted_R2.fastq -0 /dev/null -n*<br/>

# COUNT reads in R1 fastq files to check<br/>
**Raw (initial) R1 fastq file**<br/>

*cat raw_R1.fastq | echo -e "Rawfile_R1\t" $((\`wc -l\`/4)) > counts.txt*<br/>

**Enriched R1 fastq file**<br/>

*cat enriched_R1.fastq | echo -e "Enriched_R1\t" $((\`wc -l\`/4)) >> counts.txt*<br/>

**Depleted R1 fastq file**<br/>

*cat depleted_R1.fastq | echo -e "Depleted_R1\t" $((\`wc -l\`/4)) >> counts.txt*<br/>

**Sum Enriched+Depleted R1 fastq files**<br/>

*cat enriched_R1.fastq depleted_R1.fastq | echo -e "Enriched+Depleted\t" $((\`wc -l\`/4)) >> counts.txt*<br/>

**Display counts on screen**<br/>

*head counts.txt*<br/>

The sum of enriched and depleted reads counts should equal those of the raw R1 fastq file<br/>

# Extra: Single-end data<br/>

In case you are working with single end data such as as long reads<br/>
You may use the following command:<br/>

*bwa mem -t 40 query.fasta raw_single-end.fastq | samtools sort > query.sorted.bam*<br/>
**Enrich**<br/>
*samtools view -b -F 4 query.sorted.bam > depleted.sorted.bam*<br/>
**Deplete**<br/>
*samtools view -b -f 4 query.sorted.bam > depleted.sorted.bam*<br/>




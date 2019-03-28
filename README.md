# Workflow to sort Paired-ends fastq<br/>

Sometimes, we need to ENRICH and/or DEPLETE paired-ends Illumina data set for a particular target<br/>

Below are example command lines using BWA-mem + Samtools to do so<br/>

The Raw Illumina data files are ```raw_R1.fastq``` and ```raw_R2.fastq```<br/>
The file with scaffolds is ```query.fasta```<br/>

BWA-mem v0.7.15-r1140 and Samtools v1.6 are in the path on a Linux platform<br/>

**1) Index the file containing scaffolds**<br/>

```
bwa index query.fasta
```

**2) Run the actual read mapping and pipe results to samtools (adjust -t parameters to your number of cores)**<br/>

```
bwa mem -t 40 query.fasta raw_R1.fastq raw_R2.fastq | samtools sort > query.sorted.bam
```

**3) Index the sorted bam file**<br/>

```
samtools index query.sorted.bam
```

# ENRICH target (skip unmapped reads)<br/>
**Enrich by skipping unmapped reads ```-F 4``` (thus keeping mapped reads)**<br/>

```
samtools view -b -F 4 query.sorted.bam > enriched.sorted.bam
```

(Use flags ```-f 2 -F 2052``` if you want to only keep proper pairs and exclude unmapped reads + secondary alignments)<br/>

**Extract the reads to fastq file**<br/>

```
samtools fastq enriched.sorted.bam -1 enriched_R1.fastq -2 enriched_R2.fastq -0 /dev/null -n
```

# and/or DEPLETE target (keep unmapped reads)<br/>
**Deplete by keeping unmapped reads ```-f 4``` (thus skipping mapped reads)**<br/>

```samtools view -b -f 4 query.sorted.bam > depleted.sorted.bam```<br/>

(Unlike above, we cannot sort proper pairs since such information is gained via mapping and we are segregating unmapped reads...)

**Extract reads to fastq file**<br/>

```samtools fastq depleted.sorted.bam -1 depleted_R1.fastq -2 depleted_R2.fastq -0 /dev/null -n```<br/>

# COUNT reads in R1 fastq files to check<br/>
**Raw (initial) R1 fastq file**<br/>

```cat raw_R1.fastq | echo -e "Rawfile_R1\t" $((\`wc -l\`/4)) > counts.txt```<br/>

**Enriched R1 fastq file**<br/>

```cat enriched_R1.fastq | echo -e "Enriched_R1\t" $((\`wc -l\`/4)) >> counts.txt```<br/>

**Depleted R1 fastq file**<br/>

```cat depleted_R1.fastq | echo -e "Depleted_R1\t" $((\`wc -l\`/4)) >> counts.txt```<br/>

**Sum Enriched + Depleted R1 fastq files**<br/>

```cat enriched_R1.fastq depleted_R1.fastq | echo -e "Enriched+Depleted\t" $((\`wc -l\`/4)) >> counts.txt```<br/>

**Display counts on screen**<br/>

```head counts.txt```<br/>

The sum of enriched and depleted R1 reads counts should equal those of the raw R1 fastq file<br/>
(if you used flags ```-f 2 -F 2052``` to enrich, the sum will not be equal but less)<br/>

# Extra: Single-end data<br/>

```bwa index query.fasta```<br/>
```bwa mem -t 40 query.fasta raw_single_end.fastq | samtools sort > query.sorted.bam```<br/>
```samtools index query.sorted.bam```<br/>
**Enrich** (use ```-F 4``` to exclude unmapped or ```-F 2052``` to exclude unmapped + secondary alignments)<br/> 
```samtools view -b -F 4 query.sorted.bam > enriched.sorted.bam```<br/>
```samtools fastq  enriched.sorted.bam -F 4 -0 enriched_SE.fastq```<br/> 
**Deplete**<br/>
```samtools view -b -f 4 query.sorted.bam > depleted.sorted.bam```<br/>
```samtools fastq  depleted.sorted.bam -f 4 -0 depleted_SE.fastq```<br/>
**Count** (as above replacing with the appropriate fastq file name)<br/>



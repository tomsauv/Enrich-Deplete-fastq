# Mapping-tid-bits

Sometimes, we need to ENRICH and/or DEPLETE next generation sequencing data sets for particular scaffolds<br/>

Below are some command lines using BWA-mem + Samtools to do so for paired end Illumina data sets<br/>

Note that analyses were run on a linux platform with 40 cores<br/>
BWA-mem v0.7.15-r1140 and Samtools v1.6 are in the path<br/>
All files are placed and output in the current directory<br/>

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

*samtools fastq -F 4 enriched.sorted.bam > enriched.fastq*

# and/or DEPLETE (keep non-matching reads)
**Deplete by pulling out reads not mapping to query (with -f flag and integer 4)**

*samtools view -b -f 4 query.sorted.bam > depleted.sorted.bam*

**Extract reads to fastq file**

*samtools fastq -f 4 depleted.sorted.bam > depleted.fastq*

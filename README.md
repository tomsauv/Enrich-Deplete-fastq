# Mapping-tid-bits

Sometimes, we need to deplete/enrich next generation sequencing data sets for a particular genomic compartment
Below are some command lines using BWA-mem + Samtools to do so for paired end Illumina data sets

Note that analyses were run on a linux platform with 40 cores and BWA-mem v0.7.15-r1140 and Samtools v1.6 in the path

Raw zipped data (as .fastq.gz) are in the current directory. Files will be output in the current directory

We want to enrich/deplete based on scaffolds/sequences found in the "query.fasta" file that is also placed in the current directory

**1) Index the sequence file**

*bwa index query.fasta*

**2) Run the actual read mapping and pipe results to samtools (adjust -t parameters to your number of cores)**

*bwa mem -t 40 query.fasta file_R1.fastq.gz file_R2.fastq.gz | samtools sort > query.sorted.bam*

**3) Index the sorted bam file**

*samtools index query.sorted.bam*

# ENRICH
**4) Enrich by pulling out reads mapping to the query (with -F flag and integer 4)

*samtools view -b -F 4 query.sorted.bam > enriched.sorted.bam*

**5) Extract the reads mapping to query in a fastq file**

*samtools fastq -F 4 enriched.sorted.bam > enriched.fastq*

# DEPLETE
# Deplete by pulling out reads not mapping to query (with -f flag and integer 4)
samtools view -b -f 4 query.sorted.bam > depleted.sorted.bam
# extract reads not mapping to query in a fastq file
samtools fastq -f 4 depleted.sorted.bam > depleted.fastq


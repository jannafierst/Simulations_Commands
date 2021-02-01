### Simulations_Commands
###commands used for various analyses used for Sim paper 

###ART Commands

###NanoSim commands

###Subsampling Illumina data (fixed number of reads)
#-s100 is the seed this MUST BE THE SAME for both commands for the reads to match

seqtk sample -s100 read1.fq 10000 > sub1.fq
seqtk sample -s100 read2.fq 10000 > sub2.fq

###Subsampling ONT data (proportion of reads)

seqtk sample -s100 read.fq 0.5 > sub.fq

###Canu Assembly

canu -p [sample] -d [out_directory_name] -nanopore reads.fastq genomeSize=XXXM

###Canu Correction

canu -correct -p [sample] -d [out_directory_name] -nanopore reads.fastq genomeSize=XXXM

###MaSuRCA Assemblies 

###Pilon Polishing (this script can be looped to do 4 rounds of polishing)  

GENOME=/path/to/genome.fasta

FORWARD=/path/to/reads1.fq

REVERSE=/path/to/reads2.fq

#index genome
bwa index ${GENOME}
    
bwa mem -M -t 48 \
         	${GENOME} ${FORWARD} ${REVERSE} > bwa.sam

samtools view -Sb bwa.sam  > bwa.bam         	

#Sort and index the BAM
samtools sort bwa.bam -o bwa.sort
samtools index bwa.sort

#Pilon it 
java -Xmx300G -jar /data/jmsutton/anaconda3/share/pilon-1.23-0/pilon-1.23.jar --genome ${GENOME}  --frags  bwa.sort --output r1 --outdir r1

###Quast Analysis

###BUSCO Analysis

###Generating Blob plots 


###Filtering Genomes with Whitelist from Blobtools 

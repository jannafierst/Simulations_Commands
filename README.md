### Simulations_Commands
###commands used for various analyses used for Sim paper 

###ART Commands

art_illumina -ss HS25 -sam -i ./$file -l 150 -f 20 -p -m 270 -s 30 -o PE

###NanoSim commands

/jlf/janna/bama_new/NanoSim-master/NanoSim-master/src/read_analysis.py -i ERR2092776.fasta -r caenorhabditis_elegans.PRJNA13758.WBPS12.genomic.fa -o celegans -t 4

/jlf/janna/bama_new/NanoSim-master/NanoSim-master/src/simulator.py linear -n 5000000 -r caenorhabditis_elegans.PRJNA13758.WBPS12.genomic.fa -c celegans

#Subsampling Illumina data (fixed number of reads)
#-s100 is the seed this MUST BE THE SAME for both commands for the reads to match

seqtk sample -s100 read1.fq 10000 > sub1.fq
seqtk sample -s100 read2.fq 10000 > sub2.fq

#Subsampling ONT data (proportion of reads)

seqtk sample -s100 read.fq 0.5 > sub.fq

#Canu Assembly

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

# Running Quast without a reference

python /data/jdmillwood/anaconda3/pkgs/quast-5.0.2/quast.py -t 12 --plots-format pdf  /data/jmsutton/356/356_flye/assembly.fasta -o ./flye

# Running Quast with a reference

python /data/jdmillwood/anaconda3/pkgs/quast-5.0.2/quast.py -t 12 --plots-format pdf -r Caenorhabditis_elegans.WBcel235.dna.toplevel.fa /data/jmsutton/simulations/celegans/2M/2M.contigs.fasta -o ./2M

###BUSCO Analysis

# A script to run BUSCO with the nematoda_odb10 dataset

#!/bin/bash

export PATH="/data/jlfierst/anaconda3/bin:$PATH"
export AUGUSTUS_CONFIG_PATH="/data/jlfierst/anaconda3/config/"

export BUSCO_CONFIG_FILE="/data/jlfierst/anaconda3/config/config.ini"

busco -c 12 -m genome -i /data/jmsutton/simulations/celegans/1M/1M.contigs.fasta -o 1M --lineage_dataset nematoda_odb10 

###SIDR

https://github.com/damurdock/SIDR

###Jellyfish
jellyfish count -m 21 -s 8G -t 10 -C -o 21mer /Path/To/Reads/reads.fasta

###Generating Blob plots 
=======
# Blobtoolkit
### Create Blob Directory

/path/to/blobtoolkit/blobtools2/blobtools create --fasta /path/to/genome.fasta --taxdump /path/to/blobtoolkit/taxdump AssemblyName

### Run blastn

/path/to/blastn -db nt -query genome.fasta -outfmt '6 qseqid staxids bitscore std'  -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 16 -out blast.out

### Run blastx

/path/to/diamond blastx --query genome.fasta --db reference_proteomes.fasta.gz --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 > diamond.out

### Map Illumina reads 

minimap2 -ax map-ont -t 16 genome.fasta /path/to/raw/reads/reads.fastq | samtools sort -@16 -O BAM -o assembly.reads.bam 

### Add all data files to Blob Directory

path/to/blobtools2/blobtools add --hits blast.out --hits diamond.out --taxrule bestsumorder --taxdump ~/taxdump --cov assembly.reads.bam --busco /path/to/busco/run/full_table.tsv AssemblyName

### Open BlobToolKit Viewer to view BlobPlot

path/to/blobtools2/blobtools host /path/to/working/blob/info/directory

###Filtering Genomes with Whitelist

## awk to split .fasta 

#!/bin/bash

awk '/^>scaffold/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}' final.assembly.fasta

cat scaffolds_to_keep.txt | while read line; do cat "scaffold_"$line".fasta" ; done > kept_scaffolds.fasta

rm "scaffold"*".fasta"

## checking the output

cat kept_scaffolds.fasta | grep '^>' | tr -d '\>scaffold\_' | sort -n | diff - scaffolds_to_keep.txt 

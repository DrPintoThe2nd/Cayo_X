#!bin/bash

#wget ftp://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
#wget ftp://ftp.ensembl.org/pub/release-104/gff3/macaca_mulatta/Macaca_mulatta.Mmul_10.104.gff3.gz
#pigz -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
#pigz -d Macaca_mulatta.Mmul_10.104.gff3.gz

#$mamba env create --name xyalign -f xyalign_environment.yml
#$conda activate xyalign
#run "$bash xyalign_index.sh Macaca_mulatta.Mmul_10.dna.toplevel.fa" to generate masked Y (plus PAR, if known add as "--reference_mask mask.bed") and bwa and hisat2 indices

cd cleaned_reads/; while read a b; do mv "$a" "$b"; done < names.txt;
mv PT101159_R*.fastq.gz dna/;
mv PT101210_R*.fastq.gz dna/;
mv PT102485_R*.fastq.gz dna/;
mv PT102842_R*.fastq.gz dna/;
mv PT103046_R*.fastq.gz dna/;
mv PT103352_R*.fastq.gz dna/;
mv PT103760_R*.fastq.gz dna/;
mv PT103811_R*.fastq.gz dna/;
mv PT104219_R*.fastq.gz dna/;
mv PT104270_R*.fastq.gz dna/;
mv PT104733_R*.fastq.gz dna/;
mv PT104883_R*.fastq.gz dna/;
mv PT105239_R*.fastq.gz dna/;
mv PT105290_R*.fastq.gz dna/;
mv *.fastq.gz rna/;
cd ..;

xyalign --PREPARE_REFERENCE --ref $1 --output_dir Mmul_10 --sample_id Mmul_10 --cpus 4 --x_chromosome X --y_chromosome Y #--reference_mask mask.bed 

bwa index Mmul_10/reference/xyalign_noY.masked.fa

hisat2-build -p 4 Mmul_10/reference/xyalign_noY.masked.fa Mmul_10/reference/xyalign_noY.masked

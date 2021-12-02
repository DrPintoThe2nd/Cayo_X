#!bin/bash

#wget ftp://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
#wget ftp://ftp.ensembl.org/pub/release-104/gff3/macaca_mulatta/Macaca_mulatta.Mmul_10.104.gff3.gz
#pigz -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
#pigz -d Macaca_mulatta.Mmul_10.104.gff3.gz

#$mamba env create --name xyalign -f xyalign_environment.yml
#$conda activate xyalign
#run "$xyalign.sh Macaca_mulatta.Mmul_10.dna.toplevel.fa" to generate masked Y (plus PAR, if known add as "--reference_mask mask.bed") and bwa and hisat2 indices

xyalign --PREPARE_REFERENCE --ref $1 --output_dir Mmul_10 --sample_id Mmul_10 --cpus 4 --x_chromosome X --y_chromosome Y #--reference_mask mask.bed 

bwa index Mmul_10/reference/xyalign_noY.masked.fa

hisat2-build -p 4 Mmul_10/reference/xyalign_noY.masked.fa xyalign_noY.masked


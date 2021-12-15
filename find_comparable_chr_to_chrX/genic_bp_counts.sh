#!bin/bash

echo "ensuring conda environment is active."
source activate Cayo

echo "downloading gff file."
wget http://ftp.ensembl.org/pub/release-105/gff3/macaca_mulatta/Macaca_mulatta.Mmul_10.105.gff3.gz;
echo "unzipping gff using pigz."
pigz -d Macaca_mulatta.Mmul_10.105.gff3.gz;

echo "converting to bed file using bedops."
gff2bed < Macaca_mulatta.Mmul_10.105.gff3 | grep 'gene' > Macaca_mulatta.Mmul_10.105.bed;
bedextract --list-chr Macaca_mulatta.Mmul_10.105.bed > chr.txt;

echo "extracting bp counts per chromosome..."
cat chr.txt | while read line
do
bedextract $line Macaca_mulatta.Mmul_10.105.bed | awk '/gene/ { len +=$2-$3} END {print len}' >> counts.txt
done
paste chr.txt counts.txt > size.txt
sed -i 's/\-//g' size.txt
echo "sorting output."
sort -n -r -k 2 size.txt > genic_chr_size.tsv
echo "Done. Output available in 'genic_chr_size.tsv'."
head -20 genic_chr_size.tsv

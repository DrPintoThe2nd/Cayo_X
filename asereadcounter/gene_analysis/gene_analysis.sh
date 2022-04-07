#!bin/bash

conda activate Cayo

#get gene names file from Ensembl BioMart, etc
cut -f 2 Mmul10_GRCh38.tsv | sort | uniq > gene_list.txt
cp gene_list.txt GTEx/

#add bed files to directory *_chrX_transcripts_rmdups_allelebalance.bed
ls *.bed > beds.txt 
sed -i 's/_chrX_transcripts_rmdups_allelebalance.bed$//g' beds.txt

cat beds.txt | while read line
do
mkdir $line
done

#nested for loop that extracts each gene from each experiment into its own file
var1=$(cat beds.txt)
var2=$(cat gene_list.txt)
for i in ${var1} ; do
  for x in ${var2} ; do
	grep "Name=${x}-" ${i}\_chrX_transcripts_rmdups_allelebalance.bed >> ${i}/${i}\_${x}\.txt
  done
done

#delete all empty files from all experiments
find . -type f -empty -print -delete 

#make genes file for each experiment
cat beds.txt | while read line
do 
ls $line\/$line\_*.txt > $line\/column1.txt 
sed -i "s/^$line\/$line\_//g" $line\/column1.txt  
sed -i 's/.txt$//g' $line\/column1.txt  
done

#makes column 2 for each experiment
cat beds.txt | while read line
do
cd $line
	cat column1.txt | while read line1 
	do
	cat *_$line1\.txt | datamash mean 19 >> column2.txt
	done
cd ..
done

#combine columns one and two for each experiment
cat beds.txt | while read line
do
paste $line\/column1.txt $line\/column2.txt > $line\/$line\_allele_balance.tsv 
done

#identify all genes expressed across all experiments
#monkeys
cat PT10*/column1.txt | sort | uniq > test.txt
sed -i 's/$/\$/' test.txt
egrep -v -f test.txt gene_list.txt > Macaque_missing.txt







#repeat for human data
cd GTEx;

#add bed files to directory *_chrX_transcripts_rmdups_allelebalance.bed
ls *.bed > beds.txt 
sed -i 's/_chrX_transcripts_rmdups_allelebalance.bed$//g' beds.txt

cat beds.txt | while read line
do
mkdir $line
done

#nested for loop that extracts each gene from each experiment into its own file
var1=$(cat beds.txt)
var2=$(cat gene_list.txt)
for i in ${var1} ; do
  for x in ${var2} ; do
	grep "Name=${x}-" ${i}\_chrX_transcripts_rmdups_allelebalance.bed >> ${i}/${i}\_${x}\.txt
  done
done

#delete all empty files from all experiments
find . -type f -empty -print -delete 

#make genes file for each experiment
cat beds.txt | while read line
do 
ls $line\/$line\_*.txt > $line\/column1.txt 
sed -i "s/^$line\/$line\_//g" $line\/column1.txt  
sed -i 's/.txt$//g' $line\/column1.txt  
done

#makes column 2 for each experiment
cat beds.txt | while read line
do
cd $line
	cat column1.txt | while read line1 
	do
	cat *_$line1\.txt | datamash mean 19 >> column2.txt
	done
cd ..
done

#combine columns one and two for each experiment
cat beds.txt | while read line
do
paste $line\/column1.txt $line\/column2.txt > $line\/$line\_allele_balance.tsv 
done

#humans
cat GTEX*/column1.txt | sort | uniq > test.txt
sed -i 's/$/\$/' test.txt
egrep -v -f test.txt gene_list.txt > ../Human_missing.txt


#calculate gene activation status per species, per tissue, etc
cd ..

mkdir activation_status

find . -type f -name '*_allele_balance.tsv' -exec cp {} ./activation_status/ \;

cd activation_status

wc -l PT10* > Cayo_wc.txt
wc -l GTEX_1* > GTEx_wc.txt

touch join.tmp
for file in PT10*_adrenal_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp PT_adrenal_allele_balance.tsv
awk '{print NF}' PT_adrenal_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in PT10*_gonads_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp PT_gonads_allele_balance.tsv
awk '{print NF}' PT_gonads_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in PT10*_heart_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp PT_heart_allele_balance.tsv
awk '{print NF}' PT_heart_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in PT10*_liver_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp PT_liver_allele_balance.tsv
awk '{print NF}' PT_liver_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in PT10*_lung_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp PT_lung_allele_balance.tsv
awk '{print NF}' PT_lung_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in PT10*_pituitary_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp PT_pituitary_allele_balance.tsv
awk '{print NF}' PT_pituitary_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in GTEX_*_adrenal_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp GTEX_adrenal_allele_balance.tsv
awk '{print NF}' GTEX_adrenal_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in GTEX_*_gonad_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp GTEX_gonads_allele_balance.tsv
awk '{print NF}' GTEX_gonads_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in GTEX_*_heart_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp GTEX_heart_allele_balance.tsv
awk '{print NF}' GTEX_heart_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in GTEX_*_liver_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp GTEX_liver_allele_balance.tsv
awk '{print NF}' GTEX_liver_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in GTEX_*_lung_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp GTEX_lung_allele_balance.tsv
awk '{print NF}' GTEX_lung_allele_balance.tsv | sort -nu | tail -n 1

touch join.tmp
for file in GTEX_*_pituitary_allele_balance.tsv; do
    join -a1 -a2 -t $'\t' -o auto -e NA join.tmp "$file" > join.tmp.1
    mv join.tmp.1 join.tmp
done
mv join.tmp GTEX_pituitary_allele_balance.tsv
awk '{print NF}' GTEX_pituitary_allele_balance.tsv | sort -nu | tail -n 1


echo "GTEX" > dataset.txt
echo "PT" >> dataset.txt

ls PT_*.tsv > samples.txt
sed -i 's/PT_//g' samples.txt
sed -i 's/.tsv//g' samples.txt


for VARIABLE in GTEX PT
do
	cat samples.txt | while read line
	do
	awk -F\NA '{print NF-1}' $VARIABLE\_$line\.tsv > tmp
	paste $VARIABLE\_$line\.tsv tmp > $VARIABLE\_$line\_counts.tsv
	done
done

mkdir R_output
Rscript activation_status.R


cd R_ouput/
ls *.tsv > files.txt
sed -i 's/.tsv//g' files.txt

cat files.txt | while read line
do
	egrep 'TRUE[[:space:]]FALSE' $line\.tsv > $line\_inactive.tsv
	egrep 'FALSE[[:space:]]FALSE' $line\.tsv > $line\_variable.tsv
	egrep 'FALSE[[:space:]]TRUE' $line\.tsv > $line\_escape.tsv
done

cut -f 1 Cayo_*_escape.tsv | sort | uniq > Cayo_escape_genes.txt
cut -f 1 GTEX_*_escape.tsv | sort | uniq > GTEX_escape_genes.txt

wc -l Cayo_escape_genes.txt
wc -l GTEX_escape_genes.txt

cut -f 1 Cayo_*_inactive.tsv | sort | uniq > Cayo_inactive_genes.txt
cut -f 1 GTEX_*_inactive.tsv | sort | uniq > GTEX_inactive_genes.txt

wc -l Cayo_inactive_genes.txt
wc -l GTEX_inactive_genes.txt

cut -f 1 Cayo_*_variable.tsv | sort | uniq > Cayo_variable_genes.txt
cut -f 1 GTEX_*_variable.tsv | sort | uniq > GTEX_variable_genes.txt

wc -l Cayo_variable_genes.txt
wc -l GTEX_variable_genes.txt


egrep -f Cayo_escape_genes.txt GTEX_escape_genes.txt > 		shared_escape_genes.txt
egrep -f Cayo_inactive_genes.txt GTEX_inactive_genes.txt > 	shared_inactive_genes.txt
egrep -f Cayo_variable_genes.txt GTEX_variable_genes.txt > 	shared_variable_genes.txt

wc -l shared_escape_genes.txt
wc -l shared_inactive_genes.txt
wc -l shared_variable_genes.txt














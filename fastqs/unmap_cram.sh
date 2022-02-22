#!bin/bash

ls *.cram > samples.txt
sed -i "s/.cram//g" samples.txt

cat samples.txt | while read line 
do
echo "cram to bam"
samtools view -@6 -b -T ../GCF_000001405.39_GRCh38.p13_genomic.fna -o $line\.bam $line\.cram

echo "unmapping them reads and pairing them"
samtools fastq -c9 -@6 -n -o $line\.fq.gz  $line\.bam;
repair.sh in=$line\.fq.gz out1=$line\_R1.fq.gz out2=$line\_R2.fq.gz outs=$line\_singletons.fq.gz repair;
done


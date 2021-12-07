import os
configfile: "Cayo_config.json"
dna_samples = config["dna"]
rna_samples = config["rna"]
genomes_to_use = config["Mmul10_XX"]
cleaned_reads = config["cleaned_reads"]
mapped_reads = config["mapped_reads"]
called_reads = config["called_reads"]
reference_dir = config["reference_dir"]
genome = config["Mmul10_XX"]
bams = config["bams"]

rule all:
	input:
#hisat2_map_reads rule
		expand("mapped_reads/rna/{sample}.bam", sample=rna_samples),
#bwa map_reads rule
		expand("mapped_reads/dna/{sample}.bam", sample=dna_samples),
#index_stat rule
		expand("mapped_reads/rna/{sample}.bam.bai", sample=rna_samples),
		expand("mapped_reads/rna/{sample}.bam.flagstat", sample=rna_samples),
		expand("mapped_reads/rna/{sample}.bam.stats", sample=rna_samples),
		expand("mapped_reads/dna/{sample}.bam.bai", sample=dna_samples),
		expand("mapped_reads/dna/{sample}.bam.flagstat", sample=dna_samples),
		expand("mapped_reads/dna/{sample}.bam.stats", sample=dna_samples),

# -----------------------------------------------------------------------------------
# -- note that this pipeline is currently configured assuming all samples are female
# -----------------------------------------------------------------------------------
# ------------------
# Read mapping rules
# ------------------

rule hisat2_map_reads:
	input:
		idx_check = expand(os.path.join(reference_dir + "xyalign_noY.masked.{suffix}.ht2"), suffix=["1", "2", "3", "4", "5", "6", "7", "8"]),
		fq1 = "cleaned_reads/rna/{sample}_R1.fastq.gz",
		fq2 = "cleaned_reads/rna/{sample}_R2.fastq.gz",
	output:
		"mapped_reads/rna/{sample}.bam",
	params:
		threads = 40,
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
	shell:
		"""
		hisat2 -p {params.threads} --dta --rg-id {params.id} --rg SM:{params.sm} --rg LB:{params.lb} --rg PU:{params.pu} --rg PL:{params.pl} -x Mmul_10/reference/xyalign_noY.masked -1 {input.fq1} -2 {input.fq2} | samtools sort -@4 -O bam - -o {output};
		"""

rule bwa_map_reads:
	input:
		idx_check = os.path.join(reference_dir + "xyalign_noY.masked.fa.bwt"),
		fq1 = "cleaned_reads/dna/{sample}_R1.fastq.gz", 
		fq2 = "cleaned_reads/dna/{sample}_R2.fastq.gz", 
	output:
		"mapped_reads/dna/{sample}.bam",
	params:
		threads = 40,
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		genome = genome,
	shell:
		"""
		bwa mem -t {params.threads} -R "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' {params.genome} {input.fq1} {input.fq2} | samblaster | samtools fixmate -@ 4 - | samtools sort -@ 4 -O bam - -o {output};
		"""

rule index_stat:
	input:
		rna_bam = expand("mapped_reads/rna/{sample}.bam", sample=rna_samples),
		dna_bam = expand("mapped_reads/dna/{sample}.bam", sample=dna_samples),
	output:
		rna_bai = expand("mapped_reads/rna/{sample}.bam.bai", sample=rna_samples),
		rna_flagstat = expand("mapped_reads/rna/{sample}.bam.flagstat", sample=rna_samples),
		rna_stats = expand("mapped_reads/rna/{sample}.bam.stats", sample=rna_samples),
		dna_bai = expand("mapped_reads/dna/{sample}.bam.bai", sample=dna_samples),
		dna_flagstat = expand("mapped_reads/dna/{sample}.bam.flagstat", sample=dna_samples),
		dna_stats = expand("mapped_reads/dna/{sample}.bam.stats", sample=dna_samples),
	shell:
		"""
		sambamba index -t4 {input.rna_bam};
		sambamba flagstat -t4 {input.rna_bam} > {output.rna_flagstat};
		samtools stats {input.rna_bam} | grep ^SN | cut -f 2- > {output.rna_stats};
		sambamba index -t4 {input.dna_bam};
		sambamba flagstat -t4 {input.dna_bam} > {output.dna_flagstat};
		samtools stats {input.dna_bam} | grep ^SN | cut -f 2- > {output.dna_stats};
		"""

#PERL5LIB=/home/bpinto2/miniconda3/envs/Cayo/lib/perl5/5.32/ hisat2 -p 1 --dta -x Mmul_10/reference/xyalign_noY.masked -1 cleaned_reads/rna/PT102441_R1.fastq.gz -2 cleaned_reads/rna/PT102441_R2.fastq.gz | samtools sort -@2 -O bam - > mapped/reads/rna/PT102441.bam;

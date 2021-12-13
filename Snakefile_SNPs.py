import os

configfile: "AGAVE_config.json"

dna_samples = config["dna"]
rna_samples = config["rna"]
genomes_to_use = config["Mmul10_XX"]
cleaned_reads = config["cleaned_reads"]
mapped_reads = config["mapped_reads"]
called_reads = config["called_reads"]
reference_dir = config["reference_dir"]
genome = config["Mmul10_XX"]

rule all:
	input:
#hisat2_map_reads rule
		expand("mapped_reads/rna/{sample}.bam", sample=rna_samples),
#bwa map_reads rule
		expand("mapped_reads/dna/{sample}.bam", sample=dna_samples),
#index_stat rule
		expand("mapped_reads/rna/{sample}.bam.bai", sample=rna_samples),
		expand("mapped_reads/rna/{sample}.bam.stats", sample=rna_samples),
		expand("mapped_reads/dna/{sample}.bam.bai", sample=dna_samples),
		expand("mapped_reads/dna/{sample}.bam.stats", sample=dna_samples),
#gatk_gvcf_chrx rule
		expand("called_reads/{sample}.chrx.g.vcf.gz", sample=dna_samples),
#gatk_combinegvcfs_chrx rule
		"merged_vcfs/Mmul10_XX.chrx.gatk.combineVCF.vcf.gz",
#gatk_genotypegvcf_chrx rule 
		"merged_vcfs/Mmul10_XX.chrx.gatk.called.raw.vcf.gz",		
#hard_filter rule
		"merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.vcf.gz",
#subset_individuals_hard_filter rule
                expand("merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.{sample}.vcf.gz", sample=dna_samples),
		expand("merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.{sample}.vcf.gz.tbi", sample=dna_samples),
#gatk_selectheterozygous rule
		expand("Xchr_vcfs/Mmul_10_XX.chrx.gatk.called.hard.filter.het.{sample}.vcf.gz", sample=dna_samples),
#gatk_gvcf_autosomes rule
		expand("called_autosomes/{sample}.chr{chr_n}.g.vcf.gz", sample=dna_samples, chr_n=config["autosomes"]),
#gatk_combinegvcfs_autosomes rule
		expand("merged_vcfs/Mmul_10.chr{chr_n}.gatk.combinegvcf.g.vcf.gz", chr_n=config["autosomes"]),
#gatk_genotypegvcf_autosomes rule
		expand("merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.raw.vcf.gz", chr_n=config["autosomes"]),
#hard_filter_autosomes rule
		expand("merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.vcf.gz", chr_n=config["autosomes"]),
                expand("merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.vcf.gz.tbi", chr_n=config["autosomes"]),
#subset_individuals_hard_filter_autosomes rule
		expand("merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi", sample=dna_samples, chr_n=config["autosomes"]),
#gatk_selectheterozygous_autosomes rule
		expand("autosomes_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz", sample=dna_samples, chr_n=config["autosomes"]),

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
		threads = 4,
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
		threads = 4,
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		genome = genome,
	shell:
		"""
		bwa mem -t {params.threads} -R '@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' {params.genome} {input.fq1} {input.fq2} | samblaster | samtools fixmate -@ 4 - - | samtools sort -@ 4 -O bam - -o {output};
		"""

rule index_rna:
	input:
		"mapped_reads/rna/{sample}.bam",
	output:
		"mapped_reads/rna/{sample}.bam.bai",
	shell:
		"""
		samtools index {input};
		"""

rule index_dna:
	input:
		"mapped_reads/dna/{sample}.bam",
	output:
		"mapped_reads/dna/{sample}.bam.bai",
	shell:
		"""
		samtools index {input};
		"""

rule stats_rna:
	input:
		"mapped_reads/rna/{sample}.bam",
	output:
		"mapped_reads/rna/{sample}.bam.stats",
	shell:
		"""
		samtools stats {input} | grep '^SN' | cut -f 2- > {output};
		"""

rule stats_dna:
	input:
		"mapped_reads/dna/{sample}.bam",
	output:
		"mapped_reads/dna/{sample}.bam.stats",
	shell:
		"""
		samtools stats {input} | grep '^SN' | cut -f 2- > {output};
		"""

# -----------------
# Call variants dna
# -----------------

rule gatk_gvcf_chrx:
	input:
		bam = "mapped_reads/dna/{sample}.bam",
		bai = "mapped_reads/dna/{sample}.bam.bai",
	output:
		"called_reads/{sample}.chrx.g.vcf.gz",
	params:
		x_chrom = config["chrx"],
		genome = genome,
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.genome} -I {input.bam} -L {params.x_chrom} --emit-ref-confidence GVCF -O {output};
		"""

rule gatk_combinegvcfs_chrx:
	input:
		gvcfs = expand("called_reads/{sample}.chrx.g.vcf.gz", sample=dna_samples),
	output:
		"merged_vcfs/Mmul10_XX.chrx.gatk.combineVCF.vcf.gz",
	params:
		genome = genome,
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.genome} {variant_files} -O {output}
			""")

rule gatk_genotypegvcf_chrx:
	input:
		gvcf = "merged_vcfs/Mmul10_XX.chrx.gatk.combineVCF.vcf.gz",
	output:
		vcf = "merged_vcfs/Mmul10_XX.chrx.gatk.called.raw.vcf.gz",
	params:
		genome = genome,
	threads:
		4
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.genome} -V {input.gvcf} -O {output.vcf};
		"""

# -----------
# Hard filter X
# -----------

rule hard_filter:
	input:
		"merged_vcfs/Mmul10_XX.chrx.gatk.called.raw.vcf.gz",
	output:
		vcf = "merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.vcf.gz",
		idx = "merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.vcf.gz.tbi",
	params:
		genome = genome,
	shell:
		"""
		gatk SelectVariants -R {params.genome} -V {input} -O {output} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -select "AN >= 4 && MQ > 40.0 && QD > 5.0 && DP >= 10.0 && DP <= 1000.0";
		touch -c {output.idx};
		"""

# ----------------------------------------
# Further processing VCF
# 1. Subset VCF files for each individual
# 2. Keep only the heterozygous sites
# ----------------------------------------

rule subset_individuals_hard_filter:
	input:
		"merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.vcf.gz",
	output:
		vcf = "merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.{sample}.vcf.gz",
		index = "merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.{sample}.vcf.gz.tbi",
	params:
		sample = lambda wildcards: config[wildcards.sample]["SM"],
	shell:
		"""
		bcftools view -Oz -s {params.sample} {input} > {output.vcf};
		tabix -p vcf {output.vcf};
		"""

#After subsetting for each individual. In some individuals, the genotypes could be homozygous for the reference. This next rule is to remove these sites.

rule gatk_selectheterozygous:
	input:
		"merged_vcfs/Mmul10_XX.chrx.gatk.called.hard.filter.{sample}.vcf.gz",
	output:
		"Xchr_vcfs/Mmul_10_XX.chrx.gatk.called.hard.filter.het.{sample}.vcf.gz",
	params:
		genome = genome,
	shell:
		"""
		gatk SelectVariants -R {params.genome} -V {input} -O {output} -select "AC == 1"
		"""

# ---------------------------------------
# Processing of autosomes (chr1 to chr20)
# ---------------------------------------

rule gatk_gvcf_autosomes:
	input:
		bam = "mapped_reads/dna/{sample}.bam",
		bai = "mapped_reads/dna/{sample}.bam.bai",
	output:
		"called_autosomes/{sample}.chr{chr_n}.g.vcf.gz",
	params:
		genome = genome,
		chr_n = "{chr_n}",
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.genome} -I {input.bam} -L {params.chr_n} --emit-ref-confidence GVCF --output {output};
		"""

rule gatk_combinegvcfs_autosomes:
	input:
		gvcfs = expand("called_autosomes/{sample}.chr{{chr_n}}.g.vcf.gz", zip, sample=dna_samples, chr_n=config["autosomes"]),
	output:
		"merged_vcfs/Mmul_10.chr{chr_n}.gatk.combinegvcf.g.vcf.gz",
	params:
		genome = genome,
		chr_n = "{chr_n}",
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.genome} {variant_files} --intervals {params.chr_n} -O {output}
			""")

rule gatk_genotypegvcf_autosomes:
	input:
		expand("merged_vcfs/Mmul_10.chr{{chr_n}}.gatk.combinegvcf.g.vcf.gz", chr_n=config["autosomes"]),
	output:
		"merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.raw.vcf.gz",
	params:
		genome = genome,
	threads:
		4
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.genome} -V {input} -O {output}
		"""

# ------------------------
# Hard Filter on autosomes
# ------------------------

rule hard_filter_autosomes:
	input:
		"merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.raw.vcf.gz"
	output:
		vcf = "merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.vcf.gz",
		idx = "merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.vcf.gz.tbi",
	params:
		genome = genome,
	shell:
		"""
		gatk SelectVariants -R {params.genome} -V {input} -O {output} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -select "AN >= 4 && MQ > 40.0 && QD > 7.0 && DP >= 10.0 && DP <= 1000.0";
		touch -c {output.idx};
		"""

rule subset_individuals_hard_filter_autosomes:
	input:
		"merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.vcf.gz",
	output:
		vcf = "merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz",
		index = "merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi",
	params:
		sample = lambda wildcards: config[wildcards.sample]["SM"],
	shell:
		"""
		bcftools view -Oz -s {params.sample} {input} > {output.vcf};
		tabix -p vcf {output.vcf};
		"""

#After subsetting for each individual. In some individuals, the genotypes could be homozygous for the reference. This next rule is to remove these sites.

rule gatk_selectheterozygous_autosomes:
	input:
		"merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz", 
	output:
		"autosomes_vcfs/Mmul_10.chr{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz",
	params:
		genome = genome,
	shell:
		"""
		gatk SelectVariants -R {params.genome} -V {input} -O {output} -select "AC == 1";
		"""

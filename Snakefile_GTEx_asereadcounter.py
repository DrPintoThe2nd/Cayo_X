import os

configfile: "AGAVE_GTEx_config.json"

dna_samples = config["dna"]
rna_samples = config["rna"]
genomes_to_use = config["GRCh38_XX"]
cleaned_reads = config["cleaned_reads"]
mapped_reads = config["mapped_reads"]
called_reads = config["called_reads"]
reference_dir = config["reference_dir"]
genome = config["GRCh38_XX"]

combiC = []
for key in config["dna_rna"]:
	for item in config["dna_rna"][key]:
		combiC.append((key, item))
		
import itertools
combiList=list()
for c in combiC:
	combiList.append(c[0]+"-"+c[1])
#	print(combiList)

rule all:
	input:
		expand("asereadcounter/{combo}_Xchr.tsv", combo = combiList), #run asereadcounter
		expand("asereadcounter/{combo}_chr{chr_n}.tsv", combo = combiList, chr_n=config["autosomes"])

# ----
# chrX
# ----
# Run ASEReadCounter X
rule gatk_asereadcounter_minDepth10_X:
	input:
		bam = "mapped_reads/rna/{rna}.bam",
		sites = "Xchr_vcfs/GRCh38_XX.chrx.gatk.called.hard.filter.het.{dna}.vcf.gz",
	output:
		"asereadcounter/{dna}-{rna}_Xchr.tsv",
	params:
		genome = genome,
	shell:
		"""
		gatk ASEReadCounter -R {params.genome} -O {output} -I {input.bam} -V {input.sites} -min-depth 10 --min-mapping-quality 10 --min-base-quality 10 
		"""

# ----
# autosomes
# ----# Run ASEReadCounter autosomes
rule gatk_asereadcounter_minDepth10_autosomes:
	input:
		bam = "mapped_reads/rna/{rna}.bam",
		sites = "autosomes_vcfs/GRCh38.chr{chr_n}.gatk.called.hard.filter.het.{dna}.vcf.gz",
	output:
		"asereadcounter/{dna}-{rna}_chr{chr_n}.tsv"
	params:
		genome = genome,
	shell:
		"""
		gatk ASEReadCounter -R {params.genome} -O {output} -I {input.bam} -V {input.sites} -min-depth 10 --min-mapping-quality 10 --min-base-quality 10 
		"""

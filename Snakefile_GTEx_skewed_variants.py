import os

configfile: "AGAVE_GTEx_config.json"

rule all:
	input:
		expand("asereadcounter/{sample}/{sample}_{tissue}_chrX.bed", sample=config["sample_ids"], tissue=config["all_tissues"]),
		expand("asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["autosomes"]),
		expand("asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcript_stats.txt", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),

rule convert_chrX_data_to_bed:
	input:
		"asereadcounter/{sample}/{sample}_{tissue}_Xchr.tsv",
	output:
		"asereadcounter/{sample}/{sample}_{tissue}_chrX.bed",
	params:
		script = config["convert_unphased_data_to_bed_script"],
	shell:
		"""
		python {params.script} --unphased_data {input} --outfile {output};
		"""

rule convert_autosomes_data_to_bed:
	input:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.tsv",
	output:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed",
	params:
		script = config["convert_unphased_data_to_bed_script"],
	shell:
		"""
		python {params.script} --unphased_data {input} --outfile {output};
		"""

rule bedtools_intersect_data:
	input:
		transcripts = "find_comparable_chr_to_X/Homo_sapiens.GRCh38.105.transcripts.bed",
		asereadcounter = "asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed",
	output:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed",
	shell:
		"""
		bedtools intersect -a {input.transcripts} -b {input.asereadcounter} -wa -wb > {output};
		"""

rule remove_duplicates:
	input:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed",
	output:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed",
	params:
		script = config["remove_duplicates_script"],
	shell:
		"""
		python {params.script} {input} {output};
		"""

rule compute_allele_balance:
	input:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed",
	output:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed",
	params:
		script = config["compute_allele_balance_script"],
	shell:
		"""
		python {params.script} {input} {output};
		"""

rule tabulate_statistics:
	input:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed",
	output:
		"asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcript_stats.txt",
	params:
		script = config["tabulate_statistics_script"],
		sample = "{sample}",
		tissue = "{tissue}",
		chrom = "chr{chr_n}",
	shell:
		"""
		python {params.script} --sampleID {params.sample} --tissue {params.tissue} --chrom {params.chrom} --infile {input} --outfile {output};
		"""


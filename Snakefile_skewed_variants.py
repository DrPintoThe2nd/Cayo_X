import os
configfile: "AGAVE_config.json"

rule all:
	input:
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chrX.bed", sample=config["sample_ids"], tissue=config["all_tissues"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chrX.bed", sample=config["missing_lung_samples"], tissue=config["missing_lung_tissues"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chrX.bed", sample=config["missing_gonads_samples"], tissue=config["missing_gonads_tissues"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["autosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed", sample=config["missing_lung_samples"], tissue=config["missing_lung_tissues"], chr_n=config["autosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed", sample=config["missing_gonads_samples"], tissue=config["missing_gonads_tissues"], chr_n=config["autosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed", sample=config["missing_lung_samples"], tissue=config["missing_lung_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed", sample=config["missing_gonads_samples"], tissue=config["missing_gonads_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed", sample=config["missing_lung_samples"], tissue=config["missing_lung_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed", sample=config["missing_gonads_samples"], tissue=config["missing_gonads_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed", sample=config["missing_lung_samples"], tissue=config["missing_lung_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed", sample=config["missing_gonads_samples"], tissue=config["missing_gonads_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcript_stats.txt", sample=config["sample_ids"], tissue=config["all_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcript_stats.txt", sample=config["missing_lung_samples"], tissue=config["missing_lung_tissues"], chr_n=config["all_chromosomes"]),
		expand("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcript_stats.txt", sample=config["missing_gonads_samples"], tissue=config["missing_gonads_tissues"], chr_n=config["all_chromosomes"]),

rule convert_chrX_data_to_bed:
	input:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_Xchr.tsv",
	output:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chrX.bed",
	params:
		script = config["convert_unphased_data_to_bed_script"],
	shell:
		"""
		python {params.script} --unphased_data {input} --outfile {output};
		"""

rule convert_autosomes_data_to_bed:
	input:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.tsv",
	output:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed",
	params:
		script = config["convert_unphased_data_to_bed_script"],
	shell:
		"""
		python {params.script} --unphased_data {input} --outfile {output};
		"""

rule bedtools_intersect_data:
	input:
		transcripts = "/data/CEM/wilsonlab/lab_generated/cayo/find_comparable_chr_to_chrX/Macaca_mulatta.Mmul_10.105.transcripts.bed",
		asereadcounter = "/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}.bed",
	output:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed",
	shell:
		"""
		bedtools intersect -a {input.transcripts} -b {input.asereadcounter} -wa -wb > {output};
		"""

rule remove_duplicates:
	input:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts.bed",
	output:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed",
	params:
		script = config["remove_duplicates_script"],
	shell:
		"""
		python {params.script} {input} {output};
		"""

rule compute_allele_balance:
	input:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups.bed",
	output:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed",
	params:
		script = config["compute_allele_balance_script"],
	shell:
		"""
		python {params.script} {input} {output};
		"""

rule tabulate_statistics:
	input:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcripts_rmdups_allelebalance.bed",
	output:
		"/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/{sample}/{sample}_{tissue}_chr{chr_n}_transcript_stats.txt",
	params:
		script = config["tabulate_statistics_script"],
		sample = "{sample}",
		tissue = "{tissue}",
		chrom = "chr{chr_n}",
	shell:
		"""
		python {params.script} --sampleID {params.sample} --tissue {params.tissue} --chrom {params.chrom} --infile {input} --outfile {output};
		"""


#Environment: CustomFilterVCF
import os

configfile: "AGAVE_config.json"

rule all:
	input:
#extract_Xstats:
		"vcf_stats/Xchr_stats.txt",
#extract_stats:
		expand("vcf_stats/chr{chr_n}_stats.txt", chr_n=config["autosomes"]),
#plot_Xstats:
		"vcf_stats/Xchr_AN_plot.png",
		"vcf_stats/Xchr_QD_plot.png",
		"vcf_stats/Xchr_DP_plot.png",
		"vcf_stats/Xchr_MQ_plot.png",
		"vcf_stats/Xchr_DP_statistics.txt", 
#plot_stats:
		expand("vcf_stats/chr{chr_n}_AN_plot.png", chr_n=config["autosomes"]),
		expand("vcf_stats/chr{chr_n}_QD_plot.png", chr_n=config["autosomes"]),
		expand("vcf_stats/chr{chr_n}_DP_plot.png", chr_n=config["autosomes"]),
		expand("vcf_stats/chr{chr_n}_MQ_plot.png", chr_n=config["autosomes"]),
		expand("vcf_stats/chr{chr_n}_DP_statistics.txt", chr_n=config["autosomes"]),

rule extract_Xstats:
	input:
		xvcf = "merged_vcfs/Mmul10_XX.chrx.gatk.called.raw.vcf.gz",
	output:
		xstats = "vcf_stats/Xchr_stats.txt",
	params:
		AN = "AN",
		QD = "QD",
		MQ = "MQ",
		DP = "DP"
	shell:
		"""
		python extract_stats_from_vcf.py {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.xvcf} --outfile {output.xstats};
		"""

rule extract_stats:
	input:
		"merged_vcfs/Mmul_10.chr{chr_n}.gatk.called.raw.vcf.gz",
	output:
		"vcf_stats/chr{chr_n}_stats.txt",
	params:
		AN = "AN",
		QD = "QD",
		MQ = "MQ",
		DP = "DP"
	shell:
		"""
		python extract_stats_from_vcf.py {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input} --outfile {output};
		"""

rule plot_Xstats:
	input:
		"vcf_stats/Xchr_stats.txt",
	output:
		xAN_plot = "vcf_stats/Xchr_AN_plot.png",
		xQD_plot = "vcf_stats/Xchr_QD_plot.png",
		xDP_plot = "vcf_stats/Xchr_DP_plot.png",
		xMQ_plot = "vcf_stats/Xchr_MQ_plot.png",
		xDP_statistics = "vcf_stats/Xchr_DP_statistics.txt",
	shell:
		"""
		Rscript plot_stats.R {input} {output.xAN_plot} {output.xQD_plot} {output.xDP_plot} {output.xMQ_plot} {output.xDP_statistics};
		"""

rule plot_stats:
	input:
		"vcf_stats/chr{chr_n}_stats.txt",
	output:
		AN_plot = "vcf_stats/chr{chr_n}_AN_plot.png",
		QD_plot = "vcf_stats/chr{chr_n}_QD_plot.png",
		DP_plot = "vcf_stats/chr{chr_n}_DP_plot.png",
		MQ_plot = "vcf_stats/chr{chr_n}_MQ_plot.png",
		DP_statistics = "vcf_stats/chr{chr_n}_DP_statistics.txt",
	shell:
		"""
		Rscript plot_stats.R {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics};
		"""

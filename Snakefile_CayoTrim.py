import os

configfile: "Cayo_config.json"

sample = config["samples_ALL"]
raw_reads = config["raw_reads"]
raw_fastqc = config["raw_fastqc"]
trimmed_reads = config["trimmed_reads"]
cleaned_reads = config["cleaned_reads"]
cleaned_fastqc = config["cleaned_fastqc"]

#Misc notes/code (to be ignored):
#mamba create -n Cayo gatk bbmap fastqc multiqc hisat2 samtools bcftools vcftools r-base trim-galore snakemake bwa
#conda env export > Cayo_environment.yml
#PERL5LIB=/home/bpinto2/miniconda3/envs/Cayo/lib/perl5/5.32/

#To run this pipeline on any machine running linux, run
#git clone https://github.com/DrPintoThe2nd/Cayo_X.git
#where the raw fastq.gz files from SRA are placed in the "fastqs" directory

#Then, edit this line of code to include the path to your current directory with the above directories nested within as <path-to-working-directory>, note that
# forward slashed are needed to induce literal interpretations of backslashes in the path
#$pwd
#$sed 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' AGAVE_CayoTrim_config.json > CayoTrim_config.json

#Lastly, create (assuming you have mamba installed) and activate the conda env and conduct a dry run of the pipeline
#$mamba env create --name Cayo -f Cayo_environment.yml
#$conda activate Cayo
#$snakemake --use-conda -np -s Snakefile_CayoTrim.py

rule all:
	input:
#fastqc_analysis rule
		expand(os.path.join(raw_fastqc, "{sample}_R1_001_fastqc.html"), sample=sample),
		expand(os.path.join(raw_fastqc, "{sample}_R2_001_fastqc.html"), sample=sample),
#trim_galore_pe rule
		expand(os.path.join(trimmed_reads + "{sample}_R1_001_val_1.fq.gz"), sample=sample),
		expand(os.path.join(trimmed_reads + "{sample}_R2_001_val_2.fq.gz"), sample=sample),
#remove_PCR_dups rule
		expand(os.path.join(cleaned_reads + "{sample}_R1.fastq.gz"), sample=sample),
		expand(os.path.join(cleaned_reads + "{sample}_R2.fastq.gz"), sample=sample),
#fastqc_cleaned rule
		expand(os.path.join(cleaned_fastqc, "{sample}_R1_fastqc.html"), sample=sample),
		expand(os.path.join(cleaned_fastqc, "{sample}_R2_fastqc.html"), sample=sample),
#multiqc_all rule
		expand(os.path.join(raw_fastqc + "multiqc_raw/multiqc_report.html"), sample=sample),
		expand(os.path.join(trimmed_reads + "multiqc_trimmed/multiqc_report.html"), sample=sample),
		expand(os.path.join(cleaned_fastqc + "multiqc_cleaned/multiqc_report.html"), sample=sample),

#runs fastqc on all samples listed in JSON file with file names ending in "_R1[/2]_001.fastq.gz"
rule fastqc_analysis:	
	input:
		fqc_in1 = os.path.join(raw_reads + "{sample}_R1_001.fastq.gz"),
		fqc_in2 = os.path.join(raw_reads + "{sample}_R2_001.fastq.gz"),
	output:
		os.path.join(raw_fastqc + "{sample}_R1_001_fastqc.html"),
		os.path.join(raw_fastqc + "{sample}_R2_001_fastqc.html"),
	params:
		raw_fastqc = raw_fastqc
	shell:
		"""
		fastqc --threads 4 {input.fqc_in1} -o {params.raw_fastqc};
		fastqc --threads 4 {input.fqc_in2} -o {params.raw_fastqc}
		"""

#trims reads using default Trim Galore! parameters inclusing adapters and low-quality regions, simultaniously runs fastqc on each trimmed file
rule trim_galore_pe:
	input:
		trim = [os.path.join(raw_reads + "{sample}_R1_001.fastq.gz"), os.path.join(raw_reads + "{sample}_R2_001.fastq.gz")],
	output:
		os.path.join(trimmed_reads + "{sample}_R1_001_val_1.fq.gz"),
		os.path.join(trimmed_reads + "{sample}_R1_001_val_1_fastqc.html"),
		os.path.join(trimmed_reads + "{sample}_R2_001_val_2.fq.gz"),
		os.path.join(trimmed_reads + "{sample}_R2_001_val_2_fastqc.html"),
	params:
		trimmed_reads = trimmed_reads
	shell:
		"""
		trim_galore --paired --cores 4 --fastqc {input.trim} -o {params.trimmed_reads};
		"""

#removes duplicated read pairs using bbmap and compresses read files by re-ordering them
rule remove_PCR_dups:
	input:
		in1 = os.path.join(trimmed_reads + "{sample}_R1_001_val_1.fq.gz"),
		in2 = os.path.join(trimmed_reads + "{sample}_R2_001_val_2.fq.gz"),
	output:
		out1 = os.path.join(cleaned_reads + "{sample}_R1.fastq.gz"),
		out2 = os.path.join(cleaned_reads + "{sample}_R2.fastq.gz"),
	params:
		cleaned_reads = cleaned_reads
	shell:
		"""
		clumpify.sh -Xmx4g dedupe=t tmpdir={params.cleaned_reads} in={input.in1} in2={input.in2} out={output.out1} out2={output.out2}
		"""

#runs a final fastqc on the unduplicated read pairs from bbmap
rule fastqc_cleaned:	
	input:
		fqc_in1 = os.path.join(cleaned_reads + "{sample}_R1.fastq.gz"),
		fqc_in2 = os.path.join(cleaned_reads + "{sample}_R2.fastq.gz"),
	output:
		os.path.join(cleaned_fastqc + "{sample}_R1_fastqc.html"),
		os.path.join(cleaned_fastqc + "{sample}_R2_fastqc.html"),
	params:
		cleaned_fastqc = cleaned_fastqc
	shell:
		"""
		fastqc --threads 4 {input.fqc_in1} -o {params.cleaned_fastqc};
		fastqc --threads 4 {input.fqc_in2} -o {params.cleaned_fastqc}
		"""

#runs multiqc for all types of reads used and generated above
rule multiqc_all:
	input:
		expand(os.path.join(raw_fastqc + "{sample}_R1_001_fastqc.html"), sample=sample),
		expand(os.path.join(trimmed_reads + "{sample}_R1_001_val_1_fastqc.html"), sample=sample),
		expand(os.path.join(cleaned_fastqc + "{sample}_R1_fastqc.html"), sample=sample),
	output:
		os.path.join(raw_fastqc + "multiqc_raw/multiqc_report.html"),
		os.path.join(trimmed_reads + "multiqc_trimmed/multiqc_report.html"),
		os.path.join(cleaned_fastqc + "multiqc_cleaned/multiqc_report.html"),
	params:
		raw_fastqc = raw_fastqc ,
		trimmed_reads = trimmed_reads ,
		cleaned_fastqc = cleaned_fastqc ,
	shell:
		"""
		multiqc -f -o {params.raw_fastqc}/multiqc_raw {params.raw_fastqc};
		multiqc -f -o {params.trimmed_reads}/multiqc_trimmed {params.trimmed_reads};
		multiqc -f -o {params.cleaned_fastqc}/multiqc_cleaned {params.cleaned_fastqc}
		"""

# Cayo Macaque XCI and GTEx XCI

To run this pipeline on any machine running linux, run 

(Step 1)

```

git clone https://github.com/DrPintoThe2nd/Cayo_X.git

```

where the raw fastq.gz files from SRA are placed in the "fastqs" directory. 

Then, edit this line of code to include the path to your current directory with the above directories nested within as <path-to-working-directory>, note that forward slashes are needed to induce literal interpretations of backslashes in the path:

(Step 2)
 
```
 
pwd #get current path for next step
 
sed 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' AGAVE_config.json > Cayo_config.json

```

Then, create (assuming you have mamba installed) and activate the conda env and conduct a dry run of the pipeline (full pipeline assumes 4 cores):

(Step 3)

```

 mamba env create --name Cayo -f Cayo_environment.yml

 conda activate Cayo

 snakemake --use-conda -np -s Snakefile_CayoTrim.py #dry-run, see 'Snakemake_CayoTrim.sbatch' for example script for your machine

 conda deactivate 
 
```

To prepare the reference genome for mapping, use these commands to download and index the Y-masked genome assembly (assumes 4 available cores and pigz installed). Of note, xyalign and snakemake are incompatible packages, hence swapping conda environments and escaping the snakemake pipeline. Importantly, this step also generates a folder and subfolders within a directory called "Mmul10" which is called for in later steps.
 
(Step 4)
 
```
 
 wget ftp://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
 
 pigz -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz

 mamba env create --name xyalign -f xyalign_environment.yml
 
 conda activate xyalign
 
 bash xyalign_index_rename.sh Macaca_mulatta.Mmul_10.dna.toplevel.fa
 
 conda deactivate

```
 
To map reads and call SNPs using this new genome, we'll use snakemake again within our original Cayo conda env:

(Step 5)
 
```
 
conda activate Cayo
 
snakemake --use-conda -np -s Snakefile_SNPs.py #dry-run, see 'Snakemake_SNPs.sbatch' for example script for your machine
 
conda deactivate
 
```
 
To run GATK's ASEReadcounter, then rename the output into directly interpretable outputs:
 
(Step 6)
 
 ```
 
conda activate Cayo
 
snakemake --use-conda -np -s Snakefile_asereadcounter.py #dry-run, see 'Snakemake_asereadcounter.sbatch' for example script for your machine
 
bash asereadcounter/rename_asereadcounter.sh

snakemake --use-conda -np -s Snakefile_skewed_variants.py #dry-run, see 'Snakemake_skewed_variants.sbatch' for example script for your machine


#edit this path in-line for your machine
sed -i 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' asereadcounter/compute_mean_median_ratio_all_chr_macaques.R
 
bash asereadcounter/run_compute_mean_median_ratio_all_chr_macaques.sh
 
conda deactivate
 
 ```
 
 (Step 7)
 Plotting mean-median ratios across autosomes and the X chromosome.
 
 ```
 sed -i 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' asereadcounter/plot_violin_plot_allindividuals_alltissues_autosomes.R;
 sed -i 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' asereadcounter/plot_violin_plot_allindividuals_alltissues_chrX.R;
 sed -i 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' asereadcounter/plot_violin_plot_allindividuals_alltissues_chrX_compared_autosomes.R;
 
 Rscript asereadcounter/plot_violin_plot_allindividuals_alltissues_autosomes.R;
 Rscript asereadcounter/plot_violin_plot_allindividuals_alltissues_chrX.R;
 Rscript asereadcounter/plot_violin_plot_allindividuals_alltissues_chrX_compared_autosomes.R;
 
 ```
 
 (Step 8)
 Replicate analysis in humans using GTEx data (controlled access data).
 
 First, you need to download the read files. If downloading CRAM files for WGS, you'll need to convert to BAM (using GRCh38), strip the reads, and repair the reads before beginning. Simple scripts for this are available in the "fastqs" folder.
 
 Next step is to replicate the Cayo Macaque pipeline for the GTEx data. Starting with trimming reads:
 ```
 pwd #get current path for next step
 
 sed -i 's/\/scratch\/bpinto2\/Cayo_GTEx/<path-to-working-directory>/g' AGAVE_GTEx_config.json
 
 conda activate Cayo
 
 snakemake --use-conda -np -s Snakefile_GTExTrim.py #dry-run, see 'Snakemake_CayoTrim.sbatch' for example script for your machine

 conda deactivate 
 
 ```
 Next, prepare the reference genome:
 
 ```
 conda activate Cayo
 
 wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
 
 pigz -d Homo_sapiens.GRCh38.dna.toplevel.fa.gz
 
 grep 'REF' Homo_sapiens.GRCh38.dna.toplevel.fa | head -29 > name.lst
 
 seqtk subseq Homo_sapiens.GRCh38.dna.toplevel.fa name.lst > Homo_sapiens.GRCh38.dna.chromosomes.fa

 conda deactivate 
 
 conda activate xyalign
 
 bash xyalign_index.sh Homo_sapiens.GRCh38.dna.chromosomes.fa
 
 conda deactivate
 
 ```
 

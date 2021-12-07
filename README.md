# Cayo_X

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

 snakemake --use-conda -np -s Snakefile_CayoTrim.py

 conda deactivate 
 
```

To prepare the reference genome for mapping, use these commands to download and index the Y-masked genome assembly (assumes 4 available cores and pigz installed). Of note, xyalign and snakemake are incompatible packages, hence swapping conda environments and escaping the snakemake pipeline. Importantly, this step also generates a folder and subfolders within a directory called "Mmul10" which is called for in later steps.
 
(Step 4)
 
```
 
 wget ftp://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
 
 pigz -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz

 mamba env create --name xyalign -f xyalign_environment.yml
 
 conda activate xyalign
 
 bash xyalign_index.sh Macaca_mulatta.Mmul_10.dna.toplevel.fa
 
 conda deactivate

```
 
To map reads and call SNPs using this new genome, we'll use snakemake again within our original Cayo conda env:

(Step 5)
 
```
 
conda activate Cayo
 

 
```
 

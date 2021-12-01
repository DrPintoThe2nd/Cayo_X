# Cayo_X

To run this pipeline on any machine running linux, run

git clone https://github.com/DrPintoThe2nd/Cayo_X.git

where the raw fastq.gz files from SRA are placed in the "fastqs" directory. Then, edit this line of code to include the path to your current directory with the above directories nested within as <path-to-working-directory>, note that forward slashed are needed to induce literal interpretations of backslashes in the path


 """
 
 pwd
 
 sed 's/\/data\/CEM\/wilsonlab\/lab_generated\/cayo/<path-to-working-directory>/g' AGAVE_CayoTrim_config.json > CayoTrim_config.json

 """

Lastly, create (assuming you have mamba installed) and activate the conda env and conduct a dry run of the pipeline:

"""

 mamba env create --name Cayo -f Cayo_environment.yml

 conda activate Cayo

 snakemake --use-conda -np -s Snakefile_CayoTrim.py

 """

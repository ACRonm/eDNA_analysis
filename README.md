# eDNA_analysis
eDNA analysis repo for ETT2024

# Setup
Download your fasta.gz files and put them into /data/ETT_12_S_fasta

# Running the scripts
## R
The R script needs to be moved into the main directory for most of the functions to work, I will work on fixing this.

You can then run the "eDNA_processing_stage2.r script, which will output the species heat tree and the meta.csv and species.csv files.

## Python
Running the main.py script will do some analysis based on the raw data that has been processed from the fasta files. So far functionality is limited to outputting correlation pairs for genetic diversity and Heavy metal concentrations, which will be saved as a .csv to the /data directory.

You'll need to have python installed and be familiar with installing libraries. If you would like, you're welcome to reach out and I can help.

I am working on doing the same for nutrients.

This should give interesting information.

Aiden

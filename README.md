# Scripts employed to analyse associations between rDNA intragenomic frequencies and traits in the UK Biobank

The scripts compiled here have been employed to analyse the data and generate the figures 
for the manuscript: 
"Germline sequence variation within the ribosomal DNA is associated with human complex traits"

The folder `scripts_analysis` includes shell scripts employed to obtain
the intragenomic variant frequency (IGF) estimates and their association
with traits from the UK Biobank (UKB) data, plus a bed file with the coordinates
of the regions of the Hg38 assembly identified as analogue to the rDNA 
transcriptional unit.
The scripts in this folder are sequential (`s01` to `s06`), and they include:
- `s01_extract_rdna_rap.sh`: extraction of reads mapping to the rDNA analogue 
regions from UKB `cram` files. Intended to be executed in the UKB Research Analysis Platform (RAP).
- `s02_convert.sh`: conversion of `bam` file constructed in the previous step to `fastq`.
- `s03_trim.sh`: fixing and trimming of the paired-end `fastq` files obtained on the previous step.
- `s04_realign.sh`: alignment with `bowtie2` to the rDNA reference unit of the 
trimmed `fastq` files from the previous step.
- `s05_lofreq.sh`: estimation of IGFs from the `bam` files obtained from `bowtie2`.
- `s06_phesant.sh`: association between IGFs and UKB traits.

The data generated from executing the scripts above was processed and visualised using
the `R` notebooks included on the `scripts_figures` folder. 

Finally, the `var_extract` folder includes the in-house `python` script used to obtain 
rDNA variant frequencies from `bam` files on samples unsuitable for `lofreq`, 
such as RNA-seq data.

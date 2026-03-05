phesant_fd=$1
ifd=$2
ofd=$3
variant=$4
num_parts=$5
part=$6

wd=$(pwd)

phenotypes="${ifd}pheno_unrel.csv"
afs="${ifd}phesant_afs.csv"
covariates="${ifd}phesant_cov_unrel.csv"

echo "Running: ${variant} part ${part}"

mkdir -p ${ofd}

Rscript ${phesant_fd}/WAS/phenomeScan.r \
    --phenofile="${wd}/${phenotypes}" \
    --variablelistfile="${phesant_fd}variable-info/outcome-info.tsv" \
    --datacodingfile="${phesant_fd}variable-info/data-coding-ordinal-info.txt" \
    --traitofinterestfile="${wd}/${afs}" \
    --traitofinterest=${variant} \
    --confounderfile="${wd}/${covariates}" \
    --genetic="FALSE" \
    --resDir="${wd}${ofd}" \
    --userId="eid" \
    --numParts="${num_parts}" \
    --partIdx="${part}"

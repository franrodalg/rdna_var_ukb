ifd=$1
ofd=$2

ifn=$(ls ${ifd}/*.bam)
sample=`basename ${ifn} .bam`

echo "${sample}: Started"

mkdir -p $ofd
ofd=$(realpath ${ofd})

echo "${sample}: Converting"

cd $ifd

samtools fastq -@ $NSLOTS ${ifn} \
    -1 ${ofd}/${sample}_1.fq.gz -2 ${ofd}/${sample}_2.fq.gz -n

echo "${sample}: Finished"

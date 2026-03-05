ifd=$1
ofd=$2
gfd=$3

sample=$( basename $(ls ${ifd}/*_1.fq.gz) _1_val_1_.fz.gz)

echo "${sample}: Started"

mkdir -p $ofd
ofd=$(realpath ${ofd})
ofn="${ofd}/${sample}.bam"

cd $ifd

R1="${sample}_1_val_1.fq.gz"
R2="${sample}_2_val_2.fq.gz"

echo "${sample}: Aligning"

bowtie2 -x ${gfd}KY962518_looped_2120 \
    --no-unal \
    -p $NSLOTS \
    -1 $R1 -2 $R2 | \
    samtools view -bS \
    -o ${ofn} \
    -@ $NSLOTS -

echo "${sample}: Sorting"

( samtools sort -o ${ofn}.sorted -@ $NSLOTS ${ofn} &&
    echo "${sample}: Sorted" &&
    mv ${ofn}.sorted ${ofn} &&
    samtools index ${ofn} ) ||
    ( echo "${sample}: Error" && exit )

echo "${sample}: Finished"

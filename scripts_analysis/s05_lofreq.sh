ifd=$1
ofd=$2
gfd=$3

sample=$( basename $(ls ${ifd}/*.bam) .bam)

echo "${sample}: Started"

rfn="${gfd}KY962518_looped_2120.fasta"

bam=${ifd}${sample}.bam

rm -rf $ofd && mkdir -p $ofd
qbam="${ofd}${sample}.bam"

echo "${sample}: Adding qualities to alignments"

lofreq indelqual --dindel --ref $rfn -o ${qbam}.tmp ${bam}
lofreq alnqual -b ${qbam}.tmp $rfn > ${qbam}
rm ${qbam}.tmp

lofreq index $qbam

echo "${sample}: Calling rDNA variants"

lofreq call --call-indels \
    -f $rfn \
    -o ${ofd}${sample}.many.vcf \
    -r KY962518.1_looped_2120:250-15452 \
    -a 1 \
    -b 1 \
    -B \
    --no-default-filter \
    --use-orphan \
    --force-overwrite \
    $qbam

lofreq call --call-indels \
    -f $rfn \
    -o ${ofd}${sample}.unfilt.vcf \
    -r KY962518.1_looped_2120:250-15452 \
    --use-orphan \
    --force-overwrite \
    $qbam

rm -f ${ofd}${sample}.hc.vcf
lofreq filter \
    -i ${ofd}${sample}.unfilt.vcf \
    -o ${ofd}${sample}.hc.vcf \
    -v 2000 \
    -a 0.1

rm -r ${qbam}*

echo "${sample}: Finished"

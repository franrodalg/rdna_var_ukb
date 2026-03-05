ifd=$1
ofd=$2

R1=$(ls ${ifd}/*_1.fq.gz)
R2=$(ls ${ifd}/*_2.fq.gz)

sample=`basename ${R1} _1.fq.gz`

echo "${sample}: Started"

echo "${sample}: Fixing"

mkdir -p $ofd

t1="${ofd}${sample}_1.fq"
t2="${ofd}${sample}_2.fq"
ts="${ofd}${sample}_s.fq"

~/bbmap/repair.sh -Xmx20g in=$R1 in2=$R2 out=$t1 out2=$t2 outs=$ts

echo "${sample}: Trimming"

( trim_galore \
    --paired \
    --cores $(( $NSLOTS / 4 )) \
    --gzip \
    -o $ofd \
    $t1 $t2 ) ||
( echo "${sample}: Trimming error" && exit )

rm ${ofd}*.fq
rm $R1
rm $R2

echo "${sample}: Finished"

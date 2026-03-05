fn="/mnt/project/data/crams/crams.txt"

while read f; do

  eid=$( basename "$f" .cram | cut -d_ -f1 )
  echo "Sample: ${eid}"
  cram="/mnt/project/${f}"
  crai="${cram}.crai"
  bed="hg38_rdna_transcriptional_minimal.bed"

  samtools view -@ 14 -L "${bed}" -b -M -X "${cram}" "${crai}" > ${eid}.bam

done < ${fn}

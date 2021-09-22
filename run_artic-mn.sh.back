#!/bin/bash -l

set -e

export NXF_OFFLINE='TRUE'

run_name=$(basename $(pwd))

nextflow /home/mdh/shared/software_modules/artic-mn/artic-mn/artic-mn.nf \
  -c /home/mdh/shared/software_modules/artic-mn/artic-mn/artic-mn.config \
  --fastq_dir fastq_pass/ \
  -with-report ${run_name}/nf_report.html \
  -with-trace ${run_name}/nf_trace.tsv \
  -resume


# Look for failed samples
awk -F '\t' '$5 == "FAILED" && $4 !~ /\(NTC/ {print $5,$4,$6,$7,$2}' ${run_name}/nf_trace.tsv | sort > ${run_name}/${run_name}_errors.tsv

# Remove pangolin
rm /panfs/roc/groups/7/mdh/shared/software_modules/artic-mn/containers/${USER}/staphb-pangolin-latest.img

# Adjust permissions
chmod -R 775 ../${run_name}

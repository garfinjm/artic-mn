#!/usr/bin/env bash

set -e

primerset=$1
extra_config=""

# Make sure a valid primer set is given
if [ "${primerset}" = "V1" ] || [ "${primerset}" = "V2" ] || [ "${primerset}" = "V3" ] || [ "${primerset}" = "V4" ] || [ "${primerset}" = "V1200" ]
then
  echo "Running analysis with ${primerset} primers"
else
  echo -e "\nInvalid Primerset:${primerset}\n\nPick one:\nV1\nV2\nV3\nV4\nV1200"
  exit 1
fi

# Additional processing for V1200
if [ "${primerset}" = "V1200" ]
then
  extra_config="-c /home/mdh/shared/software_modules/artic-mn/artic-mn/midnight.config"
fi

export NXF_OFFLINE='TRUE'

run_name=$(basename $(pwd))

nextflow /home/mdh/shared/software_modules/artic-mn/artic-mn/artic-mn.nf \
  ${extra_config} \
  -c /home/mdh/shared/software_modules/artic-mn/artic-mn/artic-mn.config \
  --primers ${primerset} \
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
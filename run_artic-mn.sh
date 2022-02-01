#!/usr/bin/env bash

set -e

primerset=$1
extra_config=""

# Make sure a valid primer set is given
if [ "${primerset}" = "V4" ] || [ "${primerset}" = "V4.1" ] || [ "${primerset}" = "V1200" ]
then
  echo "Running analysis with ${primerset} primers"
else
  echo -e "\nInvalid Primerset:${primerset}\n\nPick one:V4,V4.1,V1200"
  exit 1
fi

# Untested primer sets
# ${primerset}" = "V1" ] || [ "${primerset}" = "V2" ] || [ "${primerset}" = "V3" ]


# Additional processing for V1200
if [ "${primerset}" = "V1200" ]
then
  extra_config="-c /home/mdh/shared/software_modules/artic-mn/artic-mn/midnight.config"
fi

export NXF_OFFLINE='TRUE'

run_name=$(basename $(pwd))

nextflow /home/mdh/shared/software_modules/artic-mn/artic-mn/artic-mn.nf \
  -c /home/mdh/shared/software_modules/artic-mn/artic-mn/artic-mn.config \
  ${extra_config} \
  --primers ${primerset} \
  --fastq_dir fastq_pass/ \
  -with-report ${run_name}/nf_report.html \
  -with-trace ${run_name}/nf_trace.tsv \
  -with-dag ${run_name}/dag.html \
  -resume

# Look for failed samples
awk -F '\t' '$5 == "FAILED" && $4 !~ /\(NTC/ {print $5,$4,$6,$7,$2}' ${run_name}/nf_trace.tsv | sort > ${run_name}/errors_${run_name}.tsv
repeats=$(cut -d' ' -f 3 ${run_name}/errors_${run_name}.tsv | uniq | tr -d '(),')
if [ "$repeats" ];
then
  echo -e "\n\033[0;31mCareful, it looks like the samples below had errors during analysis:";
  echo -e "$repeats";
  echo -e "\n\033[0m";
fi

echo "Careful, the read pileups for the samples below need review near position 19210:"
for file in ./${run_name}/vcfs/pass/*pass.vcf ./${run_name}/vcfs/review/*pass.vcf
do
  awk '(($2 == 19209 || $2 == 19212 || $2 == 19214)&&($5 == "G" || $5 == "GG"))||($2 == 19217 && $5 == "A") {print FILENAME; exit}' $file | cut -d"/" -f5 | cut -d"_" -f1
done

# Remove pangolin
rm ./work/singularity/staphb-pangolin-latest.img

# Adjust permissions
echo "Just adjusting permissions on this directory now, you can start transfering to artic_output";
chmod -R 775 ../${run_name};
echo "OK, all done now.";

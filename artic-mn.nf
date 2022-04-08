#!/usr/bin/env nextflow

//Description: Adaptation of staphb-wf monroe ont_assembly, which is a wrapper/extension of ARTIC Network nCoV-2019 Bioinformatics SOP
//Author of artic-mn.nf: Jake Garfin (jake.garfin@state.mn.us)
//Link: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Authors of staphb-wf monroe ont_assembly: Kelsey Florek and Abigail Shockey
//Link: https://github.com/StaPH-B/staphb_toolkit/blob/main/staphb_toolkit/workflows/monroe/monroe_ont_assembly.nf
//Email: kelsey.florek@slh.wisc.edu

// default params
params.fastq_dir = ""
params.outdir = workflow.launchDir.getName()
params.primers = "V4"
params.medaka_model = "r941_min_high_g360"
params.run_name = workflow.launchDir.getName()






Channel
    .fromPath( "${params.fastq_dir}/*.fastq.gz")
    .ifEmpty { exit 1, "Cannot find any fastq files in: ${params.fastq_dir}" }
    .set { demultiplexed_reads }


// run guppyplex to filer on read length, assume quality was checked during basecalling
process artic_guppyplex {
  errorStrategy 'retry'
  maxRetries 1

  input:
    path(reads) from demultiplexed_reads

  output:
    tuple env(name), env(subrun), env(guppyplex_count), path("guppyplex_*.tsv") into guppyplex_results
    tuple env(name), env(subrun), path("*.fastq") into guppyplexed_reads

  shell:
'''
name=$(basename -s '.fastq.gz' !{reads})
subrun=$(echo $name | cut -d"_" -f2-)

artic guppyplex \
--skip-quality-check \
--min-length !{params.min_length} \
--max-length !{params.max_length} \
--directory . \
--output ${name}.fastq \
| tee guppyplex_${name}.tsv

guppyplex_count=$(cut -f2 guppyplex_${name}.tsv)
'''
}

// run artic minion with medaka to create assemblies, vcfs, and bams
process artic_medaka_pipeline {
  tag "$name"

  errorStrategy { name.substring(0,3).equals('NTC') ? 'ignore' : task.attempt <= maxRetries ? 'retry' : 'ignore' }
  maxRetries 1

  input:
    tuple val(name), val(subrun), path(fastq) from guppyplexed_reads

  output:
    tuple env(name), path("*.fail.vcf"), path("*.pass.vcf"), path("*.primertrimmed.rg.sorted.bam") into artic_medaka_results
    tuple env(name), path("*.consensus.fasta") into coverage_queue

  shell:
'''
name=!{name}
subrun=!{subrun}
artic minion \
--min-depth 50 \
--medaka --medaka-model !{params.medaka_model} \
--normalise !{params.normalise} \
--threads !{task.cpus} \
--scheme-directory /primer-schemes \
--read-file !{fastq} \
SARS-CoV-2/!{params.primers} \
!{name}
'''
}

// Trim ends and check breadth of coverage on assemblies
process coverage_check {
  tag "$name"
  
  stageInMode 'copy'

  input:
    tuple val(name), path(fasta)from coverage_queue

  output:
    tuple env(name), path("${name}.fasta"), env(cov_perc), env(cov_qc), env(cov_conf) into coverage_results
    tuple env(name), path("${name}.fasta") into vadr_queue
    tuple env(name), path("${name}.fasta") into pangolin_queue
    env(NTCA_status) into NTCA_queue optional true
    env(NTCB_status) into NTCB_queue optional true


  shell:
'''
name=!{name}
sid=$(echo "!{name}" | cut -d"_" -f 1)
min_bases=26913

pretrim_char_count=$(seqkit stats -G N -a -T !{fasta} | cut -f 5 | tail -1)
pretrim_n_count=$(seqkit stats -G N -a -T !{fasta} | cut -f 12 | tail -1)

# Check trimming to make sure there is a sequence in the output fasta
if [[ 1 == $(echo "($pretrim_char_count-$pretrim_n_count)==0" | bc) ]]
then
  echo "This fasta is empty! Lets add some filler."
  seqkit -is replace -p "^n+|n+\$" -r "NNN" !{fasta} > !{name}.fasta
else
  echo "This fasta has actual sequence in it. Lets trim the Ns off the ends."
  seqkit -is replace -p "^n+|n+\$" -r "" !{fasta} > !{name}.fasta
fi

# Trim terminal Ns and calculate coverage
char_count=$(seqkit stats -G N -a -T !{name}.fasta | cut -f 5 | tail -1)
n_count=$(seqkit stats -G N -a -T !{name}.fasta | cut -f 12 | tail -1)
cov_perc=$(echo "scale=4; (($char_count-$n_count)/29903)*100" | bc -l)

if [[ 1 == $(echo "($char_count-$n_count)>=$min_bases" | bc) ]]
then
  cov_qc="Pass";
  cov_conf="Confirmed"
else
  cov_qc="Fail";
  cov_conf="Probable"
fi

# Handle NTCA status
if [ "${sid}" = "NTC-A" ]; then
  if [[ 1 == $(echo "($cov_perc) >= 10" | bc) ]]; then
    NTCA_status="Fail";
  else
    NTCA_status="Pass"
  fi
fi

# Handle NTCB status
if [ "${sid}" = "NTC-B" ]; then
  if [[ 1 == $(echo "($cov_perc) >= 10" | bc) ]]; then
    NTCB_status="Fail";
  else
    NTCB_status="Pass"
  fi
fi
'''
}

// Convert NTCA and NTCB queues into value channels, set as not found as needed
NTCA_value = NTCA_queue.first().ifEmpty("NTC-A not found")
NTCB_value = NTCB_queue.first().ifEmpty("NTC-B not found")

// run vadr, flag genomes needing review
process vadr {
  tag "$name"

  input:
    tuple val(name), path(fasta) from vadr_queue

  output:
    tuple env(name), env(vadr_status), path("*_vadr.tar.gz"), path("${name}/${name}.vadr.alt.list") into vadr_results

  shell:
'''
name=!{name}

/opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
--minlen 50 \
--maxlen 30000 \
!{fasta} 1> trimmed.fasta

if [ -s trimmed.fasta ]
then
  v-annotate.pl \
  --split \
  --cpu !{task.cpus} \
  --glsearch -s -r \
  --nomisc \
  --mkey sarscov2 \
  --lowsim5seq 6 \
  --lowsim3seq 6 \
  --alt_fail lowscore,insertnn,deletinn \
  --mdir /opt/vadr/vadr-models/ \
  trimmed.fasta !{name}

else
  mkdir !{name}
  echo "Vadr not run, zero coverage" > !{name}/no_coverage.txt
  echo "Vadr not run, zero coverage" > !{name}/!{name}.vadr.alt.list
fi


if [ -s !{name}/!{name}.vadr.pass.list ]
then
  vadr_status="Pass"
else
  vadr_status="Fail"
fi

tar -zcvf !{name}_vadr.tar.gz !{name}
'''
}

// run pangolin
process pangolin_typing {
  tag "$name"

  publishDir "${params.outdir}/pangolin", mode: 'copy'

  input:
  tuple val(name), path(fasta) from pangolin_queue

  output:
  tuple env(name), path ("pusher_*.csv") into pangolin_files

  tuple env(name), \
  env(pusher_lineage), \
  env(pusher_version), \
  env(pusher_status), \
  env(pusher_note) into pangolin_results

  shell:
'''
name=!{name}

pangolin -t 1 --analysis-mode accurate --outfile pusher_!{fasta}.csv !{fasta}

pusher_lineage=$(tail -n 1 pusher_!{fasta}.csv | cut -d, -f 2)
pusher_version=$(tail -n 1 pusher_!{fasta}.csv | cut -d, -f 9)
pusher_status=$(tail -n 1 pusher_!{fasta}.csv | cut -d, -f 14)
pusher_note=$(tail -n 1 pusher_!{fasta}.csv | cut -d, -f 15)


# Fix sometimes blank outputs
#if [ -z "${plearn_scorpio}" ]
#then
#  plearn_scorpio="null"
#fi
#if [ -z" ${plearn_note}" ]; then
#  plearn_note="null"
#fi
#if [ -z "${pusher_note}" ]; then
#  pusher_note="null"
#fi
'''
}


sample_summary_queue = guppyplex_results.join(artic_medaka_results).join(coverage_results).join(vadr_results).join(pangolin_results)


process sample_summary {
  tag "$name"
  
  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'pass/*.fasta*'
  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'fail/*.fasta*'
  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'review/*.fasta*'
  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'NTC_fail/*.fasta*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'pass/*.vcf*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'fail/*.vcf*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'review/*.vcf*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'NTC_fail/*.vcf*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'pass/*.bam*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'fail/*.bam*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'review/*.bam*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'NTC_fail/*.bam*'
  publishDir "${params.outdir}/guppyplex/", mode: 'copy', pattern: 'pass/guppyplex_*'
  publishDir "${params.outdir}/guppyplex/", mode: 'copy', pattern: 'fail/guppyplex_*'
  publishDir "${params.outdir}/guppyplex/", mode: 'copy', pattern: 'review/guppyplex_*'
  publishDir "${params.outdir}/guppyplex/", mode: 'copy', pattern: 'NTC_fail/guppyplex_*'
  publishDir "${params.outdir}/vadr/zips/", mode: 'copy', pattern: 'pass/*_vadr.tar.gz'
  publishDir "${params.outdir}/vadr/zips/", mode: 'copy', pattern: 'fail/*_vadr.tar.gz'
  publishDir "${params.outdir}/vadr/zips/", mode: 'copy', pattern: 'review/*_vadr.tar.gz'
  publishDir "${params.outdir}/vadr/zips/", mode: 'copy', pattern: 'NTC_fail/*_vadr.tar.gz'
  publishDir "${params.outdir}/vadr/", mode: 'copy', pattern: 'pass/*.vadr.alt.list'
  publishDir "${params.outdir}/vadr/", mode: 'copy', pattern: 'fail/*.vadr.alt.list'
  publishDir "${params.outdir}/vadr/", mode: 'copy', pattern: 'review/*.vadr.alt.list'
  publishDir "${params.outdir}/vadr/", mode: 'copy', pattern: 'NTC_fail/*.vadr.alt.list'

  input:
  val(NTCA_status) from NTCA_value
  val(NTCB_status) from NTCB_value
  tuple val(name), \
  val(subrun), \
  val(guppyplex_count), \
  path(guppyplex_file), \
  path(fail_vcf), \
  path(pass_vcf), \
  path(bam_file), \
  path(fasta_file), \
  val(cov_perc), \
  val(cov_qc), \
  val(cov_conf), \
  val(vadr_status), \
  path(vadr_zip), \
  path(vadr_list), \
  val(pusher_lineage), \
  val(pusher_version), \
  val(pusher_status), \
  val(pusher_note) from sample_summary_queue

  output:
  path ("*_summary.csv") into run_summary
  path ("*_harvest.csv") into harvest_summary
  tuple val(name), path("pass/*.fasta"), path("pass/*.fail.vcf"), path("pass/*.pass.vcf"), path("pass/*.primertrimmed.rg.sorted.bam"), path("pass/*_vadr.tar.gz"), path("pass/*vadr.alt.list"), path("pass/guppyplex_*") optional true
  tuple val(name), path("fail/*.fasta"), path("fail/*.fail.vcf"), path("fail/*.pass.vcf"), path("fail/*.primertrimmed.rg.sorted.bam"), path("fail/*_vadr.tar.gz"), path("fail/*vadr.alt.list"), path("fail/guppyplex_*") optional true
  tuple val(name), path("review/*.fasta"), path("review/*.fail.vcf"), path("review/*.pass.vcf"), path("review/*.primertrimmed.rg.sorted.bam"), path("review/*_vadr.tar.gz"), path("review/*vadr.alt.list"), path("review/guppyplex_*") optional true
  tuple val(name), path("NTC_fail/*.fasta"), path("NTC_fail/*.fail.vcf"), path("NTC_fail/*.pass.vcf"), path("NTC_fail/*.primertrimmed.rg.sorted.bam"), path("NTC_fail/*_vadr.tar.gz"), path("NTC_fail/*vadr.alt.list"), path("NTC_fail/guppyplex_*") optional true

  shell:
'''
# Create environment variables for each input
name="!{name}"
sid=$(echo "!{name}" | cut -d"_" -f 1)
subrun="!{subrun}"
guppyplex_count="!{guppyplex_count}"
guppyplex_file="!{guppyplex_file}"
vcf_fail="!{fail_vcf}"
vcf_pass="!{pass_vcf}"
bam_file="!{bam_file}"
fasta_file="!{fasta_file}"
cov_perc="!{cov_perc}"
cov_qc="!{cov_qc}"
cov_conf="!{cov_conf}"
vadr_status="!{vadr_status}"
vadr_zip="!{vadr_zip}"
vadr_list="!{vadr_list}"
pusher_lineage="!{pusher_lineage}"
pusher_version="!{pusher_version}"
pusher_status="!{pusher_status}"
pusher_note="!{pusher_note}"

control_set=${subrun: -1}

echo "Control Set: ${control_set}"
echo "NTC-A: !{NTCA_status}"
echo "NTC-B: !{NTCB_status}"

echo "My files:"
echo "name: !{name}"
echo "sid: ${sid}"
echo "subrun: !{subrun}"
echo "guppyplexed_reads: !{guppyplex_count}"
echo "guppyplex_file: !{guppyplex_file}"
echo "vcf_fail: !{fail_vcf}"
echo "vcf_pass: !{pass_vcf}"
echo "bam_file: !{bam_file}"
echo "fasta_file: !{fasta_file}"
echo "cov_perc: !{cov_perc}"
echo "cov_qc: !{cov_qc}"
echo "cov_conf: !{cov_conf}"
echo "vadr_status: !{vadr_status}"
echo "vadr_zip: !{vadr_zip}"
echo "vadr_list: !{vadr_list}"
echo "pusher_lineage: !{pusher_lineage}"
echo "pusher_version: !{pusher_version}"
echo "pusher_status: !{pusher_status}"
echo "pusher_note: !{pusher_note}"

# Make result sorting dirs
mkdir pass
mkdir fail
mkdir NTC_fail
mkdir review


if ([ ${control_set} = "A" ] && [ !{NTCA_status} = "Fail" ]) || ([ ${control_set} = "B" ] && [ !{NTCB_status} = "Fail" ]); then
  echo "NTC control failure, adjusting results summary files"
  cov_qc="Coverage NTC Failure"
  cov_conf="Coverage NTC Failure"
  pusher_lineage="Coverage NTC Failure"
  pusher_version="Coverage NTC Failure"
  pusher_status="Coverage NTC Failure"
  pusher_note="Coverage NTC Failure"
  mv !{guppyplex_file} !{fail_vcf} !{pass_vcf} !{bam_file} !{fasta_file} !{vadr_zip} !{vadr_list} ./NTC_fail
elif [ "!{cov_qc}" = "Pass" ]; then
  if [ "!{vadr_status}" = "Pass" ]; then
    mv !{guppyplex_file} !{fail_vcf} !{pass_vcf} !{bam_file} !{fasta_file} !{vadr_zip} !{vadr_list} ./pass
  else
    mv !{guppyplex_file} !{fail_vcf} !{pass_vcf} !{bam_file} !{fasta_file} !{vadr_zip} !{vadr_list} ./review
  fi
else
  mv !{guppyplex_file} !{fail_vcf} !{pass_vcf} !{bam_file} !{fasta_file} !{vadr_zip} !{vadr_list} ./fail
fi

# Create summary strings
echo "${sid},${subrun},${guppyplex_count},${cov_perc},${cov_qc},${cov_conf},${vadr_status},${pusher_lineage},${pusher_version},${pusher_status},${pusher_note}" > !{name}_summary.csv

if  [ "${cov_qc}" = "Pass" ]; then
  echo "${sid},${pusher_lineage},${pusher_version},${pusher_status},,," > !{name}_harvest.csv
elif [ "${cov_qc}" = "Coverage NTC Failure" ]; then
  echo "${sid},${pusher_lineage},${pusher_version},${pusher_status},none,none,none" > !{name}_harvest.csv
else
  echo "${sid},${pusher_lineage},${pusher_version},${pusher_status},none,none,none" > !{name}_harvest.csv
fi

'''
}


// summarize run
process run_summary {
  tag "run_summary"

  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path (summary_files) from run_summary.collect()

  output:
  path("${params.run_name}_summary.csv")

  shell:
'''
echo "id,run,guppyplexed_reads,cov_perc,cov_status,cov_confidence,vadr_status,pusher_lineage,pusher_version,pusher_status,pusher_note" > temp.txt

cat *_summary.csv | sort >> temp.txt

mv temp.txt !{params.run_name}_summary.csv
'''
}

// summarize run for harvest
process harvest_summary {
  tag "harvest_summary"

  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path (harvest_files) from harvest_summary.collect()

  output:
  path("${params.run_name}_harvest.csv")

  shell:
'''
echo "key,lineage,pangoLEARN,quality,Virus_ID,GENBANK,GISAID" > temp.txt

cat *_harvest.csv | sort >> temp.txt

mv temp.txt !{params.run_name}_harvest.csv
'''
}

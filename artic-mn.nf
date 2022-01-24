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
    tuple env(name), path("*.fastq") into polish_files
    path("guppyplex_*.tsv") into guppyplex_summary_queue

  shell:
    '''
    name=$(basename -s '.fastq.gz' !{reads})

    artic guppyplex \
    --skip-quality-check \
    --min-length !{params.min_length} \
    --max-length !{params.max_length} \
    --directory . \
    --output ${name}.fastq \
    | tee guppyplex_${name}.tsv
    '''
}

// run artic minion with medaka to create assemblies, vcfs, and bams
process artic_medaka_pipeline {
  tag "$name"

  errorStrategy { name.substring(0,3).equals('NTC') ? 'ignore' : task.attempt <= maxRetries ? 'retry' : 'ignore' }
  maxRetries 1

  input:
    tuple val(name), path(fastq) from polish_files

  output:
    path "*{.primertrimmed.rg.sorted.bam,.vcf.gz,.fail.vcf}"
    path("*.consensus.fasta") into typing_queue
    tuple val(name), path("*.consensus.fasta"), path("*.fail.vcf"), path("*.pass.vcf"), path("*.primertrimmed.rg.sorted.bam") into coverage_queue

  shell:
    '''
    artic minion \
    --medaka --medaka-model !{params.medaka_model} \
    --normalise !{params.normalise} \
    --threads !{task.cpus} \
    --scheme-directory /primer-schemes \
    --read-file !{fastq} \
    SARS-CoV-2/!{params.primers} \
    !{name}
    '''
}

// Check breadth of coverage on assemblies, only send >~90% complete to vadr
process coverage_check {
  tag "$name"

  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'fail/*.consensus.fasta*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'fail/*.vcf*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'fail/*.bam*'

  input:
    tuple val(name), path(fasta), path(fail_vcf), path(pass_vcf), path(bam) from coverage_queue

  output:
    path "single_covQC_*" into coverage_summary_queue
    tuple val(name), path("pass/*.consensus.fasta"), path("pass/*.fail.vcf"), path("pass/*.pass.vcf"), path("pass/*.primertrimmed.rg.sorted.bam") optional true into vadr_queue
    tuple val(name), path("fail/*.consensus.fasta"), path("fail/*.fail.vcf"), path("fail/*.pass.vcf"), path("fail/*.primertrimmed.rg.sorted.bam") optional true

  shell:
    '''
    mkdir pass
    mkdir fail
    min_bases=26913
    
    seqkit -is replace -p "^n+|n+\$" -r "" !{fasta} > temp_assembly.fa
    mv temp_assembly.fa !{fasta}
    char_count=$(seqkit stats -G N -a -T !{fasta} | cut -f 5 | tail -1)
    n_count=$(seqkit stats -G N -a -T !{fasta} | cut -f 12 | tail -1)
    percent=$(echo "scale=4; (($char_count-$n_count)/29903)*100" | bc -l)
    
    if [[ 1 == $(echo "($char_count-$n_count)>=$min_bases" | bc) ]]
    then
      echo -e "!{name},Pass,${percent},Confirmed" 2>&1 > single_covQC_!{name}.csv
      mv !{fasta} !{fail_vcf} !{pass_vcf} !{bam} ./pass
    else
      echo -e "!{name},Fail,${percent},Probable" 2>&1 > single_covQC_!{name}.csv
      mv !{fasta} !{fail_vcf} !{pass_vcf} !{bam} ./fail
    fi
    '''
}

// Summarize guppyplex post-filtering read counts
process guppyplex_collect {
  tag "$params.run_name"

  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(guppyplex_log) from guppyplex_summary_queue.collect()

  output:
    path "guppyplex_*.csv"

  shell:
    '''
    cat guppyplex_* > guppyplex_!{params.run_name}.csv
    sed -i 's/\t/,/g' guppyplex_!{params.run_name}.csv
    '''
}

// summarize breadth of coverage
process coverage_collect {
  tag "$params.run_name"

  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(cov_log) from coverage_summary_queue.collect()

  output:
    path "covQC_${params.run_name}.csv"

  shell:
    '''
       echo -e "accession,status,perc_genome,prob/conf" > covQC_!{params.run_name}.csv
       cat single_covQC_* | sed  "s/_!{params.run_name}//g" >> covQC_!{params.run_name}.csv
    '''
}

// run vadr and sort genomes that need review
process vadr {
  tag "$name"

  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'pass/*.consensus.fasta*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'pass/*.vcf*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'pass/*.bam*'
  publishDir "${params.outdir}/vadr/", mode: 'copy', pattern: 'pass/*_vadr', type: 'dir'
  publishDir "${params.outdir}/assemblies/", mode: 'copy', pattern: 'review/*.consensus.fasta*'
  publishDir "${params.outdir}/vcfs/", mode: 'copy', pattern: 'review/*.vcf*'
  publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: 'review/*.bam*'
  publishDir "${params.outdir}/vadr/", mode: 'copy', pattern: 'review/*_vadr', type: 'dir'

  input:
    tuple val(name), path(fasta), path(fail_vcf), path(pass_vcf), path(bam) from coverage_queue from vadr_queue

  output:
    tuple val(name), path("pass/*.consensus.fasta"), path("pass/*.fail.vcf"), path("pass/*.pass.vcf"), path("pass/*.primertrimmed.rg.sorted.bam"), path("pass/*_vadr") optional true
    tuple val(name), path("review/*.consensus.fasta"), path("review/*.fail.vcf"), path("review/*.pass.vcf"), path("review/*.primertrimmed.rg.sorted.bam"), path("review/*_vadr") optional true

  shell:
    '''
    mkdir pass
    mkdir review
    
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
    --minlen 50 \
    --maxlen 30000 \
    !{fasta} 1> trimmed.fasta

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
    trimmed.fasta !{name}_vadr

    #/opt/vadr/vadr/miniscripts/vadr-map-model-coords.pl\
    #!{name}_vadr/!{name}_vadr".vadr.alt.list" \
    #/opt/vadr/vadr-models/sarscov2.mmap  NC_045512 \
    #1> !{name}_vadr/!{name}_vadr".vadr.alt.coords"

    if [ -s !{name}_vadr/!{name}_vadr.vadr.pass.list ]
    then
      mv !{fasta} !{fail_vcf} !{pass_vcf} !{bam} !{name}_vadr ./pass
    else
      mv !{fasta} !{fail_vcf} !{pass_vcf} !{bam} !{name}_vadr ./review
    fi
    '''
}

// run pangolin, and create a data-drop file for Harvest LIMS
process pangolin_typing {
  tag "$params.run_name"

  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(assembly) from typing_queue.collect()

  output:
  path "lineage_report_*.csv"
  path "usher_lineage_report_*.csv"
  path "harvest_*.csv"

  shell:
    '''
    cat *.fasta > all_assemblies.fasta
    pangolin -t 8 all_assemblies.fasta --outfile lineage_report_!{params.run_name}.csv

    # Again, this time using usher
    pangolin -t 8 --usher all_assemblies.fasta --outfile usher_lineage_report_!{params.run_name}.csv

    # Make a Harvest LIMS compatiable lineage file
    echo "taxon,Key,lineage,pangoLEARN,quality" > harvest_!{params.run_name}.csv
    awk 'BEGIN {FS=",";OFS=","} NR>1 {print $1,substr($1,1,9),$2,substr($10,6,2)"/"substr($10,9,2)"/"substr($10,1,4),$12}' lineage_report_!{params.run_name}.csv >> harvest_!{params.run_name}.csv
    perl -pi -e "s/\\n/\\r\\n/" harvest_!{params.run_name}.csv

    rm all_assemblies.fasta
  '''
}
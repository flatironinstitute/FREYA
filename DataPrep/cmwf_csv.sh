#!/bin/bash

function log {
    echo $(date +%F_%T) $$ $1
}

# Arguments:
#   - CSV file describing phenotypes. It should have three columns, e.g.:
#
#       SampleID,Histology,Patient
#       1A,N,1
#       1C,M,1
#       1E,B,1
#       1F,M,1
#
#     The first line is assumed to be a header line and is
#     discarded. This file is used to generate a list of fastq
#     files. For example, '1C' indicates we should process a fastq file
#     with the name:
#
#       ${FastqDir}/1C.fq.gz
#
#     The third column provides the id of the dog. Consolidated
#     results for a dog use this id as the basename for files.
#
#   - Directory containing fq.gz files (FastqDir above)
#   - Directory for storing output
#   - One or more stages, see list below (optional, if omitted all stages are run)
phenoCsv=$1
fastqDir=$(readlink -f $2)
outDir=$(readlink -f $3)

log "Runnning on $(hostname)."

Stages="hisat2 fastqc dexcount aorrg markdup splitncr hapcall varfilt snpeff mpicker"
declare -A aaStages
doMe=1
for s in ${Stages}; do aaStages[$s]=${doMe}; done

shift 3
for s in $@
do
    if [[ ${aaStages[$s]+_} ]]
    then
	doMe=2
	aaStages[$s]=${doMe}
    else
	echo "Unknown processing stage: ""$s"", should be one of ${Stages}."
	exit 1
    fi
done


# The singularity image has many of the dependencies installed in
# standard locations, so we don't need to do much to PATH. See the
# Dockerfile used to create the docker image/container converted into
# the singularity image.  A few you other bits and pieces are located
# in this directory:
CMWF_ROOT=/usr/cmwf

export PATH=${CMWF_ROOT}/disBatch:$PATH

GATKJar=${CMWF_ROOT}/jars/GenomeAnalysisTK-3.8-0.jar
PicardJar=${CMWF_ROOT}/jars/picard.jar
SnpEffJar=${CMWF_ROOT}/jars/snpEff/snpEff.jar

# References. Paths mentioned here have to be accessible from within
# singularity (using -B /some/path/prefix).
#
# Annotation gff (?)
DC_GFF=/mnt/ceph/users/carriero/KileyGraim/LanceParsons/canine_breast_cancer/data/genomes/canFam3.1.91/snpEff/CanFam3.1.91/genes.gff
# Genome (have to have .fai file for this too).
CFFA=/mnt/ceph/users/carriero/KileyGraim/LanceParsons/canine_breast_cancer/data/genomes/canFam3.1.91/dna/Canis_familiaris.CanFam3.1.dna.toplevel.fa
# Genome prepped for use with hisat2
HSX=/mnt/ceph/users/carriero/KileyGraim/LanceParsons/canine_breast_cancer/data/genomes/canFam3.1.91/hisat2/Canis_familiaris.CanFam3.1.dna.toplevel

# Tunable parameters
#
# AddOrReplaceReadGroups: RGPL and RGPU options depend on the
# sequencer use.
RGPLParam="HISEQ2000"
RGPUParam="Princeton"
#
# SnpEff: Genome should reflect reference build and snpEff.config
# should have an entry for it. As above, config path must be
# accessible.
SnpEffConfigParam="/mnt/ceph/users/carriero/KileyGraim/LanceParsons/canine_breast_cancer/data/genomes/canFam3.1.91/snpEff/snpEff.config"
SnpEffGenomeParam="CanFam3.1.91"

# ** IMPORTANT **  ** IMPORTANT **  ** IMPORTANT **  ** IMPORTANT ** 
# We use the following to distribute subtasks over nodes:
#
# https://github.com/flatironinstitute/disBatch
#
# Each stage of the processing generates a task file (written to the
# "tasks" subdirectory of the output directory). Each of these files
# contains one or more lines of commands to be invoked. disBatch.py is
# given a list of nodes to use and a task file. It then assigns a
# certain number of tasks to each node, waits for one task to complete
# then assigns another, repeating this until all the lines of the task
# file have been executed. If run in a SLURM allocation, disBatch.py
# will automatically detect the list of nodes, otherwise you can use
# "-s host:cores" (multiple times if you wish) to describe the nodes
# to be used. In the simplest case, "-s localhost:$(nproc)" would
# direct disBatch.py to use all of your current node. You can specify
# additional nodes as long as you can ssh to them without needing a
# password. NB: as used in this script, disBatch.py doesn't make use
# of the per-node core count.
#
# To give a sense of scale: To process the dataset reported in a test
# run using three 28-core nodes with a fair amount of memory (512 GB)
# and with disBatch.py directed to run up to three tasks per node (up
# to nine tasks in flight at any given time) took about 40 hours.
# 
# If your nodes have less memory or fewer processors, you may want to
# reduce the number of concurrent tasks. Adjust the following
# accordingly.
TasksPerNode=6

# Set DB_TASK_PREFIX to something like "singularity run -B /file/system/path "
# if tasks need to be run via a singularity container.

function dbRetry () {
    # On some platforms, we've encountered sporadic SIGSEGV errors
    # with java. This wrapper function is used to invoke disBatch.py,
    # check for errors (reported in the status file produced by
    # disBatch.py) and rerun any tasks that failed the first time. If
    # the new run reports any failed tasks, we stop execution.

    bn=$1
    tpn=""
    [[ $# == 2 ]] && tpn="-t $2"
    
    log "Launching ${bn} tasks."
    
    p0="${bn}_db_0"
    disBatch.py -K -p ${p0} ${tpn} "${outDir}/tasks/${bn}Tasks"
    sn0="${p0}_status.txt"
    if $(egrep -q '^R' ${sn0})
    then
	log "Detected errors in ${sn0}, rerunning."
	p1="${bn}_db_1"
	disBatch.py -K -p ${p1} ${tpn} -r ${sn0} -R "${outDir}/tasks/${bn}Tasks"
	sn1="${p1}_status.txt"
	if $(egrep -q '^R' ${sn1})
	then
	    echo "Still errors in ${sn1}, giving up."
	    echo "Review and then consider running something like:"
	    echo "  disBatch.py -K -p ${p1/_1/_2} ${tpn} -r -R ${sn1} \"${outDir}/tasks/${bn}Tasks\""
	    exit 1
	fi
    fi
}

# Not all of these are currently used.
HapCallCores=14
Hisat2Cores=14
SplitNCRCores=14
VarFiltCores=14

mkdir -p ${outDir} && ( cd ${outDir} ; mkdir -p bams dexseq_count fastqc hisat2 logs tasks vcfs )

# We append to these files, so make sure they start clean (i.e., with just DB_TASK_PREFIX if specified).
pushd ${outDir}/tasks
for tf in hisat2Tasks fastqcTasks dexTasks aorrgTasks mdTasks sncrTasks hcTasks vfTasks seTasks
do
    echo "#DISBATCH PREFIX ${DB_TASK_PREFIX}" > ${tf}
done
popd

shopt -s nullglob

declare -A fqs missingfqs
if [[ ${aaStages["hisat2"]} == ${doMe} || ${aaStages["fastqc"]} == ${doMe} ]]
then
    # Sanity check the CSV file.
    for sample in $(awk -F, 'NR > 1{print $1}' ${phenoCsv})
    do
	fn="${fastqDir}/${sample}.fq.gz"
	if [[ -e ${fn} ]]
	then
	    fqs[${sample}]=${fn}
	else
	    missingfqs[${sample}]=${fn}
	fi
    done
fi
if [[ ${#missingfqs[@]} -gt 0 ]]
then
    echo "Missing fastqs:"
    for fq in "${!missingfqs[@]}"
    do
	echo "${fq}	${missingfqs[${fq}]}"
    done
    exit 1
fi

# Run hisat2 and fastqc
for bn in ${!fqs[@]}
do
    fn=${fqs[${bn}]}
    bam=${outDir}/hisat2/${bn}.bam

    echo "bash -c \"( ${CMWF_ROOT}/hisat2-2.0.4/hisat2  --threads ${Hisat2Cores} -x ${HSX} -U ${fn} | samtools view -Sbh -o ${bam} - )\" &> ${outDir}/logs/${bn}_hisat2.log" >> ${outDir}/tasks/hisat2Tasks
    echo "${CMWF_ROOT}/FastQC/fastqc  --quiet --outdir ${outDir}/fastqc ${fn}" >> ${outDir}/tasks/fastqcTasks
done
[[ ${aaStages["hisat2"]} == ${doMe} ]] && { log "hisat2 tasks launched." ; dbRetry hisat2 ${TasksPerNode} ; }
[[ ${aaStages["fastqc"]} == ${doMe} ]] && { log "fastqc tasks launched." ; dbRetry fastqc ${TasksPerNode} ; }

# Run dexseq_count and vc prep.
for bam in ${outDir}/hisat2/*.bam
do
    bn=$(basename ${bam} .bam)

    echo "python ${CMWF_ROOT}/DEXSeq/inst/python_scripts/dexseq_count.py --format bam ${DC_GFF} ${bam} ${outDir}/dexseq_count/${bn}.txt &> ${outDir}/logs/${bn}_dexseq_count.log" >> ${outDir}/tasks/dexTasks
    echo "java -jar ${PicardJar} AddOrReplaceReadGroups USE_JDK_DEFLATER=true USE_JDK_INFLATER=true I=${bam} O=${outDir}/bams/${bn}_sorted.bam SO=coordinate RGID=${bn} RGLB=${bn} RGPL=${RGPLParam} RGPU=${RGPUParam} RGSM=${bn} &> ${outDir}/logs/${bn}_aorrg.log" >> ${outDir}/tasks/aorrgTasks
    echo "java -jar ${PicardJar} MarkDuplicates USE_JDK_DEFLATER=true USE_JDK_INFLATER=true I=${outDir}/bams/${bn}_sorted.bam O=${outDir}/bams/${bn}_dedupped.bam QUIET=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${outDir}/bams/${bn}_output.metrics &> ${outDir}/logs/${bn}_md.log" >> ${outDir}/tasks/mdTasks
    echo "java -jar ${GATKJar} -jdk_deflater -jdk_inflater -T SplitNCigarReads -R ${CFFA} -I ${outDir}/bams/${bn}_dedupped.bam -o ${outDir}/bams/${bn}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS &> ${outDir}/logs/${bn}_sncr.log" >> ${outDir}/tasks/sncrTasks
done 
[[ ${aaStages["dexcount"]} == ${doMe} ]]  && { log "dexcount tasks launched." ; dbRetry dex ${TasksPerNode} ; }
[[ ${aaStages["aorrg"]} == ${doMe} ]]    && { log "aorrg tasks launched." ; dbRetry aorrg ${TasksPerNode} ; }
[[ ${aaStages["markdup"]} == ${doMe} ]]  && { log "markdup tasks launched." ; dbRetry md ${TasksPerNode} ; }
[[ ${aaStages["splitncr"]} == ${doMe} ]] && { log "splitncr tasks launched." ; dbRetry sncr ${TasksPerNode} ; }

if [[ ${aaStages["hapcall"]} == ${doMe} ]]
then
    # Generate locations used to parallelize HaplotypeCaller by "location" (chromosomes/contigs).
    rm -f  *_loc.list
    awk '$2 < 1000000{print $1 >> "small_loc.list" ; next}{print $1 >> $1"_loc.list"}' ${CFFA}.fai

    # Generate list of input bams for vc.
    ls ${outDir}/bams/*_split.bam > split_bams.list

    # Run HaplotypeCaller
    for l in *_loc.list
    do
	bn=$(basename $l .list)
	echo "java -jar ${GATKJar} -jdk_deflater -jdk_inflater -T HaplotypeCaller -R ${CFFA} -I split_bams.list -L $l -dontUseSoftClippedBases  -stand_call_conf 20.0 -o ${outDir}/vcfs/${bn}.vcf &> ${outDir}/logs/${bn}_hc.log"
    done >> ${outDir}/tasks/hcTasks
    log "hapcall tasks launched."
    dbRetry hc ${TasksPerNode}
fi

if [[ ${aaStages["varfilt"]} == ${doMe} ]]
then
    # Run VariantFiltration
    for vcf in ${outDir}/vcfs/*.vcf
    do
	bn=$(basename ${vcf} .vcf)
	echo "java -jar ${GATKJar} -jdk_deflater -jdk_inflater -T VariantFiltration -R ${CFFA} -V ${vcf} -o ${vcf/.vcf/_vf.vcf}  -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' &> ${outDir}/logs/${bn}_vf.log"
    done >> ${outDir}/tasks/vfTasks
    log "varfilt tasks launched."
    dbRetry vf ${TasksPerNode}
fi

if [[ ${aaStages["snpeff"]} == ${doMe} ]]
then
    log "snpeff started."
    echo "bash -c \"cat ${outDir}/vcfs/*_vf.vcf | java -jar ${SnpEffJar} -c ${SnpEffConfigParam} -v ${SnpEffGenomeParam} > ${outDir}/vcfs/CanineAll.ann.vcf 2> ${outDir}/logs/snpeff.log\"" >> ${outDir}/tasks/seTasks
    # Only one task. We execute it this way to offload the work to an allocated compute node.
    dbRetry se ${TasksPerNode}
fi

if [[ ${aaStages["mpicker"]} == ${doMe} ]]
then
    log "mutation picker started."
    # This shouldn't be too resource intensive, so we may run it on the submission node.
    ${DB_TASK_PREFIX} python ${CMWF_ROOT}/VCF_mutation_picker.0.5.py ${phenoCsv} ${outDir}/vcfs/CanineAll.ann.vcf ${outDir}/vcfs/mutations_genesOnly.csv > ${outDir}/vcfs/CanineAll_mp.txt 2> ${outDir}/logs/mpicker.log
fi

# Remind the user to clean up after reviewing the results. Provide
# text of commands that can be cut-and-pasted to do so.
cat <<EOF
Check output. If it looks good, you can run these clean up commands:

/bin/rm -f ${outDir}/bams/*_sorted.bam ${outDir}/bams/*_dedupped.bam
/bin/rm -f ${outDir}/hisat2/*.bam
/bin/rm -f ${outDir}/vcfs/*_loc.vcf ${outDir}/vcfs/*_loc.vcf.idx

EOF
log "Done."

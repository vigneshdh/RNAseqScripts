#!/usr/bin/env nextflow

/*
* USAGE: nextflow run SalminRNAseq_workflow.nf -qs 8
* Note that "-qs" is similar to "#PBS -t" and will only run a specified # of jobs at a time.
* This script creates hard links to data that exists in nextflow's work directory for everything.
* @Author - Vignesh Dhandapani(UoB)
*/


/*
 * Set parameter values here. Only change the values within quotations
 * Make sure the paths are correct for your dataset.
*/

/*Path to project folder*/
projectPath = "/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi"


/*Path to transcriptome files IF SALMON INDEX HAS NOT BEEN RUN (note is hashed out) */
/*reference_transcriptome = "INSERT REFERENCE TRANSCRIPTOME.FASTA" */



/*Path to indexed transcriptome files IF SALMON INDEX HAS ALREADY BEEN RUN*/
index_reference_dir_all = "/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/Dmagna_total_transcripts_index"
index_reference_dir_longest = "/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/Dmagna_total_transcripts_index"


/*Path to raw fastq files.*/
rawDataPath = "/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/mergedFastq"
params.reads = "$rawDataPath/*_{1,2}.fq" /*Provide pattern that will match all FASTQ files*/
Channel
    .fromFilePairs(params.reads, flat: true)
    .ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
    .set {read_pairs}

/*CompGen options. List memory in gigabytes like in suggestions below*/
myQueue = '--qos colboujk'
salmonMemory = '10'
salmonCPU = '4'
trimMemory = '10'
trimCPU = '2'

/*Module versions*/
salmonMod = 'salmon/0.8.2'
trimVersion = '0.33' /*Put the version here only*/
fastqcMod = 'fastqc/0.11.5'

/*
* Trimming options. Change trimming options here and note that $trimVersion is used to make sure
* the version is called consistently.
*/
trimOptions = "ILLUMINACLIP:/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/adapter_seqs.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:20"

/*Output paths*/
trimPath = "$projectPath/trimmomatic"
fastqcPath = "$projectPath/fastqc_afterTrim"
salmonquantallPath = "$projectPath/Salmon_quant_all"
salmonquantlongestPath = "$projectPath/Salmon_quant_longest"

/*
* Step 1. Trimming
* WARNING: considers '1' a valid exit status to get around wrapper error
*/

process trimmomatic {
    executor 'slurm'
    cpus trimCPU
    clusterOptions myQueue
    memory "$trimMemory GB"
    time '120h'
   /* module "trimmomatic/$trimVersion"
    */
    publishDir trimPath, mode: 'link'
    validExitStatus 0,1

    input:
    set val(id), file(read1), file(read2) from read_pairs

    output:
    set val(id), "${read1.baseName}.qualtrim.paired.fastq", "${read2.baseName}.qualtrim.paired.fastq" into fastqChannel
    set val(id), "${read1.baseName}.qualtrim.paired.fastq", "${read2.baseName}.qualtrim.paired.fastq" into salmonChannel_1
    set val(id), "${read1.baseName}.qualtrim.paired.fastq", "${read2.baseName}.qualtrim.paired.fastq" into salmonChannel_2
    file "*.trimmomatic.log"
    file "*.qualtrim.unpaired.fastq"

    """
    java -jar ${TMATIC_ROOT}/trimmomatic-0.32.jar PE \
    -threads $trimCPU -phred33 -trimlog ${id}.trimmomatic.log ${read1} ${read2} \
    ${read1.baseName}.qualtrim.paired.fastq ${read1.baseName}.qualtrim.unpaired.fastq \
    ${read2.baseName}.qualtrim.paired.fastq ${read2.baseName}.qualtrim.unpaired.fastq $trimOptions
    """
}



/* Step 2. FASTQC of trimmed reads
*/
process runFASTQC {
    executor 'slurm'
    cpus 1
    clusterOptions myQueue
    memory '2 GB'
    time '120h'
    /* module fastqcMod
    */
    publishDir fastqcPath, mode: 'copy'

    input:
    set pair_id, file(read1), file(read2) from fastqChannel

    output:
    file "*.html"
    file "*.zip"

    """
    fastqc -o ./ --noextract ${read1}
    fastqc -o ./ --noextract ${read2}
    """
}

/* Step 3.b.2 Salmon alignment and counting of trimmed reads on d. magna all
*/
process runSalmonquantall {
    executor 'slurm'
    cpus salmonCPU
    clusterOptions myQueue
    memory '32 GB'
    time '120h'
  /*  module salmonMod
  */
    publishDir salmonquantallPath, mode: 'copy'
    tag "$pair_id"

    input:
    set pair_id, file(read1), file(read2) from salmonChannel_2

    output:
    file "$pair_id"

    """
    /gpfs/bb/dhandapv/local/bin/salmon quant -p $salmonCPU -i $index_reference_dir_all -l A --incompatPrior 0.0 -1 $read1 -2 $read2 -o $pair_id --seqBias --gcBias
    """
}

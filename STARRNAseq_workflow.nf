#!/usr/bin/env nextflow

/*
* USAGE: nextflow run STARRNAseq_workflow.nf -qs 8
* Note that "-qs" is similar to "#PBS -t" and will only run a specified # of jobs at a time.
* This script creates hard links to data that exists in nextflow's work directory for everything.
* @Author - Vignesh Dhandapani (UoB)
*/


/*
 * Set parameter values here. Only change the values within quotations
 * Make sure the paths are correct for your dataset.
*/

/*Path to project folder*/
projectPath = "/rds/projects/2018/colboujk-testing-the-waters"



params.genome = "$projectPath/data/genome/STAR_2.5.1b-noAnno-dmagset7finloc9c"
genome_file = file(params.genome)
params.gtf = "$projectPath/data/genome/rev_dmagset7finloc9c.puban.gtf"
gtf_file = file(params.gtf)



/*Path to transcriptome files IF SALMON INDEX HAS NOT BEEN RUN (note is hashed out) */
/*reference_transcriptome = "INSERT REFERENCE TRANSCRIPTOME.FASTA" */



/*Path to indexed transcriptome files IF SALMON INDEX HAS ALREADY BEEN RUN*/
 * index_reference_dir_all = "/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/Dmagna_total_transcripts_index"
 * index_reference_dir_longest = "/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/Dmagna_total_transcripts_index"
*/


/*Path to raw fastq files.*/
rawDataPath = "/rds/projects/2018/colboujk-testing-the-waters/mergedFastq"
params.reads = "$rawDataPath/*_{1,2}.fq" /*Provide pattern that will match all FASTQ files*/
Channel
    .fromFilePairs(params.reads, flat: true)
    .ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
    .set {read_pairs}

/*CompGen options. List memory in gigabytes like in suggestions below*/
myQueue = '--qos colboujk'
starMemory = '10'
starCPU = '4'
countMemory = '5'
countCPU = '4'
trimMemory = '10'
trimCPU = '2'

/*Module versions*/
starMod = 'apps/star-aligner/v2.5.2b'
trimVersion = '0.33' /*Put the version here only*/
fastqcMod = 'fastqc/0.11.5'

/*
* Trimming options. Change trimming options here and note that $trimVersion is used to make sure
* the version is called consistently.
*/
trimOptions = "ILLUMINACLIP:/rds/projects/2018/colboujk-testing-the-waters/xiaojingLi/adapter_seqs.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:20"


/* Alignment & Counting options */
overhang = '99'  /* This is the value listed in STAR's sjdbOverhang parameter */
countMetaFeat = 'gene_id' /* meta-feature used to group features by. Usually "gene_id" */
countFeature = 'exon' /* Feature in annotation that is considered for counting. Usually "exon" */
strand = '2' /* Strandedness of reads: 0 = unstranded, 1 = forward stranded, 2 = reverse stranded */
orientation = 'fr' /* Read orientation of PE reads; f = forward strand, r = reverse strand. Can be 'fr', 'ff', or 'rf', but normally is 'fr'. */



/*Output paths*/
trimPath = "$projectPath/results/trimmomatic"
fastqcPath = "$projectPath/results/fastqc_afterTrim"
alignPath = "$projectPath/results/alignments"
countPath = "$projectPath/results/counts"


/* salmonquantallPath = "$projectPath/Salmon_quant_all"
 * salmonquantlongestPath = "$projectPath/Salmon_quant_longest"
*/

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
    set val(id), "${read1.baseName}.qualtrim.paired.fastq", "${read2.baseName}.qualtrim.paired.fastq" into starChannel
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


/*
* Step 3. STAR alignment
*/

limitM = starMemory + '000000000'

process runSTAR {
    executor 'slurm'
    cpus starCPU
    clusterOptions myQueue
    memory "$starMemory GB"
    /* module "$starMod"
    */
    publishDir alignPath, mode: 'link'

    input:
    file genome_file
    set pair_id, file(read1), file(read2) from starChannel

    output:
    file "${pair_id}_Aligned.sortedByCoord.out.bam" into align_channel
    file "${pair_id}_Log.*.out"
    file "${pair_id}_SJ.out.tab"
    file "${pair_id}__STARgenome/*.t*"

    """
    STAR --runThreadN $starCPU \
     --genomeDir ${genome_file} \
     --readFilesIn ${read1} ${read2}\
     --sjdbGTFfile ${gtf_file} \
     --outFileNamePrefix ${pair_id}_ \
     --sjdbGTFtagExonParentGene $countMetaFeat \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbOverhang $overhang \
     --limitGenomeGenerateRAM $limitM
    """
}



/*
* Step 4. Count data with FeatureCounts
* WARNING: considers '1' a valid exit status to get around wrapper error
*/
process featureCount {
    executor 'slurm'
    cpus countCPU
    clusterOptions myQueue
    memory "$countMemory GB"
    /* module "$featCountMod"
    */
    publishDir countPath, mode: 'link'
    validExitStatus 0,1

    input:
    file in_count from align_channel

    output:
    file "*_featCounts.txt*"

    """
    $HOME/local/miniconda2/bin/featureCounts -T $countCPU -s $strand -S $orientation -C -p -g $countMetaFeat \
    -t $countFeature -o ${in_count.baseName}_featCounts.txt -a ${gtf_file} ${in_count}
    """
}

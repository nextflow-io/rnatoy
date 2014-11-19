/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.pair1 = "$baseDir/data/ggal/*_1.fq.gz"
params.pair2 = "$baseDir/data/ggal/*_2.fq.gz"
params.suffix = 8
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.chunkSize = 1_0000_000

log.info "R N A T O Y   P I P E L I N E      "
log.info "================================="
log.info "genome             : ${params.genome}"
log.info "annotat            : ${params.annot}"
log.info "pair1              : ${params.pair1}"
log.info "pair2              : ${params.pair2}"
log.info "chunkSize          : ${params.chunkSize}"
log.info "suffix             : ${params.suffix}" 

 
/*
 * emits all reads ending with "_1" suffix and map them to pair containing the common
 * part of the name
 */

count1 = 0
reads1 = Channel
            .fromPath( params.pair1 )
            .splitFastq(by: params.chunkSize, meta: 'path') { chunk, path -> 
               def split = cacheableFile([path,params.chunkSize,count1++], "chunk1_${count1}.fastq")
               if(!split.exists()) split.text = chunk 
               log.debug "read1 split: $split"
               tuple( path.name[0..-params.suffix], split ) 
            }

/*
 * as above for "_2" read pairs
 */

count2 = 0
reads2 = Channel
            .fromPath( params.pair2 )
            .splitFastq(by: params.chunkSize, meta: 'path') { chunk, path ->
               def split = cacheableFile([path,params.chunkSize,count2++], "chunk2_${count2}.fastq")
               if(!split.exists()) split.text = chunk
               log.debug "read2 split: $split"
               tuple( path.name[0..-params.suffix], split )
            } 
     
/*
 * Match the pairs emittedb by "read1" and "read2" channels having the same 'key'
 * and emit a new pair containing the expected read-pair files
 */
read_pairs = reads1
        .phase(reads2)
        .map { pair1, pair2 -> [ pair1[0], pair1[1], pair2[1] ] }
 
/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)
 
 
/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    input:
    file genome_file
     
    output:
    file 'genome.index*' into genome_index
       
    """
    bowtie2-build ${genome_file} genome.index
    """
 
}
 
/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
     
    input:
    file 'genome.index.fa' from genome_file 
    file annotation_file
    file genome_index from genome_index.first()
    set pair_id, file(read1), file(read2) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam
 
    """
    tophat2 -p ${task.cpus} --GTF $annotation_file genome.index ${read1} ${read2}
    """
}

group_bam = bam.groupTuple()
  

process merge {
    
    input: 
    set pair_id, file('bam?') from group_bam 
    
    output: 
    set pair_id, file('*.merge.bam') into merged_bam 
    
    """
    count=`ls -l bam* | wc -l`
    if [ "\$count" -gt "1" ]; then 
    samtools merge ${pair_id}.merge.bam bam*
    else
    cp bam* ${pair_id}.merge.bam  
    fi
    """

}


/*
 * Step 4. Assemples the transcript by using the "cufflinks" 
 */
process makeTranscript {
    input:
    file 'anno.gtf' from annotation_file
    set pair_id, file(bam_file) from merged_bam
     
    output:
    set pair_id, file('transcripts.gtf') into transcripts
 
    """
    cufflinks --no-update-check -q -p ${task.cpus} -g anno.gtf ${bam_file}
    """
}
 
/*
 * Step 5. Save trabscripts files and print them
 */
transcripts
  .subscribe { tuple ->
    tuple[1].copyTo( "transcript_${tuple[0]}.gtf" )
    println "Saving: transcript_${tuple[0]}.gtf"
  }






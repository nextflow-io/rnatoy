/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.pair1 = "$baseDir/data/ggal/*_1.fq"
params.pair2 = "$baseDir/data/ggal/*_2.fq"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
 
/*
 * emits all reads ending with "_1" suffix and map them to pair containing the common
 * part of the name
 */
reads1 = Channel
    .fromPath( params.pair1 )
    .map {  path -> [ path.baseName[0..-2], path ] }
  
/*
 * as above for "_2" read pairs
 */
reads2 = Channel
    .fromPath( params.pair2 )
    .map {  path -> [ path.baseName[0..-2], path ] }
     
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
    file genome_file
    file genome_index from genome_index.first()
    set pair_id, file(read1), file(read2) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam
 
    """
    tophat2 genome.index ${read1} ${read2}
    """
}
 
 
/*
 * Step 3. Assemples the transcript by using the "cufflinks" 
 */
process makeTranscript {
    input:
    set pair_id, bam_file from bam
     
    output:
    set pair_id, 'transcripts.gtf' into transcripts
 
    """
    cufflinks ${bam_file}
    """
}
 
/*
 * Step 4. Collects the trabscripts files and print them
 */
transcripts
  .collectFile() {
     [ "${it[0]}transcript", it[1] ]
  }
  .subscribe { println it }
  
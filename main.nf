/*
 * Copyright (c) 2013, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'RNA-Toy'.
 *
 *   RNA-Toy is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNA-Toy is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNA-Toy.  If not, see <http://www.gnu.org/licenses/>.
 */
 
/* 
 * Proof of concept Nextflow based RNA-Seq pipeline
 * 
 * Authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Emilio Palumbo <emiliopalumbo@gmail.com> 
 */ 

 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.pair1 = "$baseDir/data/ggal/*_1.fq.gz"
params.pair2 = "$baseDir/data/ggal/*_2.fq.gz"
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

 
/*
 * emits all reads ending with "_1" suffix and map them to pair containing the common
 * part of the name
 */

reads1 = Channel
            .fromPath( params.pair1 )
            .ifEmpty { error "Cannot find any reads matching: ${params.pair1}" }
            .map { path -> 
                    def prefix = readPrefix(path, params.pair1)
                    log.debug "Splitting read1: [$prefix] $path"
                    tuple( prefix, path ) 
                 }
            .splitFastq(by: params.chunkSize, file: true)

/*
 * as above for "_2" read pairs
 */

reads2 = Channel
            .fromPath( params.pair2 )
            .ifEmpty { error "Cannot find any reads matching: ${params.pair2}" }
            .map { path -> 
                    def prefix = readPrefix(path, params.pair2)
                    log.debug "Splitting read2: [$prefix] $path"
                    tuple( prefix, path ) 
                 }
            .splitFastq(by: params.chunkSize, file:true)
     
/*
 * Match the pairs emittedb by "read1" and "read2" channels having the same 'key'
 * and emit a new pair containing the expected read-pair files
 */
read_pairs = reads1
        .phase(reads2)
        .ifEmpty { error "Cannot find any matching reads" }
        .map { pair1, pair2 -> tuple(pair1[0], pair1[1], pair2[1]) }
 
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
    cufflinks --no-update-check -q -p ${task.cpus} -G anno.gtf ${bam_file}
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



/* 
 * Helper function, given a file Path 
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 * 
 * For example: 
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 * 
 * Returns: 
 *   'file_alpha'
 */
 
def readPrefix( actual, template ) {

    final fileName = actual instanceof Path ? actual.name : actual.toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') ) 
        filePattern = '*' + filePattern 
  
    def regex = filePattern.replace('.','\\.').replace('*','(.*)').replace('?','(.?)')

    def matcher = (fileName =~ /$regex/  )
    if( matcher.matches() ) { 
        def end = matcher.end(matcher.groupCount() )      
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') ) 
          prefix=prefix[0..-2]
          
        return prefix
    }
    
    return null
}
  



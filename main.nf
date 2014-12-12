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
params.pair1 = "$baseDir/data/ggal/*_1.fq"
params.pair2 = "$baseDir/data/ggal/*_2.fq"
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

log.info "R N A T O Y   P I P E L I N E      "
log.info "================================="
log.info "genome             : ${params.genome}"
log.info "annotat            : ${params.annot}"
log.info "pair1              : ${params.pair1}"
log.info "pair2              : ${params.pair2}"

 
/*
 * emits all reads ending with "_1" suffix and map them to pair containing the common
 * part of the name
 */
reads1 = Channel
    .fromPath( params.pair1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.pair1}" }
    .map { path -> 
       def prefix = readPrefix(path, params.pair1)
       tuple(prefix, path) 
    }
  
/*
 * as above for "_2" read pairs
 */
reads2 = Channel
    .fromPath( params.pair2 )
    .ifEmpty { error "Cannot find any reads matching: ${params.pair2}" }
    .map { path -> 
       def prefix = readPrefix(path, params.pair2)
       tuple(prefix, path) 
    }
     
/*
 * Match the pairs emittedb by "read1" and "read2" channels having the same 'key'
 * and emit a new pair containing the expected read-pair files
 */
read_pairs = reads1
                .phase(reads2)
                .ifEmpty { error "Cannot find any matching reads" }
                .map { pair1, pair2 -> [ pair1[0], pair1[1], pair2[1] ] }
 
/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)
 

/*
 * Step 1. Prepare the transcriptome from reference genome and annotation
 */
 process prepareTranscriptome {
    input:
    file genome_file
    file annotation_file
    
    output:
    file "transcriptome.idx.fa" into transcriptome
    
    """
    rsem-prepare-reference -gtf ${annotation_file} ${genome_file} transcriptome
    """

 }

/*
 * Step 2. Builds the transcriptome index required by the mapping process
 */
process buildIndex {
    input:
    file transcriptome_file from transcriptome
     
    output:
    file 'transcriptome.index*' into transcriptome_index
       
    """
    bowtie2-build ${transcriptome_file} transcriptome.index
    """
 
}

/*
 * Step 3. Maps each read-pair by using the Bowtie2 mapper tool
 */
process mapping {
     
    input:
    file transcriptome_index from transcriptome_index.first()
    set pair_id, file(read1), file(read2) from read_pairs
 
    output:
    set pair_id, "tx_map.bam" into bam
 
    """
    bowtie2 -x transcriptome.index -p ${task.cpus} -1 ${read1} -2 ${read2} | samtools view -@ ${task.cpus} -Sb - | samtools sort -@ ${task.cpus} - tx_map
    """
}

/*
 * Step 3. Assemples the transcript by using the "cufflinks" 
 */
process makeTranscript {
    input:
    file 'anno.gtf' from annotation_file
    set pair_id, file(bam_file) from bam
     
    output:
    set pair_id, file('transcripts.gtf') into transcripts
 
    """
    cufflinks --no-update-check -q -p ${task.cpus} -G anno.gtf ${bam_file}
    """
}
 
/*
 * Step 4. Collects the trabscripts files and print them
 */
transcripts
  .subscribe { tuple ->
    def fileName = "transcript_${tuple[0]}.gtf" 
    tuple[1].copyTo(fileName)
    println "Saving: $fileName"
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
 
def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

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
  

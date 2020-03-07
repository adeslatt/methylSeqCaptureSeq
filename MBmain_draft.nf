#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/MBmethylseq
========================================================================================
 nf-core/MBmethylseq Analysis Pipeline.
----------------------------------------------------------------------------------------
*/



/**********************************
 *      PIPELINE PROCESSES        *
 **********************************/

#NotDone
	#downloadfiles --> bam2fastq
	#index_reference_hg19_withbiscuit


/*
 * Fastqc
 */
 
// main.nf
params.reads = false

reads = Channel.fromFilePairs(params.reads, size: 2)

process fastqc {

    tag "$name"
    publishDir "results", mode: 'copy'
    container 'flowcraft/fastqc:0.11.7-1'

    input:
    set val(name), file(reads) from reads

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc $reads
    """
}

/*
 * trim reads with trim_galore
 */

process trim_galore {
        tag "$name"
        publishDir "results/trim_galore", mode: 'copy',
       	container 'mskaccess/trim_galore:0.6.2'
        input:
        set val(name), file(reads) from reads

        output:
        file "*fq.gz" into trimmed_reads_for_alignment
        file "*trimming_report.txt" into trim_galore_results
        file "*_fastqc.{zip,html}"
    

        script:
            """
            trim_galore --paired --fastqc --gzip $reads
            """
        }


/*
 * Align paired-end read files to the human genome 19 with biscuit
 */
 
reads = Channel.fromFilePairs(params.reads, size: 2)

process biscuitalign {
    tag "$name"
    publishDir "results/biscuitalign", mode: 'copy'
    
    input:
    set val(name), file(reads) from trimmed_reads_for_alignment
    file index from biscuit_hg19_index_align.collect()
    
    output:
    set val(name), file("*.bam") into biscuit2bams
    
    script:
        """
            biscuit align -t 10 ${params.biscuitindex} $read1 $read2 
            		| samtools view -hb > ${name}_biscuit.bam
           
        """
    }

/*
 * Sort BAM files with SAMTOOLS
 */

process sortbam {
    tag "sortbam: $name"
    
    input:
    set val(name), file(bam) from biscuit2bams
    
    output:
    set val(name), file("${name}.sorted.bam") into sortedBams
    
    script:
    avail_mem=""
    """
    samtools sort \\
        $bam \\
        -@ ${task.cpus} ${avail_mem} \\
        -o ${name}.sorted.bam
    samtools index ${name}.sorted.bam
    """
}

/*
 * Mark duplicates and sort with picard tools
 */

process markduplicates {
    tag "markdups: $name"
    
    input:
    set val(name), file(bam) from sortedBams
    
    output:
    set val(name), file("${name}.sorted.nodup.bam") into markedDups, bamstatsSamples
    file("${name}.metric.txt") into markduplicatesMultiqc
    
    script:
    // Runtime java Mark duplicates options
    markdup_java_options="-Xmx${task.memory.toGiga()}g"
    """
    picard ${markdup_java_options} MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${name}.sorted.nodup.bam \\
        METRICS_FILE=${name}.metric.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

#NotDone 
#filterbam_bamtools filter

/*
 * Generate BAM file statistics
 */
 
process bamstats {
    tag "bamstats: $name"

    publishDir = [path: "results/flagstats", mode: 'copy', overwrite: 'true' ]
    
    input:
    set val(name), file(bam) from bamstatsSamples
    
    output:
    file("${name}.sorted.nodup.bam.flagstat.txt") into flagstatMultiqc
    
    script:
    """
    samtools flagstat $bam > ${name}.sorted.nodup.bam.flagstat.txt
    samtools view $bam | cut -f 10 | perl -ne 'chomp;print length(\$_) . "\n"' | sort | uniq -c >> ${name}.sorted.nodup.bam.flagstat.txt
    """
}

#NotDone
#Biscuitpileup
#Biscuitvcf2bed
#wgbs_tools_bam2pat

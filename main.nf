//main.nf
#genomeindexing 
#biscuit index <path_to_genome_folder> hg19.fa

#alignment
#biscuit align -t 30 <path_to_genome_folder> hg19.fa A1/*1.fq.gz A1/*2.fq.gz | samtools view -Sb > A1.bam
e with accession list to fetch from SRA
      --gencodeFile                 Path to input gencode GTF file
            --keyFile                     Path to a keyfile used to fetch restricted access datasets with SRAtools
	          -profile                      Configuration profile to use. Can use multiple (comma separated)
		                                      Available: short-test, key-test, ...
						          Generic:
							        --skipTrimming                Sets if trimming has to be skipped
								      --readType                    Specifies what type of input reads are to be processed: single end or paired end
								            --readLength                  Specifies the read length
									        """.stripIndent()
}

/*********************************
 *      CHANNELS SETUP           *
 *********************************/

alternativeSplicingTypeList = ['a3ss', 'a5ss', 'mxe', 'ri', 'se']
junctionCountTypeList       = ['jc','jcec']
countingTypeList            = ['ijc','sjc','inc','inclen','skiplen']

key_file = file(params.keyFile)

if(!params.starIndexPath) {
  exit 1, "Cannot find STAR index path for rMATs"
}

if(!params.splitNumber) {
    exit 1, "Cannot find a splitNumber value set for createMatrices"
}

if(params.accessionList) {
    Channel
        .fromPath( params.accessionList )
        .splitText()
        .map{ it.trim() }
        .dump(tag: 'AccessionList content')
        .set { accessionIDs }
} else {
    exit 1, "Accession list file not provided"
}

if(params.gencodeFile) {
    Channel
        .value(file(params.gencodeFile))
        .ifEmpty { error "Cannot find any Gencode GTF annotaton for parameter --gencode: ${params.gencodeFile}" }
        .set { gencodeFile }    
} else {
    exit 1, "Gencode annotation file not provided"
}

if (!params.readLength) {
    exit 1, "Read length parameter not provided"
}
if (!params.readType) {
    exit 1, "Read type parameter not provided"
}

/**********************************
 *      PIPELINE PROCESSES        *
 **********************************/

/*
 * Get accession samples from SRA with or without keyFile
 */
if(!params.reads) {
    process getAccession {
        tag "${accession}"
        
        input:
        val accession from accessionIDs
        file keyFile from key_file
        
        output:
        set val(accession), file("*.fastq.gz") into readFiles
        
        script:
        def vdbConfigCmd = keyFile.name != 'NO_FILE' ? "vdb-config --import ${keyFile} ./" : ''
        """
	        $vdbConfigCmd
		        fasterq-dump $accession --threads ${task.cpus} --split-3
			        pigz *.fastq
				        """
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.readType == 'single' ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set { readFiles }
}

// ------------------------------------------------------------------------------   
//	QUALITY CONTROL 
// ------------------------------------------------------------------------------   

/**
	Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

	This step is OPTIONAL. De-duplication should be carried on iff you are
    using PCR amplification (in this case identical reads are technical artefacts)
	but not otherwise (identical reads will identify natural duplicates).
*/

process dedup {
	
    tag "$name"
    
	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }
		
	input:
	tuple val(name), file(reads) 

	output:
	tuple val(name), path("${name}_dedup*.fq.gz"), emit: to_synthetic_contaminants
	path("dedup_mqc.yaml"), emit: dedup_log
	
	when:
	params.mode != "characterisation" && params.dedup

	script:
	// This is to deal with single and paired end reads
	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
	def output = params.singleEnd ? "out=\"${name}_dedup.fq.gz\"" :  "out1=\"${name}_dedup_R1.fq.gz\" out2=\"${name}_dedup_R2.fq.gz\""
	
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo \"$task.memory\" | sed 's/ //g' | sed 's/B//g')
    clumpify.sh -Xmx\"\$maxmem\" $input $output qin=$params.qin dedupe subs=0 threads=${task.cpus} &> dedup_mqc.txt
	
	# MultiQC doesn't have a module for clumpify yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_dedup_log.sh > dedup_mqc.yaml
	"""
}

/**
	Quality control - STEP 2. A decontamination of synthetic sequences and artefacts 
	is performed.
*/

process remove_synthetic_contaminants {
	
	tag "$name"
	
	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	input:
	tuple file(artefacts), file(phix174ill), val(name), file(reads) 
   
	output:
	tuple val(name), path("${name}_no_synthetic_contaminants*.fq.gz"), emit: to_trim
	path("synthetic_contaminants_mqc.yaml"), emit: synthetic_contaminants_log
	
	when:
	params.mode != "characterisation"

   	script:
	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
	def output = params.singleEnd ? "out=\"${name}_no_synthetic_contaminants.fq.gz\"" :  "out=\"${name}_no_synthetic_contaminants_R1.fq.gz\" out2=\"${name}_no_synthetic_contaminants_R2.fq.gz\""
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	bbduk.sh -Xmx\"\$maxmem\" $input $output k=31 ref=$phix174ill,$artefacts qin=$params.qin threads=${task.cpus} ow &> synthetic_contaminants_mqc.txt
	
	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_remove_synthetic_contaminants_log.sh > synthetic_contaminants_mqc.yaml
	"""
}


/**
	Quality control - STEP 3. Trimming of low quality bases and of adapter sequences. 
	Short reads are discarded. 
	
	If dealing with paired-end reads, when either forward or reverse of a paired-read
	are discarded, the surviving read is saved on a file of singleton reads.
*/


process trim {

	tag "$name"
	
	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }
	
	input:
	tuple file(adapters), val(name), file(reads)
	
	output:
	tuple val(name), path("${name}_trimmed*.fq.gz"), emit: to_decontaminate
	path("trimming_mqc.yaml"), emit: trimming_log
	
	when:
	params.mode != "characterisation"

   	script:
	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
	def output = params.singleEnd ? "out=\"${name}_trimmed.fq.gz\"" :  "out=\"${name}_trimmed_R1.fq.gz\" out2=\"${name}_trimmed_R2.fq.gz\" outs=\"${name}_trimmed_singletons.fq.gz\""
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	bbduk.sh -Xmx\"\$maxmem\" $input $output ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe ow &> trimming_mqc.txt

	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_trimming_log.sh > trimming_mqc.yaml
	"""
}


/**
	Quality control - STEP 4. Decontamination. Removes external organisms' contamination, 
	using given genomes. 

	When an indexed contaminant (pan)genome is not provided, the index_foreign_genome process is run 
	before the decontamination process. This process require the FASTA file of the contaminant (pan)genome.
*/

process index_foreign_genome {
	
	//tag "index foreign genome"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	input:
	file(foreign_genome)

	output:
	path "ref/", emit: ref_foreign_genome

	when:
	params.mode != "characterisation" && params.foreign_genome_ref == ""

	script:
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	# This step will have a boilerplate log because the information saved by bbmap are not relevant
	bbmap.sh -Xmx\"\$maxmem\" ref=$foreign_genome &> foreign_genome_index_mqc.txt
	"""
}

process decontaminate {

    tag "$name"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*QCd.fq.gz"

	input:
	tuple path(ref_foreign_genome), val(name), file(reads)

	output:
	tuple val(name), path("*_QCd.fq.gz"), emit: qcd_reads
	path("decontamination_mqc.yaml"), emit: decontaminate_log

	when:
	params.mode != "characterisation"

	script:
	// When paired-end are used, decontamination is carried on independently on paired reads
	// and on singleton reads thanks to BBwrap, that calls BBmap once on the paired reads
	// and once on the singleton ones, merging the results on a single output file
	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\",\"${reads[2]}\" in2=\"${reads[1]}\",null"
	def output = "outu=\"${name}_QCd.fq.gz\" outm=\"${name}_contamination.fq.gz\""
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	bbwrap.sh -Xmx\"\$maxmem\"  mapper=bbmap append=t $input $output minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred path="./" qin=$params.qin threads=${task.cpus} untrim quickmatch fast ow &> decontamination_mqc.txt

	# MultiQC doesn't have a module for bbwrap yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_decontamination_log.sh > decontamination_mqc.yaml
	"""
}


// ------------------------------------------------------------------------------   
//	QUALITY ASSESSMENT 
// ------------------------------------------------------------------------------   


process quality_assessment {
	
    tag "$name"
	
	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_fastqc
    } else {
        container params.docker_container_fastqc
    }
	
	publishDir "${params.outdir}/${params.prefix}/fastqc", mode: 'copy' //,

    input:
    tuple val(name), file(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_log
	
	when:
	params.mode != "characterisation"

    script:
    """
    fastqc -q $reads
    """
}


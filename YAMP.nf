#!/usr/bin/env nextflow
	
/**
Yet Another Metagenomic Pipeline (YAMP)
Copyright (C) 2017-2021	Dr Alessia Visconti 	      
	      
This script is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
	
This script is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
	
You should have received a copy of the GNU General Public License
along with this script. If not, see <http://www.gnu.org/licenses/>.
	
For any bugs or problems found, please go to:
- https://github.com/alesssia/YAMP/issues
*/


def versionMessage() 
{
	log.info"""
	 
	YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: ${workflow.manifest.version} 
	""".stripIndent()
}

def helpMessage() 
{
	log.info"""
		
	YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: ${workflow.manifest.version} 
	
	This pipeline is distributed in the hope that it will be useful
	but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.
	
	Please report comments and bugs at https://github.com/alesssia/YAMP/issues.
	Check https://github.com/alesssia/YAMP for updates, and refer to
	https://github.com/alesssia/YAMP/wiki for more details.
	
	Usage: 
	nextflow run YAMP.nf --reads1 R1 --reads2 R2 --prefix mysample --outdir path --mode MODE  
	[options] [-with-docker|-with-singularity]
	
	Mandatory arguments:
	--reads1   R1      Forward (if paired-end) OR all reads (if single-end) file path
	[--reads2] R2      Reverse reads file path (only if paired-end library layout)
	--prefix   prefix  Prefix used to name the result files
	--outdir   path    Output directory (will be outdir/prefix/)
	--mode     <QC|characterisation|complete>
	
	Options:
	--singleEnd     <true|false>   whether the layout is single-end
	--dedup         <true|false>   whether to perform de-duplication

	Please refer to nextflow.config for more options.
		
	YAMP supports FASTQ and compressed FASTQ files.
	""".stripIndent()
}

/**
Prints version when asked for
*/
if (params.version) {
	versionMessage()
	exit 0
}

/**
Prints help when asked for
*/

if (params.help) {
	helpMessage()
	exit 0
}

/**
STEP 0. 
	
Checks input parameters and (if it does not exists) creates the directory 
where the results will be stored (aka working directory). 
Initialises the log file.
	
The working directory is named after the prefix and located in the outdir 
folder. The log file, that will save summary statistics, execution time,
and warnings generated during the pipeline execution, will be saved in the 
working directory as "prefix.log".
*/


//Checking user-defined parameters	
if (params.mode != "QC" && params.mode != "characterisation" && params.mode != "complete") {
	exit 1, "Mode not available. Choose any of <QC, characterisation, complete>"
}	

if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

//--reads2 can be omitted when the library layout is "single" (indeed it specifies single-end
//sequencing)
if (params.mode != "characterisation" && !params.singleEnd && (params.reads2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
}

//--reads1 and --reads2 can be omitted (and the default from the config file used instead) 
//only when mode is "characterisation". Obviously, --reads2 should be always omitted when the
//library layout is single.
if (params.mode != "characterisation" && ( (!params.singleEnd && (params.reads1 == "null" || params.reads2 == "null")) || (params.singleEnd && params.reads1 == "null")) ) {
	exit 1, "Please set the reads1 and/or reads2 parameters"
}

//Creates working dir
workingpath = params.outdir + "/" + params.prefix
workingdir = file(workingpath)
if( !workingdir.exists() ) {
	if( !workingdir.mkdirs() ) 	{
		exit 1, "Cannot create working directory: $workingpath"
	} 
}	


// Header log info
log.info """---------------------------------------------
YET ANOTHER METAGENOMIC PIPELINE (YAMP) 
---------------------------------------------

Analysis introspection:

"""

def summary = [:]

summary['Starting time'] = new java.util.Date() 
//Environment
summary['Environment'] = ""
summary['Pipeline Name'] = 'YAMP'
summary['Pipeline Version'] = workflow.manifest.version

summary['Config Profile'] = workflow.profile
summary['Resumed'] = workflow.resume
		
summary['Nextflow version'] = nextflow.version.toString() + " build " + nextflow.build.toString() + " (" + nextflow.timestamp + ")"

summary['Java version'] = System.getProperty("java.version")
summary['Java Virtual Machine'] = System.getProperty("java.vm.name") + "(" + System.getProperty("java.vm.version") + ")"

summary['Operating system'] = System.getProperty("os.name") + " " + System.getProperty("os.arch") + " v" +  System.getProperty("os.version")
summary['User name'] = System.getProperty("user.name") //User's account name

summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container

if (params.mode != "characterisation") 
{
	if (params.enable_conda){
		summary['BBmap'] = "bioconda::bbmap=38.87-0"
		summary['FastQC'] = "bioconda::fastqc=0.11.9-0"
		summary['MultiQC'] = "bioconda::multiqc=1.9-1"
	} else if (workflow.containerEngine == 'singularity') {
	   	summary['BBmap'] = "https://depot.galaxyproject.org/singularity/bbmap:38.87--h1296035_0"
		summary['FastQC'] = "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
		summary['MultiQC'] = "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"
	} else if (workflow.containerEngine == 'docker') {
    	summary['BBmap'] = "quay.io/biocontainers/bbmap:38.87--h1296035_0"
		summary['FastQC'] = "quay.io/biocontainers/fastqc:0.11.9--0"
		summary['MultiQC'] = "quay.io/biocontainers/multiqc:1.9--py_1"
	} else {
		summary['BBmap'] = "No container information"
		summary['FastQC'] = "No container information"
		summary['MultiQC'] = "No container information"
	}
}

// TODO: complete
if (params.mode != "QC")
{
	if (params.enable_conda){
		summary['biobakery'] = "biobakery::humann"
	} else if (workflow.containerEngine == 'singularity') {
		summary['biobakery'] = "biobakery/workflows:3.0.0.a.6.metaphlanv3.0.7"
	} else if (workflow.containerEngine == 'docker') {
		summary['biobakery'] = "biobakery/workflows:3.0.0.a.6.metaphlanv3.0.7"
	} else {
		summary['biobakery'] = "No container information"
	}
}

if(workflow.profile == 'awsbatch'){
	summary['AWS Region'] = params.awsregion
	summary['AWS Queue'] = params.awsqueue
}

//General
summary['Running parameters'] = ""
summary['Reads'] = "[" + params.reads1 + ", " + params.reads2 + "]"
summary['Prefix'] = params.prefix
summary['Running mode'] = params.mode

if (params.mode != "characterisation") 
{
	summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'
	summary['Performing de-duplication'] = params.dedup

	//Trimming
	summary['Trimming parameters'] = ""
	summary['Input quality offset'] = params.qin == 33 ? 'ASCII+33' : 'ASCII+64'
	summary['Min phred score'] = params.phred
	summary['Min length'] = params.minlength
	summary['Shorter kmer'] = params.mink 
	summary['Max Hamming distance'] = params.hdist 

	//Decontamination
	summary['Decontamination parameters'] = ""
	if (params.foreign_genome_ref != "") {
		summary['Contaminant (pan)genome'] = params.foreign_genome_ref + " (indexed)"
	} else if (	params.foreign_genome_ref == "") {
		summary['Contaminant (pan)genome'] = params.foreign_genome
	}	
	summary['Min alignment identity'] = params.mind
	summary['Max indel length'] = params.maxindel
	summary['Max alignment band'] = params.bwr
}

if (params.mode != "QC")
{
	//TODO add MP/QUIIME/HM3 options
}

//Folders
summary['Folders'] = ""
summary['Output dir'] = workingpath
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Lunching dir'] = workflow.launchDir

//FIXME: adding ending time


log.info summary.collect { k,v -> "${k.padRight(27)}: $v" }.join("\n")
log.info ""


/**
	Prepare workflow introspection

	This process adds the workflow introspection (also printed at runtime) in the logs
	This is NF-CORE code.
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'workflow-summary'
    description: "This information is collected when the pipeline is started."
    section_name: 'YAMP Workflow Summary'
    section_href: 'https://github.com/alesssia/yamp'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd>$v</dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/**
	Gets software version. 

	This process ensures that software version are included in the logs.
*/

process get_software_versions {

	//Starting the biobakery container. I need to run metaphlan and Humann to get
	//their version number (due to the fact that they live in the same conrtainer)
    conda (params.enable_conda ? params.conda_biobakery : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	output:
	file "software_versions_mqc.yaml" into software_versions_yaml

	script:
	// TODO: Get all tools to print their version number here
	//I am using a multi-containers scenarios, supporting conda, docker, and singularity
	//with the software at a specific version (the same for all platforms). Therefore, I
	//will simply parse the version from there. Perhaps overkill, but who cares?  
	//This is not true for the biobakery suite (metaphlan/humann) which extract the 
	//information at runtime from the actula commands (see comment above)
	"""
	echo $workflow.manifest.version > v_pipeline.txt
	echo $workflow.nextflow.version > v_nextflow.txt

	echo $params.conda_fastqc | cut -d= -f 2 > v_fastqc.txt
	echo $params.conda_bbmap | cut -d= -f 2 > v_bbmap.txt
	echo $params.conda_multiqc | cut -d= -f 2 > v_multiqc.txt
	
	metaphlan --version > v_metaphlan.txt
	humann --version > v_humann.txt
	
	scrape_software_versions.py > software_versions_mqc.yaml
	"""
}

/**
	Creates a set of channels for input read files.
	- read_files_fastqc is used for the first QC assessment (on the raw reads)
	- read_files_dedup  is used for the deduplication step (which is optional and may skip to trimming)
	- read_files_trim   is used for the decontamination from synthetic contaminants (used only if
	  deduplication is not run)
*/

if (params.singleEnd) {
	Channel
	.from([[params.prefix, [file(params.reads1)]]])
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
} else {
	Channel
	.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] )
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
}

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
    conda (params.enable_conda ? params.conda_bbmap : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }
		
	input:
	tuple val(name), file(reads) from read_files_dedup

	output:
	tuple val(name), path("${name}_dedup*.fq.gz") into to_synthetic_contaminants
	file "dedup_mqc.yaml" into dedup_log
	
	when:
	params.mode != "characterisation" && params.dedup

	script:
	// This is to deal with single and paired end reads
	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
	def output = params.singleEnd ? "out=\"${name}_dedup.fq.gz\"" :  "out1=\"${name}_dedup_R1.fq.gz\" out2=\"${name}_dedup_R2.fq.gz\""
	
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo \"$task.memory\" | sed 's/ //g' | sed 's/B//g')
	echo \"$reads\"
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

//When the de-suplication is not done, the raw file should be pushed in the corret channel
//FIXME: make this also optional?
if (!params.dedup) {
	to_synthetic_contaminants = read_files_synthetic_contaminants
	dedup_log = Channel.from(file("$baseDir/assets/no_dedup.yaml"))
}

// Defines channels for resources file 
Channel.fromPath( "${params.artefacts}", checkIfExists: true ).set { artefacts }
Channel.fromPath( "${params.phix174ill}", checkIfExists: true ).set { phix174ill }

process remove_synthetic_contaminants {
	
	tag "$name"
	
	//Enable multicontainer settings
    conda (params.enable_conda ? params.conda_bbmap : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	input:
	tuple file(artefacts), file(phix174ill), val(name), file(reads) from artefacts.combine(phix174ill).combine(to_synthetic_contaminants)
   
	output:
	tuple val(name), path("${name}_no_synthetic_contaminants*.fq.gz") into to_trim
	file "synthetic_contaminants_mqc.yaml" into synthetic_contaminants_log
	
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

// Defines channels for resources file 
Channel.fromPath( "${params.adapters}", checkIfExists: true ).set { adapters }

process trim {

	tag "$name"
	
	//Enable multicontainer settings
    conda (params.enable_conda ? params.conda_bbmap : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }
	
	input:
	tuple file(adapters), val(name), file(reads) from adapters.combine(to_trim) 
	
	output:
	tuple val(name), path("${name}_trimmed*.fq.gz") into to_decontaminate
	file "trimming_mqc.yaml" into trimming_log
	
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

// Defines channels for foreign_genome file 
Channel.fromPath( "${params.foreign_genome}", checkIfExists: true ).set { foreign_genome }

process index_foreign_genome {

	tag "$name"

	//Enable multicontainer settings
    conda (params.enable_conda ? params.conda_bbmap : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	input:
	file(foreign_genome) from foreign_genome

	output:
	path("ref/", type: 'dir') into ref_foreign_genome
	
	when:
	params.mode != "characterisation" && params.foreign_genome_ref == ""

	script:
	"""
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	# This step will have no logging because the information saved by bbmap are not relevant
	bbmap.sh -Xmx\"\$maxmem\" ref=$foreign_genome &> foreign_genome_index_mqc.txt
	"""
}

//Channel.fromPath( "${params.foreign_genome_ref}", checkIfExists: true ).set { ref_foreign_genome }

//When the indexed contaminant (pan)genome is already available, its path should be pushed in the corret channel
if (params.foreign_genome_ref != "") {
	ref_foreign_genome = Channel.from(file(params.foreign_genome_ref))
}

process decontaminate {

    tag "$name"

	//Enable multicontainer settings
    conda (params.enable_conda ? params.conda_bbmap : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*QCd.fq.gz"

	input:
	tuple path(ref_foreign_genome), val(name), file(reads) from ref_foreign_genome.combine(to_decontaminate)

	output:
	tuple val(name), path("*_QCd.fq.gz") into qcd_reads
	tuple val(name), path("*_QCd.fq.gz") into to_profile_taxa_decontaminated
	tuple val(name), path("*_QCd.fq.gz") into to_profile_functions_decontaminated
	file "decontamination_mqc.yaml" into decontaminate_log

	when:
	params.mode != "characterisation"

	script:
	// When paired-end are used, decontamination is carried on idependently on paired reads
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
    conda (params.enable_conda ? params.conda_fastqc : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_fastqc
    } else {
        container params.docker_container_fastqc
    }
	
	publishDir "${params.outdir}/${params.prefix}/fastqc", mode: 'copy' //,

    input:
    set val(name), file(reads) from read_files_fastqc.mix(qcd_reads)

    output:
    path "*_fastqc.{zip,html}" into fastqc_log
	
	when:
	params.mode != "characterisation"

    script:
    """
    fastqc -q $reads
    """
}

// ------------------------------------------------------------------------------   
//  COMMUNITY CHARACTERISATION 
// ------------------------------------------------------------------------------   

// The user will specify the clean file either as a single clean file (that is the YAMP
// default behaviour), or as two files (foward/reverse). ]
// In the former case, the user will set singleEnd = true and only one file will be 
// selected and used direcly for taxa and community profiling.
// In the latter case, the user will set singleEnd = false and provide two files, that will
// be merged before feeding the relevant channels for profiling.
if (params.mode == "characterisation" && params.singleEnd) {
	Channel
	.from([[params.prefix, [file(params.reads1)]]])
	.into { reads_profile_taxa; reads_profile_functions }
	
	//Initialise empty channels
	reads_merge_paired_end_cleaned = Channel.empty()
	merge_paired_end_cleaned_log = Channel.empty()
} else if (params.mode == "characterisation" && !params.singleEnd) {
	Channel
	.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] )
	.set { reads_merge_paired_end_cleaned }
	
	//Stage boilerplate log
	merge_paired_end_cleaned_log = Channel.from(file("$baseDir/assets/merge_paired_end_cleaned_mqc.yaml"))
} else if (params.mode != "characterisation")
	to_merge_paired_end_cleaned = Channel.empty()

process merge_paired_end_cleaned {

	tag "$name"
		
	input:
	tuple val(name), file(reads) from reads_merge_paired_end_cleaned
	
	output:
	tuple val(name), path("*_QCd.fq.gz") into to_profile_taxa_merged
	tuple val(name), path("*_QCd.fq.gz") into to_profile_functions_merged
	
	when:
	params.mode == "characterisation" && !params.singleEnd

   	script:
	"""
	# This step will have no logging because the information are not relevant
	# I will simply use a boilerplate YAML to record that this has happened
	# If the files were not compressed, they will be at this stage
	if (file ${reads[0]} | grep -q compressed ) ; then
	    cat ${reads[0]} ${reads[1]} > ${name}_QCd.fq.gz
	else
		cat ${reads[0]} ${reads[1]} | gzip > ${name}_QCd.fq.gz
	fi
	"""
}

/**
	Community Characterisation - STEP 1. Performs taxonomic binning and estimates the 
	microbial relative abundancies using MetaPhlAn and its databases of clade-specific markers.
*/


process profile_taxa {

    tag "$name"

	//Enable multicontainer settings
    conda (params.enable_conda ? params.conda_biobakery : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*.{biom,tsv}"
	
	input:
	tuple val(name), file(reads) from to_profile_taxa_decontaminated.mix(to_profile_taxa_merged).mix(reads_profile_taxa)
	file (bowtie2db) from Channel.fromPath( params.metaphlan_databases, type: 'dir' )
	
	output:
	tuple val(name), path("*.biom") into to_alpha_diversity
	tuple val(name), path("*_metaphlan_bugs_list.tsv") into to_profile_function_bugs
	file "profile_taxa_mqc.yaml" into profile_taxa_log
	
	when:
	params.mode != "QC"
	
	script:
	"""
	#If a file with the same name is already present, Metaphlan2 used to crash, leaving this here just in case
	rm -rf ${name}_bt2out.txt
	
	metaphlan --input_type fastq --tmp_dir=. --biom ${name}.biom --bowtie2out=${name}_bt2out.txt --bowtie2db $bowtie2db --bt2_ps ${params.bt2options} --nproc ${task.cpus} $reads ${name}_metaphlan_bugs_list.tsv &> profile_taxa_mqc.txt
	
	# MultiQC doesn't have a module for Metaphlan yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_profile_taxa_log.sh > profile_taxa_mqc.yaml
	"""
}


/**
	Community Characterisation - STEP 2. Performs the functional annotation using HUMAnN.
*/

process profile_function {
	
    tag "$name"

	//Enable multicontainer settings
    conda (params.enable_conda ? params.conda_biobakery : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*.{tsv,log}"
	
	input:
	tuple val(name), file(reads) from to_profile_functions_decontaminated.mix(to_profile_functions_merged).mix(reads_profile_functions)
	tuple val(name), file(metaphlan_bug_list) from to_profile_function_bugs
	file (chocophlan) from Channel.fromPath( params.chocophlan, type: 'dir' )
	file (uniref) from Channel.fromPath( params.uniref, type: 'dir')
	
    output:
	file "*_HUMAnN.log"
	file "*_genefamilies.tsv"
	file "*_pathcoverage.tsv"
	file "*_pathabundance.tsv"
	file "profile_functions_mqc.yaml" into profile_functions_log

	when:
	params.mode != "QC"

	script:
	"""
	#HUMAnN will uses the list of species detected by the profile_taxa process
	humann --input $reads --output . --output-basename ${name} --taxonomic-profile $metaphlan_bug_list --nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads ${task.cpus} --memory-use minimum &> ${name}_HUMAnN.log 
	
	# MultiQC doesn't have a module for humann yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_profile_functions.sh ${name} ${name}_HUMAnN.log > profile_functions_mqc.yaml
 	"""
}

// ------------------------------------------------------------------------------   
//	MULTIQC LOGGING
// ------------------------------------------------------------------------------   


/**
	Generate Logs. 

	Logs generate at each analysis step are collected and processed with MultiQC 
*/

// Stage config files
multiqc_config = file(params.multiqc_config)

process log {
	
	publishDir "${params.outdir}/${params.prefix}", mode: 'copy'

    conda (params.enable_conda ? params.conda_multiqc : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_multiqc
    } else {
        container params.docker_container_multiqc
    }

	input:
	file multiqc_config
	file workflow_summary from create_workflow_summary(summary)
	file "software_versions_mqc.yaml" from software_versions_yaml
	path "fastqc/*" from fastqc_log.collect().ifEmpty([])
	file "dedup_mqc.yaml" from dedup_log.ifEmpty([])
	file "synthetic_contaminants_mqc.yaml" from synthetic_contaminants_log.ifEmpty([])
	file "trimming_mqc.yaml" from trimming_log.ifEmpty([])
	file "decontamination_mqc.yaml" from decontaminate_log.ifEmpty([])
	file "merge_paired_end_cleaned_mqc.yaml" from merge_paired_end_cleaned_log.ifEmpty([])
	file "profile_taxa_mqc.yaml" from profile_taxa_log.ifEmpty([])
	
	output:
	path "*multiqc_report*.html" into multiqc_report
	path "*_data*"

	script:
	"""
	multiqc --config $multiqc_config . -f
	mv multiqc_report.html ${params.prefix}_multiqc_report_${params.mode}.html
	mv multiqc_data ${params.prefix}_multiqc_data_${params.mode}
	"""
}


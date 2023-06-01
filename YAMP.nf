#!/usr/bin/env nextflow
	
/**
Yet Another Metagenomic Pipeline (YAMP)
Copyright (C) 2017-2023	Dr Alessia Visconti 	      
	      
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


/**
  Import processes from modules
*/

nextflow.enable.dsl = 2

include { dedup; remove_synthetic_contaminants; trim; index_foreign_genome; decontaminate; quality_assessment } from './modules/quality_control'
include { profile_taxa; profile_function; alpha_diversity } from './modules/community_characterisation'
include { merge_paired_end_cleaned; get_software_versions; log } from './modules/house_keeping'


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
  nextflow run YAMP.nf --reads1 R1 --reads2 R2 --prefix prefix --outdir path [options] 
  
  Mandatory arguments:
    --reads1   R1      Forward (if paired-end) OR all reads (if single-end) file path
    [--reads2] R2      Reverse reads file path (only if paired-end library layout)
    --prefix   prefix  Prefix used to name the result files
    --outdir   path    Output directory (will be outdir/prefix/)
  
  Main options:
    --mode       <QC|characterisation|complete>
    --singleEnd  <true|false>   whether the layout is single-end
    --dedup      <true|false>   whether to perform de-duplication
  
  Other options:
  BBduk parameters for removing synthetic contaminants and trimming:
    --qin                 <33|64> Input quality offset 
    --kcontaminants       value   kmer length used for identifying contaminants
    --phred               value   regions with average quality BELOW this will be trimmed 
    --minlength           value   reads shorter than this after trimming will be discarded
    --mink                value   shorter kmer at read tips to look for 
    --hdist               value   maximum Hamming distance for ref kmer
    --artefacts           path    FASTA file with artefacts
    --phix174ill          path    FASTA file with phix174_ill
    --adapters            path    FASTA file with adapters         
  
  BBwrap parameters for decontamination:
    --foreign_genome      path    FASTA file for contaminant (pan)genome
    --foreign_genome_ref  path    folder for for contaminant (pan)genome (pre indexed)
    --mind                value   approximate minimum alignment identity to look for
    --maxindel            value   longest indel to look for
    --bwr                 value   restrict alignment band to this
  
  MetaPhlAn parameters for taxa profiling:
    --metaphlan_databases path    folder for the MetaPhlAn database
    --bt2options          value   BowTie2 options
  
  HUMANn parameters for functional profiling:
    --chocophlan          path    folder for the ChocoPhlAn database
    --uniref              path	  folder for the UniRef database

YAMP supports FASTQ and compressed FASTQ files.
"""
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

/**
	Prepare workflow introspection

	This process adds the workflow introspection (also printed at runtime) in the logs
	This is NF-CORE code.
*/


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
	if (workflow.containerEngine == 'singularity') {
	   	summary['BBmap'] = "https://depot.galaxyproject.org/singularity/bbmap:38.87--h1296035_0"
		summary['FastQC'] = "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
	} else if (workflow.containerEngine == 'docker') {
    	summary['BBmap'] = "quay.io/biocontainers/bbmap:38.87--h1296035_0"
		summary['FastQC'] = "quay.io/biocontainers/fastqc:0.11.9--0"
	} else {
		summary['BBmap'] = "No container information"
		summary['FastQC'] = "No container information"
	}
}

if (params.mode != "QC")
{
	if (workflow.containerEngine == 'singularity') {
		summary['biobakery'] = "biobakery/workflows:3.0.0.a.6.metaphlanv3.0.7"
		summary['qiime'] = "qiime2/core:2020.8"
	} else if (workflow.containerEngine == 'docker') {
		summary['biobakery'] = "biobakery/workflows:3.0.0.a.6.metaphlanv3.0.7"
		summary['qiime'] = "qiime2/core:2020.8"
	} else {
		summary['biobakery'] = "No container information"
		summary['qiime'] = "No container information"
	}
}

if (workflow.containerEngine == 'singularity') {
	summary['MultiQC'] = "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"
} else if (workflow.containerEngine == 'docker') {
	summary['MultiQC'] = "quay.io/biocontainers/multiqc:1.9--py_1"
} else {
	summary['MultiQC'] = "No container information"
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
summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'

if (params.mode != "characterisation") 
{
	summary['Performing de-duplication'] = params.dedup

	//remove_synthetic_contaminants 
	summary['Synthetic contaminants'] = ""
	summary['Artefacts'] = params.artefacts
	summary['Phix174ill'] = params.phix174ill

	//Trimming
	summary['Adapters'] = params.adapters
	summary['Trimming parameters'] = ""
	summary['Input quality offset'] = params.qin == 33 ? 'ASCII+33' : 'ASCII+64'
	summary['Min phred score'] = params.phred
	summary['Min length'] = params.minlength
	summary['kmer lenght'] = params.kcontaminants
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
    //BowTie2 databases for metaphlan
	summary['MetaPhlAn parameters'] = ""
    summary['MetaPhlAn database'] = params.metaphlan_databases
    summary['Bowtie2 options'] = params.bt2options
  
    // ChocoPhlAn and UniRef databases
	summary['HUMAnN parameters'] = ""
	summary['Chocophlan database'] = params.chocophlan
	summary['Uniref database'] = params.uniref
}

//Folders
summary['Folders'] = ""
summary['Output dir'] = workingpath
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Lunching dir'] = workflow.launchDir

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



workflow { 
	
	/**
		The channel for the input read files is used 
		- for the first QC assessment (on the raw reads)
		- for the deduplication step (which is optional and may skip to trimming)
		- for the decontamination from synthetic contaminants (used only if deduplication is not run)
		- for the characterisation, if no QC is performed beforehand
	*/

	if (params.singleEnd) {
		Channel.from([[params.prefix, [file(params.reads1)]]]).set{ read_files }
	} else {
		Channel.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] ).set{ read_files }
	}
	
	// Introspection
	get_software_versions()
	
	// ------------------------------------------------------------------------------   
	//	QUALITY CONTROL 
	// ------------------------------------------------------------------------------   

	/**
		Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

		This step is OPTIONAL. De-duplication should be carried on iff you are
	    using PCR amplification (in this case identical reads are technical artefacts)
		but not otherwise (identical reads will identify natural duplicates).
	*/

	dedup(read_files)

	/**
		Quality control - STEP 2. A decontamination of synthetic sequences and artefacts 
		is performed.
	*/

	//When the de-duplication is not done, the raw file should be pushed in the correct channel
	//FIXME: make this also optional?
	if (!params.dedup) { 
		to_synthetic_contaminants = read_files
		dedup_log = Channel.from(file("$baseDir/assets/no_dedup.yaml"))
	} else if (params.dedup) { 
		to_synthetic_contaminants = dedup.out.to_synthetic_contaminants
	} 

	// Defines channels for resources file 
	Channel.fromPath( "${params.artefacts}", checkIfExists: true ).set { artefacts }
	Channel.fromPath( "${params.phix174ill}", checkIfExists: true ).set { phix174ill }

	remove_synthetic_contaminants(artefacts.combine(phix174ill).combine(to_synthetic_contaminants))

	/**
		Quality control - STEP 3. Trimming of low quality bases and of adapter sequences. 
		Short reads are discarded. 
	
		If dealing with paired-end reads, when either forward or reverse of a paired-read
		are discarded, the surviving read is saved on a file of singleton reads.
	*/

	// Defines channels for resources file 
	Channel.fromPath( "${params.adapters}", checkIfExists: true ).set { adapters }

	trim(adapters.combine(remove_synthetic_contaminants.out.to_trim))


	/**
		Quality control - STEP 4. Decontamination. Removes external organisms' contamination, 
		using given genomes. 

		When an indexed contaminant (pan)genome is not provided, the index_foreign_genome process is run 
		before the decontamination process. This process require the FASTA file of the contaminant (pan)genome.
	*/

	// Defines channels for foreign_genome file 
	Channel.fromPath( "${params.foreign_genome}", checkIfExists: true ).set { foreign_genome }

	//Stage boilerplate log when the contaminant (pan)genome is indexed
	//When the indexed contaminant (pan)genome is already available, its path should be pushed in the correct channel
	if (params.mode != "characterisation" && params.foreign_genome_ref == "") {
		index_foreign_genome_log = Channel.from(file("$baseDir/assets/foreign_genome_indexing_mqc.yaml"))
		index_foreign_genome(foreign_genome)
		ref_foreign_genome = index_foreign_genome.out.ref_foreign_genome
	} else if (params.mode != "characterisation" && params.foreign_genome_ref != "") {
		index_foreign_genome_log = Channel.empty()
		ref_foreign_genome = Channel.fromPath(params.foreign_genome_ref, checkIfExists: true)
	} else {
		index_foreign_genome_log = Channel.empty()
		ref_foreign_genome = Channel.empty()
	}

	decontaminate(ref_foreign_genome.combine(trim.out.to_decontaminate))

	// ------------------------------------------------------------------------------   
	//	QUALITY ASSESSMENT 
	// ------------------------------------------------------------------------------   

	quality_assessment(read_files.mix(decontaminate.out.qcd_reads))

	// ------------------------------------------------------------------------------   
	//  COMMUNITY CHARACTERISATION 
	// ------------------------------------------------------------------------------   

	// The user will specify the clean file either as a single clean file (that is the YAMP
	// default behaviour), or as two files (forward/reverse). ]
	// In the former case, the user will set singleEnd = true and only one file will be 
	// selected and used directly for taxa and community profiling.
	// In the latter case, the user will set singleEnd = false and provide two files, that will
	// be merged before feeding the relevant channels for profiling.
	if (params.mode == "characterisation" && params.singleEnd) {
		//sets the paths passed in input
		reads_profile = read_files 
	
		//Initialise empty channel for log
		merge_paired_end_cleaned_log = Channel.empty()
	} else if (params.mode == "characterisation" && !params.singleEnd) {
		//merges the paths passed in input	
		merge_paired_end_cleaned(read_files)
		reads_profile = merge_paired_end_cleaned.out.to_profile_merged
	
		//Stage boilerplate log
		merge_paired_end_cleaned_log = Channel.from(file("$baseDir/assets/merge_paired_end_cleaned_mqc.yaml"))
	} else if (params.mode != "characterisation")
	{
		//Uses the decontaminated reads from previous steps
		reads_profile = decontaminate.out.qcd_reads
		
		//Initialise empty channel for log
		merge_paired_end_cleaned_log = Channel.empty()
	}

	/**
		Community Characterisation - STEP 1. Performs taxonomic binning and estimates the 
		microbial relative abundances using MetaPhlAn and its databases of clade-specific markers.
	*/

	// Defines channels for bowtie2_metaphlan_databases file 
	Channel.fromPath( params.metaphlan_databases, type: 'dir', checkIfExists: true ).set { bowtie2_metaphlan_databases }

	profile_taxa(reads_profile, bowtie2_metaphlan_databases)

	/**
		Community Characterisation - STEP 2. Performs the functional annotation using HUMAnN.
	*/

	// Defines channels for bowtie2_metaphlan_databases file 
	Channel.fromPath( params.chocophlan, type: 'dir', checkIfExists: true ).set { chocophlan_databases }
	Channel.fromPath( params.uniref, type: 'dir', checkIfExists: true ).set { uniref_databases }

	profile_function(reads_profile, profile_taxa.out.to_profile_function_bugs, chocophlan_databases, uniref_databases)

	/**
		Community Characterisation - STEP 3. Evaluates several alpha-diversity measures. 

	*/

	alpha_diversity(profile_taxa.out.to_alpha_diversity)

	// ------------------------------------------------------------------------------   
	//	MULTIQC LOGGING
	// ------------------------------------------------------------------------------   

	/**
		Generate Logs. 

		Logs generate at each analysis step are collected and processed with MultiQC 
	*/

	// Stage config files
	multiqc_config = file(params.multiqc_config)
	
	log(multiqc_config, 
		create_workflow_summary(summary), 
		get_software_versions.out.software_versions_yaml, 
		quality_assessment.out.fastqc_log.collect().ifEmpty([]), 
		dedup.out.dedup_log.ifEmpty([]), 
		remove_synthetic_contaminants.out.synthetic_contaminants_log.ifEmpty([]), 
		trim.out.trimming_log.ifEmpty([]), 
		index_foreign_genome_log.ifEmpty([]), 
		decontaminate.out.decontaminate_log.ifEmpty([]), 
		merge_paired_end_cleaned_log.ifEmpty([]), 
		profile_taxa.out.profile_taxa_log.ifEmpty([]), 
		profile_function.out.profile_functions_log.ifEmpty([]), 
		alpha_diversity.out.alpha_diversity_log.ifEmpty([])
	   )

}


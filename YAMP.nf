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
	--singleEnd     <true|false>
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
summary['Pipeline Name'] = 'YAMP'
summary['Pipeline Version'] = workflow.manifest.version

summary['Nextflow version'] = nextflow.version.toString() + " build " + nextflow.build.toString() + " (" + nextflow.timestamp + ")"

summary['Java version'] = System.getProperty("java.version") //Java Runtime Environment version
summary['Java Virtual Machine'] = System.getProperty("java.vm.name") + "(" + System.getProperty("java.vm.version") + ")"

summary['Reads'] = "[" + params.reads1 + ", " + params.reads2 + "]"
summary['Prefix'] = params.prefix
summary['Running mode'] = params.mode
summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Performing de-duplication'] = params.dedup

summary['Output dir'] = workingpath
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Lunching dir'] = workflow.launchDir

summary['Config Profile'] = workflow.profile
summary['Resumed'] = workflow.resume
		
summary['Command Line'] = workflow.commandLine	

summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container

if(workflow.profile == 'awsbatch'){
	summary['AWS Region'] = params.awsregion
	summary['AWS Queue'] = params.awsqueue
}

//FIXME: adding ending time
summary['Operating system'] = System.getProperty("os.name") + " " + System.getProperty("os.arch") + " v" +  System.getProperty("os.version")
summary['User name'] = System.getProperty("user.name") //User's account name
summary['Starting time'] = new java.util.Date() 


log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
log.info ""


/**
	Prepare workflow introspection

	This process adds the workflow introspection (also printed at runtime) in the logs
	This is NF-CORE code.
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'YAMP-summary'
    description: "This information is collected when the pipeline is started."
    section_name: 'YAMP Workflow Summary'
    section_href: 'https://github.com/alesssia/yamp'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/**
	Gets software version. 

	This process ensures that software version are included in the logs.
*/

process get_software_versions {

	output:
	path 'software_versions_mqc.yaml' into software_versions_yaml

	script:
	// TODO nf-core: Get all tools to print their version number here
	//I am using a multi-containers scenarios, supporting conda, docker, and singularity
	//with the software at a specific version (the same for all platforms). Therefore, I
	//will simply parse the version from there. 
	//Overkill, TODO: think of something better 
	"""
	echo $workflow.manifest.version > v_pipeline.txt
	echo $workflow.nextflow.version > v_nextflow.txt

	echo $params.conda_fastqc | cut -d= -f 2 > v_fastqc.txt
	echo $params.conda_bbmap | cut -d= -f 2 > v_bbmap.txt
	echo $params.conda_multiqc | cut -d= -f 2 > v_multiqc.txt

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

if(params.singleEnd) {
	Channel
	.from([params.prefix, params.reads1, params.reads2])
	.map { row -> [ row[0], [file(row[1])]] }
	.ifEmpty { exit 1, "Internal error: params.readPaths was empty - no input files supplied" }
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
} else {
	Channel
	.from([params.prefix, params.reads1, params.reads2] )
	.map { row -> [ row[0], [file(row[1]), file(row[2])]] }
	.ifEmpty { exit 1, "Internal error: params.readPaths was empty - no input files supplied" }
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
}


/**
	Generate Logs. 

	Logs generate at each analysis step are collected and processed with MultiQC 
*/

// Stage config files
multiqc_config = file(params.multiqc_config)

process multiqc {
	publishDir workingdir, mode: 'copy'

    conda (params.enable_conda ? params.conda_multiqc : null)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_multiqc
    } else {
        container params.docker_container_multiqc
    }

	input:
	file multiqc_config
	// TODO nf-core: Add in log files from your new processes for MultiQC to find!
	file workflow_summary from create_workflow_summary(summary)
	file ('software_versions/*') from software_versions_yaml
	// file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
	// file ('dedup_mqc.yaml') from dedup_log
	// file ('synthetic_contaminants_mqc.yaml') from synthetic_contaminants_log
	// file ('trimming_mqc.yaml') from trimming_log
	
	output:
	path "*multiqc_report.html" into multiqc_report
	path "*_data"

	script:
	// TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
	"""
	multiqc --config $multiqc_config .
	"""
}


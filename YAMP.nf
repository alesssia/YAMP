#!/usr/bin/env nextflow
	
/**
	Yet Another Metagenomic Pipeline (YAMP)
	Copyright (C) 2017 	Dr Alessia Visconti 	      
	      
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
	
	For any bugs or problems found, please contact us at:
	- alessia.visconti@kcl.ac.uk ; 
	- https://github.com/alesssia/YAMP/issues
*/

version='0.9.4.1'
timestamp='20180424'

/**
	Prints version when asked for
*/
if (params.version) {
	System.out.println("")
	System.out.println("YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: $version ($timestamp)")
	exit 1
}

/**
	Prints help when asked for
*/

if (params.help) {
	System.out.println("")
	System.out.println("YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: $version ($timestamp)")
	System.out.println("This pipeline is distributed in the hope that it will be useful")
	System.out.println("but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.")
	System.out.println("")
	System.out.println("Please report comments and bugs to alessia.visconti@kcl.ac.uk")
	System.out.println("or at https://github.com/alesssia/YAMP/issues.")
	System.out.println("Check https://github.com/alesssia/YAMP for updates, and refer to")
	System.out.println("https://github.com/alesssia/YAMP/wiki for more details.")
	System.out.println("")
	System.out.println("Usage: ")
	System.out.println("   nextflow run YAMP.nf --reads1 R1 --reads2 R2 --prefix mysample --outdir path --mode MODE  ")
	System.out.println("                [options] [-with-docker|-with-singularity]")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --reads1   R1      Forward (if paired-end) OR all reads (if single-end) file path")
	System.out.println("    [--reads2] R2      Reverse reads file path (only if paired-end library layout)")
	System.out.println("    --prefix   prefix  Prefix used to name the result files")
	System.out.println("    --outdir   path    Output directory (will be outdir/prefix/)")
	System.out.println("    --mode     <QC|characterisation|complete>")
	System.out.println("Options:")
	System.out.println("    --librarylayout <single|paired>")
	System.out.println("    --dedup         <true|false>   whether to perform de-duplication")
	System.out.println("    --keepQCtmpfile <true|false>   whether to save QC temporary files")
	System.out.println("    --keepCCtmpfile <true|false>   whether to save community characterisation temporary files")
	System.out.println("Please refer to nextflow.config for more options.")
	System.out.println("")
	System.out.println("Container:")
	System.out.println("    Docker image to use with -with-docker|-with-singularity options is")
	System.out.println("    'docker://alesssia/yampdocker'")
	System.out.println("")
	System.out.println("YAMP supports FASTQ and compressed FASTQ files.")
	System.out.println("")
    exit 1
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

if (params.librarylayout != "paired" && params.librarylayout != "single") { 
	exit 1, "Library layout not available. Choose any of <single, paired>" 
}   

if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

//--reads2 can be omitted when the library layout is "single" (indeed it specifies single-end
//sequencing)
if (params.mode != "characterisation" && params.librarylayout != "single" && (params.reads2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
}

//--reads1 and --reads2 can be omitted (and the default from the config file used instead) 
//only when mode is "characterisation". Obviously, --reads2 should be always omitted when the
//library layout is single.
if (params.mode != "characterisation" && ( (params.librarylayout == "paired" && (params.reads1 == "null" || params.reads2 == "null")) ||			
 							 			   params.librarylayout == "single" && params.reads1 == "null") ) {
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

//Creates main log file
mylog = file(params.outdir + "/" + params.prefix + "/" + params.prefix + ".log")

//Logs headers
mylog <<  """---------------------------------------------
YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: $version ($timestamp)
---------------------------------------------
	
Copyright (C) 2017 Dr Alessia Visconti <alessia.visconti@kcl.ac.uk>

This pipeline is distributed in the hope that it will be useful
but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.

Please report comments and bugs to alessia.visconti@kcl.ac.uk
or at https://github.com/alesssia/YAMP/issues.						  
Check https://github.com/alesssia/YAMP for updates, and refer to
https://github.com/alesssia/YAMP/wiki for more details.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   						   	   
"""
	   
//Fetches information on OS and java versions, including user name
osname = System.getProperty("os.name") //Operating system name
osarch = System.getProperty("os.arch") //Operating system architecture
osversion = System.getProperty("os.version") //Operating system version
osuser = System.getProperty("user.name") //User's account name

javaversion = System.getProperty("java.version") //Java Runtime Environment version
javaVMname = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
javaVMVersion = System.getProperty("java.vm.version") //Java Virtual Machine implementation version

//Gets starting time		
sysdate = new java.util.Date() 
		
//Logs starting time and other information about the run		
mylog << """ 
Analysis starting at $sysdate by user: $osuser
Analysed sample(s): $params.reads1 and $params.reads2
Results will be saved at $workingdir
New files will be saved using the '$params.prefix' prefix

Analysis mode? $params.mode
Library layout? $params.librarylayout
Saving QC temporary files? $params.keepQCtmpfile
Saving community characterisation temporary files? $params.keepCCtmpfile
Performing de-duplication? $params.dedup	

------------
	
Analysis introspection:

Operating System:
	name:         $osname
	architecture: $osarch
	version:      $osversion

Java
	version: $javaversion
	Java Virtual Machine: $javaVMname ; version: $javaVMVersion

nextflow:
	version:   $nextflow.version	
	build:     $nextflow.build	
	timestamp: $nextflow.timestamp	
			
Container:
	Docker image: $workflow.container		

Repository:
	url:            $workflow.repository
	Git commit ID:  $workflow.commitId
	Git branch/tag: $workflow.revision

Analysis environment:

	projectDir: $workflow.projectDir	
	launchDir:  $workflow.launchDir
	workingDir: $workflow.workDir	
		
	command line: $workflow.commandLine	
	
	Run name:   $workflow.runName	
	Session ID: $workflow.sessionId	
	profile:    $workflow.profile
			
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
	   
""" 

/**
	Quality Assessment - STEP 1. Assessment of read quality of FASTQ file, 
	done by means of FastQC. Multiple  plots are generated to show average 
	phred quality scores and other metrics.

	This step will generate (for each end, if the layout is paired) an HTML page,
	showing a summary of the results and a set of plots offering a visual guidance
	and for assessing the quality of the sample, and a zip file, that includes the 
	HTML page, all the images and a text report. The script will save the HTML page 
	and the text report and delete the archive.
	Several information on multiple QC parameters are logged.
*/


//Defines channel with <step, readfile, label, stem> as input for fastQC script:
//- step defines which step of the analysis is performed (QA of raw reads is the 1st,
//  QA of the trimmed reads is the 4th, QA of the decontaminated reads is the 6th), and
//  it will be used to sort the step's logs into the final YAMP log file.
//- readfile is the reads FASTQ file, 
//- label indicated whethere it is forward (R1) or reverse (R2) strand -- it is empty 
//  for single-end reads.
//- stem indicates which reads are QA'd (raw/trimmed/decontaminated reads), and it is
//  used to name the report files.
//These parameters are used also for the trimmed and decontaminated reads channel

//When layout is single, the params reads2 is not used.
//This channel will be processed by the process Quality Assessment, that is reported
//at the end of the QC section of this script (that runs all the QC assessed tasks)
if (params.librarylayout == "paired") {
	rawreads = Channel.from( ['1', file(params.reads1), '_R1', '_rawreads'], ['1', file(params.reads2), '_R2', '_rawreads'] )
}
else {
	rawreads = Channel.value( ['1', file(params.reads1), '', '_rawreads'] )
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
//	QUALITY CONTROL (QC)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

/**
	Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

	If the layout is "paired", two FASTQ files are outputted, one for each paired-end.
	If "single", a single FASTQ file will be generated.
	This step is OPTIONAL. De-duplication should be carried on iff you are
    using PCR amplification (in this case identical reads are technical artefacts)
	but not otherwise (identical reads will identify natural duplicates).
*/

// Defines channel with <readfile1, readfile2> as input for de-duplicates
// When layout is single, the params.reads2 is not used, so nevermind its value
if (params.mode == "QC" && params.librarylayout == "paired") {
	todedup = Channel.value( [file(params.reads1), file(params.reads2)] )
} 
else if (params.mode == "QC") {
	todedup = Channel.value( [file(params.reads1), "null"] )
}
else {
	todedup = Channel.value( ["null1", "null2"] )
}

process dedup {
	
	input:
	set file(in1), file(in2) from todedup

	output:
	file  ".log.2" into log2
	file("${params.prefix}_dedupe*.fq") into totrim
	file("${params.prefix}_dedupe*.fq") into topublishdedupe

	when:
	(params.mode == "QC" || params.mode == "complete") && params.dedup

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 1 [De-duplication] at \$sysdate\" > .log.2
	echo \" \" >> .log.2
	
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo \"$task.memory\" | sed 's/ //g' | sed 's/B//g')
	
	#Defines command for de-duplication
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		CMD=\"clumpify.sh -Xmx\"\$maxmem\" in1=$in1 in2=$in2 out1=${params.prefix}_dedupe_R1.fq out2=${params.prefix}_dedupe_R2.fq qin=$params.qin dedupe subs=0 threads=${task.cpus}\"
	else
		CMD=\"clumpify.sh -Xmx\"\$maxmem\" in=$in1 out=${params.prefix}_dedupe.fq qin=$params.qin dedupe subs=0 threads=${task.cpus}\"
	fi
					
	#Logs version of the software and executed command (BBmap prints on stderr)
	version=\$(clumpify.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using clumpify.sh in \$version \" >> .log.2
	echo \"Executing command: \$CMD \" >> .log.2
	echo \" \" >> .log.2
	
	#De-duplicates
	exec \$CMD 2>&1 | tee tmp.log

	#Logs some figures about sequences passing de-duplication
	echo  \"Clumpify's de-duplication stats: \" >> .log.2
	echo \" \" >> .log.2
	sed -n '/Reads In:/,/Duplicates Found:/p' tmp.log >> .log.2
	echo \" \" >> .log.2
	totR=\$(grep \"Reads In:\" tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
	remR=\$(grep \"Duplicates Found:\" tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
	survivedR=\$((\$totR-\$remR))
	percentage=\$(echo \$survivedR \$totR | awk '{print \$1/\$2*100}' )
	echo \"\$survivedR out of \$totR paired reads survived de-duplication (\$percentage%, \$remR reads removed)\" >> .log.2
	echo \" \" >> .log.2

	#Measures and logs execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"STEP 1 (Quality control) terminated at \$sysdate (\$exectime seconds)\" >> .log.2
	echo \" \" >> .log.2
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.2
	echo \" \" >> .log.2
	"""
}

/**
	Quality control - STEP 2. Trimming of low quality bases and of adapter sequences. Short reads
	are discarded. A decontamination of synthetic sequences is also pefoermed.
	If dealing with paired-end reads, when either forward or reverse of a paired-read
	are discarded, the surviving read is saved on a file of singleton reads.

	If layout is "paired", three compressed FASTQ file (forward/reverse paired-end 
	and singleton reads) are outputed, if "single", only one compressed FASTQ file
	is returned. Several information are logged.
*/


//When the de-suplication is not done, the raw file should be pushed in the corret channel
if (!params.dedup) {
	if (params.librarylayout == "paired") {
		totrim = Channel.value( [file(params.reads1), file(params.reads2)] )
	} 
	else {
		totrim = Channel.value( [file(params.reads1), "null"] )
	}
}

//When single-end reads are used, the input tuple (singleton) will not match input set 
//cardinality declared by 'trim' process (pair), so I push two mock files in the channel,
//and then I take only the first two files
mocktrim = Channel.from("null")
process trim {

	input:
   	set file(reads1), file(reads2) from totrim.concat(mocktrim).flatMap().take(2).buffer(size : 2)
	file(adapters) from Channel.from( file(params.adapters) )
	file(artifacts) from Channel.from( file(params.artifacts) )
	file(phix174ill) from Channel.from( file(params.phix174ill) )

	output:
	file  ".log.3" into log3
	file("${params.prefix}_trimmed*.fq") into trimmedreads
	file("${params.prefix}_trimmed*.fq") into todecontaminate
	file("${params.prefix}_trimmed*.fq") into topublishtrim 	

	when:
	params.mode == "QC" || params.mode == "complete"

   	script:
	"""	
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 2 [Trimming] at \$sysdate\" > .log.3
	echo \" \" >> .log.3

	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	#Defines command for trimming of adapters and low quality bases
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 in2=$reads2 out=${params.prefix}_trimmed_R1_tmp.fq out2=${params.prefix}_trimmed_R2_tmp.fq outs=${params.prefix}_trimmed_singletons_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe ow\"
	else
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 out=${params.prefix}_trimmed_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe ow\"
	fi
	
	#Logs version of the software and executed command (BBMap prints on stderr)
	version=\$(bbduk.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using bbduk.sh in \$version \" >> .log.3
	echo \"Using adapters in $params.adapters \" >> .log.3
	echo \"Using synthetic contaminants in $params.phix174ill and in $params.artifacts \" >> .log.3
	echo \" \" >> .log.3
	echo \"Executing command: \$CMD \" >> .log.3
	echo \" \" >> .log.3
	
	#Trims adapters and low quality bases	
	exec \$CMD 2>&1 | tee tmp.log
	
	#Logs some figures about sequences passing trimming
	echo  \"BBduk's trimming stats (trimming adapters and low quality reads): \" >> .log.3
	sed -n '/Input:/,/Result:/p' tmp.log >> .log.3
	echo \" \" >> .log.3			
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		unpairedR=\$(wc -l ${params.prefix}_trimmed_singletons_tmp.fq | cut -d\" \" -f 1)
		unpairedR=\$((\$unpairedR/4))
		echo  \"\$unpairedR singleton reads whose mate was trimmed shorter preserved\" >> .log.3
		echo \" \" >> .log.3
	fi

	#Defines command for removing synthetic contaminants
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_R1_tmp.fq in2=${params.prefix}_trimmed_R2_tmp.fq out=${params.prefix}_trimmed_R1.fq out2=${params.prefix}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=${task.cpus} ow\"
	else
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_tmp.fq out=${params.prefix}_trimmed.fq k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=${task.cpus} ow\"
	fi

	#Logs executed command
	echo \"Executing command: \$CMD \" >> .log.3
	echo \" \" >> .log.3
	
	#Removes synthetic contaminants
	exec \$CMD 2>&1 | tee tmp.log

	#Logs some figures about sequences passing deletion of contaminants
	echo  \"BBduk's trimming stats (synthetic contaminants): \" >> .log.3
	sed -n '/Input:/,/Result:/p' tmp.log >> .log.3
	echo \" \" >> .log.3

	#Removes synthetic contaminants and logs some figures (singleton read file, 
	#that exists iif the library layout was 'paired')
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_singletons_tmp.fq out=${params.prefix}_trimmed_singletons.fq k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=${task.cpus} ow\"
		
		echo \"Executing command: \$CMD \" >> .log.3
		echo \" \" >> .log.3
	
		#Removes synthetic contaminants
		exec \$CMD 2>&1 | tee tmp.log		
							
		#Logs some figures about sequences passing deletion of contaminants
		echo  \"BBduk's trimming stats (synthetic contaminants, singleton reads): \" >> .log.3
		sed -n '/Input:/,/Result:/p' tmp.log >> .log.3
		echo \" \" >> .log.3
	fi
	
	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${params.prefix}_trimmed*_tmp.fq 

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"STEP 2 (Quality Control) terminated at \$sysdate (\$exectime seconds)\" >> .log.3
	echo \" \" >> .log.3
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.3
	echo \"\" >> .log.3
	"""
}


/**
	Quality assessed - STEP 2. Quality control after trimming.

	An output similar to that generated by STEP 1 is produced.
*/

//This comment is a place holder. The execution of this step is delegated to the process
//Quality Assessment, that is reported at the end of the QC section of this script and
//that runs all the QC assessed tasks.

/**
	Quality control - STEP 3. Decontamination. Removes external organisms' contamination, 
	using a previously created index. When paired-end are used, decontamination is 
	carried on idependently on paired reads and on singleton reads thanks to BBwrap, 
	that calls BBmap once on the paired reads and once on the singleton ones, merging
	 the results on a single output file.

	Two files are outputted: the FASTQ of the decontaminated reads (including both
	paired-reads and singletons) and that of the contaminating reads (that can be
	used for refinements/checks).
	Please note that if keepQCtmpfile is set to false, the file of the contaminating 
	reads is discarded

	TODO: use BBsplit for multiple organisms decontamination, or fix the ref to a 
	FASTA file pangenome
*/

//When single-end reads are used, the input tuple (singleton) will not match input set 
//cardinality declared by 'trim' process (triplet), so I push two mock files in the channel,
//and then I take only the first three files.
mockdecontaminate = Channel.from("null", "null")
process decontaminate {
	
	publishDir  workingdir, mode: 'move', pattern: "*_clean.fq.gz"
		
	input:
	set file(infile1), file(infile2), file(infile12) from todecontaminate.concat(mockdecontaminate).flatMap().take(3).buffer(size : 3)
	file(refForeingGenome) from Channel.from( file(params.refForeingGenome, type: 'dir') )
	
	output:
	file "*_clean.fq.gz"
	file  ".log.5" into log5
	file "${params.prefix}_clean.fq" into decontaminatedreads
	file "${params.prefix}_clean.fq" into toprofiletaxa
	file "${params.prefix}_clean.fq" into toprofilefunctionreads
	file "${params.prefix}_cont.fq" into topublishdecontaminate

	when:
	params.mode == "QC" || params.mode == "complete"
	
	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 3 [Decontamination] at \$sysdate\" > .log.5
	echo \" \" >> .log.5

	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	#Defines command for decontamination
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		CMD=\"bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=$infile1,$infile12 in2=$infile2,null outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred path=$refForeingGenome qin=$params.qin threads=${task.cpus} untrim quickmatch fast ow\"
	else
		CMD=\"bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=$infile1 outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred path=$refForeingGenome qin=$params.qin threads=${task.cpus} untrim quickmatch fast ow\"
	fi
	
	#Logs version of the software and executed command (BBmap prints on stderr)
	version=\$(bbwrap.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using bbwrap.sh in \$version \" >> .log.5
	echo \"Using contaminant (pan)genome indexed in $params.refForeingGenome \" >> .log.5
	echo \" \" >> .log.5
	echo \"Executing command: \$CMD \" >> .log.5
	echo \" \" >> .log.5
	
	#Decontaminates
	exec \$CMD 2>&1 | tee tmp.log
	
	#Logs some figures about decontaminated/contaminated reads
	echo  \"BBwrap's human decontamination stats (paired reads): \" >> .log.5
	sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | head -17 >> .log.5
	echo \" \" >> .log.5
	sed -n '/Read 2 data:/,/N Rate:/p' tmp.log >> .log.5
	echo \" \" >> .log.5
	
	if [ \"$params.librarylayout\" = \"paired\" ]; then
		echo  \"BBmap's human decontamination stats (singletons reads): \" >> .log.5
		sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | tail -17 >> .log.5
		echo \" \" >> .log.5
	fi

	gzip -c ${params.prefix}_clean.fq > ${params.prefix}_clean.fq.gz

	nClean=\$(wc -l ${params.prefix}_clean.fq | cut -d\" \" -f 1)
	nClean=\$((\$nClean/4))
	nCont=\$(wc -l ${params.prefix}_cont.fq | cut -d\" \" -f 1)
	nCont=\$((\$nCont/4))
	echo \"\$nClean reads survived decontamination (\$nCont reads removed)\" >> .log.5
	echo \" \" >> .log.5

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"STEP 3 (Quality Control) terminated at \$sysdate (\$exectime seconds)\" >> .log.5
	echo \" \" >> .log.5
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.5
	echo \"\" >> .log.5
	"""
}



/**
	Quality Assessment - STEP 3. Quality control after decontamination. The clean FASTQ 
	file producted by the decontamination step is assessed with FastQC.

	An output similar to that generated by STEP 1 is produced
*/

//This comment is a place holder. The execution of this step is delegated to the process
//Quality Assessment, that is reported below and that runs all the QC assessed tasks.


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
//	QUALITY ASSESSMENT AND VISUALISATION
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

//Creates the correct objects for the quality assessment, by merging the files derived from
//trimming and decontamination and the step number, label and step.
if (params.librarylayout == "paired") {
	trimmedreads2qc = Channel.from('4').combine(trimmedreads.flatMap().merge( Channel.from( ['_R1', '_R2'] ) ){ a, b -> [a, b] }).combine(Channel.from('_trimmedreads'))
} 
else {
	trimmedreads2qc = Channel.from('4').combine(trimmedreads.flatMap()).combine( Channel.from( '' ) ).combine(Channel.from('_trimmedreads'))
}
decontaminatedreads2qc = Channel.from('6').combine(decontaminatedreads).combine( Channel.from( '' ) ).combine(Channel.from('_decontaminatedreads'))

//Creates the channel which performs the QC
toQC = rawreads.mix(trimmedreads2qc, decontaminatedreads2qc) 

//Process performing all the Quality Assessment
process qualityAssessment {
	
	publishDir  workingdir, mode: 'move', pattern: "*.{html,txt}"
	  	
	input:
   	set val(step), file(reads), val(label), val(stem) from toQC

	output:
	file ".log.$step$label" into logQC
	file "${params.prefix}*_fastqc.html" 
	file "${params.prefix}*_fastqc_data.txt" 

	when:
	params.mode == "QC" || params.mode == "complete"

   	script:
	"""	
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. [Assessment of read quality] at \$sysdate\" > .log.$step$label
	echo \"File being analysed: $reads\" >> .log.$step$label
	echo \" \" >> .log.$step$label
	
	#Logs version of the software and executed command
	version=\$(fastqc --version) 
	CMD=\"fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} $reads\"
	
	echo \"Using \$version \" >> .log.$step$label
	echo \"Executing command \$CMD \" >> .log.$step$label
	echo \" \" >> .log.$step$label
	
	#Does QC, extracts relevant information, and removes temporary files
	bash fastQC.sh $reads ${params.prefix}${stem}${label} ${task.cpus} $reads
	
	#Logging QC statistics (number of sequences, Pass/warning/fail, basic statistics, duplication level, kmers)
	base=\$(basename $reads)
	bash logQC.sh \$base ${params.prefix}${stem}${label}_fastqc_data.txt .log.$step$label
				
	#Measures and log execution time			
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"Quality assessment on $reads terminated at \$sysdate (\$exectime seconds)\" >> .log.$step$label
	echo \" \" >> .log.$step$label	
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.$step$label
	echo \" \" >> .log.$step$label	
	"""	
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
//  COMMUNITY CHARACTERISATION 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

//When doing community characterisation from previously QC'd files, this file should be 
//pushed in the corret channels.
//Please note that the name of the file externally QC'ed should be the same of the 
//one geneared by the YAMP, that is prefix_clean.fq, and should include ALL the reads.
if (params.mode == "characterisation") {
	toprofiletaxa =  Channel.from( file("$workingdir/${params.prefix}_clean.fq") )
	toprofilefunctionreads = Channel.from( file("$workingdir/${params.prefix}_clean.fq") )
}


/**
	Community Characterisation - STEP 1. Performs taxonomic binning and estimates the 
	microbial relative abundancies. MetaPhlAn2 and its databases of clade-specific markers
	are used to infers the presence and relative abundance of the organisms (at the specie/ 
	strain level) that are present in the sample and to estimate their relative abundance.

	Two files are outputted: a tab-separated file reporting the species detected and their
	relative abundance, and a BIOM file, that will be used to evaluate alpha (STEP 8) and
	beta diversity.
*/


process profileTaxa {

	publishDir  workingdir, mode: 'copy', pattern: "*.{biom,tsv}"
	
	input:
	file(infile) from toprofiletaxa
	file(mpa_pkl) from Channel.from( file(params.mpa_pkl) )
	file(bowtie2db) from Channel.fromPath( params.bowtie2db, type: 'dir' )

    output:
	file ".log.7" into log7
	file "${params.prefix}.biom" into toalphadiversity
	file "${params.prefix}_metaphlan_bugs_list.tsv" into toprofilefunctionbugs
	file "${params.prefix}_bt2out.txt" into topublishprofiletaxa

	when:
	params.mode == "characterisation" || params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Community Characterisation. STEP 1 [Taxonomic binning and profiling] at \$sysdate\" > .log.7
	echo \" \" >> .log.7
	
	#If a file with the same name is already present, Metaphlan2 will crash
	rm -rf ${params.prefix}_bt2out.txt
	
	#Defines command for estimating abundances
	CMD=\"metaphlan2.py --input_type fastq --tmp_dir=. --biom ${params.prefix}.biom --bowtie2out=${params.prefix}_bt2out.txt --mpa_pkl $mpa_pkl  --bowtie2db $bowtie2db/$params.bowtie2dbfiles --bt2_ps $params.bt2options --nproc ${task.cpus} $infile ${params.prefix}_metaphlan_bugs_list.tsv\"

	#Logs version of the software and executed command 
	#MetaPhlAn prints on stderr
	version=\$(metaphlan2.py --version 2>&1 >/dev/null | grep \"MetaPhlAn\")
	echo \"Using \$version \" >> .log.7
	echo \"Using BowTie2 database in $params.bowtie2db \" >> .log.7
	echo \" \" >> .log.7
	echo \"Executing command: \$CMD \" >> .log.7
	
	echo \" \" >> .log.7

	#Estimates microbial abundances
	exec \$CMD 2>&1 | tee tmp.log

	#Sets the prefix in the biom file
	sed -i 's/Metaphlan2_Analysis/${params.prefix}/g' ${params.prefix}.biom
	sed -i 's/Metaphlan2_Analysis/${params.prefix}/g' ${params.prefix}_metaphlan_bugs_list.tsv

	#Logs some info
	tree=(kingdom phylum class order family genus species)
	for i in {2..7}
	do
		c=\$(sed '1d' ${params.prefix}_metaphlan_bugs_list.tsv | cut -d\"|\" -f \$i | grep -v \"k__\" | cut -f 1  | sort | uniq | sed '/^\\s*\$/d' | wc -l | cut -d\" \" -f 1)
		echo \"\$c \${tree[((\$i-1))]} found\" >> .log.7
	done

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"\" >> .log.7
	echo \"STEP 1 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> .log.7
	echo \" \" >> .log.7
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.7
	echo \"\" >> .log.7
	"""
}

/**
	Community Characterisation - STEP 2. Evaluates alpha-diversity, that is the 
	mean species diversity the given sample. Please note that the alpha diversity 
	is the only per-sample measure, so it is the only one evaluated by this module. 
 	
	If a newick tree is provided as input (see QIIME documentation for details), a 
	further and more reliable phylogenetic measure is evaluated (i.e., PD_whole_tree).

	One text file listing the alpha-diversity values, evaluated by means of
	multiple measure, is outputted.
*/

process alphaDiversity {

	publishDir  workingdir, mode: 'move', pattern: "*.{tsv}"
	
	input:
	file(infile) from toalphadiversity
	file(treepath) from Channel.from( file(params.treepath) )
	
    output:
	file ".log.8" into log8
	file "${params.prefix}_alpha_diversity.tsv"
	
	when:
	params.mode == "characterisation" || params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Community Characterisation. STEP 2 [Evaluating alpha-diversity] at \$sysdate\" > .log.8
	echo \" \" >> .log.8

	#It checks if the profiling was successful, that is if identifies at least three species
	n=\$(grep -o s__ $infile | wc -l  | cut -d\" \" -f 1)
	if (( n > 3 ))
	then
		#Defines command -- if the tree path is not specified, not all the alpha 
		#measures can be evaluated (that is, PD_whole_tree is skipped)
		if [ $params.treepath == null ]
		then
			CMD=\"alpha_diversity.py -i $infile -o ${params.prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong\"
		else
			CMD=\"alpha_diversity.py -i $infile -o ${params.prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong,PD_whole_tree -t $treepath\"
		fi
		
		#Logs version of the software and executed command
		version=\$(alpha_diversity.py --version) 
		echo \"Using \$version \" >> .log.8
		if [ $params.treepath == null ]
		then
			echo \"Newick tree not used, PD_whole_tree skipped\" >> .log.8
		else
			echo \"Using Newick tree in $params.treepath\" >> .log.8
		fi
		echo \" \" >> .log.8
		echo \"Executing command: \$CMD \" >> .log.8
		echo \" \" >> .log.8
		
		#Evaluates alpha diversities, redirect is done here because QIIME gets it as an extra parameter
		exec \$CMD 2>&1 | tee tmp.log
	else
		#Also if the alpha are not evaluated the file should be created in order to be returned
		echo \"Not enough classified species detected (N=\$n). Analysis skipped.\" >> .log.8
		touch ${params.prefix}_alpha_diversity.tsv 
	fi
	
	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"\" >> .log.8
	echo \"STEP 2 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> .log.8
	echo \" \" >> .log.8
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.8
	echo \"\" >> .log.8
	"""
}


/**
	Community Characterisation - STEP 3. Performs the functional annotation using HUMAnN2.
	HUMAnN2 will bypasses the taxomonic profiling step (since it has already
	been performed) and uses the list of specied detected on step 7.
	While the aligners are forced to be Bowtie2 and DIAMOND, the user can
	select the UniRef database to use (UniRef50, UniRef90).

	It outputs several files, and some of them will be removed if keepCCtmpfile is
	set to false. Namely, it creates:
	- three tab-separated files representing the gene families, and the pathways'
	  coverahe and abundancies
	- a SAM file representing the full alignment from Bowtie2 (saved only if keepCCtmpfile
	  is set to true)
    - two tab-separated file representing the reduced aligned reads from both
	  Bowtie2 and DIAMOND (saved only if keepCCtmpfile is set to true)
	- two FASTA file, representing the unaligned reads from both Bowtie2 and
	  DIAMOND (saved only if keepCCtmpfile is set to true)
	- a log of the execution
*/

process profileFunction {

	publishDir  workingdir, mode: 'copy', pattern: "*.{tsv,log}"
	
	input:
	file(cleanreads) from toprofilefunctionreads
	file(metaphlanbuglist) from toprofilefunctionbugs
	file(chocophlan) from Channel.fromPath( params.chocophlan, type: 'dir' )
	file(uniref) from Channel.fromPath( params.uniref, type: 'dir' )
	
    output:
	file ".log.9" into log9
	file "${params.prefix}_HUMAnN2.log"
	file "${params.prefix}_genefamilies.tsv"
	file "${params.prefix}_pathcoverage.tsv"
	file "${params.prefix}_pathabundance.tsv"
	
	//Those may or may be not kept, according to the value of the keepCCtmpfile parameter
	set ("${params.prefix}_bowtie2_aligned.sam", "${params.prefix}_bowtie2_aligned.tsv", "${params.prefix}_diamond_aligned.tsv", "${params.prefix}_bowtie2_unaligned.fa", "${params.prefix}_diamond_unaligned.fa") into topublishhumann2	

	when:
	params.mode == "characterisation" || params.mode == "complete"

	script:
	"""
	#Measures execution time
 	sysdate=\$(date)
 	starttime=\$(date +%s.%N)
 	echo \"Performing Community Characterisation. STEP 3 [Performing functional annotation] with HUMAnN2 at \$sysdate\" > .log.9
 	echo \" \" >> .log.9

	#Defines HUMAnN2 command taking advantages of the MetaPhlAn2's results
	CMD=\"humann2 --input $cleanreads --output . --output-basename ${params.prefix} --taxonomic-profile $metaphlanbuglist --nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads ${task.cpus} --memory-use minimum\"
	
	#Logs version of the software and executed command
	#HUMAnN2 prints on stderr
	version=\$(humann2 --version 2>&1 >/dev/null | grep \"humann2\") 
	echo \"Using \$version \" >> .log.9
	echo \"Using ChocoPhlAn database in $params.chocophlan \" >> .log.9
	echo \"Using UniRef database in $params.uniref \" >> .log.9
	echo \" \" >> .log.9
	echo \"Executing command: \$CMD > ${params.prefix}_HUMAnN2.log\" >> .log.9
	echo \" \" >> .log.9
	
	#Performs functional annotation, redirect is done here because HUMAnN2 freaks out
	#This is  also reported in the log.
	exec \$CMD 2>&1 | tee ${params.prefix}_HUMAnN2.log 

	#If `|| true` is not add, nextflow stops... WTF 
	grep \"Total species selected from prescreen:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Selected species explain\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Unaligned reads after nucleotide alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Total gene families after translated alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Unaligned reads after translated alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	echo \"More information on HUMAnN2 run are available in the ${params.prefix}_HUMAnN2.log file\" >> .log.9 

	#Some of temporary files (if they exist) may be moved in the working directory, 
	#according to the keepCCtmpfile parameter. Others (such as the bowties2 indexes), 
	#are always removed. Those that should be moved, but have not been created by 
	#HUMAnN2, are now created by the script (they are needed as output for the channel)
	files=(${params.prefix}_bowtie2_aligned.sam ${params.prefix}_bowtie2_aligned.tsv ${params.prefix}_diamond_aligned.tsv ${params.prefix}_bowtie2_unaligned.fa ${params.prefix}_diamond_unaligned.fa)
	for i in {1..5}
	do
		if [ -f ${params.prefix}_humann2_temp/\${files[((\$i-1))]} ]
		then
			mv ${params.prefix}_humann2_temp/\${files[((\$i-1))]} .
		else
			touch \${files[((\$i-1))]}
		fi
	done
	rm -rf ${params.prefix}_humann2_temp/

 	#Measures and log execution time
 	endtime=\$(date +%s.%N)
 	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
 	sysdate=\$(date)
 	echo \"\" >> .log.9
 	echo \"STEP 3 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> .log.9
 	echo \" \" >> .log.9
 	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.9
 	echo \"\" >> .log.9
 	"""
}


/**
	CLEANUP 1. Collapses all the logs resulting from the QC in the main one, 
	and removes them.
*/

process logQC {

	input:
	file(tolog)  from logQC.flatMap().mix(log2, log3, log5).toSortedList( { a, b -> a.name <=> b.name } )

	when:
	params.mode == "QC" || params.mode == "complete"

	script:
	"""
	cat $tolog >> $mylog
	"""
}

/**
	CLEANUP 2. Saves the temporary files generate during QC (if the users requested so)
*/
	
	
process saveQCtmpfile {

	publishDir  workingdir, mode: 'copy'
		
	input:
	file (tmpfile) from topublishdedupe.mix(topublishtrim, topublishdecontaminate).flatMap()

	output:
	file "*.fq.gz"

	when:
	(params.mode == "QC" || params.mode == "complete") && params.keepQCtmpfile
		
	script:
	"""
	gzip --force -c $tmpfile > ${tmpfile}.gz
	"""
}

/**
	CLEANUP 3. Collapses all the logs resulting from the community characterisation steps
	in the main one, and removes them.
*/

process logCC {

	input:
	file(tolog) from log7.mix(log8, log9).flatMap().toSortedList( { a, b -> a.name <=> b.name } )
	
	when:
	params.mode == "characterisation" || params.mode == "complete"
		
	script:
	"""
	cat $tolog >> $mylog
	"""
}

/**
	CLEANUP 4. Saves the temporary files generate during the community characterisation 
	(if the users requested so)
*/
	
	
process saveCCtmpfile {

	publishDir  workingdir, mode: 'copy'
		
	input:
	file (tmpfile) from topublishprofiletaxa.mix(topublishhumann2).flatMap()

	output:
	file "$tmpfile"

	when:
	(params.mode == "characterisation" || params.mode == "complete") && params.keepCCtmpfile
		
	script:
	"""
	echo $tmpfile
	"""
}

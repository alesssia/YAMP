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
	along with this script.  If not, see <http://www.gnu.org/licenses/>.
	
	For any bugs or problems found, please contact us at
	alessia.visconti@kcl.ac.uk
*/

	
/**
	STEP 0. 
	
	Checks input parameters (mode only) and creates the so-called working directory
	were the results are stored (if it does not exists) as well as the log file.
	
	The working directory is named after the prefix and located in the outdir 
	folder. The log file, that will save summary statistics, execution time,
	and warnings generated during the pipeline execution, will be saved in the 
	working directory as "prefix.log".
*/
	
if (params.mode != "QC" && params.mode != "characterisation" && params.mode != "complete") {
	exit 1, "Mode not available. Chose any of <QC, characterisation, complete>"
}

//--reads1 and --reads2 can be omitted (and the default read from the config file) only
//when mode is "characterisation"
if (params.mode != "characterisation" && (params.reads1 == "/dev/null" || params.reads2 == "/dev/null") ) {
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
YET ANOTHER METAGENOMIC PIPELINE (YAMP) v 0.9.2
---------------------------------------------
	
Copyright (C) 2017 Dr Alessia Visconti   <alessia.visconti@kcl.ac.uk>

This pipeline is distributed in the hope that it will be useful
but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.

Please report comments and bugs to: alessia.visconti@kcl.ac.uk						  

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   						   	   
"""
	   
//Logs starting time and other information about the run
sysdate = new java.util.Date() 
mylog << """ 
Analysis starting at $sysdate 
Analysed samples are: $params.reads1 and $params.reads2
Working directory set to $workingdir
New files will be saved using the '$params.prefix' prefix

Analysis mode? $params.mode
Saving QC temporary files (if performed)? $params.keepQCtmpfile
Saving community characterisation temporary files (if performed)? $params.keepCCtmpfile
	     
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
	   
""" 
	   		

/**
	STEP 1. Assessment of read quality of FASTQ file, done by FastQC. Multiple 
	plots are generated to show average phred quality scores and other metrics.

	This step will generate, for each paired-end, an HTML page, showing a summary 
	of the results and a set of plots offering a visual guidanceand for assessing 
	the quality of the sample, and a zip file, that includes the HTML page, all 
	the images and a text report. The script will save the HTML page and the text 
	report and delete the archive.
	Several information on multiple QC parameters are logged.
*/

//Defines channel with <readfile, label> as input for fastQC script
rawreads = Channel.from( [file(params.reads1), '1'], [file(params.reads2), '2'] )

//This step is not a prerequisite for any other, so it does not redirect any of its output to channels
process qualityAssessmentRaw {
	
	publishDir  workingdir, mode: 'copy', pattern: "*.{html,txt}"
	  	
	input:
   	set file(reads), val(label) from rawreads

	output:
	file ".log.1$label" into log1
	file "${params.prefix}_R${label}_fastqc.html" 
	file "${params.prefix}_R${label}_fastqc_data.txt" 

	when:
	  params.mode == "QC" || params.mode == "complete"

   	script:
	"""		
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing STEP 1 [Assessment of read quality of FASTQ file] at \$sysdate on $reads\" > .log.1$label
	echo \" \" >> .log.1$label
	
	#Does QC, extracts relevant information, and removes temporary files
	bash fastQC.sh $reads ${params.prefix}_R${label} $threads $reads
	
	#Logging QC statistics (number of sequences, Pass/warning/fail, basic statistics, duplication level, kmers)
	base=\$(basename $reads)
	bash logQC.sh \$base ${params.prefix}_R${label}_fastqc_data.txt .log.1$label
				
	#Measures and log execution time			
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"STEP 1 on $reads terminated at \$sysdate (\$exectime seconds)\" >> .log.1$label
	echo \" \" >> .log.1$label	
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.1$label
	echo \" \" >> .log.1$label	
	"""	
}



/**
	STEP 2. De-duplication of paired reads. Only exact duplicates are removed.

	Two FASTQ files are outputted, one for each paired-end.
*/

// Defines channel with <readfile1, readfile2> as input for de-duplicates
todedup = Channel.value( [file(params.reads1), file(params.reads2)] )
process dedup {
	
	input:
	set file(in1), file(in2) from todedup

	output:
	file  ".log.2" into log2
	set file("${params.prefix}_dedupe_R1.fq.gz"), file("${params.prefix}_dedupe_R2.fq.gz") into totrim
	set file("${params.prefix}_dedupe_R1.fq.gz"), file("${params.prefix}_dedupe_R2.fq.gz") into topublishdedupe

	when:
	  params.mode == "QC" || params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing STEP 2 [De-duplication] at \$sysdate\" > .log.2
	echo \" \" >> .log.2

	#De-duplicates
	clumpify.sh -Xmx$maxmem in1=$in1 in2=$in2 out1=${params.prefix}_dedupe_R1.fq.gz out2=${params.prefix}_dedupe_R2.fq.gz qin=$params.qin dedupe subs=0 threads=$threads &> tmp.log

	#Logs some figures about sequences passing de-duplication
	echo  \"BBduk's de-duplication stats: \" >> .log.2
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
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"STEP 2 terminated at \$sysdate (\$exectime seconds)\" >> .log.2
	echo \" \" >> .log.2
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.2
	echo \" \" >> .log.2
	"""
}

/**
	STEP 3. Trimming of low quality bases and of adapter sequences. Short reads
	are discarded. A decontamination of synthetic sequences is also pefoermed.
	When either forward or reverse of a paired-end read are discarded, the
	surviving reads are saved on a file of singleton reads.

	Three compressed FASTQ file (forward/reverse paired-end and singleton reads)
	are outputed. Several information are logged.
*/

//Two outputs will be produced (and pushed into two different channels): one
//included the paired trimmed reads, that will be QC'ed, the other the paired
//trimmed reads along with the singleton reads, that will be decontaminated
process trim {

	input:
   	set file(reads1), file(reads2) from totrim
	file(adapters) from Channel.from( file(params.adapters) )
	file(artifacts) from Channel.from( file(params.artifacts) )
	file(phix174ill) from Channel.from( file(params.phix174ill) )

	output:
	file  ".log.3" into log3
	set file("${params.prefix}_trimmed_R1.fq"), file("${params.prefix}_trimmed_R2.fq")  into trimmedreads
	set file("${params.prefix}_trimmed_R1.fq"), file("${params.prefix}_trimmed_R2.fq"), file("${params.prefix}_trimmed_singletons.fq") into todecontaminate
	set file("${params.prefix}_dedupe_R1.fq.gz"), file("${params.prefix}_dedupe_R2.fq.gz") into toremovetrimming
	set file("${params.prefix}_trimmed_R1.fq"), file("${params.prefix}_trimmed_R2.fq"), file("${params.prefix}_trimmed_singletons.fq") into topublishtrim 	

	when:
	  params.mode == "QC" || params.mode == "complete"

   	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing STEP 3 [Trimming] at \$sysdate\" > .log.3
	echo \" \" >> .log.3

	#Trim adapters and low quality sequences
	bbduk.sh -Xmx$maxmem in=$reads1 in2=$reads2 out=${params.prefix}_trimmed_R1_tmp.fq out2=${params.prefix}_trimmed_R2_tmp.fq outs=${params.prefix}_trimmed_singletons_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin threads=$threads tbo tpe ow &> tmp.log
	
	#Logs some figures about sequences passing trimming
	echo  \"BBduk's trimming stats (trimming adapters and low quality sequences): \" >> .log.3
	sed -n '/Input:/,/Result:/p' tmp.log >> .log.3
	echo \" \" >> .log.3
	unpairedR=\$(wc -l ${params.prefix}_trimmed_singletons_tmp.fq | cut -d\" \" -f 1)
	unpairedR=\$((\$unpairedR/4))
	echo  \"\$unpairedR singleton reads whose mate was trimmed shorter preserved\" >> .log.3
	echo \" \" >> .log.3

	#Removes synthetic contaminants (paired reads)
	bbduk.sh -Xmx$maxmem in=${params.prefix}_trimmed_R1_tmp.fq in2=${params.prefix}_trimmed_R2_tmp.fq out=${params.prefix}_trimmed_R1.fq out2=${params.prefix}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=$threads ow &> tmp.log

	#Logs some figures about sequences passing deletion of contaminants
	echo  \"BBduk's trimming stats (synthetic contaminants, paired reads): \" >> .log.3
	sed -n '/Input:/,/Result:/p' tmp.log >> .log.3
	echo \" \" >> .log.3

	#Removes synthetic contaminants (singleton reads)
	bbduk.sh -Xmx$maxmem in=${params.prefix}_trimmed_singletons_tmp.fq out=${params.prefix}_trimmed_singletons.fq k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=$threads ow &> tmp.log

	#Logs some figures about sequences passing deletion of contaminants
	echo  \"BBduk's trimming stats (synthetic contaminants, singleton reads): \" >> .log.3
	sed -n '/Input:/,/Result:/p' tmp.log >> .log.3
	echo \" \" >> .log.3
	
	#Removing tmp files
	rm ${params.prefix}_trimmed_R1_tmp.fq ${params.prefix}_trimmed_R2_tmp.fq ${params.prefix}_trimmed_singletons_tmp.fq

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"STEP 3 terminated at \$sysdate (\$exectime seconds)\" >> .log.3
	echo \" \" >> .log.3
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.3
	echo \"\" >> .log.3
	"""
}


/**
	STEP 4. Quality control after trimming (paired reads only).

	An output similar to that generated by STEP 1 is produced
*/

process qualityAssessmentTrimmed {

	publishDir  workingdir, mode: 'copy', pattern: "*.{html,txt}"
	
	//Merges the files with their labels (that is 1 for R1 and 2 for R2)
	input:
   	set file(reads), val(label) from trimmedreads.flatMap().merge( Channel.from( ['1', '2'] ) ){ a, b -> [a, b] }

	output:
	file  ".log.4$label" into log4
	file "${params.prefix}_trimmed_R${label}_fastqc.html"
	file "${params.prefix}_trimmed_R${label}_fastqc_data.txt"

	when:
	  params.mode == "QC" || params.mode == "complete"

   	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing STEP 4 [Assessment of read quality of trimmed FASTQ file] at \$sysdate on $reads\" > .log.4$label
	echo \" \" >> .log.4$label

	#Does QC, extracts relevant information, and removes temporary files
	bash fastQC.sh $reads ${params.prefix}_trimmed_R${label} $threads

	#Logging QC statistics (number of sequences, Pass/warning/fail, basic statistics, duplication level, kmers)
	base=\$(basename $reads)
	bash logQC.sh \$base ${params.prefix}_trimmed_R${label}_fastqc_data.txt .log.4$label

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"STEP 4 on $reads terminated at \$sysdate (\$exectime seconds)\" >> .log.4$label
	echo \" \" >> .log.4$label
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.4$label
	echo \" \" >> .log.4$label
	"""
}


/**
	STEP 5. Decontamination.
	Removes contamination with external organisms, that in this specific case is
	linmited to human contamination. Decontamination is carried on idependently
	on paired reads and on singleton reads thanks to BBwrap, that calls BBmap once
	on the paired reads and once on the singleton ones, merging the results on a
	single output file.

	Two files are outputted: the FASTQ of the decontaminated reads (including both
	paired-reads and singletons) and that of the contaminating reads (that can be
	used for refinements/checks).
	Please note that if keepQCtmpfile is set to false, the file of the contaminating 
	reads is discarded

	TODO: use BBsplit for multiple organisms decontamination
*/

process decontaminate {

	publishDir  workingdir, mode: 'copy', pattern: "*_clean.fq"
		
	input:
	set file(infile1), file(infile2), file(infile12) from todecontaminate
	file(refForeingGenome) from Channel.from( file(params.refForeingGenome, type: 'dir') )
	
	output:
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
	echo \"Performing STEP 5 [Decontamination] at \$sysdate\" > .log.5
	echo \" \" >> .log.5

	bbwrap.sh  -Xmx$maxmem mapper=bbmap append=t in1=$infile1,$infile12 in2=$infile2,null outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred path=$refForeingGenome qin=$params.qin threads=$threads untrim quickmatch fast ow &> tmp.log
	
	#Logs some figures about decontaminated/contaminated reads
	echo  \"BBmap's human decontamination stats (paired reads): \" >> .log.5
	sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | head -17 >> .log.5
	echo \" \" >> .log.5
	sed -n '/Read 2 data:/,/N Rate:/p' tmp.log >> .log.5
	echo \" \" >> .log.5
	echo  \"BBmap's human decontamination stats (singletons reads): \" >> .log.5
	sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | tail -17 >> .log.5
	echo \" \" >> .log.5

	nClean=\$(wc -l ${params.prefix}_clean.fq | cut -d\" \" -f 1)
	nClean=\$((\$nClean/4))
	nCont=\$(wc -l ${params.prefix}_cont.fq | cut -d\" \" -f 1)
	nCont=\$((\$nCont/4))
	echo \"\$nClean reads survived decontamination (\$nCont reads removed)\" >> .log.5
	echo \" \" >> .log.5

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"STEP 5 terminated at \$sysdate (\$exectime seconds)\" >> .log.5
	echo \" \" >> .log.5
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.5
	echo \"\" >> .log.5
	"""
}



/**
	STEP 6. Quality control after decontamination. The clean FASTQ file producted by
	the decontamination step is assessed with FastQC.

	An output similar to that generated by STEP 1 is produced
*/

process qualityAssessmentClean {

	publishDir  workingdir, mode: 'copy', pattern: "*.{html,txt}"
	
	input:
	file reads from decontaminatedreads

	output:
	file ".log.6" into log6
	file "${params.prefix}_clean_fastqc.html" 
	file "${params.prefix}_clean_fastqc_data.txt" 

	when:
	  params.mode == "QC" || params.mode == "complete"

   	script:
	"""		
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing STEP 6 [Assessment of read quality of decontaminated FASTQ file] at \$sysdate\" >  .log.6
	echo \" \" >> .log.6

	#Does QC, extracts relevant information, and removes temporary files
	bash fastQC.sh $reads ${params.prefix}_clean $threads $reads

	#Logging QC statistics (number of sequences, Pass/warning/fail, basic statistics, duplication level, kmers)
	base=\$(basename $reads)
	bash logQC.sh \$base ${params.prefix}_clean_fastqc_data.txt .log.6
		
	#Measures and log execution time			
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"STEP 6 on $reads terminated at \$sysdate (\$exectime seconds)\" >> .log.6
	echo \" \" >> .log.6	
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.6
	echo \" \" >> .log.6	
	"""	
}


//When doing community characterisation from previously QC'd files, this file should be 
//pushed in the corret channels.
//Please note that the name of the file externally QC'ed should be the same of the 
//one geneared by the YAMP, that is prefix_clean.fq, and should include ALL the reads.
if (params.mode == "characterisation") {
	toprofiletaxa =  Channel.from( file("$workingdir/${params.prefix}_clean.fq") )
	toprofilefunctionreads = Channel.from( file("$workingdir/${params.prefix}_clean.fq") )
}

/**
	STEP 7. Performs taxonomic binning and estimates the microbial relative abundancies.
	MetaPhlAn2 and its databases of clade-specific markers are used to infers the presence
	and relative abundance of the organisms (at the specie level) that are present in the
	sample and to estimate their relative abundance.

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
	echo \"Performing STEP 7 [Estimates microbial abundancies] at \$sysdate\" > .log.7
	echo \" \" >> .log.7
	
	rm -rf ${params.prefix}_bt2out.txt
	metaphlan2.py --input_type fastq --tmp_dir=. --biom ${params.prefix}.biom --bowtie2out=${params.prefix}_bt2out.txt --mpa_pkl $mpa_pkl  --bowtie2db $bowtie2db/$params.bowtie2dbfiles --bt2_ps $params.bt2options --nproc $threads $infile ${params.prefix}_metaphlan_bugs_list.tsv

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
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"\" >> .log.7
	echo \"STEP 7 terminated at \$sysdate (\$exectime seconds)\" >> .log.7
	echo \" \" >> .log.7
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.7
	echo \"\" >> .log.7
	"""
}

/**
	STEP 8. Evaluates alpha-diversity, that is the mean species diversity the
	given sample. Please note that the alpha diversity is the only per-sample
	measure, so it is the only one evaluated by this module. If a newick tree
 	is provided as input (see QIIME documentation for details), a more reliable
	phylogenetic measure is evaluated.

	One text file, listing the alpha-diversity values evaluated by means of
	multiple measure, is outputted.
*/

process alphaDiversity {

	publishDir  workingdir, mode: 'copy', pattern: "*.{tsv}"
	
	input:
	file(infile) from toalphadiversity
	
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
	echo \"Performing STEP 8 [Evaluating alpha-diversity] at \$sysdate\" > .log.8
	echo \" \" >> .log.8

	#It checks if the profiling was successful, that is if identifies at least three species
	n=\$(grep -o s__ $infile | wc -l  | cut -d\" \" -f 1)
	if (( n > 3 ))
	then
		#If the tree path is not specified, not all the alpha measures can be evaluated (that is,
		#PD_whole_tree is skipped)
		if [ $params.treepath == null ]
		then
			alpha_diversity.py -i $infile -o ${params.prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong &> /dev/null
		else
			alpha_diversity.py -i $infile -o ${params.prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong,PD_whole_tree &> /dev/null
		fi
	else
		#Also if the alpha are not evaluated the file should be created in order to be returned
		echo \"Not enough species detected (N=\$n). Analysis skipped.\" >> .log.8
		touch ${params.prefix}_alpha_diversity.tsv 
	fi
	
	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
	sysdate=\$(date)
	echo \"\" >> .log.8
	echo \"STEP 8 terminated at \$sysdate (\$exectime seconds)\" >> .log.8
	echo \" \" >> .log.8
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.8
	echo \"\" >> .log.8
	"""
}


/**
	STEP 9. Performs the functional annotation using HUMAnN2.
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
 	echo \"Performing STEP 9 [Performing functional annotation] with HUMAnN2 at \$sysdate\" > .log.9
 	echo \" \" >> .log.9

	#Runs HUMAnN2 taking advantages of the MetaPhlAn2's results
	humann2 --input $cleanreads --output . --output-basename ${params.prefix} --taxonomic-profile $metaphlanbuglist --nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads $threads --memory-use minimum > ${params.prefix}_HUMAnN2.log

	#If `|| true` is not add, nextflow stops... WTF 
	grep \"Total species selected from prescreen:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Selected species explain\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Unaligned reads after nucleotide alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Total gene families after translated alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Unaligned reads after translated alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true

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
 	exectime=\$(echo \"\$endtime - \$starttime\" | bc)
 	sysdate=\$(date)
 	echo \"\" >> .log.9
 	echo \"STEP 9 terminated at \$sysdate (\$exectime seconds)\" >> .log.9
 	echo \" \" >> .log.9
 	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.9
 	echo \"\" >> .log.9
 	"""
}


/**
	CLEANUP 1. Collapses all the logs resulting from the QC in the main one, and removed them.
*/

process logQC {

	input:
	set file(log11), file(log12), file(log2), file(log3), file(log41), file(log42), file(log5), file(log6) from log1.flatMap().mix(log2, log3, log4.flatMap(), log5, log6).toSortedList( { a, b -> a.name <=> b.name } )

	when:
	params.mode == "QC" || params.mode == "complete"

	script:
	"""
	cat $log11 $log12 $log2 $log3 $log41 $log42 $log5 $log6 >> $mylog
	"""
}

/**
	CLEANUP 2. Saves the temporary files generate during QC (if the user requested so)
*/

	
process saveQCtmpfile {

	publishDir  workingdir, mode: 'copy'
		
	input:
	file (tmpfile) from topublishdedupe.mix(topublishtrim, topublishdecontaminate).flatMap()

	output:
	file "$tmpfile"

	when:
	(params.mode == "QC" || params.mode == "complete") && params.keepQCtmpfile
		
	script:
	"""
	echo $tmpfile
	"""
}

/**
	CLEANUP 3. Collapses all the logs resulting from the community characterisation steps
	in the main one, and removed them.
*/

process logCC {

	input:
	set file(log7), file(log8), file(log9) from log7.mix(log8, log9).flatMap().toSortedList( { a, b -> a.name <=> b.name } )
	
	when:
	params.mode == "characterisation" || params.mode == "complete"
		
	script:
	"""
	cat $log7 $log8 $log9 >> $mylog
	"""
}

/**
	CLEANUP 4. Saves the temporary files generate during the community characterisation 
	(if the user requested so)
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

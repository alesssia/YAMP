/**
	Gets software version. 

	This process ensures that software version are included in the logs.
*/

process get_software_versions {

	//Starting the biobakery container. I need to run metaphlan and Humann to get
	//their version number (due to the fact that they live in the same container)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	output:
	path("software_versions_mqc.yaml"), emit: software_versions_yaml

	script:
	//I am using a multi-containers scenarios, supporting docker and singularity
	//with the software at a specific version (the same for all platforms). Therefore, I
	//will simply parse the version from there. Perhaps overkill, but who cares?  
	//This is not true for the biobakery suite (metaphlan/humann) which extract the 
	//information at runtime from the actual commands (see comment above)
	"""
	echo $workflow.manifest.version > v_pipeline.txt
	echo $workflow.nextflow.version > v_nextflow.txt

	echo $params.docker_container_fastqc | cut -d: -f 2 > v_fastqc.txt
	echo $params.docker_container_bbmap | cut -d: -f 2 > v_bbmap.txt
	
	metaphlan --version > v_metaphlan.txt
	humann --version > v_humann.txt
	echo $params.docker_container_qiime2 | cut -d: -f 2 > v_qiime.txt
	
	echo $params.docker_container_multiqc | cut -d: -f 2 > v_multiqc.txt
	
	scrape_software_versions.py > software_versions_mqc.yaml
	"""
}

// The user will specify the clean file either as a single clean file (that is the YAMP
// default behaviour), or as two files (forward/reverse). ]
// In the former case, the user will set singleEnd = true and only one file will be 
// selected and used directly for taxa and community profiling.
// In the latter case, the user will set singleEnd = false and provide two files, that will
// be merged before feeding the relevant channels for profiling.
process merge_paired_end_cleaned {

	tag "$name"

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}_QCd.fq.gz"), emit: to_profile_merged

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

// ------------------------------------------------------------------------------   
//	MULTIQC LOGGING
// ------------------------------------------------------------------------------   


/**
	Generate Logs. 

	Logs generate at each analysis step are collected and processed with MultiQC 
*/

process log {
	
	publishDir "${params.outdir}/${params.prefix}", mode: 'copy'

    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_multiqc
    } else {
        container params.docker_container_multiqc
    }

	input:
	file multiqc_config
	file workflow_summary
	file "software_versions_mqc.yaml"
	path "fastqc/*" 
	file "dedup_mqc.yaml"
	file "synthetic_contaminants_mqc.yaml"
	file "trimming_mqc.yaml"
	file "foreign_genome_indexing_mqc.yaml"
	file "decontamination_mqc.yaml"
	file "merge_paired_end_cleaned_mqc.yaml"
	file "profile_taxa_mqc.yaml"
	file "profile_functions_mqc.yaml"
	file "alpha_diversity_mqc.yaml"
	
	output:
	path("*multiqc_report*.html"), emit: multiqc_report
	path("*multiqc_data*")

	script:
	"""
	multiqc --config $multiqc_config . -f
	mv multiqc_report.html ${params.prefix}_multiqc_report_${params.mode}.html
	mv multiqc_data ${params.prefix}_multiqc_data_${params.mode}
	"""
}




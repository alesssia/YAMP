

/**
	Community Characterisation - STEP 1. Performs taxonomic binning and estimates the 
	microbial relative abundances using MetaPhlAn and its databases of clade-specific markers.
*/

process profile_taxa {

    tag "$name"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*.{biom,tsv}"
	
	input:
	tuple val(name), file(reads)
	file (bowtie2db)
	
	output:
	tuple val(name), path("*.biom"), emit: to_alpha_diversity
	tuple val(name), path("*_metaphlan_bugs_list.tsv"), emit: to_profile_function_bugs
	path("profile_taxa_mqc.yaml"), emit: profile_taxa_log
	
	when:
	params.mode != "QC"
	
	script:
	"""
	#If a file with the same name is already present, Metaphlan2 used to crash, leaving this here just in case
	rm -rf ${name}_bt2out.txt
	
	metaphlan --input_type fastq --tmp_dir=. --biom ${name}.biom --bowtie2out=${name}_bt2out.txt --bowtie2db $bowtie2db --bt2_ps ${params.bt2options} --add_viruses --sample_id ${name} --nproc ${task.cpus} $reads ${name}_metaphlan_bugs_list.tsv &> profile_taxa_mqc.txt
	
	# MultiQC doesn't have a module for Metaphlan yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_profile_taxa_log.sh ${name}_metaphlan_bugs_list.tsv > profile_taxa_mqc.yaml
	"""
}


/**
	Community Characterisation - STEP 2. Performs the functional annotation using HUMAnN.
*/

process profile_function {
	
    tag "$name"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*.{tsv,log}"
	
	input:
	tuple val(name), file(reads)
	tuple val(name), file(metaphlan_bug_list)
	file (chocophlan)
	file (uniref)
	
    output:
	path("*_HUMAnN.log")
	path("*_genefamilies.tsv")
	path("*_pathcoverage.tsv")
	path("*_pathabundance.tsv")
	path("profile_functions_mqc.yaml"), emit: profile_functions_log

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


/**
	Community Characterisation - STEP 3. Evaluates several alpha-diversity measures. 

*/

process alpha_diversity {

    tag "$name"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_qiime2
    } else {
        container params.docker_container_qiime2
    }

	publishDir "${params.outdir}/${params.prefix}", mode: 'copy', pattern: "*.{tsv}"
	
	input:
	tuple val(name), file(metaphlan_bug_list)
		
    output:
	path("*_alpha_diversity.tsv")
	path("alpha_diversity_mqc.yaml"), emit: alpha_diversity_log
	
	when:
	params.mode != "QC"

	script:
	"""
	#It checks if the profiling was successful, that is if identifies at least three species
	n=\$(grep -o s__ $metaphlan_bug_list | wc -l  | cut -d\" \" -f 1)
	if (( n <= 3 )); then
		#The file should be created in order to be returned
		touch ${name}_alpha_diversity.tsv 
	else
		echo $name > ${name}_alpha_diversity.tsv
		qiime tools import --input-path $metaphlan_bug_list --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ${name}_abundance_table.qza &> /dev/null
		for alpha in ace berger_parker_d brillouin_d chao1 chao1_ci dominance doubles enspie esty_ci fisher_alpha gini_index goods_coverage heip_e kempton_taylor_q lladser_pe margalef mcintosh_d mcintosh_e menhinick michaelis_menten_fit osd pielou_e robbins shannon simpson simpson_e singles strong
		do
			qiime diversity alpha --i-table ${name}_abundance_table.qza --p-metric \$alpha --output-dir \$alpha &> /dev/null
			qiime tools export --input-path \$alpha/alpha_diversity.qza --output-path \${alpha} &> /dev/null
			value=\$(sed -n '2p' \${alpha}/alpha-diversity.tsv | cut -f 2)
		    echo -e  \$alpha'\t'\$value 
		done >> ${name}_alpha_diversity.tsv  
	fi

	# MultiQC doesn't have a module for qiime yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash generate_alpha_diversity_log.sh \${n} > alpha_diversity_mqc.yaml	
	"""
}



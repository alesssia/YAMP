# MultiQC doesn't have a module for human yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

tot_species_prescreeing=$(grep "Total species selected from prescreen:" $2 | cut -d: -f 2 | sed 's/ //g')
selected_species_explain=$(grep "Selected species explain" $2 | cut -d" " -f 4-)
unaligned_reads_nucleotide=$(grep "Unaligned reads after nucleotide alignment:" $2 | cut -d: -f 2 | sed 's/ //g')
tot_gene_family=$(grep "Total gene families after translated alignment:" $2 | cut -d: -f 2 | sed 's/ //g')
unaligned_reads_translated=$(grep "Unaligned reads after translated alignment:" $2 | cut -d: -f 2 | sed 's/ //g')


# Dump to YAML (header)
echo "id: 'humann'"
echo "section_name: 'HUMAnN'" 
echo "section_href: 'https://github.com/alesssia/yamp'" 
echo "plot_type: 'html'" 
echo "description: 'This information is collected at run time from the software output.'" 
echo "data: |" 
echo "    <dl class="dl-horizontal">" 
echo  "        <dt>Selected from prescreen</dt><dd>"${tot_species_prescreeing}" species</dd>" 
echo  "        <dt>Selected species explain</dt><dd>"${selected_species_explain}"</dd>" 
echo  "        <dt>Unaligned</dt><dd>"${unaligned_reads_nucleotide}" reads unaligned after nucleotide alignment</dd>" 
echo  "        <dt>Unaligned</dt><dd>"${unaligned_reads_translated}" reads unaligned after translated alignment</dd>" 
echo  "        <dt>Total gene families</dt><dd>"${tot_gene_family}" (after translated alignment)</dd>" 
echo  "        <dt>Complete log</dt><dd>"$1"_HUMAnN.log</dd>" 
echo "    </dl>" 


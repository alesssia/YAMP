# MultiQC doesn't have a module for Metaphlan yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

# Dump to YAML (header)
echo "id: 'metaphlan'"
echo "section_name: 'MetaPhlAn'" 
echo "section_href: 'https://github.com/alesssia/yamp'" 
echo "plot_type: 'html'" 
echo "description: 'This information is collected at run time from the software output.'" 
echo "data: |" 
echo "    <dl class="dl-horizontal">" 

#Looks for warnings
warning=$(grep WARNING profile_taxa_mqc.txt | wc -l | cut -d" " -f 1)
if (( warning != 0 )); then
	echo  "        <dt>Multiple merged species:</dt><dd>Present, additional column listing added</dd>" 
fi

#Logs some info regarding the taxomic tree
echo "        <dt>"Detected"</dt><dd></dd>" 
tree=(Kingdom Phylum Class Order Family Genus Species)
for i in {1..7}
do
	c=$(sed '1d' $1 | cut -d"|" -f $i | grep -v "k__" | cut -f 1  | sort | uniq | sed '/^\\s*\$/d' | wc -l | cut -d" " -f 1)
	echo  "        <dt>"${tree[(($i-1))]}"</dt><dd>"$c"</dd>" 
done

echo "    </dl>" 

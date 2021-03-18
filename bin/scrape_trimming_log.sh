# MultiQC doesn't have a module for bbduk yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

# Dump to YAML (header)
echo "id: 'trimming'"
echo "section_name: 'YAMP Trimming'" 
echo "section_href: 'https://github.com/alesssia/yamp'" 
echo "plot_type: 'html'" 
echo "description: 'This information is collected at run time from the software output.'" 
echo "data: |" 
echo "    <dl class="dl-horizontal">" 

# Dump some figures about sequences passing trimming 
while IFS= read -r line
do
	echo $line | sed 's/  */ /g' | awk -F: '{print "        <dt>"$1":</dt><dd>"$2"</dd>" }' 
done < <(sed -n '/Input:/,/Result:/p'  trimming_mqc.txt) 

# Logs some more (to standardize with other outputs)
totR=$(grep "Input:" trimming_mqc.txt | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
remR=$(grep "Total Removed:" trimming_mqc.txt | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
survivedR=$(($totR-$remR))
percentage=$(echo $survivedR $totR | awk '{print $1/$2*100}' )
time=$(grep "Time:" trimming_mqc.txt | cut -d: -f 2 | cut -f 2 | sed 's/s\./s/g')

# If singleEnd, there could be some reads without mate.
# The file is generated also if empty.
if ls *_trimmed_singletons.fq.gz 1> /dev/null 2>&1; then
	file=$(ls *_trimmed_singletons.fq.gz)
	unpairedR=$(zcat $file | wc -l | cut -d" " -f 1)
	unpairedR=$(($unpairedR/4))
	echo "        <dt>Unpaired reads:</dt><dd>"$unpairedR"</dd>" 
fi

echo "        <dt>Survived de-duplication:</dt><dd>"$survivedR"</dd>" 
echo "        <dt>Survived (percentage)</dt><dd>"$percentage"</dd>" 
echo "        <dt>Total time:</dt><dd>"$time"</dd>" 
echo "    </dl>" 


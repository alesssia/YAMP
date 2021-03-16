# MultiQC doesn't have a module for bbduk yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

# Dump to YAML (header)
echo "id: 'synthetic-contaminants'" > synthetic_contaminants_mqc.yaml
echo "section_name: 'YAMP Remove synthetic contaminants'" >> synthetic_contaminants_mqc.yaml
echo "section_href: 'https://github.com/alesssia/yamp'" >> synthetic_contaminants_mqc.yaml
echo "plot_type: 'html'" >> synthetic_contaminants_mqc.yaml
echo "description: 'This information is collected at run time from the software output.'" >> synthetic_contaminants_mqc.yaml
echo "data: |" >> synthetic_contaminants_mqc.yaml
echo "    <dl class="dl-horizontal">" >> synthetic_contaminants_mqc.yaml

# Dump some figures about sequences passing trimming
while IFS= read -r line
do
	echo $line | sed 's/  */ /g' | awk -F: '{print "        <dt>"$1":</dt><dd>"$2"</dd>" }'
done < <(sed -n '/Input:/,/Result:/p' synthetic_contaminants_mqc.txt) >> synthetic_contaminants_mqc.yaml

# Log some figures about sequences passing deletion of contaminants
totR=$(grep "Input:" synthetic_contaminants_mqc.txt | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
remR=$(grep "Contaminants:" synthetic_contaminants_mqc.txt | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
survivedR=$(($totR-$remR))
percentage=$(echo $survivedR $totR | awk '{print $1/$2*100}' )
percentage=`printf "%.2f" $percentage`
time=$(grep "Time:" synthetic_contaminants_mqc.txt | cut -d: -f 2 | cut -f 2 | sed 's/s\./s/g')

# Logs some more (to standardize with other outputs)
echo "        <dt>Surviving:</dt><dd>"$survivedR" ("$percentage"%)</dd>" >> synthetic_contaminants_mqc.yaml
echo "        <dt>Total time:</dt><dd>"$time"</dd>" >> synthetic_contaminants_mqc.yaml
echo "    </dl>" >> synthetic_contaminants_mqc.yaml



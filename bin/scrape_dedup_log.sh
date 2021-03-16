# MultiQC doesn't have a module for clumpify yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

# Log some figures about sequences passing de-duplication
totR=$(grep "Reads In:" dedup_mqc.txt | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
remR=$(grep "Duplicates Found:" dedup_mqc.txt | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
survivedR=$(($totR-$remR))
percentage=$(echo $survivedR $totR | awk '{print $1/$2*100}' )
percentage=`printf "%.2f" $percentage`
time=$(grep "Total time:" dedup_mqc.txt | cut -d: -f 2 | cut -f 2 | sed 's/s\.//g')

# Dump to YAML
echo "id: 'deduplication'" > dedup_mqc.yaml
echo "section_name: 'YAMP Deduplication'" >> dedup_mqc.yaml
echo "section_href: 'https://github.com/alesssia/yamp'" >> dedup_mqc.yaml
echo "plot_type: 'html'" >> dedup_mqc.yaml
echo "description: 'This information is collected at run time from the software output.'" >> dedup_mqc.yaml
echo "data: |" >> dedup_mqc.yaml
echo "    <dl class="dl-horizontal">" >> dedup_mqc.yaml
echo "        <dt>Reads In:</dt><dd>"$totR"</dd>" >> dedup_mqc.yaml
echo "        <dt>Duplicated Found:</dt><dd>"$remR"</dd>" >> dedup_mqc.yaml
echo "        <dt>Surviving:</dt><dd>"$survivedR" ("$percentage"%)</dd>" >> dedup_mqc.yaml
echo "        <dt>Total time:</dt><dd>"$time"</dd>" >> dedup_mqc.yaml
echo "    </dl>" >> dedup_mqc.yaml



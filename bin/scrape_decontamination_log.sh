# MultiQC doesn't have a module for bbwrap yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

nClean=$(gunzip -c *_QCd.fq.gz | wc -l | cut -d" " -f 1)
nClean=$(($nClean/4))
nCont=$(gunzip -c *_contamination.fq.gz | wc -l | cut -d" " -f 1)
nCont=$(($nCont/4))
totR=$(($nCont+$nClean))
percentage=$(echo $nClean $totR | awk '{print $1/$2*100}' )
percentage=`printf "%.2f" $percentage`
time1=$(grep "Total time:" decontamination_mqc.txt | head -1 | cut -d: -f 2 | cut -f 2 | sed 's/s\./s/g')



# Dump to YAML
echo "id: 'decontamination'" > decontamination_mqc.yaml
echo "section_name: 'YAMP Decontamination'" >> decontamination_mqc.yaml
echo "section_href: 'https://github.com/alesssia/yamp'" >> decontamination_mqc.yaml
echo "plot_type: 'html'" >> decontamination_mqc.yaml
echo "description: 'This information is collected at run time from the software output.'" >> decontamination_mqc.yaml
echo "data: |" >> decontamination_mqc.yaml
echo "    <dl class="dl-horizontal">" >> decontamination_mqc.yaml
echo "        <dt>Reads In:</dt><dd>"$totR"</dd>" >> decontamination_mqc.yaml
echo "        <dt>Contaminant Found:</dt><dd>"$nCont"</dd>" >> decontamination_mqc.yaml
echo "        <dt>Survived decontamination:</dt><dd>"$nClean"</dd>" >> decontamination_mqc.yaml
echo "        <dt>Survived (percentage):</dt><dd>"$percentage"</dd>" >> decontamination_mqc.yaml
echo "        <dt>Time (paired reads):</dt><dd>"$time1"</dd>" >> decontamination_mqc.yaml

if ls *_trimmed_singletons.fq.gz 1> /dev/null 2>&1; then
	time2=$(grep "Total time:" decontamination_mqc.txt | tail -1 | cut -d: -f 2 | cut -f 2 | sed 's/s\./s/g')
	echo "        <dt>Time (singletons reads):</dt><dd>"$time2"</dd>" >> decontamination_mqc.yaml
fi

echo "    </dl>" >> decontamination_mqc.yaml



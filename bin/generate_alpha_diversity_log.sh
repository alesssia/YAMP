# MultiQC doesn't have a module for qiime alpha diversity yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

# Dump to YAML
echo "id: 'alpha-diversity'" 
echo "section_name: 'QIIME'"
echo "section_href: 'https://github.com/alesssia/yamp'"
echo "plot_type: 'html'" 
echo "description: 'This information is collected at run time from the software output.'"
echo "data: |"
echo "    <dl class="dl-horizontal">"
echo "        <dt>Number of species:</dt><dd>"$1"</dd>"

if (( $1 < 3 )); then
	echo "        <dt>qiime diversity alpha</dt><dd>non performed</dd>"
else
	echo "        <dt>qiime diversity alpha</dt><dd>performed</dd>"
fi 

echo "    </dl>"


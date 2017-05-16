#  Yet Another Metagenomic Pipeline (YAMP)
#  Copyright (C) 2017 	Dr Alessia Visconti
#  	      
#  This script is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this script.  If not, see <http://www.gnu.org/licenses/>.
#  
#  For any bugs or problems found, please contact us at
#  alessia.visconti@kcl.ac.uk


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This scripts executes fastQC, extracts relevant information from the output 
# (that is the html file and the textual summary stored in the zip archive) and
# removes the temporary files (that is, the zip archive).
# It requires three parameters, namely:
# $1, the fastq file to assess (e.g., mysample.fq.gz)
# $2, the label used to rename the summary files (e.g., "mysample")
# $3, the number of threads
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#This has been added to deal with the Docker image (fastQC requires X11 to be set)
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

#Performs QC
fastqc --quiet --noextract --format fastq --outdir=. --threads $3 $1

mv *_fastqc.html ${2}_fastqc.html
unzip -p *_fastqc.zip ${base}*/fastqc_data.txt > ${2}_fastqc_data.txt  
rm -rf *.zip	
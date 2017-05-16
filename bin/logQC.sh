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
# This scripts extracts QC metrics from fastQC summary output
# It requires three parameters, namely:
# $1, the file that went through QC (e.g., mysample.fq.gz or "original file")
# $2, the fastQC summary file (e.g.,  mysample.fq.1_data.txt), stored in the
#	  archive created by fastQC
# $3, the path to the log file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

echo "Summary of $1's basic statistic" >> $3

#Number of sequences
nR=$(grep "Total Sequences" $2 | cut -f 2)
echo "$1's total reads: $nR" >> $3

#Pass/warning/fail
grep ">>" $2 | grep -v MODULE | grep -v "Basic Statistics" | sed 's/>>//g' >> $3
echo "" >> $3
#Basic statistic
sed -n '/Basic Statistics/,/END_MODULE/p' $2 | grep -v MODULE | grep -v "Basic Statistics" | sed 's/>>//g'  >> $3
echo "" >> $3
#Duplication level
sed -n '/Sequence Duplication Levels/,/END_MODULE/p' $2 | grep -v MODULE | grep -v "Basic Statistics" | sed 's/>>//g'  >> $3
echo "" >> $3
#Kmers
sed -n '/Kmer Content/,/END_MODULE/p' $2 | grep -v MODULE | grep -v "Basic Statistics" | sed 's/>>//g' >> $3
echo "" >> $3





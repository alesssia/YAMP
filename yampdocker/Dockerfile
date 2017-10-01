FROM ubuntu:xenial
FROM continuumio/anaconda

MAINTAINER Alessia Visconti <alessia.visconti@gmail.com>

#Updates Ubuntu and intalls some packages that will be used either
#in for building other software of directly by YAMP
RUN apt-get update --fix-missing && apt-get install -q -y build-essential 
RUN apt-get install -q -y wget curl unzip gzip openjdk-7-jdk bc python perl samtools 

#Updateds anaconda and uses it to install software used by YAMP
RUN conda update conda -y
RUN conda install -c bioconda  metaphlan2=2.6.0 qiime=1.9.1 humann2=0.9.9 diamond=0.8.36

#Install the last two software, namely FastQC and BBmap
RUN wget -q -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    rm fastqc.zip

RUN wget -q -O BBMap_37.09.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_37.09.tar.gz && \
	tar -xzf BBMap_37.09.tar.gz && \
	mv bbmap/ /opt/ && \
	rm BBMap_37.09.tar.gz

RUN curl -fsSL get.nextflow.io | bash && \
	mv nextflow /opt/

ENV PATH $PATH:/opt/::/opt/FastQC/:/opt/bbmap


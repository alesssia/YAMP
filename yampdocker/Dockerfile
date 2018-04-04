FROM frolvlad/alpine-glibc:alpine-3.7

MAINTAINER Alessia Visconti <alessia.visconti@gmail.com>

#Installs packages that will be used either for building other software 
#or directly by YAMP
RUN apk --update add --no-cache bash procps wget curl gzip perl mesa-gl 

#Installs miniconda 
RUN wget -q -O Miniconda2-4.4.10-Linux-x86_64.sh http://repo.continuum.io/miniconda/Miniconda2-4.4.10-Linux-x86_64.sh && \
	bash Miniconda2-4.4.10-Linux-x86_64.sh -f -b -p /opt/conda && \
	rm Miniconda2-4.4.10-Linux-x86_64.sh

#Exports conda path
ENV PATH $PATH:/opt/conda/bin/

#Updateds conda and uses it to install software used by YAMP, as well AWS CLI
#that is required to use YAMP on AWS Batch
RUN conda update conda -y
RUN conda install -c bioconda -y bbmap=37.10 fastqc=0.11.5 metaphlan2=2.6.0 qiime=1.9.1 humann2=0.9.9 
RUN conda install -c conda-forge -y awscli
RUN conda clean --tarballs -y

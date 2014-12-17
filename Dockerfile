FROM debian:wheezy
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

RUN apt-get update --fix-missing && \
  apt-get install -q -y bc wget curl vim nano unzip make gcc g++ gfortran && \
  apt-get clean 

#
# Create the home folder 
#
RUN mkdir -p /root
ENV HOME /root

#
# Install pre-requistes
#
RUN apt-get install -q -y samtools
  
#
# RNA-Seq tools 
# 

RUN wget -q -O bowtie.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download && \
  unzip bowtie.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.3/ /opt/bowtie && \
  rm bowtie.zip 
  
RUN wget -q http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
  tar xf cufflinks-2.2.1.Linux_x86_64.tar.gz -C /opt/ && \
  ln -s /opt/cufflinks-2.2.1.Linux_x86_64/ /opt/cufflinks && \
  rm cufflinks-2.2.1.Linux_x86_64.tar.gz 
  
  
RUN wget -q http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.12.Linux_x86_64.tar.gz && \
  tar xf tophat-2.0.12.Linux_x86_64.tar.gz -C /opt/ && \
  ln -s /opt/tophat-2.0.12.Linux_x86_64/ /opt/tophat && \
  rm tophat-2.0.12.Linux_x86_64.tar.gz  
  
#
# Finalize environment
#
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bowtie:/opt/tophat:/opt/cufflinks
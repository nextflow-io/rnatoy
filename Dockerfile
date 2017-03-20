FROM pditommaso/dkrbase:1.2

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

#
# Install pre-requistes
#
RUN apt-get update --fix-missing && \
  apt-get install -q -y samtools python 
  
#
# RNA-Seq tools 
# 

RUN wget -q -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/bowtie2-2.2.7-linux-x86_64.zip/download && \
  unzip bowtie.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.7/ /opt/bowtie && \
  rm bowtie.zip 
  
RUN \
  wget -q http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz -O- \
  | tar xz -C /opt/ && \
  ln -s /opt/cufflinks-2.2.1.Linux_x86_64/ /opt/cufflinks 
  
  
RUN \
  wget -q https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz -O- \
  | tar xz -C /opt/ && \
  ln -s /opt/tophat-2.1.0.Linux_x86_64/ /opt/tophat 
  
#
# Finalize environment
#
ENV PATH=$PATH:/opt/bowtie:/opt/tophat:/opt/cufflinks

FROM pditommaso/dkrbase:1.1

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

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
  
RUN wget -q http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
  tar xf cufflinks-2.2.1.Linux_x86_64.tar.gz -C /opt/ && \
  ln -s /opt/cufflinks-2.2.1.Linux_x86_64/ /opt/cufflinks && \
  rm cufflinks-2.2.1.Linux_x86_64.tar.gz 

RUN apt-get install -q -y libncurses5-dev libncursesw5-dev && \
    wget -q http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.19.tar.gz && \
    tar xf rsem-1.2.19.tar.gz -C /opt/ && \
    make -C /opt/rsem-1.2.19 && \
    ln -s /opt/rsem-1.2.19/ /opt/rsem && \
    rm rsem-1.2.19.tar.gz  

  
#
# Finalize environment
#
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bowtie2:/opt/rsem:/opt/cufflinks

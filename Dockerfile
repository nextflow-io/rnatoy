FROM pditommaso/dkrbase:1.1

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

#
# Install pre-requistes
#
RUN apt-get install -q -y python bzip2 libncurses5-dev libncursesw5-dev
  
#
# Samtools
#  
RUN wget -q -O samtools-1.1.tar.bz2 'http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2F&ts=1420911249&use_mirror=heanet' && \
  bzip2 -d samtools-1.1.tar.bz2 && \
  tar xf samtools-1.1.tar && \
  rm samtools-1.1.tar && \
  cd samtools-1.1 && make && make install 

  
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

RUN wget -q http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.19.tar.gz && \
    tar xf rsem-1.2.19.tar.gz -C /opt/ && \
    make -C /opt/rsem-1.2.19 && \
    ln -s /opt/rsem-1.2.19/ /opt/rsem && \
    rm rsem-1.2.19.tar.gz  

  
#
# Finalize environment
#
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bowtie:/opt/rsem:/opt/cufflinks

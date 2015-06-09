FROM pditommaso/sl65base
MAINTAINER paolo.ditommaso@gmail.com


RUN wget -q -O- 'http://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F0.1.18%2F&ts=1430746369&use_mirror=kent' \
  | tar xj && \
  cd samtools-* && \ 
  make && \ 
  find . -maxdepth 1 -perm /a+x -exec cp {} /usr/local/bin \; && \
  cd - && rm -rf samtools-*

RUN wget -q -O bowtie.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download && \
  unzip bowtie.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.3/ /opt/bowtie && \
  rm bowtie.zip 
  
RUN \
  wget -q http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz -O- \
  | tar xz -C /opt/ && \
  ln -s /opt/cufflinks-2.2.1.Linux_x86_64/ /opt/cufflinks 
  
  
RUN \
  wget -q http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.12.Linux_x86_64.tar.gz -O- \
  | tar xz -C /opt/ && \
  ln -s /opt/tophat-2.0.12.Linux_x86_64/ /opt/tophat 
  
#
# Finalize environment
#
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bowtie:/opt/tophat:/opt/cufflinks 
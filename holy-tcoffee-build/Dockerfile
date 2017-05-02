FROM cbcrg/dkrbase:1.2
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

ENV HOME /root 
WORKDIR /root

#
# Install deps
# 
RUN apt-get update --fix-missing && \
  apt-get install -y antiword curl wget build-essential automake openjdk-7-jre-headless git ghostscript unzip

#
# Perl modules
#
RUN cpanm -n Net::SSLeay XML::Simple SOAP::Lite

#
# Installbuilder 
# 
RUN curl -s  https://s3-eu-west-1.amazonaws.com/cbcrg-eu/tcoffee-ci/installbuilder-8.6.0.tar.gz | tar xz &&\
  ln -s $HOME/installbuilder-8.6.0 $HOME/installbuilder

#
# Blast+
# 
RUN curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz | tar xz &&\
  ln -s $HOME/ncbi-blast-2.2.28+ $HOME/ncbi-blast 
  
ENV PATH=$HOME/ncbi-blast/bin:$PATH  


#
# Install argtable2
# 
RUN curl -fsSL http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz | tar xz &&\
  cd argtable2-13 &&\
  ./configure &&\
  make &&\
  make install

ENV OSNAME=linux OSARCH=x64 WORKSPACE=$HOME PERL_MM_USE_DEFAULT=1 PERL_EXTUTILS_AUTOINSTALL=--defaultdeps



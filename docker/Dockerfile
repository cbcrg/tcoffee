FROM debian:jessie
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

RUN apt-get update --fix-missing && \
  apt-get install -q -y bc vim procps ghostscript unzip curl wget && \
  apt-get install -q -y perl-modules libnet-ssleay-perl libcrypt-ssleay-perl libssl-dev libexpat1 libexpat1-dev liblwp-protocol-https-perl && \
  apt-get install -q -y libgfortran3 libglib2.0-0 libgomp1 && \
  apt-get clean  
  
RUN mkdir -p /root
ENV HOME /root  
WORKDIR /root

#
# Blast+
# 
RUN curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz | tar xz -C /opt &&\
  ln -s /opt/ncbi-blast-2.2.28+ /opt/ncbi-blast 

#
# Add T-Coffee
# 
ADD tcoffee /opt/tcoffee

#
# Environment 
#
ENV PATH=$PATH:/opt/ncbi-blast/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/ TEMP=/tmp PERL5LIB=/opt/tcoffee/perl/lib/perl5 DIR_4_TCOFFEE=/opt/tcoffee EMAIL_4_TCOFFEE=tcoffee.msa@gmail.com CACHE_4_TCOFFEE=/tmp/cache/ LOCKDIR_4_TCOFFEE=/tmp/lck/ TMP_4_TCOFFEE=/tmp/  




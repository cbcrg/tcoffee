set -x
set -e
if [ ! -e $HOME/ncbi-blast-2.2.28+ ]; then
  cd $HOME
  wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz
  tar xf ncbi-blast-2.2.28+-x64-linux.tar.gz
  rm -rf ncbi-blast-2.2.28+-x64-linux.tar.gz
fi


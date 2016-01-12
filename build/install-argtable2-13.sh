set -x
set -e

curl -fsSL http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz | tar xz 
cd argtable2-13 
./configure 
make 
make install
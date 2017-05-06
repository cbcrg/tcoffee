set -x
set -e

sudo wget -q --no-check-certificate cpanmin.us -O /usr/local/bin/cpanm && sudo chmod +x /usr/local/bin/cpanm  
sudo cpanm -n Net::SSLeay XML::Simple SOAP::Lite
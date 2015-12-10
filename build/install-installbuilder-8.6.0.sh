set -x
set -e
if [ ! -e $HOME/installbuilder-8.6.0 ]; then
  cd $HOME
  wget -q https://s3-eu-west-1.amazonaws.com/cbcrg-eu/tcoffee-ci/installbuilder-8.6.0.tar.gz
  tar xf installbuilder-8.6.0.tar.gz
  rm -rf installbuilder-8.6.0.tar.gz
fi
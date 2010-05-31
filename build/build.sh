#
# common vars 
#
_DIST_FILE=~/Dropbox/latest_distributions/T-COFFEE_distribution_$(OSNAME).tar.gz
_G_USR=cbcrg.lab
_G_PWD=Zf4Na7vf8SX8

#
# clean current sandbox content 
#
if [ -d ~/Dropbox ] 
then
rm -rf ~/Dropbox
fi

cd $WORKSPACE/main/t_coffee/src

#
# rename temporary makefile and run it
#
rm makefile
mv makefile_2 makefile
make distribution

#
# check that the distribution file has been  created
#
if [ ! -f $_DIST_FILE ] 
then 
echo "Destination file has not been created: $_DIST_FILE"
exit 1
if

#
# upload to google code
#
cd $WORKSPACE/main/build
./googlecode_upload.py -p tcoffee -s 'T-Coffee latest build distribution' -l tcoffee -u $_G_USR -w $_G_PWD $_DIST_FILE 
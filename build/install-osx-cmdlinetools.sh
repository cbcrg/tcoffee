set -x
set -e

#
# See http://apple.stackexchange.com/a/107314
# 
sudo rm -rf /Applications/Xcode.app/
wget -q --no-check-certificate https://s3-eu-west-1.amazonaws.com/cbcrg-eu/tcoffee-ci/commandline_tools_os_x_mountain_lion_for_xcode__march_2014.dmg
hdiutil attach commandline_tools_os_x_mountain_lion_for_xcode__march_2014.dmg
cd "/Volumes/Command Line Tools (Mountain Lion)"/
sudo installer -verbose -pkg "Command Line Tools (Mountain Lion).mpkg" -target /

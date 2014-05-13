#!/bin/bash
#
# Headless building script 
#
# Used environment variables 
# 
# $WORKSPACE	: root path containing all build artifacts 
# $OSNAME	: the target platform name i.e. linux, macosx, windows		
# $OSARCH	: the target platform architecture i.e. i386, x64
# $BUILD_REPO   : (optional) shared path to cache thirdy parts tools compiled binaries
# $GIT_REVISION : (optional) the svn revision used to mark this buld
# $VERSION 	: (optional) the version used to tag this build
# $BUILD_REPO	: (optional) the path where compiled binaries will be cached
# $USER_BIN	: (deprecated) path where store binary to be used 
#
# Required dependencies: 
# - wget 
# - antiword
# - installbuild (http://installbuilder.bitrock.com)
#

#
# Release flag  
# 
export RELEASE=${RELEASE:-0}


#
# default SVN revision number 
#
if [ -z $GIT_REVISION ]; then 
GIT_REVISION=`( cd $WORKSPACE/tcoffee; git rev-parse --short HEAD; )`
fi

if [ "$GIT_REVISION" == "" ]; then 
  echo 'Missing $GIT_REVISION value. Cannot continue the build process.' 
  exit 1
fi

#
# The build number is generated automatically by the build containg (Hudson/Jenkins)
#
if [ -z $BUILD_NUMBER ]; then 
	BUILD_NUMBER=0
fi

#
# Define the VERSION number 
#
if [[ (-z $VERSION) || ($VERSION == auto) ]]; then 
	export VERSION="`cat $WORKSPACE/tcoffee/lib/version_number.version`.$GIT_REVISION"
fi

#
# The date timestamp string contains also the svn revision number
#
if [ -z $DATE ]; then 
export DATE=`date +"%Y-%m-%d %H:%M:%S"`
fi

#
# BUILD_INFO
#
if [ -z $BUILD_INFO ]; then 
export BUILD_INFO="`date +"%Y-%m-%d %H:%M:%S"` - Revision $GIT_REVISION - Build $BUILD_NUMBER"
fi


# default bin path 
if [ -z $USER_BIN ]; then 
export USER_BIN=$WORKSPACE/bin/
fi

#
# default third party binaries cache location 
#
if [ -z $BUILD_REPO ]; then 
export BUILD_REPO=$WORKSPACE/repo
fi

#
# default install builder location 
#
if [ -z $INSTALLER ]; then 
	INSTALLER=~/installbuilder/bin/builder
	if [ $OSNAME == "macosx" ]
	then
	INSTALLER=~/installbuilder/bin/Builder.app/Contents/MacOS/installbuilder.sh
	fi
fi

# Flag DO_TEST, if true test are executed (default: false)
if [ -z $DO_TEST ]; then 
DO_TEST=0
fi 

#
# The Dropbox folder where to copy produces binaries 	
#
if [ -z $DROPBOX ]; then 
DROPBOX=$HOME/Dropbox/upstream; 
fi

#
# script directives
#

set -e
set -u
set -o nounset
set -o errexit
#set -x

#
# other common variables 
#
SANDBOX=$WORKSPACE/sandbox
PERLM=$SANDBOX/perl

_SRC=$WORKSPACE/tcoffee/t_coffee/src
_G_USR=cbcrg.lab
_G_PWD=Zf4Na7vf8SX8


#
# Define the distribution directory that will contain produced artifacts
#
DIST_ROOT=$SANDBOX/distributions

if [ $RELEASE == 1 ]; then
DIST_BASE=$SANDBOX/distributions/Stable/
else
DIST_BASE=$SANDBOX/distributions/Beta/
fi

# Distribution package file name
DIST_DIR=$DIST_BASE/$VERSION/$OSNAME
DIST_NAME=T-COFFEE_distribution_$VERSION.tar.gz
DIST_HOST='tcoffeeo@tcoffee.org:~/public_html/Packages/'

# Installer package file name 
INST_NAME=T-COFFEE_installer_"$VERSION"_"$OSNAME"_"$OSARCH"

UNTARED=$SANDBOX/untared_distributions/T-COFFEE_distribution_"$VERSION"
TCDIR=$SANDBOX/build

#
# exported variabled (required by 'generic_makefile')
#
export HOME2=$SANDBOX


#
# Display the current environment
#
function env() 
{
  echo "[ env ]"

  echo "- WORKSPACE   : $WORKSPACE"
  echo "- OSNAME      : $OSNAME"
  echo "- OSARCH      : $OSARCH"
  echo "- VERSION     : $VERSION"
  echo "- DATE        : $DATE"
  echo "- RELEASE     : $RELEASE"
  echo "- BUILD_INFO  : $BUILD_INFO"
  echo "- BUILD_REPO  : $BUILD_REPO"
  echo "- GIT_REVISION: $GIT_REVISION"
  echo "- USER_BIN    : $USER_BIN"
  echo ". SANDBOX     : $SANDBOX"
  echo ". _SRC        : $_SRC"
  echo ". UNTARED     : $UNTARED"  
  echo ". INSTALLER   : $INSTALLER"
  echo ". TCDIR       : $TCDIR"
  echo ". PERLM       : $PERLM"
  echo ". DIST_BASE   : $DIST_BASE"
  echo ". DIST_DIR    : $DIST_DIR"
  echo ". DIST_NAME   : $DIST_NAME"
  echo ". DIST_HOST   : $DIST_HOST"
  echo ". INST_NAME   : $INST_NAME"
  echo ". DO_TEST     : $DO_TEST"

}

#
# clean current sandbox content 
#
function clean() 
{
	echo "[ clean ]"
	rm -rf $SANDBOX

}

#
# Execute legacy doc_test target 
#
function doc_test() { 
	echo "[ doc_test ]"

	set +u
	if [ -z $TEST_HTML_PREFIX ]; 	then TEST_HTML_PREFIX="all"; fi
	if [ -z $TEST_STOP ]; 			then TEST_STOP="error"; fi
	if [ -z $TEST_DELETE ]; 		then TEST_DELETE="never"; fi
	if [ -z $TEST_SANDBOX ];		then TEST_SANDBOX="$WORKSPACE/test-results/all"; fi
	if [ -z $TEST_OUTPUT ]; 		then TEST_OUTPUT="$WORKSPACE/test-results/index.html"; fi
	if [ -z $TEST_FILES ];			then TEST_FILES="-R ./all"; fi

	TEST_CMDLINE="--var tcoffee.home=$TCDIR --stop $TEST_STOP --delete $TEST_DELETE --sandbox-dir \"$TEST_SANDBOX\" --html-path-prefix \"$TEST_HTML_PREFIX\" -o \"$TEST_OUTPUT\" $TEST_ARGS $TEST_FILES"
	echo Test parameters: $TEST_CMDLINE
	set -u
	
	# remove previous result (if any)
	rm -rf $TEST_SANDBOX
	
	# run tests 
	cd $WORKSPACE/tcoffee/testsuite/
	
	set +e
	java -jar black-coffee.jar $TEST_CMDLINE

	if [ $? != 0 ]; then
		echo "Some test FAILED. Check result file: $TEST_OUTPUT "
		exit 2
	fi


	if [ $? == 0 ]; then
		echo "All tests PASSED. Check result file: $TEST_OUTPUT "
	fi
	set -e	
}


#
# rename temporary makefile and run it
#
function build_dist() 
{
	echo "[ build_dist ]"
	cd $_SRC
	make distribution || true

	# check that the distribution file has been  created
	DIST_FILE=$SANDBOX/distributions/$DIST_NAME
	if [ ! -f $DIST_FILE ] 
	then 
		echo "Destination file has not been created: $DIST_FILE"
		exit 1
	fi

	# Move created package to distribution directory define by $DIST_DIR
	mkdir -p $DIST_BASE/$VERSION
	mv $DIST_FILE $DIST_BASE/$VERSION
	
	# Create a file containing the latest version number 
	echo $VERSION > $DIST_BASE/.version

}

#
# Upload all packages 
#
function upload() 
{
	echo "[ upload_distribution ]"

	scp -B -2 -r -i $WORKSPACE/tcoffee/build/tcoffee_org_id $DIST_BASE $DIST_HOST
}


#
# Compile T-Coffee distribution
# - distribution sources are located at $UNTARED path
# - target binaries will be located at $TCDIR 
#
function build_binaries()
{
	echo "[ build_binaries ]"

	rm -rf $TCDIR
	mkdir -p $TCDIR
	mkdir -p $TCDIR/bin

	# Make sure that does not exist already binaries
	rm -rf $UNTARED/bin
	mkdir -p $UNTARED/bin

	# create t-coffee binaries installation
	cd $UNTARED
	./install all -tclinkdb=./tclinkdb.txt -repo=$BUILD_REPO -tcdir=$TCDIR -exec=$TCDIR/bin || true
    
	# Check that the binary has successfully compiled 
	if [ ! -f $TCDIR/bin/t_coffee ] 
	then 
		echo "Target 't_coffee' binary has not been compiled"
		exit 1
	fi    
    
    
	# add perl modules 
	#cp -r $PERLM/lib/perl5/ $TCDIR/perl
	mkdir -p $TCDIR/perl
	cp $WORKSPACE/tcoffee/build/cpanm  $TCDIR/perl	
        chmod +x $TCDIR/perl/cpanm

	# add gfortran libraries
	if [ $OSNAME == "macosx" ] 
	then 
		mkdir -p $TCDIR/gfortran
		cp /usr/local/gfortran/lib/libgfortran.3.dylib $TCDIR/gfortran
		cp /usr/local/gfortran/lib/libgfortran.a $TCDIR/gfortran
		cp /usr/local/gfortran/lib/libgfortran.la $TCDIR/gfortran
	fi
	
	#
	# add extra pack packages 
	#
	# + hmmtop
	cp $BUILD_REPO/hmmtop/2.1/$OSNAME-$OSARCH/* $TCDIR/plugins/$OSNAME

	# + secondary_struc.py (by Carsten)
	cp $WORKSPACE/tcoffee/build/extra/secondary_struc.py $TCDIR/plugins/$OSNAME/secondary_struc.py
	chmod +x $TCDIR/plugins/$OSNAME/secondary_struc.py
	
	# HH_seach
	cp $WORKSPACE/tcoffee/build/extra/HHsearch_1.5.1/$OSNAME-$OSARCH/* $TCDIR/plugins/$OSNAME/	
}



#
# download perl packages 
#
function build_perlm() {
	echo "[ build_perlm ]"

	chmod +x $WORKSPACE/tcoffee/build/cpanm
	$WORKSPACE/tcoffee/build/cpanm -q -n -l $PERLM -L $PERLM SOAP::Lite XML::Simple LWP

}

#
# create platform independent paltform distribution
#
function pack_binaries() {
	echo "[ pack_binaries ]"
	echo Package name: $INST_NAME 

	# remove the t_coffee binaries from 'plugins' folder 
	# it have to exist in 'bin' folder  
	rm -rf $TCDIR/plugins/$OSNAME/t_coffee

    # copy the sources 
    rm -rf $TCDIR/src
    cp -r $UNTARED/t_coffee_source $TCDIR/src
    find $TCDIR/src -name '*.o_*' | xargs rm -rf {} 
    find $TCDIR/src -name 't_coffee' | xargs rm -rf {} 
    find $TCDIR/src/ -type f -exec chmod 644 '{}' \;

	# invoke the install builder 
	mkdir -p $DIST_DIR
	"$INSTALLER" build $WORKSPACE/tcoffee/build/tcoffee-installer.xml --setvars product_version=$VERSION untared=$UNTARED osname=$OSNAME tcdir=$TCDIR outdir=$DIST_DIR outname=$INST_NAME
	
	# mac osx specific step 
	if [ $OSNAME == "macosx" ]
	then
	$WORKSPACE/tcoffee/build/mkdmg.sh $DIST_DIR/$INST_NAME.app
	mv $DIST_DIR/$INST_NAME.app.dmg $DIST_DIR/$INST_NAME.dmg 
	rm -rf $DIST_DIR/$INST_NAME.app
	fi
	
	# add execution attribute to the generated binary 
	if [ $OSNAME == "linux" ]
	then
	mv $DIST_DIR/$INST_NAME.run $DIST_DIR/$INST_NAME.bin
	chmod u+x $DIST_DIR/$INST_NAME.bin
	fi
	
}

function pack_tarball() {
	
    # creates a tar with all precompiled binaries
    if [[ $OSARCH == "x64" && $OSNAME == 'linux' ]]
    then
      echo "[ pack_tarball ]"

      cd $SANDBOX 
      $DIST_DIR/$INST_NAME.bin --mode unattended --user_email tcoffee.msa@gmail.com --prefix $PWD/$INST_NAME
      tar -zcvf $INST_NAME.tar.gz $INST_NAME/
      mv $INST_NAME.tar.gz $DIST_DIR/
      rm -rf $PWD/$INST_NAME
      cd -
    fi
}



function build_and_pack_stable() {
	echo "[ build_and_pack_stable ]"
		
	export CFLAGS="-O3"
	export INST_NAME=T-COFFEE_installer_"$VERSION"_"$OSNAME"_"$OSARCH"
	
	build_binaries	
	pack_binaries
	pack_tarball
} 

function build_and_pack_debug() {
	echo "[ build_and_pack_debug ]"
	
	export CFLAGS="-g -O1"
	export INST_NAME=T-COFFEE_installer_"$VERSION"_"$OSNAME"_"$OSARCH"_debug
	export DATE="`date +"%Y-%m-%d %H:%M:%S"` - Revision $GIT_REVISION - DEBUG"
	
	build_binaries	
	pack_binaries
} 

#
# Copy to dropbox for publishing
#

function copy_to_dropbox() {
	echo "[ copy_to_dropbox ]"

	mkdir -p $DROPBOX

	if [ $RELEASE == 1 ]; then 
	DROPBOX_DIR=$DROPBOX/stable/"$OSNAME"_"$OSARCH"
	else 
	DROPBOX_DIR=$DROPBOX/beta/"$OSNAME"_"$OSARCH"
	fi
	
	rm -rf $DROPBOX_DIR
	mkdir -p $DROPBOX_DIR 
	cp -r $TCDIR/* $DROPBOX_DIR

} 


#
#
# Execute all T-coffee core tasks (no server related)
#
function tcoffee() {
	echo "[ tcoffee ]"

	env
	clean
	build_dist
	
	build_and_pack_stable
	build_and_pack_debug

	if [ $DO_TEST == 1 ]; then
	doc_test 
	fi
	
	copy_to_dropbox
} 


#
# when at least a parameter is specified they are invoked as function call
#
if [ $# -gt 0 ] 
then
	while [ "$*" != "" ]
	do
		echo "Target: $1"
		$1
		shift
	done
else
    echo "Usage: build <target>"
    exit 1
fi



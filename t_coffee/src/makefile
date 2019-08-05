#
# T-Coffee legacy make file -- used mainly by the distibution and packaging procedure 
#
PROGRAM=T-COFFEE
EXECUTABLE=t_coffee

AUTHOR=Cedric Notredame 
DATE=$(shell date '+%d %b %Y - %H:%M')
DISTRIBUTION_ADDRESS=www.tcoffee.org/Packages
INSTITUTE=Centro de Regulacio Genomica
EMAIL=cedric.notredame@gmail.com
SHELL := /bin/bash

#
# make it possibile to define HOME2 externally using environment
# variable named 'SANDBOX', otherwise fallback on the defult mechanism
#
ifndef HOME2
HOME2=$(HOME)/Dropbox
endif
#### These paths related to the publishing on the production machine: they can stay relative
WEB_DIR=$(HOME2)/dropBoxShared/Work/tcoffee.org/public_html/
WEB_DOC=$(WEB_DIR)/Projects/tcoffee/documentation
WEB_BIN_DIR=$(HOME2)/dropBoxShared/Work/tcoffee.org/public_html/Packages/Binaries/tcoffee/$(OSNAME)/
WEB_ARCHIVE_SOURCES=$(WEB_DIR)/Packages/Archives/

###########################################################
#                                                         #
#       COMPILERS AND FLAGS                               #
#                                                         #
#                                                         #
#                                                         #
# #########################################################

CC=g++
ifndef CFLAGS
#CFLAGS=-O3 -Wno-write-strings
CFLAGS=-ggdb -w
endif

ifdef OPENMP
    KM_FLAGS=-fopenmp
endif



##### These paths must be defined relatively to the src directory
ifndef LIB_DIR
#LIB_DIR=/Users/cnotredame/Dropbox/projects/git/tcoffee/lib/
LIB_DIR=$(CURDIR)/../../lib
endif



ifndef GIT_DIR
GIT_DIR=$(CURDIR)/../../
endif

ifndef GIT_BIN
GIT_BIN=$(GIT_DIR)/binaries
endif

ifndef EXECUTABLE_DIR
EXECUTABLE_DIR=$(CURDIR)
endif

ifndef RST_DOC
RST_DOC=$(GIT_DIR)/docs
endif

ifndef RST_DATA
RST_DATA=$(GIT_DIR)/docs/.data/
endif

ifndef HTML_DOC
HTML_DOC=$(RST_DOC)/.html
endif


ifndef BUILD_DIR
BUILD_DIR=$(GIT_DIR)/build
endif

ifndef LIB_DIR
LIB_DIR=$(GIT_DIR)/lib
endif

ifndef VER_DIR
VER_DIR=$(LIB_DIR)/version
endif

ifndef CIRCLECI_CONFIG
CIRCLECI_CONFIG=$(GIT_DIR)/.circleci/config.yml
endif

PERL_LIB_DIR=$(LIB_DIR)/perl/lib
KO_LIB=$(LIB_DIR)/data_headers

#LIST OF LICENSES
LICENSES = $(LIB_DIR)/licenses

#
# Base path for library folders
#

ifdef NORM
NORM_2_USE=util_dp_gotoh_nw_2.c
endif

ifndef NORM
NORM_2_USE=util_dp_gotoh_nw.c
endif


ifndef USER_BIN
USER_BIN=$(HOME)/bin/
endif

CODE_LIST=$(LIB_DIR)/compilation/code_list


#LIST OF PERL MODULES
SCRIPT_LIB=$(PERL_LIB_DIR)/scripts
PERL4MAKEFILE_LIB=$(PERL_LIB_DIR)/perl4makefile

##### End of PATH definition


#
# READ VERSION NUMBER as minatined and updated 
#
#




#ifndef VERSION
VERSION=$(shell $(PERL4MAKEFILE_LIB)/read_program_version.pl -path=$(VER_DIR) -print)
#endif


DISTRIBUTION_FILE=$(PROGRAM)_distribution_$(VERSION)



LATEST_DISTRIBUTION_FILE=$(PROGRAM)_distribution
DISTRIBUTIONS=$(HOME2)/distributions/
UNTARED_DISTRIBUTIONS=$(HOME2)/untared_distributions/
LATEST_DISTRIBUTIONS=$(HOME2)/latest_distributions/
DISTRIB_DIR=$(HOME2)/distributions/$(PROGRAM)_distribution_$(VERSION)
DISTRIB_DOC=$(DISTRIB_DIR)/doc


#########################################################
#                                                       #
#       DISTRIBUTION                                    #
#                                                       #
#                                                       #
#                                                       #
#########################################################

distribution:	
	echo $(VERSION)
	make  VERSION=$(VERSION) CC=cc general_distribution
	make  VERSION=$(VERSION) CC=cc examples
	make  VERSION=$(VERSION) CC=cc mcoffee
	make  doc	
	make  VERSION=$(VERSION) CC=cc packing

doc:
	cp $(RST_DOC)/README  $(DISTRIB_DIR)
	rm -f $(DISTRIB_DIR)/doc/html
	mkdir -p $(DISTRIB_DIR)/doc/html/	
	cp -r $(HTML_DOC)/* $(DISTRIB_DIR)/doc/html/
	mkdir -p $(DISTRIB_DIR)/doc/misc/
	cp $(RST_DOC)/*.doc $(DISTRIB_DIR)/doc/misc/
	cp $(RST_DOC)/*.pdf $(DISTRIB_DIR)/doc/misc/

mcoffee:
	-rm -rf $(DISTRIB_DIR)/mcoffee
	mkdir -p $(DISTRIB_DIR)/mcoffee
	cp -r $(LIB_DIR)/mcoffee/* $(DISTRIB_DIR)/mcoffee

examples:
	mkdir -p $(DISTRIB_DIR)/example/
	mkdir -p $(DISTRIB_DIR)/example/
	cp -r $(RST_DATA)/*  $(DISTRIB_DIR)/example/

test:
	mkdir -p $(DISTRIB_DIR)/test
	cp $(PWD)/../test/test.pep $(DISTRIB_DIR)/test/test.pep
	$(EXECUTABLE)  $(DISTRIB_DIR)/test/test.pep -in=Mfast_pair -outfile=$(DISTRIB_DIR)/test/reference_test.aln -outorder=input
	echo '#!/usr/bin/env perl'    > $(DISTRIB_DIR)/bin/test.pl
	echo 'chdir "./test";' >> $(DISTRIB_DIR)/bin/test.pl
	echo '$$r=system "../bin/t_coffee test.pep -method fast_pair -outfile=new_reference_test.aln -outorder=input>/dev/null";' >> $(DISTRIB_DIR)/bin/test.pl
	echo 'if ($$r==0){print "\\nInstallation of $(EXECUTABLE) Successful\\n";}'>>$(DISTRIB_DIR)/bin/test.pl
	echo 'else {print "\\nInstallation of $(EXECUTABLE) Not Successful\\n";}'>>$(DISTRIB_DIR)/bin/test.pl
	echo 'chdir "/..";'
	chmod u+x $(DISTRIB_DIR)/bin/test.pl

debug: t_coffee
	rm -fr debug
	mkdir debug
	cp *.c debug
	cp *.h debug
	mv debug $(HOME2)/debug

env:
	@echo 'Home             :' $(HOME)
	@echo 'Pwd              :' $(PWD)
	@echo 'VERSION          :' $(VERSION)
	@echo 'DATE             :' $(DATE)
	@echo 'OSNAME           :' $(OSNAME)
	@echo 'OSARCH           :' $(OSARCH)
	@echo 'HOME2            :' $(HOME2)
	@echo 'CFLAGS           :' $(CFLAGS)
	@echo 'KM_FLAGS         :' $(KM_FLAGS)
	@echo 'DISTRIB_DIR      :' $(DISTRIB_DIR)
	@echo 'DISTRIB_DOC      :' $(DISTRIB_DOC)
	@echo 'EXECUTABLE       :' $(EXECUTABLE)
	@echo 'USER_BIN         :' $(USER_BIN)
	@echo 'LIB_DIR          :' $(LIB_DIR)
	@echo 'DISTRIBUTIONS    :' $(DISTRIBUTIONS)


general_distribution:
	-rm -f $(CODE_LIST)
	mkdir -p $(DISTRIB_DIR)
	mkdir -p $(DISTRIB_DIR)/bin
	mkdir -p $(DISTRIB_DIR)/bin/$(OSNAME)
	mkdir -p $(DISTRIB_DIR)/bin/macosx	
	mkdir -p $(DISTRIB_DOC)
	mkdir -p $(DISTRIB_DIR)/$(EXECUTABLE)_source
	mkdir -p $(DISTRIB_DIR)/reference_test
	echo '' > $(DISTRIB_DIR)/install
	make -i VERSION=$(VERSION) CC=cc code
	cp $(EXECUTABLE) $(DISTRIB_DIR)/bin/$(OSNAME)
	cp $(EXECUTABLE) $(GIT_BIN)/$(OSNAME)/$(EXECUTABLE).$(VERSION)
	make -i VERSION=$(VERSION) CC=cc license

library_distribution:
	cd $(LIB_DIR)/..
	tar -cvf lib.tar $(LIB_DIR)
	gzip lib.tar
	mkdir $(DISTRIBUTIONS)
	mv lib.tar.gz $(DISTRIBUTIONS)/lib_$(VERSION).tar.gz

empty:
	-rm -f *.o_$(OSNAME)
	tar -cf old_c.tar *.c *.h
	rm -f *.c
	rm -f *.h
	echo x > core
	rm -f core

distribute: all
	echo x > core
	rm -f core
	-rm -f *.o_$(OSNAME)
code:
	echo $(OSNAME) Compilation
	make clean
	make CC="g++" $(EXECUTABLE)
	mkdir -p $(DISTRIB_DIR)/$(EXECUTABLE)_source
	make copy_sources
	$(PERL4MAKEFILE_LIB)/tclinkdb2header.pl -infile=$(SCRIPT_LIB)/install.pl -mode=Perl -db=$(KO_LIB)/tclinkdb.txt >$(DISTRIB_DIR)/install
	cp $(KO_LIB)/tclinkdb.txt $(DISTRIB_DIR)
	chmod u+x $(DISTRIB_DIR)/install
	make install2version VERSION=$(VERSION)

install2version:
	echo $(DISTRIB_DIR)/install
	$(shell $(PERL4MAKEFILE_LIB)/edit_version.pl -in $(DISTRIB_DIR)/install -tag "#_#UPDATE_VERSION" -string 'our $VERSION=\"$(VERSION))\";'
	eeetete
	$(PERL4MAKEFILE_LIB)/edit_version.pl -in $(BUILD_DIR)/web-readme.txt -tag "#_#UPDATE_VERSION" -string "$(VERSION)"	
 
# This target copies T-Coffee sources files from the library folder 
# to the distribution folder. 
# This process is a bit tricky to maintaning compatibiltiy with the legacy 
# build process. It uses two temporary files 'code_list' and 'foldes'.
# 
# The file "code_list" is generated during the code compilation and contains
# the list of all C sources requied to compile T-Coffee. The target below 
# copies all the contained files to the target distribution directory. 
# 
# During this copy it is created a file "folders" with all the 
# sub-directories containing the used C sources. 
# 
# Only the header files "*.h" in that directories will be copied to the distrubtion folder�
# 
copy_sources: 
	-rm -rf folders
	for it in `cat code_list`; do \
	  dir=$(LIB_DIR); \
	  entry=$${it/$$dir/}; \
	  path=`dirname $$entry`; \
	  file=`basename $$entry`; \
	  mkdir -p $(DISTRIB_DIR)/$(EXECUTABLE)_source/$$path; \
	  cp $$it $(DISTRIB_DIR)/$(EXECUTABLE)_source/$$path; \
	  echo "$$path" >> folders; \
	done

# Always add the current and 'data_headers' folders
	echo "." >> folders
	echo "data_headers" >> folders
# Add the 'header' files
	for path in `cat folders | uniq`; do \
	  mkdir -p $(DISTRIB_DIR)/$(EXECUTABLE)_source/$$path; \
	  cp $(LIB_DIR)/$$path/*.h $(DISTRIB_DIR)/$(EXECUTABLE)_source/$$path; \
	done
# Copy these headers that are generated/upodated automatically by the build process
	cp ./perl_header_lib.h $(DISTRIB_DIR)/$(EXECUTABLE)_source
	cp ./programmes_define.h $(DISTRIB_DIR)/$(EXECUTABLE)_source
	cp ./define_header.h $(DISTRIB_DIR)/$(EXECUTABLE)_source

# Some missing stuff  	
	cp $(LIB_DIR)/fortran/TMalign.f $(DISTRIB_DIR)/$(EXECUTABLE)_source
	cp $(LIB_DIR)/c_make_lib/dist_makefile $(DISTRIB_DIR)/$(EXECUTABLE)_source/makefile


license:
	$(PERL4MAKEFILE_LIB)/make_license.pl -file $(LICENSES)/academic_license.txt -email "$(EMAIL)" -type txt -institute "$(INSTITUTE)" -author "$(AUTHOR)" -date "$(DATE)" -program "$(PROGRAM)" > $(DISTRIB_DIR)/license.txt
	$(PERL4MAKEFILE_LIB)/make_license.pl -file $(LICENSES)/short_license_utf8.txt -email "$(EMAIL)" -type c -institute "$(INSTITUTE)" -author "$(AUTHOR)" -date "$(DATE)" -program "$(PROGRAM)" > $(DISTRIB_DIR)/short_license.c
#$(PERL4MAKEFILE_LIB)/make_license.pl -file $(LICENSES)/short_license.txt -email "$(EMAIL)" -type perl -institute "$(INSTITUTE)" -author "$(AUTHOR)" -date "$(DATE)" -program "$(PROGRAM)" > $(DISTRIB_DIR)/short_license.pl
	for file in `find $(DISTRIB_DIR)/$(EXECUTABLE)_source \( -name '*.h' -o -name '*.c' \)` ; do \
	  cat $(DISTRIB_DIR)/short_license.c $$file > $$file.new; \
	  mv $$file.new $$file; \
	done
	rm -f $(DISTRIB_DIR)/short_license.c

#install:
#   $(PERL4MAKEFILE_LIB)/tclinkdb2header.pl -infile=$(SCRIPT_LIB)/install.pl -mode=Perl -db=$(KO_LIB)/tclinkdb.txt >  install

packing:
	mkdir -p $(UNTARED_DISTRIBUTIONS)
	rm -rf $(UNTARED_DISTRIBUTIONS)
	mkdir -p $(UNTARED_DISTRIBUTIONS)	
	mkdir -p $(DISTRIBUTIONS)
	mkdir -p $(LATEST_DISTRIBUTIONS)
	mv $(DISTRIB_DIR) .
	tar -cvf $(DISTRIBUTION_FILE).tar $(DISTRIBUTION_FILE)  
	gzip $(DISTRIBUTION_FILE).tar
	mv $(DISTRIBUTION_FILE)        $(UNTARED_DISTRIBUTIONS)
	cp $(DISTRIBUTION_FILE).tar.gz $(DISTRIBUTIONS)
	mv $(DISTRIBUTION_FILE).tar.gz $(LATEST_DISTRIBUTIONS)/$(LATEST_DISTRIBUTION_FILE).tar.gz




version:
	$(shell $(PERL4MAKEFILE_LIB)/read_program_version.pl -path=$(VER_DIR) -$(release_type) -circleci=$(CIRCLECI_CONFIG))


release:
	make version release_type=$(type)
	make release2

release2:
	VERSION=$(shell $(PERL4MAKEFILE_LIB)/read_program_version.pl -path=$(VER_DIR) -print)	
	$(PERL4MAKEFILE_LIB)/edit_version.pl -in $(RST_DOC)/conf.py -tag "#_#UPDATE_VERSION_TAG1" -string "version = u'$(VERSION)'"
	$(PERL4MAKEFILE_LIB)/edit_version.pl -in $(RST_DOC)/conf.py -tag "#_#UPDATE_VERSION_TAG2" -string "release = u'$(VERSION)'"	
	sphinx-build -b singlehtml $(RST_DOC) $(HTML_DOC)
	make distribution VERSION=$(VERSION)
	make web
	git commit -a -m '$(m)'	
#	git push

web:
	rm -f $(WEB_DOC)/*.html
	cp -r $(HTML_DOC)/* $(WEB_DOC)
	cp $(EXECUTABLE_DIR)/$(EXECUTABLE) $(WEB_BIN_DIR)/$(EXECUTABLE).$(VERSION);
	mkdir -p $(WEB_ARCHIVE_SOURCES)
	cp $(LATEST_DISTRIBUTIONS)/$(LATEST_DISTRIBUTION_FILE).tar.gz $(WEB_ARCHIVE_SOURCES)/$(PROGRAM)_distribution_$(shell $(PERL4MAKEFILE_LIB)/read_program_version.pl -path=$(VER_DIR) -print).tar.gz
	cp $(LIB_DIR)/version/version_number.version  $(WEB_DIR)/Packages/T-COFFEE.version

#
# Include the main T-Coffee makefile
#
include ../../compile/makefile 

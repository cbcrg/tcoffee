PROGRAM="T-COFFEE"
EXECUTABLE=t_coffee
AUTHOR= "Cedric Notredame "

#
# make it possibile to define HOME2 externally using environment
# variable named 'SANDBOX', otherwise fallback on the defult mechanism
#
ifndef HOME2
HOME2=$(HOME)/Dropbox
endif


#########################################################
#							#
#		COMPILERS AND FLAGS			#
#							#
#							#
#							#
#########################################################

ifndef CFLAGS
CFLAGS=-O3
endif


CC=gcc


SOURCE_DIR=$(PWD)/../../$(EXECUTABLE)/lib
SCRIPT_DIR=$(PWD)/../../$(EXECUTABLE)/scripts
#########################################################
#							#
#		LIBRARY          			#
#							#
#							#
#							#
#########################################################
include $(PWD)/../../lib/c_make_lib/lib_makefile



#########################################################
#							#
#		EXECUTABLES			#
#							#
#							#
#########################################################

HEADERS	   =$(UT_LIB_HEADER) $(IO_LIB_HEADER)   $(DP_LIB_HEADER) $(DEV1_LIB_HEADER) $(DEV2_LIB_HEADER) $(DEV3_LIB_HEADER) $(DEV4_LIB_HEADER) $(FASTAL_LIB_HEADER) $(T_COFFEE_LIB_HEADER)
OBJECTS	   =$(UT_LIB_OBJECTS) $(IO_LIB_OBJECTS) $(DP_LIB_OBJECTS) $(DEV1_LIB_OBJECTS) $(DEV2_LIB_OBJECTS) $(DEV3_LIB_OBJECTS) $(DEV4_LIB_OBJECTS) $(FASTAL_LIB_OBJECTS) $(T_COFFEE_LIB_OBJECTS)


$(EXECUTABLE):  $(OBJECTS) $(USE_MEM_LIB) $(FORTRAN_LIB_OBJECTS)
	echo Platform: $(OSNAME)
	$(CC) $(CFLAGS) -o $(EXECUTABLE)  $(OBJECTS) $(USE_MEM_LIB) -lm
# 	cp $(EXECUTABLE) $(USER_BIN)/t_coffee

$(EXECUTABLE).test:  $(OBJECTS) $(USE_MEM_LIB)

all: $(EXECUTABLE)
#########################################################
#							#
#		DISTRIBUTION		#
#							#
#							#
#							#
#########################################################
distribution:
	-rm *.c
	-rm *.h
	make -i VERSION=$(VERSION) CC=cc general_distribution
	make -i VERSION=$(VERSION) CC=cc examples
	make -i VERSION=$(VERSION) CC=cc mcoffee
	make -i VERSION=$(VERSION) CC=cc packing

code:
	make CC=cc c_code
	chmod u+x $(DISTRIB_DIR)/install

doc:
	cp $(PWD)/../../$(EXECUTABLE)/doc/$(EXECUTABLE)_tutorial.* $(DISTRIB_DOC)
	antiword -m  "UTF-8" $(PWD)/../../$(EXECUTABLE)/doc/$(EXECUTABLE)_tutorial.doc >$(DISTRIB_DOC)/$(EXECUTABLE)_tutorial.txt
	cp $(PWD)/../../$(EXECUTABLE)/doc/$(EXECUTABLE)_technical.* $(DISTRIB_DOC)
	antiword -m "UTF-8" $(PWD)/../../$(EXECUTABLE)/doc/$(EXECUTABLE)_technical.doc >$(DISTRIB_DOC)/$(EXECUTABLE)_technical.txt
	cp $(PWD)/../../$(EXECUTABLE)/doc/README $(DISTRIB_DIR)
	cp $(PWD)/../../$(EXECUTABLE)/doc/README $(DISTRIB_DIR)/doc/
	cp $(PWD)/../../$(EXECUTABLE)/doc/aln_compare.doc.txt $(DISTRIB_DOC)
	cp $(PWD)/../../$(EXECUTABLE)/doc/seq_reformat.doc.txt $(DISTRIB_DOC)
	cp $(PWD)/../../$(EXECUTABLE)/doc/$(EXECUTABLE).pdf $(DISTRIB_DOC)

mcoffee:
	-rm -rf $(DISTRIB_DIR)/mcoffee
	mkdir -p $(DISTRIB_DIR)/mcoffee
	cp $(LIBRARY)/mcoffee/* $(DISTRIB_DIR)/mcoffee

examples:
	mkdir -p $(DISTRIB_DIR)/example/
	cp $(PWD)/../../$(EXECUTABLE)/doc_test/data/*  $(DISTRIB_DIR)/example/

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
	#echo "./bin/test.pl" >>$(DISTRIB_DIR)/install
	#echo "./bin/t_coffee -update" >>$(DISTRIB_DIR)/install

debug:  t_coffee
	rm -r debug
	mkdir debug
	cp *.c debug
	cp *.h debug
	mv debug $(HOME2)/debug

env:
	@echo 'Home             :' $(HOME)
	@echo 'Pwd              :' $(PWD)
	@echo 'Version          :' $(VERSION)
	@echo 'OSNAME           :' $(OSNAME)
	@echo 'OSARCH           :' $(OSARC)
	@echo 'HOME2            :' $(HOME2)
	@echo 'CFLAGS           :' $(CFLAGS)
	@echo 'DISTRIB_DIR      :' $(DISTRIB_DIR)
	@echo 'DISTRIB_DOC      :' $(DISTRIB_DOC)
	@echo 'EXECUTABLE       :' $(EXECUTABLE)
	@echo 'USER_BIN         :' $(USER_BIN)
	@echo 'DOCUMENTATION    :' $(DOCUMENTATION)
	@echo 'LIBRARY          :' $(LIBRARY)
	@echo 'DISTRIBUTIONS    :' $(DISTRIBUTIONS)


#########################################################
#							#
#		GENERIC  LIBRARY         		#
#							#
#							#
#########################################################
include $(PWD)/../../lib/c_make_lib/generic_makefile


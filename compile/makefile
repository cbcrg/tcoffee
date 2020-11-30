#
# Basic T-Coffee makefile
#
# This makefile depends on the following environment variables
# - LIB_DIR: the base path containing  T-Coffee library folders
# - CFLAGS: flags used to control the compilation process (i.e. optimization, debug symbols, etc)
# - DATE: the date to the be stamped in T-Coffee binary
# - VERSION: the version to be stamped in T-Coffee binary
#

#
# About Mac compilation flags see http://goo.gl/HLtS4
#

#
# Full Version number
#
ifndef BUILD_INFO
BUILD_INFO=$(VERSION)
endif

#
# Use current date as default if nothing is specified
#
ifndef DATE
DATE=$(shell date '+%d %b %Y - %H:%M')
endif

#
# Default program name
#
ifndef PROGRAM
PROGRAM=T-COFFEE
endif

#
# (not so much to explaing here .. :))
#
ifndef AUTHOR
AUTHOR=Cedric Notredame
endif

#
# Web path where T-Coffee distribution are published
#
ifndef DISTRIBUTION_ADDRESS
DISTRIBUTION_ADDRESS=www.tcoffee.org/Packages
endif



#LIST OF C LIBRARIES

PERL_LIB_DIR=$(LIB_DIR)/perl/lib
KO_LIB=$(LIB_DIR)/data_headers
IO_LIB=$(LIB_DIR)/io_lib
UT_LIB=$(LIB_DIR)/util_lib
DP_LIB=$(LIB_DIR)/dp_lib
MEM_LIB=$(LIB_DIR)/mem_lib
GA_LIB=$(LIB_DIR)/ga_lib
DEV1_LIB=$(LIB_DIR)/dev_lib
FASTAL_LIB=$(LIB_DIR)/fastal
KM_COFFEE_LIB=$(LIB_DIR)/km_coffee
T_COFFEE_LIB=$(LIB_DIR)/t_coffee_lib


#LIST OF PERL MODULES

FILTER_LIB=$(PERL_LIB_DIR)/filters
SCRIPT_LIB=$(PERL_LIB_DIR)/scripts
WEBSERVICES_LIB=$(PERL_LIB_DIR)/webservices
PERL4MAKEFILE_LIB=$(PERL_LIB_DIR)/perl4makefile
PERL4DATA_PROCESS_LIB=$(PERL_LIB_DIR)/perl4data_process


# ===================================================================
# 	MAKE sources search path. Include here all folder where C source
#   have to be located
# ===================================================================
VPATH=$(IO_LIB):$(UT_LIB):$(DP_LIB):$(MEM_LIB):$(GA_LIB):$(DEV1_LIB):$(FASTAL_LIB):$(KM_COFFEE_LIB):$(T_COFFEE_LIB)


#########################################################################
#				PERL LIB												#
#																		#
#																		#
#																		#
#########################################################################

#########################################################################
########################      	PERL4DATA_PROCESS_LIB   #################
#########################################################################
#1
PERL4DATA_PROCESS_OBJECTS=\
	rec_sum.pl\
	count.pl
#2
rec_sum.pl:\
	$(PERL4DATA_PROCESS_LIB)/rec_sum.pl
	cp $(PERL4DATA_PROCESS_LIB)/rec_sum.pl .

count.pl:\
	$(PERL4DATA_PROCESS_LIB)/count.pl
	cp $(PERL4DATA_PROCESS_LIB)/count.pl .


#########################################################################
######################   PERL4MAKEFILE_LIB   ############################
#########################################################################
#1
PERL4MAKEFILE_OBJECTS=\
	process_list.pl\
	make_license.pl\
	CCsed.script
#2

process_list.pl:\
	$(PERL4MAKEFILE_LIB)/process_list.pl
	cp $(PERL4MAKEFILE_LIB)/process_list.pl .
make_license.pl:\
	$(PERL4MAKEFILE_LIB)/make_license.pl
	cp $(PERL4MAKEFILE_LIB)/make_license.pl .


read_program_version:\
	$(PERL4MAKEFILE_LIB)/read_program_version.pl
	cp $(PERL4MAKEFILE_LIB)/read_program_version.pl .

append2file.pl:\
	$(PERL4MAKEFILE_LIB)/append2file.pl
	cp $(PERL4MAKEFILE_LIB)/append2file.pl .

CCsed.script:\
	$(PERL4MAKEFILE_LIB)/CCsed.script
	cp $(PERL4MAKEFILE_LIB)/CCsed.script .
#########################################################################
#########################  WEBSERVICES_LIB   ############################
#########################################################################
#1
dollar="'$$'"
quote="'"'"

WEBSERVICES_LIB_OBJECTS=\
	dalilite.pl\
	wublast.pl\
	blastpgp.pl\
	ncbiblast_lwp.pl\
	wublast_lwp.pl


#2
dalilite.pl:\
	$(WEBSERVICES_LIB)/dalilite.pl
	cp $(WEBSERVICES_LIB)/dalilite.pl .

wublast.pl:\
	$(WEBSERVICES_LIB)/wublast.pl
	cp $(WEBSERVICES_LIB)/wublast.pl .

blastpgp.pl:\
	$(WEBSERVICES_LIB)/blastpgp.pl
	cp $(WEBSERVICES_LIB)/blastpgp.pl .

wublast_lwp.pl:\
	$(WEBSERVICES_LIB)/wublast_lwp.pl
	cp $(WEBSERVICES_LIB)/wublast_lwp.pl .

ncbiblast_lwp.pl:\
	$(WEBSERVICES_LIB)/ncbiblast_lwp.pl
	cp $(WEBSERVICES_LIB)/ncbiblast_lwp.pl .

#########################################################################
########################      	SCRIPT_LIB   ############################
#########################################################################
#1
dollar="'$$'"
quote="'"'"

SCRIPT_LIB_OBJECTS=\
	msa2bootstrap.pl\
	tc_generic_method.pl\
	mmseqs2prf.pl\
	dynamic.pl\
	rnapdb2protpdb.pl\
	generic_method.tc_method\
	clustalw_method.tc_method\
	extract_from_pdb\
	install.pl\
	clean_cache.pl\
	nature_protocol.pl\
	mocca
#2
extract_from_pdb:\
	$(SCRIPT_LIB)/extract_from_pdb
	cp $(SCRIPT_LIB)/extract_from_pdb .
	echo "# `date`" >> extract_from_pdb
install.pl:\
	$(SCRIPT_LIB)/install.pl
	$(PERL4MAKEFILE_LIB)/tclinkdb2header.pl -infile=$(SCRIPT_LIB)/install.pl -mode=Perl -db=$(KO_LIB)/tclinkdb.txt >install.pl
	echo "# `date`" >> install.pl
msa2bootstrap.pl:\
	$(SCRIPT_LIB)/msa2bootstrap.pl
	cp $(SCRIPT_LIB)/msa2bootstrap.pl .
clean_cache.pl:\
	$(SCRIPT_LIB)/clean_cache.pl
	cp $(SCRIPT_LIB)/clean_cache.pl .
nature_protocol.pl:\
	$(SCRIPT_LIB)/nature_protocol.pl
	cp $(SCRIPT_LIB)/nature_protocol.pl .
clustalw_method.tc_method:\
	$(SCRIPT_LIB)/clustalw_method.tc_method
	cp $(SCRIPT_LIB)/clustalw_method.tc_method .
generic_method.tc_method:\
	$(SCRIPT_LIB)/generic_method.tc_method
	cp $(SCRIPT_LIB)/generic_method.tc_method .
tc_generic_method.pl:\
	$(SCRIPT_LIB)/tc_generic_method.pl
	cp $(SCRIPT_LIB)/tc_generic_method.pl .
	echo "my $(dollar)program=$(quote)$(PROGRAM) ($(VERSION))$(quote);">> tc_generic_method.pl
dynamic.pl:\
	$(SCRIPT_LIB)/dynamic.pl
	cp $(SCRIPT_LIB)/dynamic.pl .
	echo "my $(dollar)program=$(quote)$(PROGRAM) ($(VERSION))$(quote);">> dynamic.pl
mmseqs2prf.pl:\
	$(SCRIPT_LIB)/mmseqs2prf.pl
	cp $(SCRIPT_LIB)/mmseqs2prf.pl .
	echo "my $(dollar)program=$(quote)$(PROGRAM) ($(VERSION))$(quote);">> mmseqs2prf.pl
rnapdb2protpdb.pl:\
	$(SCRIPT_LIB)/rnapdb2protpdb.pl
	cp $(SCRIPT_LIB)/rnapdb2protpdb.pl .
t_coffee_dpa:\
	$(SCRIPT_LIB)/t_coffee_dpa
	cp $(SCRIPT_LIB)/t_coffee_dpa .
	echo "my $(dollar)program=$(quote)$(PROGRAM) ($(VERSION))$(quote);\n">> t_coffee_dpa
t_coffee_dpa2:\
	$(SCRIPT_LIB)/t_coffee_dpa2
	cp $(SCRIPT_LIB)/t_coffee_dpa2 .
mocca:\
	$(SCRIPT_LIB)/mocca
	cp $(SCRIPT_LIB)/mocca .

#########################################################################
########################      	FILTER_LIB   ############################
#########################################################################
#1
FILTER_LIB_OBJECTS=\
	RNAplfold2tclib.pl\
	fasta_seq2RNAplfold_templatefile.pl\
	fasta_seq2hmmtop_fasta.pl\
	fasta_seq2consan_aln.pl\
	clustalw_aln2fasta_aln.pl\
	seq2name_seq.pl\
	seq2intersection.pl\
	msf_aln2fasta_aln.pl\
	msa.pl\
	upp.pl\
	clustalo.pl\
	dca.pl\
	blast_aln2fasta_aln.pl\
	blast_xml2fasta_aln.pl\
	fasta_aln2fasta_aln_unique_name.pl\
	newick2name_list.pl\
	excel2fasta.pl\
	nameseq2fasta.pl\
	any_file2unix_file.pl

nameseq2fasta.pl:\
	$(FILTER_LIB)/nameseq2fasta.pl
	cp $(FILTER_LIB)/nameseq2fasta.pl .
	echo "# `date`" >> nameseq2fasta.pl

fasta_seq2hmmtop_fasta.pl:\
	$(FILTER_LIB)/fasta_seq2hmmtop_fasta.pl
	cp $(FILTER_LIB)/fasta_seq2hmmtop_fasta.pl .
	echo "# `date`" >> fasta_seq2hmmtop_fasta.pl

fasta_seq2consan_aln.pl:\
	$(FILTER_LIB)/fasta_seq2consan_aln.pl
	cp $(FILTER_LIB)/fasta_seq2consan_aln.pl .
	echo "# `date`" >> fasta_seq2consan_aln.pl
upp.pl:\
	$(FILTER_LIB)/upp.pl
	cp $(FILTER_LIB)/upp.pl .
	echo "# `date`" >> upp.pl
clustalo.pl:\
	$(FILTER_LIB)/clustalo.pl
	cp $(FILTER_LIB)/clustalo.pl .
	echo "# `date`" >> clustalo.pl
msa.pl:\
	$(FILTER_LIB)/msa.pl
	cp $(FILTER_LIB)/msa.pl .
	echo "# `date`" >> msa.pl
dca.pl:\
	$(FILTER_LIB)/dca.pl
	cp $(FILTER_LIB)/dca.pl .
	echo "# `date`" >> dca.pl

fasta_seq2RNAplfold_templatefile.pl:\
	$(FILTER_LIB)/fasta_seq2RNAplfold_templatefile.pl
	cp $(FILTER_LIB)/fasta_seq2RNAplfold_templatefile.pl .
	echo "# `date`" >> fasta_seq2RNAplfold_templatefile.pl

RNAplfold2tclib.pl:\
	$(FILTER_LIB)/RNAplfold2tclib.pl
	cp $(FILTER_LIB)/RNAplfold2tclib.pl .
	echo "# `date`" >> RNAplfold2tclib.pl

blast_aln2fasta_aln.pl:\
	$(FILTER_LIB)/blast_aln2fasta_aln.pl
	cp $(FILTER_LIB)/blast_aln2fasta_aln.pl .
	echo "# `date`" >> blast_aln2fasta_aln.pl
blast_xml2fasta_aln.pl:\
	$(FILTER_LIB)/blast_xml2fasta_aln.pl
	cp $(FILTER_LIB)/blast_xml2fasta_aln.pl .
	echo "# `date`" >> blast_xml2fasta_aln.pl
fasta_aln2fasta_aln_unique_name.pl:\
	$(FILTER_LIB)/fasta_aln2fasta_aln_unique_name.pl
	cp $(FILTER_LIB)/fasta_aln2fasta_aln_unique_name.pl .
	echo "# `date`" >> fasta_aln2fasta_aln_unique_name.pl
seq2name_seq.pl:\
	$(FILTER_LIB)/seq2name_seq.pl
	cp $(FILTER_LIB)/seq2name_seq.pl .
	echo "# `date`" >> seq2name_seq.pl
seq2intersection.pl:\
	$(FILTER_LIB)/seq2intersection.pl
	cp $(FILTER_LIB)/seq2intersection.pl .
	echo "# `date`" >> seq2intersection.pl
clustalw_aln2fasta_aln.pl:\
	$(FILTER_LIB)/clustalw_aln2fasta_aln.pl
	cp $(FILTER_LIB)/clustalw_aln2fasta_aln.pl .
	echo "# `date`" >> clustalw_aln2fasta_aln.pl
msf_aln2fasta_aln.pl:\
	$(FILTER_LIB)/msf_aln2fasta_aln.pl
	cp $(FILTER_LIB)/msf_aln2fasta_aln.pl .
	echo "# `date`" >> msf_aln2fasta_aln.pl
newick2name_list.pl:\
	$(FILTER_LIB)/newick2name_list.pl
	cp $(FILTER_LIB)/newick2name_list.pl .
	echo "# `date`" >> newick2name_list.pl
excel2fasta.pl:\
	$(FILTER_LIB)/excel2fasta.pl
	cp $(FILTER_LIB)/excel2fasta.pl .
	echo "# `date`" >> excel2fasta.pl
any_file2unix_file.pl:\
	$(FILTER_LIB)/any_file2unix_file.pl
	cp $(FILTER_LIB)/any_file2unix_file.pl .
	echo "# `date`" >> any_file2unix_file.pl

#########################################################################
########################      	PERL4C   ################################
#########################################################################
COMPLETE_PERL_OBJECTS=\
	$(PERL4DATA_PROCESS_OBJECTS)\
	$(PERL4MAKEFILE_OBJECTS)\
	$(SCRIPT_LIB_OBJECTS)\
	$(WEBSERVICES_LIB_OBJECTS)\
	$(FILTER_LIB_OBJECTS)
perl_header_lib.h:\
	$(COMPLETE_PERL_OBJECTS)
	$(FILTER_LIB)/perlfile2cheader.pl $(COMPLETE_PERL_OBJECTS) > perl_header_lib.h

#########################################################################
#				C LIB					#
#									#
#									#
#									#
#########################################################################

#########################################################################
########################      	KO_LIB      #############################
#########################################################################

define_header.h:\
	$(LIB_DIR)/coffee_defines.h \
	$(KO_LIB)/programmes_define.h \
	$(KO_LIB)/tclinkdb.txt
#
# Copy the 'define_header.h' file and udpate the compile info
#
	cp $(LIB_DIR)/define_header.h .
	sed -i -e 's/#define DATE ".*"/#define DATE "$(DATE)"/' define_header.h
	sed -i -e 's/#define AUTHOR ".*"/#define AUTHOR "$(AUTHOR)"/' define_header.h
	sed -i -e 's/#define VERSION ".*"/#define VERSION "$(VERSION)"/' define_header.h
	sed -i -e 's/#define PROGRAM ".*"/#define PROGRAM "$(PROGRAM)"/' define_header.h
	sed -i -e 's/#define BUILD_INFO ".*"/#define BUILD_INFO "$(BUILD_INFO)"/' define_header.h
	sed -i -e 's/#define DISTRIBUTION_ADDRESS ".*"/#define DISTRIBUTION_ADDRESS "$(subst /,\/,$(DISTRIBUTION_ADDRESS))"/' define_header.h


programmes_define.h:\
	$(KO_LIB)/tclinkdb.txt \
	$(PERL4MAKEFILE_LIB)/tclinkdb2header.pl \
	$(KO_LIB)/programmes_define.h
#
# Create in the current folder the file 'programmes_define.h'
#
	$(PERL4MAKEFILE_LIB)/tclinkdb2header.pl -infile=$(KO_LIB)/programmes_define.h -mode=C -db=$(KO_LIB)/tclinkdb.txt > programmes_define.h



#########################################################################
########################      	T-Coffee_LIB      #############################
#########################################################################

T_COFFEE_LIB_OBJECTS = t_coffee.o_$(OSNAME)


#########################################################################
########################      	DEV1_LIB      #############################
#########################################################################

DEV1_LIB_OBJECTS = dev1.o_$(OSNAME)



#########################################################################
########################      	FASTAL_LIB      #############################
#########################################################################


FASTAL_LIB_OBJECTS = \
		fastal.o_$(OSNAME)\
		parttree.o_$(OSNAME)\
		tree.o_$(OSNAME)\
		diagonal.o_$(OSNAME)\
		fastal_opt_parsing.o_$(OSNAME)\
		scoring.o_$(OSNAME)\
		iteration.o_$(OSNAME)

#########################################################################
########################      	IO_LIB      #############################
#########################################################################


IO_LIB_OBJECTS = \
	    io_func.o_$(OSNAME)\
		pb_util_read_sequence.o_$(OSNAME)\
		pb_util_read_seq_util.o_$(OSNAME)\
		tree_util.o_$(OSNAME)




#########################################################################
########################      	UTIL_LIB    #############################
#########################################################################

UT_LIB_OBJECTS = \
		util_constraints_list.o_$(OSNAME)\
		util_job_handling.o_$(OSNAME)\
		util_dps.o_$(OSNAME)\
		util_domain_constraints_list.o_$(OSNAME)\
		util_analyse_constraints_list.o_$(OSNAME)\
		util_aln_analyze.o_$(OSNAME)\
		aln_convertion_util.o_$(OSNAME)\
		util_declare.o_$(OSNAME)\
		hsearch.o_$(OSNAME)\
		random.o_$(OSNAME)\
		util_make_tree.o_$(OSNAME)\
		util.o_$(OSNAME)\
		reformat_struc.o_$(OSNAME)\
		reformat.o_$(OSNAME)\
		aln_compare.o_$(OSNAME)



#########################################################################
#####################      	KM-COFFEE_LIB      ##########################
#########################################################################

KM_COFFEE_OBJECTS = \
		Stack.o_$(OSNAME)\
		Vector.o_$(OSNAME)\
		classes.o_$(OSNAME)\
		km_util.o_$(OSNAME)\
		kmeans.o_$(OSNAME)\
		km_coffee.o_$(OSNAME)



#########################################################################
########################      	DP_LIB      #############################
#########################################################################

DP_LIB_OBJECTS = \
        util_graph_maln.o_$(OSNAME)\
		util_dp_clean_maln.o_$(OSNAME)\
		util_dp_ssec_pwaln.o_$(OSNAME)\
		util_dp_sim.o_$(OSNAME)\
		util_dp_mm_nw.o_$(OSNAME)\
		util_dp_gotoh_nw.o_$(OSNAME)\
		util_dp_suboptimal_nw.o_$(OSNAME)\
		util_dp_cdna_fasta_nw.o_$(OSNAME)\
		util_dp_generic_fasta_nw.o_$(OSNAME)\
		util_dp_fasta_nw.o_$(OSNAME)\
		util_dp_fasta_sw.o_$(OSNAME)\
		util_dp_gotoh_sw.o_$(OSNAME)\
		util_dp_est.o_$(OSNAME)\
		util_domain_dp_drivers.o_$(OSNAME)\
		util_dp_drivers.o_$(OSNAME)\
		util_domain_dp.o_$(OSNAME)\
		CUSTOM_evaluate_for_struc.o_$(OSNAME)\
		evaluate_for_struc.o_$(OSNAME)\
		evaluate_for_domain.o_$(OSNAME)\
		evaluate_dirichlet.o_$(OSNAME)\
		evaluate.o_$(OSNAME)\
		phylo3d.o_$(OSNAME)\
		showpair.o_$(OSNAME)\
		fsa_dp.o_$(OSNAME)\
		pavie_dp.o_$(OSNAME)


#
# Note: The following three header files are updated/created on-fly
# in the current path by the makefile
#
HEADERS	= define_header.h programmes_define.h perl_header_lib.h

OBJECTS = $(UT_LIB_OBJECTS) $(IO_LIB_OBJECTS) $(DP_LIB_OBJECTS) $(DEV1_LIB_OBJECTS) $(FASTAL_LIB_OBJECTS) $(KM_COFFEE_OBJECTS) $(T_COFFEE_LIB_OBJECTS)

# ===================================================================
# 	ALL target: create the main application file
# ===================================================================
all: t_coffee

# ===================================================================
# 	T_COFFEE target: build T-Coffee binary artifact
# ===================================================================
t_coffee: $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS) $(KM_FLAGS) -lm -o t_coffee $(OBJECTS)
# 	cp t_coffee $(HOME)/bin

# ===================================================================
# Automatically manage headers dependencies
#
# Please note: the current path incusion '-I.' must be the fist to give
# so that the header files created in the current folder will
# override the default definitions
# ===================================================================

%.d: %.h
	$(CC) $(CFLAGS) $(KM_FLAGS) -MM -I. -I$(LIB_DIR) $< -MF $@

# ===================================================================
# 	Generic compilation rule
# ===================================================================

%.o_$(OSNAME): %.c %.d
	echo "$<" >> code_list
	$(CC) $(CFLAGS) $(KM_FLAGS) -I. -I$(LIB_DIR) -c $< -o $@

# ===================================================================
# Clean stuff
# ===================================================================

.PHONY: clean

clean:
	rm -f *.o_*
	rm -f *.script
	rm -f *.pl
	rm -f *.h
	rm -f *.tc_method
	rm -f code_list
	rm -f extract_from_pdb
	rm -rf folders




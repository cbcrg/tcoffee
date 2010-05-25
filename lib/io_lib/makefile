CC=cc
#CC= /home/cnotreda/pure/purify-4.2-hpux/purify gcc -g
OSF_FLAGS=
SGI_FLAGS=


CFLAGS= -O2
#CFLAGS=  

UTIL=../io_lib/io_func.o\
../io_lib/pb_util_read_seq_util.o\
../io_lib/pb_util_read_sequence.o\
../util_lib/util.o\
../util_lib/util_declare.o\
../util_lib/random.o 

SEQ2MAT=../aln_lib/seq2mat.o 
PWCW2MAT= ../util_lib/util_constraints_list.o ../aln_lib/pwcw2mat.o\
../util_lib/aln_convertion_util.o 
MAT2EMAT=../util_lib/util_constraints_list.o ../aln_lib/mat2emat.o
MAT2ALN=../util_lib/util_constraints_list.o ../util_lib/util_make_tree.o \
../io_lib/tree_util.o ../aln_lib/mat2aln.o 
MAT2NDP_ALN=../util_lib/util_constraints_list.o ../aln_lib/mat2ndp_aln.o 
ALN2MAT=../util_lib/util_constraints_list.o \
../util_lib/aln_convertion_util.o ../aln_lib/aln2mat.o 

SEQ2MAT_OBJECTS    =$(UTIL)$(SEQ2MAT)
MAT2EMAT_OBJECTS   =$(UTIL)$(MAT2EMAT)
MAT2ALN_OBJECTS    =$(UTIL)$(MAT2ALN)
PWCW2MAT_OBJECTS   =$(UTIL)$(PWCW2MAT)
MAT2NDP_ALN_OBJECTS=$(UTIL)$(MAT2NDP_ALN)
ALN2MAT_OBJECTS    =$(UTIL)$(ALN2MAT)
EXEC= mat2emat mat2aln pwcw2mat mat2ndp_aln aln2mat
mat2emat:  $(MAT2EMAT_OBJECTS)
	   $(CC)  $(CFLAGS) -o mat2emat $(MAT2EMAT_OBJECTS) -lm

mat2aln:  $(MAT2ALN_OBJECTS)
	    $(CC)  $(CFLAGS) -o mat2aln $(MAT2ALN_OBJECTS) -lm

mat2ndp_aln: $(MAT2NDP_ALN_OBJECTS)
	    $(CC)  $(CFLAGS) -o mat2ndp_aln $(MAT2NDP_ALN_OBJECTS) -lm

pwcw2mat:  $(PWCW2MAT_OBJECTS)
	    $(CC)  $(CFLAGS) -o pwcw2mat $(PWCW2MAT_OBJECTS) -lm

aln2mat:  $(ALN2MAT_OBJECTS)
	    $(CC)  $(CFLAGS) -o aln2mat $(ALN2MAT_OBJECTS) -lm
all: $(EXEC)

clean:  $(OBJECTS)
	    touch ../*/*.c 


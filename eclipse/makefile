#
# Basic T-Coffee makefile
#
# This makefile depends on the following environment variables 
# - LIB_DIR: the base path containing  T-Coffee library folders 
# - CFLAGS: flags used to control the compilation process (i.e. optimization, debug symbols, etc) 
# - DATE: the date to the be stamped in T-Coffee binary  
# - VERSION: the version to be stamped in T-Coffee binary 
# 
 
CC=gcc

ifndef LIB_DIR 
LIB_DIR = ../lib
endif 

ifndef VERSION 
VERSION = "development@`date +%y%m%d_%H:%M`"
endif 

# include the main T-Coffee build file 
include $(LIB_DIR)/c_make_lib/lib_makefile

HEADERS = $(UT_LIB_HEADER) $(IO_LIB_HEADER)   $(DP_LIB_HEADER) $(DEV1_LIB_HEADER) $(DEV2_LIB_HEADER) $(DEV3_LIB_HEADER) $(DEV4_LIB_HEADER) $(FASTAL_LIB_HEADER) $(T_COFFEE_LIB_HEADER)
OBJECTS = $(UT_LIB_OBJECTS) $(IO_LIB_OBJECTS) $(DP_LIB_OBJECTS) $(DEV1_LIB_OBJECTS) $(DEV2_LIB_OBJECTS) $(DEV3_LIB_OBJECTS) $(DEV4_LIB_OBJECTS) $(FASTAL_LIB_OBJECTS) $(T_COFFEE_LIB_OBJECTS)


all: t_coffee
	
	
clean:
	rm `ls * | grep -v makefile`


t_coffee: $(OBJECTS) $(USE_MEM_LIB) $(FORTRAN_LIB_OBJECTS)
	$(CC) $(CFLAGS) -o t_coffee $(OBJECTS) $(USE_MEM_LIB) -lm

.SUFFIXES: .o_$(OSNAME) .pl

.c.o_$(OSNAME):
	$(CC) $(CFLAGS) -c $*.c -o $*.o_$(OSNAME)


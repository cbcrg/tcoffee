#ifndef KM_UTIL_H_
#define KM_UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>



FILE *
my_fopen(char *name_f, char *mode);

void *
my_malloc(size_t size);

void *
my_realloc(void *p, size_t size);

void *
my_calloc ( size_t num, size_t size );


char *
my_make_temp_dir(char *template, char *function, char *file);





#endif


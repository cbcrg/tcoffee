// #include "km_util.h"
#include "km_coffee_header.h"
FILE *
my_fopen(char *name_f, char *mode)
{
	FILE *name_F = fopen(name_f, mode);
	if (name_F == NULL)
	{
		fprintf(stderr, "ERROR: Cannot open file %s!\n", name_f);
		exit(1);
	}
	else
		return name_F;
}


void *
my_malloc(size_t size)
{
	void *p = malloc(size);
	if (p == NULL)
	{
		fprintf(stderr, "ERROR: Could not allocate space of size %li\n", size);
		exit(1);
	}
	return p;
}

void *
my_calloc ( size_t num, size_t size )
{
	void *p = calloc(num, size);
	if (p == NULL)
	{
		fprintf(stderr, "ERROR: Could not allocate space of size %li\n", size*num);
		exit(1);
	}
	return p;
}


void *
my_realloc(void *p, size_t size)
{
	p = realloc(p, size);
	if (p == NULL)
	{
		fprintf(stderr, "ERROR: Could not allocate space of size %li\n", size);
		exit(1);
	}
	return p;
}




char *
my_make_temp_dir(char *template, char *function)
{

	char *temp_dir_name = malloc(20 * sizeof(char*));
	sprintf(temp_dir_name, "%s", template);
// 	char hostname[50];
// 	gethostname(hostname, 50);
// 	printf("%s %s\n",hostname, temp_dir_name);
	if ((temp_dir_name = mkdtemp(temp_dir_name))==NULL)
	{
		int errsv = errno;
		fprintf(stderr, "ERROR! A temporary directory could not be created: %s\n",strerror(errsv));
		fprintf(stderr, "This error was caused by '%s' in '%s'\n", function);

		exit(-1);
	}
	return temp_dir_name;
}










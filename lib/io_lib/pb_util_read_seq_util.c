#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <ctype.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"


/*
*	Prototypes
*/


void		fatal(char *,...);
void		error(char *,...);
void		warning(char *,...);
char 		*rtrim(char *);
char 		*blank_to_(char *);	/* DES change blanks to _ */
char 		*upstr(char *);
char 		*lowstr(char *);
void 		getstr(char *,char *);
double		getreal(char *,double,double,double);
int		getint(char *,int,int,int);
void		do_system(void);
Boolean		linetype(char *,char *);
Boolean 	blankline(char *);
void		get_path(char *,char *);

/*
*	ckalloc()
*
*	Tries to allocate "bytes" bytes of memory. Exits program if failed.
*	Return value:
*		Generic pointer to the newly allocated memory.
*/

void *ckalloc(size_t bytes)
{
	register void *ret;
        extern void *vcalloc (size_t nelem, size_t elsize);
	
	if( (ret = vcalloc(bytes, sizeof(char))) == NULL)
/*
	if( (ret = vmalloc(bytes)) == NULL)
*/
		fatal("Out of memory\n");
	else
		return ret;	
	return ret;
}

/*
*	ckvrealloc()
*
*	Tries to vreallocate "bytes" bytes of memory. Exits program if failed.
*	Return value:
*		Generic pointer to the re-allocated memory.
*/

void *ckvrealloc(void *ptr, size_t bytes)
{
	register void *ret;
        extern void *vrealloc (void *ptr, size_t size);
	
	if( (ret = vrealloc(ptr, bytes)) == NULL)
		fatal("Out of memory\n");
	else
		return ret;	
	return ret;
}

/*
*	ckfree()
*
*	Tries to free memory allocated by ckalloc.
*	Return value:
*		None.
*/

void ckfree(void *ptr)
{
	vfree(ptr);
}


/*
*	fatal()
*
*	Prints error msg to stderr and exits.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void fatal( char *msg,...)
{
	va_list ap;
	
	va_start(ap,msg);
	fprintf(stderr,"\n\nFATAL ERROR: ");
	vfprintf(stderr,msg,ap);
	fprintf(stderr,"\n\n");
	va_end(ap);
	myexit(EXIT_FAILURE);
}

/*
*	error()
*
*	Prints error msg to stderr.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void error( char *msg,...)
{
	va_list ap;
	
	va_start(ap,msg);
	fprintf(stderr,"\n\nERROR: ");
	vfprintf(stderr,msg,ap);
	fprintf(stderr,"\n\n");
	va_end(ap);
}

/*
*	warning()
*
*	Prints warning msg to stderr.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void warning( char *msg,...)
{
	va_list ap;
	
	va_start(ap,msg);
	fprintf(stderr,"\n\nWARNING: ");
	vfprintf(stderr,msg,ap);
	fprintf(stderr,"\n\n");
	va_end(ap);
}


/*
*	rtrim()
*
*	Removes trailing blanks from a string
*
*	Return values:
*		Pointer to the processed string
*/

char * rtrim(char *str)
{
	register int p;

	p = strlen(str) - 1;
	
	while ( isspace(str[p]) )
		p--;
		
	str[p + 1] = EOS;
	
	return str;
}


/*
*	blank_to_()
*
*	Replace blanks in a string with underscores
*
*       Also replaces , ; : ( or ) with _
*
*	Return value:
*		Pointer to the processed string
*/

char * blank_to_(char *str)
{
	int i,p;

	
	p = strlen(str) - 1;
	for(i=0;i<=p;i++) 
	  {
	    if( strrchr(";,():",str[i]))str[i]='_';
	    else if (isspace(str[i]));
	  }
	return str;
}


/*
*	upstr()
*
*	Converts string str to uppercase.
*	Return values:
*		Pointer to the converted string.
*/

char * upstr(char *str)
{
	register char *s = str;
	
	while( (*s = toupper(*s)) )
		s++;
		
	return str;
}

/*
*	lowstr()
*
*	Converts string str to lower case.
*	Return values:
*		Pointer to the converted string.
*/

char * lowstr(char *str)
{
	register char *s = str;
	
	while( (*s = tolower(*s)) )
		s++;
		
	return str;
}

void getstr(char *instr,char *outstr)
{	
	fprintf(stdout,"%s: ",instr);
	fgets(outstr, 100, stdin);
}

double getreal(char *instr,double minx,double maxx,double def)
{
	int status;
	double ret;
	char line[MAXLINE];	
	
	while(TRUE) {
		fprintf(stdout,"%s (%.1lf-%.1lf)   [%.1lf]: ",instr,minx,maxx,def);
		fgets(line, MAXLINE, stdin);
		status=sscanf(line,"%lf",&ret);
		if(status == EOF) return def;
		if(ret>maxx) {
			fprintf(stderr,"ERROR: Max. value=%.1lf\n\n",maxx);
			continue;
		}
		if(ret<minx) {
			fprintf(stderr,"ERROR: Min. value=%.1lf\n\n",minx);
			continue;
		}
		break;
	}
	return ret;
}


int getint(char *instr,int minx,int maxx, int def)
{
	int ret,status;
	char line[MAXLINE];	

	while(TRUE) {
		fprintf(stdout,"%s (%d..%d)    [%d]: ",
		instr,(pint)minx,(pint)maxx,(pint)def);
		fgets(line, MAXLINE, stdin);
		status=sscanf(line,"%d",&ret);
		if(status == EOF) return def;
		if(ret>maxx) {
			fprintf(stderr,"ERROR: Max. value=%d\n\n",(pint)maxx);
			continue;
		}
		if(ret<minx) {
			fprintf(stderr,"ERROR: Min. value=%d\n\n",(pint)minx);
			continue;
		}
		break;
	}
	return ret;
}

void do_system()
{
	char line[MAXLINE];
	
	getstr("\n\nEnter system command",line);
	if(*line != EOS)
		system(line);
	fprintf(stdout,"\n\n");
}


Boolean linetype(char *line,char *code)
{
	return( strncmp(line,code,strlen(code)) == 0 );
}

Boolean blankline(char *line)
{
	int i;

	for(i=0;line[i]!='\n' && line[i]!=EOS;i++) {
		if( isdigit(line[i]) ||
		    isspace(line[i]) ||
		    (line[i] == '*') ||
                    (line[i] == '.')) 
			;
		else
			return FALSE;
	}
	return TRUE;
}


void get_path(char *str,char *path)
{
	register int i;
	
	strcpy(path,str);
	for(i=strlen(path)-1;i>-1;--i) {
		if(str[i]==DIRDELIM) {
			i = -1;
			break;
		}
		if(str[i]=='.') break;
	}
	if(i<0)
		strcat(path,".");
	else
		path[i+1]=EOS;
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"


void print_list(Constraint_list *CL)
     {

     fprintf ( stderr, "\nPRINT LIST");
     output_list ( CL, stderr);
     }
void save_full_list (Constraint_list *CL, char*fname)
     {
     FILE *fp;
     fp=vfopen ( fname, "w");
     fp=output_list ( CL, fp);
     fclose ( fp);
     }
FILE * output_list ( Constraint_list *CL, FILE *fp)
     {
     int a;
     
     fprintf ( fp, "\nPRINT LIST: %d Elements\n", CL->ne);
     for ( a=0; a<CL->ne; a++)fp=output_pair(CL, a, fp);
     fprintf (fp, "\n");
     return fp;
     }
FILE * output_pair (Constraint_list *CL,int p, FILE *fp)
	{	
	int a;
	fprintf (fp, "\n");
	for ( a=0; a<CL->entry_len; a++)
		{
		fprintf (fp, "%4d ", vread_clist(CL,p,a));
		}
	return fp;
	}
void print_pair (Constraint_list *CL,int p)
	{	
	int a;
	fprintf ( stderr, "\n");
	for ( a=0; a<CL->entry_len; a++)
		{
		fprintf ( stderr, "%d ", vread_clist(CL,p,a));
		}
	fprintf ( stderr, "\n");
	}

int** bin_list (Constraint_list *CL,int field, int Threshold)
        {
	int  a, c;
	int max;
	int **bin_list;
	CLIST_TYPE x;
	
	max=return_max_constraint_list (CL, CONS);
	
	bin_list=declare_int (max+1, 5);
	for (c=0,a=0; a<(CL->ne); a++)
	    if ( vread_clist(CL,a,field)!=UNDEFINED && vread_clist(CL,a,field)>Threshold )
		{
		x=vread_clist(CL,a,CONS);
		bin_list[x][0]=x;
		bin_list[x][1]++;
	        bin_list[x][2]+=vread_clist(CL, a, field);
		c++;
		}
	
	for ( a=0; a<= max; a++)
	    {
	    if (bin_list[a][0]>0)
		{
		bin_list[a][3]=bin_list[a][2]/bin_list[a][1];
		bin_list[a][4]=(a==0)?0:(bin_list[a][3]/a);
		}
	    }

	return bin_list;
	}
  

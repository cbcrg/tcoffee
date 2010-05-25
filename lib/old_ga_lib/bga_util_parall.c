#define MAIN 1
#include "ga.h"
#include <curses.h>

void save_population ( Population *POP, int g,Parameter *PARAM)
	{
	int a, b, c;
	char sname[500];
	char cname[500];
	FILE *fp;
	static Decoded_chromosome *A;
	static Chromosome *C;
	
	if ( A==NULL)
		{
		C=declare_chromosome(PARAM);
		A=declare_decoded_chromosome(PARAM);
		}
		 
	for ( a=0; a< (PARAM->BGA)->MAXPOP; a++)
		{
		extract_chrom (POP,g,a,C,PARAM);
		decode_chromosome ( C, A, PARAM);
		sprintf ( cname, "dumped_i_%s_%d.ga", (PARAM->BGA)->experience_name, a);
		output_full_decoded_chromosome ( A, cname,PARAM);
		sprintf ( sname, "i_score_%s_%d.ga", (PARAM->BGA)->experience_name, a);
		fp=vfopen ( sname, "w");
		fprintf ( fp, "%d", (int) A->score);
		fclose ( fp);
		}
 	set_flag_file ( (PARAM->BGA)->POP_DUMP_FLAG_FILE, 1);
 	}

int wait_flag (char *fname,Parameter * PARAM)
	{
	char stop;
	if ( strcmp ( fname, "/dev/null")==0)("\nWaiting for a flag with No name");
	if ( fname==NULL)crash ("\nWaiting for a flag with No name");
	while ( read_flag_file( fname)!=1)
		{
		if ( read_flag_file ((PARAM->BGA)->FINISH_FLAG_FILE)==1)
			{
			set_flag_file ( (PARAM->BGA)->NORMAL_END_FLAG_FILE, 1); 
			exit (0);
			}
		system ( "sleep 1");
		}
	}
	 	
void read_population ( Population *POP,int g, Parameter *PARAM) 	
	{
	static Decoded_chromosome *A;
	static Chromosome *C;
	int a;
	char cname[500];
	char stop;
	
	if ( A==NULL)
		{
		C= declare_chromosome ( PARAM);
		A=declare_decoded_chromosome(PARAM);
 		}
 		
 	for ( a=0; a< (PARAM->BGA)->MAXPOP; a++)
 		{
 		sprintf ( cname, "new_i_%s_%d.ga", (PARAM->BGA)->experience_name, a);
 		input_full_decoded_chromosome ( A, cname, PARAM); 
 		code_chromosome ( A, C, PARAM);
		insert_chrom ( C,POP, g, a, PARAM);
		}
	set_flag_file ((PARAM->BGA)->POP_READ_FLAG_FILE, 0);
	}
					
void set_flag_file ( char *fname, int val)
 	{
 	FILE * fp;
 	if ( fname==NULL) return;
 	if ( strcmp ( fname, "/dev/null")==0)return;
 	fp=vfopen (fname, "w");
 	fprintf ( fp, "%d", val);
 	fclose ( fp);
 	} 
int read_flag_file ( char *fname)	 
 	{
 	FILE *fp;
 	int a;
 	char stop;
 	
 	if ( strcmp ( fname, "/dev/null")==0)return;
 	if ( fname==NULL)return;
 	while ((fp= fopen ( fname, "r"))==NULL);
 	fscanf ( fp, "%d", &a);
 	fclose ( fp);
 	return a;
 	}
 	

 	
 	

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
void print_atom ( Atom*A);

float drmsd ( Alignment *A,  float max_distance, float delta);

float **** quantile_apdb_filtration ( Alignment *A, float ****residues, Constraint_list *CL,Pdb_param *PP, FILE *fp);
float **** irmsdmin_apdb_filtration ( Alignment *A, float ****residues, Constraint_list *CL,Pdb_param *PP, FILE *fp);
float get_default_pdb_probe()
{
  static float probe;
  static int set;

  if (!set)
    {
      if (getenv ("PDB_PROBE"))
	{
	  probe=atof (getenv ("PDB_PROBE"));
	}
      else
	{
	  probe=1.2;
	}
      set=1;
    }
  return probe;
}
float get_default_pdb_max()
{
  static float probe;
  static int set;

  if (!set)
    {
      if (getenv ("PDB_MAXD"))
	{
	  probe=atof (getenv ("PDB_MAXD"));
	}
      else
	{
	  probe=15;
	}
      set=1;
    }
  return probe;
}


char * pdb2contacts2lib (Sequence *S, char *mode,float max, char *name, char * iscope)
{
  //if max==0 returns contacts
  //if max >0 returns all pairs of CA within Max
  Ca_trace **T;
  FILE *fp;
  int keptD=0, totD=0;
  int intra=1;
  int inter=2;
  int all=3;
  int scope;
  int s1;
  int ***strikeS;
  int ***strikeB;
  static int **matrix;
  
  if (!iscope)scope=all;
  else if (strm (iscope, "inter"))scope=inter;
  else if (strm (iscope, "intra"))scope=intra;
  else if (strm (iscope, "all"  ))scope=all;
  else scope=all;
  
  

  
  if (!max)
    {
      if (strm (mode, "distances"))max=get_default_pdb_max();
      else max=get_default_pdb_probe();
    }
  
  
  if (!name)name=vtmpnam (NULL);
  fp=vfopen (name, "w");
  fprintf ( fp, "! TC_LIB_FORMAT_01\n");
  fprintf ( fp, "!CMT: [SOURCE] +seq2contacts %s %.2f (Angstrom)\n",mode,max);
    
  fprintf ( fp, "%d\n", S->nseq);
  for (s1=0; s1<S->nseq; s1++)fprintf ( fp, "%s %d %s\n",S->name[s1], S->len[s1], S->seq[s1]);
  

  //Collect All the Traces
  T=(Ca_trace**)vcalloc (S->nseq,sizeof (Ca_trace*));
  for (s1=0; s1<S->nseq; s1++)
    {
      char *p=seq2T_value (S, s1, "template_file", "_P_");
      if ( !p || !is_pdb_file (p));
      else
	{
	  
	   T[s1]=read_ca_trace    (p, "ATOM");
	   T[s1]=trim_ca_trace    (T[s1], S->seq[s1]);
	}
    }
  
  
  strikeS=(int***)declare_arrayN(3,sizeof (int), S->nseq, S->nseq, 2);
  strikeB=(int***)declare_arrayN(3,sizeof (int), S->nseq, S->nseq, 2);
  S=fast_get_sequence_type(S);
  
  if (!matrix)
    {
      if (strm (S->type, "RNA"))matrix=read_matrice ("strikeR");
      else matrix=read_matrice ("strikeP");;
    }
  
  //Measure the Strike Background
  if (!strm (mode, "distances")) 
    {
      int s1, s2;
      for (s1=0; s1<S->nseq; s1++)
	{
	  for (s2=0; s2<S->nseq; s2++)
	    {
	      int i, j;
	      int l1=S->len[s1];
	      int l2=S->len[s1];
	      
	      for (i=0; i<l1; i++)
		{
		  for (j=0; j<l2; j++)
		    {
		      strikeB[s1][s2][0]+=matrix[tolower(S->seq[s1][i])][tolower(S->seq [s2][j])];
		      strikeB[s1][s2][1]++;
		    }
		}
	    }
	}
    }
  
  //Run Self and Intra depending on Scope
  for (s1=0; s1<S->nseq; s1++)
    {
      if ( T[s1])
	{
	  
	  int l1=S->len[s1];
	  if (scope == intra || scope == all)
	    {
	      int i, j;
	      float **d;
	      
	      if (strm (mode, "distances"))     d=trace2ca_distances      (T[s1], max);
	      else if ( strm (mode, "closest")) d=trace2closest_contacts  (T[s1], max);
	      else if ( strm (mode, "best"))    d=trace2best_contacts     (T[s1], max);
	      else if ( strm (mode, "count"))   d=trace2count_contacts    (T[s1], max);
	      else if ( strm (mode, "contacts"))d=trace2contacts          (T[s1], max);
	      else 
		{
		  printf_exit (EXIT_FAILURE,stderr, "pdb2contacts2lib: %s is an unknown mode", mode);
		}
	      fprintf ( fp, "#%d %d\n", s1+1, s1+1);
	      for (i=0; i<l1-1; i++)
		for (j=i+1; j<l1; j++)
		  {
		    totD++;
		    if ((d[i][j])>=UNDEFINED && d[i][j]>0)
		      {
			keptD++;
			fprintf (fp, "%d %d %d\n", i+1, j+1,(int)(d[i][j]*100)); 
			strikeS[s1][s1][0]+=matrix[tolower(S->seq[s1][i])][tolower(S->seq [s1][j])];
			strikeS[s1][s1][1]++;
		      }
		  }
	      free_float (d,-1);
	    }
	  if ( scope == inter || scope == all)
	    {
	      int s2;
	      
	      for (s2=s1+1;s2<S->nseq; s2++)
		{
		  if ( T[s2])
		    {
		      int i, j;
		      int l2=S->len[s2];
		      float **d;
		      
		      if (strm (mode, "distances"))     d=traces2ca_distances      (T[s1],T[s2], max);
		      else if ( strm (mode, "closest")) d=traces2closest_contacts  (T[s1],T[s2], max);
		      else if ( strm (mode, "best"))    d=traces2best_contacts     (T[s1],T[s2], max);
		      else if ( strm (mode, "count"))   d=traces2count_contacts    (T[s1],T[s2], max);
		      else if ( strm (mode, "contacts"))d=traces2contacts          (T[s1],T[s2], max);
		      else 
			{
			  printf_exit (EXIT_FAILURE,stderr, "pdb2contacts2lib: %s is an unknown mode", mode);
			}
		      
		      fprintf ( fp, "#%d %d\n", s1+1, s2+1);
		      for (i=0; i<l1; i++)
			for (j=0; j<l2; j++)
			  {
			    totD++;
			    
			    if ((d[i][j])>=UNDEFINED && d[i][j]>0)
			      {
				//HERE ("%d %d %.3f",i,j, d[i][j]);
				keptD++;
				fprintf (fp, "%d %d %d\n", i+1, j+1,(int)(d[i][j]*100)); 
				strikeS[s1][s2][0]+=matrix[tolower(S->seq[s1][i])][tolower(S->seq [s2][j])];
				strikeS[s1][s2][1]++;
			      }
			  }
		      free_float (d,-1);
		    }
		}
	    }
	}
    }
  if (!strm (mode, "distances"))
    {
      int s1, s2;
      if (scope== intra || scope == all)
	for (s1=0; s1<S->nseq; s1++)
	  {
	    float rs=(strikeS[s1][s1][1]>0)?strikeS[s1][s1][0]/strikeS[s1][s1][1]:0;
	    float bg=(strikeB[s1][s1][1]>0)?strikeB[s1][s1][0]/strikeB[s1][s1][1]:0;
	    fprintf ( fp, "!CMT: [INFO] %s Strike Contact Score: RS: %6.2f BG: %6.2f Ratio: %4.2f NContacts: %d\n",S->name[s1], rs, bg, (bg>0)?rs/bg:0,strikeS[s1][s1][1]);
	  }
      
      if (scope == inter || scope == all)
	for (s1=0; s1<S->nseq; s1++)
	  for (s2=s1+1; s2<S->nseq; s2++)
	    {
	      float rs=(strikeS[s1][s2][1]>0)?strikeS[s1][s2][0]/strikeS[s1][s2][1]:0;
	      float bg=(strikeB[s1][s2][1]>0)?strikeB[s1][s2][0]/strikeB[s1][s2][1]:0;
	      fprintf ( fp, "!CMT: [INFO] %s vs %s Strike Contact Score: RS: %6.2f BG: %6.2f Ratio: %4.2f NContacts: %d\n",S->name[s1], S->name[s2],rs, bg, (bg>0)?rs/bg:0, strikeS[s1][s2][1]);
	    }
    }
  
  free_arrayN((void**)strikeS, 3);
  free_arrayN((void**)strikeB, 3);
  

  fprintf ( fp, "!CMT: [INFO] Size: %d\n", keptD);
  fprintf ( fp, "!CMT: [INFO] KEPT %d pairs out of %d -- %.2f%%\n", keptD, totD, (float)(keptD*100)/(float)totD);
  fprintf ( fp, "! SEQ_1_TO_N\n");
  vfclose (fp);
  return name;
}
int apdb ( int argc, char *argv[])
    {

	Constraint_list *CL=NULL;
	Sequence  *S=NULL;
	Alignment *A=NULL;
	Alignment *EA=NULL;
	Pdb_param *pdb_param;

	Fname *F=NULL;
	char *file_name;
	int a,c;

	int n_pdb;

/*PARAMETERS VARIABLES*/
	int garbage;
	char *parameters;
	FILE *fp_parameters;

	int quiet;
	char *se_name;
	FILE *le=NULL;

	char **list_file;
	int n_list;
	char **struc_to_use;
	int n_struc_to_use;

	char *aln;
	char *repeat_seq;
	char *repeat_pdb;

	char *color_mode;
	char *comparison_io;

	int  n_excluded_nb;

	float maximum_distance;
	float maximum_ratio;
	
	float similarity_threshold;
	float   md_threshold;


	int print_rapdb;

	char  *outfile;
	char  *run_name;

	char  *apdb_outfile;
	char *cache;

	char **out_aln_format;
	int    n_out_aln_format;

	char *output_res_num;
	char *local_mode;
	float filter;
	int filter_aln;
	int irmsd_graph;
	int nirmsd_graph;
	int n_template_file;
	char **template_file_list;
	char *mode;
		int prot_min_sim;
	int prot_max_sim;
	int prot_min_cov;
	int pdb_min_sim;
	int pdb_max_sim;
	int pdb_min_cov;
	int gapped;


	char *prot_blast_server;
	char *pdb_blast_server;

	char *strike;

	char *pdb_db;
	char *pdm;
	char *prot_db;
	int min_ncol;

	argv=standard_initialisation (argv, &argc);

/*PARAMETER PROTOTYPE:    READ PARAMETER FILE     */
	       declare_name (parameters);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-parameters"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "R_F"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Read the files in the parameter file" ,\
			    /*Parameter*/ &parameters          ,\
			    /*Def 1*/     "NULL"             ,\
			    /*Def 2*/     "stdin"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       if ( parameters && parameters[0])
	          {
		  argv[argc]=(char*)vcalloc ( VERY_LONG_STRING, sizeof(char));
		  a=0;
		  fp_parameters=vfopen (parameters, "r");
		  while ((c=fgetc (fp_parameters))!=EOF)argv[1][a++]=c;
		  vfclose (fp_parameters);
		  argv[argc][a]='\0';
		  argc++;
		  argv=break_list ( argv, &argc, "=:;, \n");
		  }
/*PARAMETER PROTOTYPE*/
     declare_name (se_name);
     get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-quiet"      ,\
			    /*Flag*/      &quiet        ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &se_name      ,\
			    /*Def 1*/     "stderr"   ,\
			    /*Def 2*/     "/dev/null"   ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

     le=vfopen ( se_name, "w");
     fprintf ( le, "\nPROGRAM: %s\n",argv[0]);

/*PARAMETER PROTOTYPE:        IN */
	list_file=declare_char ( 200, STRING);
	n_list=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-in"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ list_file     ,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:        IN */
	struc_to_use=declare_char ( 200, STRING);
	n_struc_to_use=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-struc_to_use"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ struc_to_use     ,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:        COMPARISON IO */
	declare_name (comparison_io);
	get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-io_format"  ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/  &comparison_io,\
			    /*Def 1*/    "hsgd0123456",\
			    /*Def 2*/     ""       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:        ALN */
	declare_name (aln);
	get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-aln"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/  &aln,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:        ALN */

	declare_name (repeat_seq);
	get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-repeat_seq"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/  &repeat_seq,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:        ALN */
	declare_name (repeat_pdb);
	get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-repeat_pdb"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/  &repeat_pdb,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    Nb to exclude     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-n_excluded_nb"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Exclude the N Nb on each side of the central residue. -1 triggers an automatic setting equal to the window size corresponding to the sphere"          ,\
			    /*Parameter*/ &n_excluded_nb          ,\
			    /*Def 1*/    "-1"         ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:   diatances to count    */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-similarity_threshold"    ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "F"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &similarity_threshold,\
			    /*Def 1*/    "70"             ,\
			    /*Def 2*/    "70"             ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:   diatances to count    */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-filter"    ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "F"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Filter by only keeping the best quantile"           ,\
			    /*Parameter*/ &filter,\
			    /*Def 1*/    "1.00"             ,\
			    /*Def 2*/    "1.00"             ,\
			    /*Min_value*/ "-1.00"          ,\
			    /*Max Value*/ "1.00"           \
		   );
/*PARAMETER PROTOTYPE:   diatances to count    */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-filter_aln"    ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Lower Case For Residues Filtered Out"           ,\
			    /*Parameter*/ &filter_aln,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "1"           \
		   );
/*PARAMETER PROTOTYPE:   diatances to count    */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-irmsd_graph"    ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Outputs the irmsd, position/position"           ,\
			    /*Parameter*/ &irmsd_graph,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "1"           \
		   );
/*PARAMETER PROTOTYPE:   diatances to count    */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-nirmsd_graph"    ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Outputs the NIRMSD VS N Removed Residues Curve"           ,\
			    /*Parameter*/ &nirmsd_graph,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "1"           \
		   );
/*PARAMETER PROTOTYPE:    -rmsd_threshold     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-md_threshold"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "F"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &md_threshold ,\
			    /*Def 1*/    "1"            ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:    -maximum distances     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maximum_distance"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "F"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &maximum_distance          ,\
			    /*Def 1*/    "10"            ,\
			    /*Def 2*/    "10"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );


/*PARAMETER PROTOTYPE:    -maximum ratio     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maximum_ratio"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "F"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &maximum_ratio          ,\
			    /*Def 1*/    "0.1"            ,\
			    /*Def 2*/    "0.1"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );


/*PARAMETER PROTOTYPE:    -print_rapdb     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-print_rapdb"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Prints the neighborhood of each pair of aligned residues, along with the associated local score"          ,\
			    /*Parameter*/ &print_rapdb          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:    RUN_NAME     */
	       declare_name (run_name);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-run_name"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &run_name      ,\
			    /*Def 1*/    "default"      ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "default"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUTFILE     */
/*PARAMETER PROTOTYPE:    OUTFILE     */
	       declare_name ( outfile);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-outfile"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &outfile      ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/     "default"             ,\
			    /*Min_value*/ "default"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUTFILE     */
	       declare_name ( apdb_outfile);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-apdb_outfile"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &apdb_outfile      ,\
			    /*Def 1*/    "stdout"      ,\
			    /*Def 2*/    "default"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    OUTPUT_FORMAT    */
	       out_aln_format=declare_char ( 200, STRING);
	       n_out_aln_format=get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-output"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  200              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ out_aln_format,\
			    /*Def 1*/    "score_html"           ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );



/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (color_mode);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-color_mode"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &color_mode ,\
			    /*Def 1*/    "apdb"            ,\
			    /*Def 2*/    "irmsd"             ,\
			    /*Min_value*/ "any"           ,\
			    /*Max Value*/ "any"            \
		   );
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (output_res_num);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-seqnos"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &output_res_num ,\
			    /*Def 1*/    "off"            ,\
			    /*Def 2*/    "on"             ,\
			    /*Min_value*/ "any"           ,\
			    /*Max Value*/ "any"            \
		   );
	       declare_name (cache);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-cache"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "use,ignore,update,local, directory name"          ,\
			    /*Parameter*/ &cache       ,\
			    /*Def 1*/    "use"      ,\
			    /*Def 2*/    "update"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

	        declare_name (local_mode);
		get_cl_param(					\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-local_mode"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Mode for choosing the Neighborhood (bubble or window)\nWhen selecting window, maximum distance becomes the window 1/2 size, in residues\nWhen using sphere, maximum_distance is the sphere radius in Angstrom"          ,\
			    /*Parameter*/ &local_mode       ,\
			    /*Def 1*/    "sphere"      ,\
			    /*Def 2*/    "window"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );

/*PARAMETER PROTOTYPE:        IN */
		template_file_list=declare_char (100, STRING);
		n_template_file=get_cl_param(			\
			    /*argc*/      argc          , \
			    /*argv*/      argv          , \
			    /*output*/    &le           ,\
			    /*Name*/      "-template_file"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1000           ,\
			    /*DOC*/       "List of templates file for the sequences",\
			    /*Parameter*/ template_file_list     ,	\
			    /*Def 1*/    "_SELF_P_",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );
		/*PARAMETER PROTOTYPE:        MODE */
		declare_name (mode);
		get_cl_param(				  \
			    /*argc*/      argc          , \
			    /*argv*/      argv          , \
			    /*output*/    &le           ,\
			    /*Name*/      "-mode"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1           ,\
			    /*DOC*/       "Mode: irmsd, ",\
			    /*Parameter*/ &mode    ,	\
			    /*Def 1*/    "irmsd",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );



	get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-prot_min_sim"        ,\
			    /*Flag*/      &prot_min_sim        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Minimum similarity between a sequence and its PDB target" ,\
			    /*Parameter*/ &prot_min_sim          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "20"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
 set_int_variable ("prot_min_sim", prot_min_sim);

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-prot_max_sim"        ,\
			    /*Flag*/      &prot_max_sim        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Maximum similarity between a sequence and its BLAST relatives" ,\
			    /*Parameter*/ &prot_max_sim          ,\
			    /*Def 1*/     "90"             ,\
			    /*Def 2*/     "100"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
 set_int_variable ("prot_max_sim", prot_max_sim);

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-prot_min_cov"        ,\
			    /*Flag*/      &prot_min_cov        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Minimum coverage of a sequence by its BLAST relatives" ,\
			    /*Parameter*/ &prot_min_cov          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "0"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
set_int_variable ("prot_min_cov", prot_min_cov);

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-pdb_min_sim"        ,\
			    /*Flag*/      &pdb_min_sim        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Minimum similarity between a sequence and its PDB target" ,\
			    /*Parameter*/ &pdb_min_sim          ,\
			    /*Def 1*/     "35"             ,\
			    /*Def 2*/     "35"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

 set_int_variable ("pdb_min_sim", pdb_min_sim);
 get_cl_param(							\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-pdb_max_sim"        ,\
			    /*Flag*/      &pdb_max_sim        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Maximum similarity between a sequence and its PDB target" ,\
			    /*Parameter*/ &pdb_max_sim          ,\
			    /*Def 1*/     "100"             ,\
			    /*Def 2*/     "0"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
 set_int_variable ("pdb_max_sim", pdb_max_sim);
 get_cl_param(							\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-pdb_min_cov"        ,\
			    /*Flag*/      &pdb_min_cov        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Minimum coverage of a sequence by its PDB target" ,\
			    /*Parameter*/ &pdb_min_cov          ,\
			    /*Def 1*/     "50"             ,\
			    /*Def 2*/     "25"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
set_int_variable ("pdb_min_cov", pdb_min_cov);



declare_name (pdb_blast_server);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-pdb_blast_server"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&pdb_blast_server       ,\
			    /*Def 1*/    "EBI"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

declare_name (prot_blast_server);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-blast"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&prot_blast_server       ,\
			    /*Def 1*/    ""      ,\
			    /*Def 2*/    ""      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 //make sure that -blast and -blast_server are both supported blast>blast_server
 if ( !prot_blast_server[0])
   {
     get_cl_param(						\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-blast_server"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&prot_blast_server       ,\
			    /*Def 1*/    "EBI"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
   }
 set_string_variable ("blast_server", prot_blast_server);



 declare_name (pdm);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-dm"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "print distance matrix to file Column_<N>.dm"          ,\
			    /*Parameter*/&pdm       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "yes"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 if (strm (pdm, "yes"))set_string_variable ("dm", pdm);

 declare_name (pdb_db);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-pdb_db"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Non Redundant PDB database"          ,\
			    /*Parameter*/&pdb_db       ,\
			    /*Def 1*/    "pdb"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 set_string_variable ("pdb_db", pdb_db);


declare_name (prot_db);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-protein_db"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
 			    /*Parameter*/&prot_db       ,\
	 		    /*Def 1*/    "uniprot"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

 get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-gapped"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
 			    /*Parameter*/&gapped       ,\
	 		    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-min_ncol"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "minimum number of columns (negative: fraction)"          ,\
 			    /*Parameter*/&min_ncol       ,\
	 		    /*Def 1*/    "4"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 // set the correct mode:
 if ( strm (argv[0], "trmsd"))sprintf (mode, "trmsd");
 if ( strm (argv[0], "strike"))sprintf (mode, "strike");


 set_string_variable ("prot_db", prot_db);


		if (argc==1){myexit (EXIT_SUCCESS);}

		if ( strm (outfile,"no"))n_out_aln_format=0;

		get_cl_param( argc, argv,&le, NULL,NULL,NULL,0,0,NULL);
		prepare_cache (cache);


		if (strm ( aln, ""))
		  sprintf ( aln, "%s", argv[1]);

		//if (!is_aln (aln))
		// {
		//  printf_exit (EXIT_FAILURE, stderr, "\n\n---- ERROR: File %s must be a valid alignment [FATAL:%s-%s]\n\n",aln,argv[0], PROGRAM);
		//}

		pdb_param=(Pdb_param*)vcalloc ( 1, sizeof(Pdb_param));

		pdb_param->similarity_threshold=similarity_threshold;

		pdb_param->md_threshold=md_threshold;
		pdb_param->maximum_distance=maximum_distance;

		if ( n_excluded_nb>0)
		  pdb_param->n_excluded_nb=n_excluded_nb;
		else if ( n_excluded_nb==-1)
		  pdb_param->n_excluded_nb=(int)((float)maximum_distance/(float)1.57);
		/* Exclude all the nb within the bubble at +1, +2, +n*/
		pdb_param->print_rapdb=print_rapdb;
		pdb_param->comparison_io=comparison_io;

		pdb_param->local_mode=local_mode;
		pdb_param->color_mode=lower_string (color_mode);
		pdb_param->filter=filter;
		pdb_param->filter_aln=filter_aln;
		pdb_param->irmsd_graph=irmsd_graph;
		pdb_param->nirmsd_graph=nirmsd_graph;

		sprintf ( list_file[n_list++], "S%s", aln);

		if (!strm (repeat_seq, ""))
		  {

		    sprintf ( template_file_list[0], "%s", process_repeat (list_file[0], repeat_seq, repeat_pdb));
		    fprintf ( le, "\n##Turn a repeat List into a Template File\n");
		    le=display_file_content (le,template_file_list[0]);
		    fprintf ( le, "\n\n");
		  }
		S=read_seq_in_n_list (list_file, n_list, NULL, NULL);

		le=display_sequences_names ( S,le,0, 0);

		if ( n_template_file)
		  {
		    fprintf ( le, "\nLooking For Sequence Templates:\n");
		    for ( a=0; a< n_template_file; a++)
		      {
			fprintf ( le, "\n\tTemplate Type: [%s] Mode Or File: [%s] [Start", template_type2type_name(template_file_list[a]), template_file_list[a]);
			S=seq2template_seq(S, template_file_list[a], F);
			fprintf ( le, "]");
		      }
		  }

		if ( !strcmp(mode, "strike"))
		{
			char *cache=get_cache_4_tcoffee();
			char * x=vtmpnam(NULL);
			FILE *strike_tmp = fopen(x, "w");

			unsigned int n_seq = S->nseq;
			unsigned int i = 0;
			for (; i < n_seq; ++i)
			{
				if (S->T[i]->P->template_file != NULL)
					fprintf(strike_tmp, "%s _P_ %s%s\n", S->name[i], cache,S->T[i]->P->template_file);
			}
			fclose(strike_tmp);
			printf("\n\nSTRIKE out:\n");
			printf_system("strike -a %s -c %s",aln, x);
			exit(EXIT_SUCCESS);
		}

		if ( !strm (run_name, "default"))
		  {
		    F=parse_fname(run_name);
		    sprintf (F->name, "%s", F->full);
		  }
		else
		  {
		    F=parse_fname (aln);
		  }
		
		if ( get_string_variable ("dm"))set_string_variable ("dm", F->name);
		

		for ( a=0; a< S->nseq; a++)
		  {
		    char *p;

		    p=seq2T_value (S, a, "template_file", "_P_");

		    if (p)sprintf (S->file[a], "%s",p);
		  }

		CL=declare_constraint_list ( S,NULL, NULL, 0,NULL, NULL);
		CL->T=(Ca_trace**)vcalloc (S->nseq,sizeof (Ca_trace*));


		for ( n_pdb=0,a=0; a<S->nseq; a++)
		  {
		    if ( !is_pdb_file ( S->file[a])){CL->T[a]=NULL;continue;}
		    CL->T[a]=read_ca_trace    (S->file[a], "ATOM");
		    CL->T[a]=trim_ca_trace    (CL->T[a], S->seq[a]);
		    (CL->T[a])->pdb_param=pdb_param;
		    n_pdb++;
		  }

		A=declare_aln (S);
		A->residue_case=KEEP_CASE;
		A=main_read_aln(aln, A);
		A->CL=CL;
				
			
		EA=copy_aln (A, EA);
		
		

		if ( strm (apdb_outfile, "default"))
		  sprintf ( apdb_outfile, "%s.apdb_result", F->name);




		if ( n_pdb<2)
		  {
		    FILE *fp;
		    fp=vfopen (apdb_outfile, "w");
		    fprintf (fp, "\nYour Alignment Does Not Contain Enough Sequences With a known Structure\n");
		    fprintf (fp, "To Use APDB, your alignment must include at least TWO sequences with a known structure.\n");
		    fprintf (fp, "These sequences must be named according to their PDB identifier, followed by the chain index (if any) ex: 1fnkA\n");
		    fprintf (fp, "[FATAL:%s]\n", PROGRAM);
		    vfclose (fp);
		  }
		else if ( strm (mode, "drmsd"))
		  {
		    float score;

		    //t_coffee -other_pg irmsd  -template_file ../TEMPLATEFILE -mode drmsd -maximum_ratio 0.3 -maximum_distance 10 -aln $.tc' list
		    score=drmsd ( A, maximum_distance, maximum_ratio);
		    
		    fprintf (stdout, "DRMSD: %.2f%% MaxDistance: %.2f +/- %.2f\n", (float)(score*(float)100), maximum_distance, maximum_ratio);
		    exit (0);
		  }
		
		else if ( strm (mode, "irmsd"))
		  {
		      EA=analyse_pdb ( A, EA, apdb_outfile);
		  }
		else if ( strm (mode, "msa2tree") || strm (mode, "trmsd"))
		  {
		    EA=msa2struc_dist ( A, EA,F->name,outfile, gapped, min_ncol);
		  }
		le=display_output_filename ( le, "APDB_RESULT", "APDB_RESULT_FORMAT_01", apdb_outfile, CHECK);

		if ( n_pdb>=2)
		  {
		    declare_name (file_name);
		    for ( a=0; a< n_out_aln_format; a++)
		      {
			if ( strm2( outfile, "stdout", "stderr"))sprintf (file_name, "%s", outfile);
			else if ( strm (outfile, "default"))
			  sprintf (file_name, "%s.%s",F->name, out_aln_format[a]);
			else
			  sprintf (file_name, "%s.%s",outfile,out_aln_format[a]);

			output_format_aln (out_aln_format[a],A,EA,file_name);
			le=display_output_filename ( le, "MSA", out_aln_format[a], file_name, CHECK);
		      }
		  }
		return EXIT_SUCCESS;
    }



Constraint_list * set_constraint_list4align_pdb (Constraint_list *CL,int seq, char *dp_mode, char *local_mode, char *param_file)
{
  static Constraint_list *PWCL;
  static Pdb_param *pdb_param;
  char **x;
  int n;

  if ( !CL)
    {
      free_constraint_list (PWCL);
      return NULL;
    }
  else if ( !PWCL)
    {
      PWCL=declare_constraint_list ( CL->S,NULL, NULL, 0,NULL, NULL);

      pdb_param=(Pdb_param*)vcalloc ( 1, sizeof(Pdb_param));
      pdb_param->N_ca=0;
      pdb_param->max_delta=2.0;
      pdb_param->maximum_distance=14;
      declare_name (pdb_param->local_mode);
      sprintf (pdb_param->local_mode, "%s", local_mode);
      pdb_param->scale=50;

      PWCL->pw_parameters_set=1;
      PWCL->S=CL->S;
      PWCL->lalign_n_top=10;
      PWCL->sw_min_dist=10;

      PWCL->T=(Ca_trace**)vcalloc ( (PWCL->S)->nseq, sizeof (Ca_trace*));

      PWCL->extend_jit=0;
      PWCL->maximise=1;
      /*PWCL->gop=-40;*/
      PWCL->gop=-50;
      PWCL->gep=-20;
      sprintf (CL->matrix_for_aa_group, "vasiliky");
      PWCL->use_fragments=0;
      PWCL->ktup=0;
      PWCL->TG_MODE=1;
    }


   if ( param_file && check_file_exists ( param_file) )
	{
	  if ( (x=get_parameter ( "-nca",              &n, param_file))!=NULL){pdb_param->N_ca=atoi(x[0]);free_char (x, -1);}
	  if ( (x=get_parameter ( "-max_delta",        &n, param_file))!=NULL){pdb_param->max_delta=atof(x[0]);free_char (x, -1);}
	  if ( (x=get_parameter ( "-maximum_distance", &n, param_file))!=NULL){pdb_param->maximum_distance=atoi(x[0]);free_char (x, -1);}
	  if ( (x=get_parameter ( "-local_mode",       &n, param_file))!=NULL){sprintf (pdb_param->local_mode, "%s",x[0]);free_char (x, -1);}
	  if ( (x=get_parameter ( "-scale",            &n, param_file))!=NULL){pdb_param->scale=atoi(x[0]);free_char (x, -1);}
	  if ( (x=get_parameter ( "-gapopen", &n, param_file))!=NULL){PWCL->gop=atoi(x[0]);free_char (x, -1);}
	  if ( (x=get_parameter ( "-gapext" , &n, param_file))!=NULL){PWCL->gep=atof(x[0]);free_char (x, -1);}

	}




   sprintf ( PWCL->dp_mode, "%s", dp_mode);

   if (strm (PWCL->dp_mode, "lalign"))sprintf (PWCL->dp_mode,"sim_pair_wise_lalign");
   else if (strm (PWCL->dp_mode, "sw"))sprintf (PWCL->dp_mode,"gotoh_pair_wise_sw");

   local_mode=pdb_param->local_mode;
   if ( strm ( local_mode, "hasch_ca_trace_nb"))      PWCL->evaluate_residue_pair=evaluate_ca_trace_nb;
   else if ( strm ( local_mode, "hasch_ca_trace_bubble")) PWCL->evaluate_residue_pair=evaluate_ca_trace_bubble;
   else if ( strm ( local_mode, "hasch_ca_trace_sap1_bubble")) PWCL->evaluate_residue_pair=evaluate_ca_trace_sap1_bubble;
   else if ( strm ( local_mode, "hasch_ca_trace_sap2_bubble")) PWCL->evaluate_residue_pair=evaluate_ca_trace_sap2_bubble;

   else if ( strm ( local_mode, "hasch_ca_trace_transversal")) PWCL->evaluate_residue_pair=evaluate_ca_trace_transversal;
   else if ( strm ( local_mode, "hasch_ca_trace_bubble_2")) PWCL->evaluate_residue_pair=evaluate_ca_trace_bubble_2;
   else if ( strm ( local_mode, "hasch_ca_trace_bubble_3")) PWCL->evaluate_residue_pair=evaluate_ca_trace_bubble_3;
   else if ( strm ( local_mode, "custom_pair_score_function1"))  PWCL->evaluate_residue_pair=custom_pair_score_function1;
   else if ( strm ( local_mode, "custom_pair_score_function2"))  PWCL->evaluate_residue_pair=custom_pair_score_function2;
   else if ( strm ( local_mode, "custom_pair_score_function3"))  PWCL->evaluate_residue_pair=custom_pair_score_function3;
   else if ( strm ( local_mode, "custom_pair_score_function4"))  PWCL->evaluate_residue_pair=custom_pair_score_function4;
   else if ( strm ( local_mode, "custom_pair_score_function5"))  PWCL->evaluate_residue_pair=custom_pair_score_function5;
   else if ( strm ( local_mode, "custom_pair_score_function6"))  PWCL->evaluate_residue_pair=custom_pair_score_function6;
   else if ( strm ( local_mode, "custom_pair_score_function7"))  PWCL->evaluate_residue_pair=custom_pair_score_function7;
   else if ( strm ( local_mode, "custom_pair_score_function8"))  PWCL->evaluate_residue_pair=custom_pair_score_function8;
   else if ( strm ( local_mode, "custom_pair_score_function9"))  PWCL->evaluate_residue_pair=custom_pair_score_function9;
   else if ( strm ( local_mode, "custom_pair_score_function10")) PWCL->evaluate_residue_pair=custom_pair_score_function10;


   else
     {
       fprintf ( stderr, "\n%s is an unknown hasch mode, [FATAL]\n", local_mode);
       myexit (EXIT_FAILURE);
     }

   if ( PWCL->T[seq]);
   else
     {
       PWCL->T[seq]=read_ca_trace (is_pdb_struc((CL->S)->name[seq]), "ATOM");
       (PWCL->T[seq])->pdb_param=pdb_param;
       PWCL->T[seq]=trim_ca_trace (PWCL->T[seq], (CL->S)->seq[seq]);
       PWCL->T[seq]=hasch_ca_trace(PWCL->T[seq]);

     }


  return PWCL;
}



int evaluate_ca_trace_nb (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {

     return (int)(neighborhood_match(CL, s1,r1, s2, r2, (CL->T[s1])->Chain,(CL->T[s2])->Chain ));
   }
int evaluate_ca_trace_sap2_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {



       return sap2_neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble );

     }
int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*
	 Function documentation: start

	 int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2)
	 This function evaluates the cost for matching two residues:

	 a1 is the cost for matching the two neighborood ( bubble type), using sap
	 a1: [0,+100], +100 is the best possible match.
	 a2 is the residue type weight:
	    min=worst substitution value
	    best=best of r1/r1, r2/r2-min

	    a2=(r1/r2 -min)/best --> a1:[0, 100]

	 score=a1*a2-->[-inf, +10000];
       */



       float a1;


       a1=(int) sap1_neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble );

       return (int)a1;


     }
int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*
	 Function documentation: start

	 int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2)
	 This function evaluates the cost for matching two residues:

	 a1 is the cost for matching the two neighborood ( bubble type)
	 a1: [-inf,+100-scale], +100-scale is the best possible match.

       */



       float a1;



       a1=(int) neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble )-((CL->T[s1])->pdb_param)->scale;

              return a1;


     }
int evaluate_ca_trace_transversal (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       return (int)(transversal_match (CL, s1, r1, s2, r2, (CL->T[s1])->Transversal,(CL->T[s2])->Transversal ));
     }

int evaluate_ca_trace_bubble_3 (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*This Mode evaluates :

	 1-The Bubble
	 2-The Match of the transversal residues
       */

       int a1, l1;
       int a2, l2;
       int a;

       l1=MAX(( (CL->T[s1])->Chain )->nb[r1][0] ,((CL->T[s2])->Chain )->nb[r2][0]);
       l2=MAX(( (CL->T[s1])->Bubble)->nb[r1][0], ((CL->T[s2])->Bubble)->nb[r2][0]);

       a1=(int)(neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble ));
       a2=(int)(transversal_match (CL, s1, r1, s2, r2, (CL->T[s1])->Transversal,(CL->T[s2])->Transversal ));

       if ( !l1 && !l2)return 0;
       a=(a1+a2)/2;
       return a;
     }
int evaluate_ca_trace_bubble_2 (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*This Mode evaluates :
	 1-The Ca neighborhood
	 2-The Bubble
       */


       return (int)((neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Chain,(CL->T[s2])->Chain )));
     }


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:START                                   */
/*                                                                                           */
/*********************************************************************************************/
float matrix_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)

     {
       /*
	 Function documentation: start

	 float matrix_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
	 This function evaluates the matrix for matching two residues:

	    min=worst substitution value
	    best=best of r1/r1, r2/r2-min

	    a2=(r1/r2 -min)/best --> a1:[0, 100]

	 score=a1*a2-->[-inf, +10000];
       */



       float a2;
       float m1, m2, m;
       static float min=0;
       int a, b;

       if ( !CL->M)
	 {
	   CL->M=read_matrice ( "pam250mt");
	   min=CL->M[0][0];
	   for ( a=0; a< 26; a++)
	     for ( b=0; b< 26; b++)min=MIN(min, CL->M[a][b]);
	 }

       if ( r1<=0 || r2<=0)return 0;
       m1=CL->M[(CL->S)->seq[s1][r1-1]-'A'][(CL->S)->seq[s1][r1-1]-'A']-min;
       m2=CL->M[(CL->S)->seq[s2][r2-1]-'A'][(CL->S)->seq[s2][r2-1]-'A']-min;
       m=MAX(m1, m2);
       a2=(CL->M[(CL->S)->seq[s1][r1-1]-'A'][(CL->S)->seq[s2][r2-1]-'A']-min)/m;

       return a2;
     }


float transversal_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	int a, l1, l2;
	float score=0;
	float delta, max_delta;
	float max;
	Pdb_param*PP;

	PP=(CL->T[s1])->pdb_param;
	max_delta=PP->max_delta;

	l1=nbs1->nb[r1][0];
	l2=nbs2->nb[r2][0];

	if ( l1!=l2 || l1<(PP->N_ca)) return 0;


	max=MAX(l1, l2)*max_delta;
	for ( delta=0,a=0; a< l2 ; a++)
	       {

		   delta+=max_delta-FABS((nbs1->d_nb[r1][a]-nbs2->d_nb[r2][a]));
	       }
	score=(delta*100)/max;



       return score;
      }

float neighborhood_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	  static float **table;
	  static int table_size;
	  int a, b, l1, l2;
	  float score=0;
	  float ins, del, sub;
	  float delta, max_delta;
	  float max;
	  Pdb_param*PP;


       PP=(CL->T[s1])->pdb_param;
       max_delta=PP->max_delta;


       if ( r1> 0 && r2 >0) {r1--; r2--;}
       else return 0;

       l1=nbs1->nb[r1][0];
       l2=nbs2->nb[r2][0];

       if (table_size< (MAX(l1, l2)+1))
	 {
	 table_size=MAX(l1, l2)+1;
	 if ( table)free_float (table, -1);
	 table=NULL;
	 }
       if ( !table) table=declare_float (table_size, table_size);


       max=MAX(l1, l2)*max_delta;
       if ( max==0)return 0;


       table[0][0]=0;
       for ( b=1; b<=l2; b++)
           {
	   table[0][b]=0;
	   }
       for ( a=1; a<=l1; a++)
           {
	   table[a][0]=0;
	   for ( b=1; b<=l2 ; b++)
	       {

		   delta=max_delta-FABS((nbs1->d_nb[r1][a]-nbs2->d_nb[r2][b]));

		   del=table[a-1][b];
		   ins=table[a][b-1];
		   sub= table[a-1][b-1]+delta;

		   if ( del >= ins && del >= sub)score=del;
		   else if ( ins >= del && ins >= sub) score=ins;
		   else score=sub;
		   table[a][b]=score;
	       }
	   }


       score=((((score)*100)/max));


       return score;
      }

float sap1_neighborhood_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	/*
	  Function documentation: start

	  float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
	  This function is adapted from Taylor, Orengo, Protein Structure Alignment JMB 1989, (208)1-22
	  It is the first function where
	  score= A/(|dra-drb|+b)

	  Function documentation: end
	*/

	  static float **table;
	  static int table_size;
	  int a, b, l1, l2;
	  float score=0;
	  float ins, del, sub;
	  float delta;
	  float max;

	  int A=50;
	  int B=5;






       if ( r1> 0 && r2 >0) {r1--; r2--;}
       else return 0;

       l1=nbs1->nb[r1][0];
       l2=nbs2->nb[r2][0];

       if (table_size< (MAX(l1, l2)+1))
	 {
	 table_size=MAX(l1, l2)+1;
	 if ( table)free_float (table, -1);
	 table=NULL;
	 }
       if ( !table) table=declare_float (table_size, table_size);


       max=MAX(l1, l2)*(A/B);
       if ( max==0)return 0;


       table[0][0]=0;
       for ( b=1; b<=l2; b++)
           {
	   table[0][b]=0;
	   }
       for ( a=1; a<=l1; a++)
           {
	   table[a][0]=0;
	   for ( b=1; b<=l2 ; b++)
	       {

		   delta=A/(FABS((nbs1->d_nb[r1][a]-nbs2->d_nb[r2][b]))+B);

		   del=table[a-1][b];
		   ins=table[a][b-1];
		   sub= table[a-1][b-1]+delta;
		   if ( del >= ins && del >= sub)score=del;
		   else if ( ins >= del && ins >= sub) score=ins;
		   else score=sub;
		   table[a][b]=score;
	       }
	   }


       score=((score*100))/(max);


       return score;
      }

float sap2_neighborhood_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	/*
	  Function documentation: start

	  float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
	  This function is adapted from Taylor, Orengo, Protein Structure Alignment JMB 1989, (208)1-22
	  It is the first function where
	  score= A/(|dra-drb|+b)

	  Function documentation: end
	*/

	  static float **table;
	  static int table_size;
	  int a, b, l1, l2;
	  float score=0;
	  float ins, del, sub;
	  float delta;
	  float max;

	  Amino_acid **pep1;
	  Amino_acid **pep2;
	  static Atom *vX_1, *vY_1, *vZ_1;
	  static Atom *vX_2, *vY_2, *vZ_2;
	  static Atom *ca1, *ca2;
	  float val;

	  int A=50;
	  int B=2;




       if ( r1> 0 && r2 >0) {r1--; r2--;}
       else return 0;

       /*Make up the referencial*/
       pep1=(CL->T[s1])->peptide_chain;
       pep2=(CL->T[s2])->peptide_chain;

       /*Get Referencial for CA1*/
       if ( (pep1[r1])->C)vX_1 =diff_atom(pep1[r1]->C,pep1[r1]->CA, vX_1);
       if ( (pep1[r1])->N)vY_1 =diff_atom(pep1[r1]->N,pep1[r1]->CA, vY_1);
       if ( (pep1[r1])->CB)vZ_1=diff_atom(pep1[r1]->CB,(pep1[r1])->CA,vZ_1);
       else vZ_1=add_atom (vX_1, vY_1, vZ_1);





       /*Get Referencial for CA2*/
       if ( (pep2[r2])->C)vX_2 =diff_atom((pep2[r2])->C,(pep2[r2])->CA, vX_2);
       if ( (pep2[r2])->N)vY_2 =diff_atom((pep2[r2])->N,(pep2[r2])->CA, vY_2);
       if ( (pep2[r2])->CB)vZ_2=diff_atom((pep2[r2])->CB,(pep2[r2])->CA, vZ_2);
       else vZ_2=add_atom (vX_2, vY_2, vZ_2);




       /*END OF GETTING REFERENCIAL*/

       /*Test
       if ( r1>1 && r2>1)
	 {
	 fprintf (stdout,"\n\t*******");

	 fprintf (stdout, "RESIDUE %d %c", r1, (CL->S)->seq[s1][r1]);
	 if ( (pep1[r1])->CA)fprintf (stdout,"\n\tCA ");print_atom (pep1[r1]->CA );
	 if ( (pep1[r1])->C)fprintf (stdout,"\n\tC  ");print_atom (pep1[r1]->C );
	 if ( (pep1[r1])->N)fprintf (stdout,"\n\tN  ");print_atom (pep1[r1]->N );
	 if ( (pep1[r1])->CB)fprintf (stdout,"\n\tCB ");print_atom (pep1[r1]->CB );
	 fprintf (stdout,"\n\t*******");
	 fprintf (stdout,"\n\tvX ");print_atom ( vX_1);
	 fprintf (stdout,"\n\tvY ");print_atom ( vY_1);
	 fprintf (stdout,"\n\tvZ ");print_atom ( vZ_1);

	 ca1= copy_atom ((pep1[r1-1])->CA, ca1);
	 ca1 =diff_atom(ca1,(pep1[r1])->CA, ca1);
	 fprintf (stdout,"\n\tca ");print_atom ( ca1);
	 fprintf ( stdout, "\n\tSQ1=%d ", (int)square_atom(ca1));
	 ca1=reframe_atom(vX_1, vY_1, vZ_1, ca1, ca1);
	 fprintf ( stdout, "\n\tSQ2=%d ", (int)square_atom(ca1));
	 fprintf (stdout,"\n\tca ");print_atom ( ca1);
	 fprintf (stdout,"\n\n");
	 }
       */

       l1=nbs1->nb[r1][0];
       l2=nbs2->nb[r2][0];

       if (table_size< (MAX(l1, l2)+1))
	 {
	 table_size=MAX(l1, l2)+1;
	 if ( table)free_float (table, -1);
	 table=NULL;
	 }
       if ( !table) table=declare_float (table_size, table_size);


       max=MAX(l1, l2)*(A/B);

       if ( max==0)return 0;


       table[0][0]=0;
       for ( b=1; b<=l2; b++)
           {
	   table[0][b]=0;
	   }

       for ( a=1; a<=l1; a++)
           {
	   ca1=copy_atom ((CL->T[s1])->structure[nbs1->nb[r1][a]], ca1);
	   ca1=diff_atom(ca1,(pep1[r1])->CA, ca1);
	   ca1=reframe_atom(vX_1, vY_1, vZ_1, ca1, ca1);

	   table[a][0]=0;
	   for ( b=1; b<=l2 ; b++)
	       {
		  ca2  =copy_atom((CL->T[s2])->structure[nbs2->nb[r2][b]], ca2);
		  ca2  =diff_atom(ca2,(pep2[r2])->CA, ca2);
		  ca2  =reframe_atom(vX_2, vY_2, vZ_2, ca2, ca2);

		  ca2=diff_atom(ca2,ca1,ca2);
		  val=square_atom (ca2);

		  val=(float)sqrt ((double)val);

		  delta=A/(val+B);


		  del=table[a-1][b];
		  ins=table[a][b-1];
		  sub= table[a-1][b-1]+delta;

		  if ( del >= ins && del >= sub)score=del;
		  else if ( ins >= del && ins >= sub) score=ins;
		  else score=sub;
		  table[a][b]=score;
	       }
	   }


       score=(((score*100))/(max)-50);


       return score;
      }

/*********************************************************************************************/
/*                                                                                           */
/*         APDB                                                                              */
/*                                                                                           */
/*********************************************************************************************/
float **** irmsdmin_apdb_filtration ( Alignment *A, float ****residues, Constraint_list *CL, Pdb_param *PP, FILE *fp)
{
  int s1, s2, a,col1, n,n2=0, t,flag;
  int **pos, **list;
  float nirmsd, min_nirmsd,max_nirmsd,ref_sum, sum, sum2;
  float **normalized_len;

  normalized_len=declare_float (A->nseq+1, A->nseq+1);
  for (s1=0; s1<A->nseq; s1++)
    {
      int l1, l2, r1, r2, p;
      for (s2=0; s2<A->nseq; s2++)
	{
	  for ( l1=l2=p=0; p< A->len_aln; p++)
	    {
	      r1=A->seq_al[s1][p];
	      r2=A->seq_al[s2][p];
	      if (!is_gap(r1) && isupper(r1))l1++;
	      if (!is_gap(r2) && isupper(r2))l2++;
	    }
	  normalized_len[s1][s2]=MIN(l1,l2);
	}
    }

  pos=aln2pos_simple (A, A->nseq);
  for ( s1=0; s1< A->nseq; s1++)
    for ( s2=0; s2<A->nseq; s2++)
      {
	if ( s1==s2) continue;
	else if (!(CL->T[A->order[s1][0]]) || !(CL->T[A->order[s2][0]]))continue;

	list=declare_int (A->len_aln, 2);

	for ( sum=0,n=0,col1=0; col1< A->len_aln; col1++)
	  {
	    if ( islower (A->seq_al[s1][col1]) || islower ( A->seq_al[s2][col1]))continue;
	    else if ( pos[s1][col1]<=0 || pos[s2][col1]<=0 ) continue;
	    else if ( residues[s1][s2][pos[s1][col1]-1][0]==0)continue;

	    list[n][0]=pos[s1][col1]-1;
	    list[n][1]=(int)100000*residues[s1][s2][pos[s1][col1]-1][4];
	    sum2+=residues[s1][s2][pos[s1][col1]-1][4];
	    n++;
	  }

	if (n==0)return residues;

	sort_int_inv (list, 2, 1,0, n-1);
	for (sum=0,a=0; a<n; a++)
	  {
	    sum+=list[a][1];
	  }
	ref_sum=sum;
	nirmsd=min_nirmsd=max_nirmsd=sum/(n*n);
	t=0;


	/*1 Find the maximum*/
	sum=ref_sum;
	for (flag=0,a=0; a< n-1; a++)
	  {
	    sum-=list[a][1];
	    nirmsd=sum/((n-(a+1))*(n-(a+1)));
	    if (nirmsd<max_nirmsd)flag=1;
	    if ((nirmsd>=max_nirmsd) && flag==1)break;
	    n2=a;
	  }

	sum=ref_sum;
	for (a=0; a<n2-1; a++)
	  {
	    sum-=list[a][1];
	    nirmsd=sum/((n-(a+1))*(n-(a+1)));


       	    if ( nirmsd<min_nirmsd)
	      {
		min_nirmsd=nirmsd;
		t=a;
		if ( PP->nirmsd_graph)
		  {
		    fprintf ( stdout, "\n_NIRMSD_GRAPH %s %s POS: %4d Removed: %4d NiRMSD: %.2f", A->name[s1], A->name[s2], list[a][0],a,(nirmsd/100000)*normalized_len[s1][s2]);
		  }
	      }
	  }

	if ( PP->print_rapdb)
	  {
	    for ( a=0; a<n; a++)
	      {
		if      ( list[a][1]>0 && a<=t)fprintf ( stdout, "\nRAPDB QUANTILE REMOVE S1: %3d S2: %3d COL: %3d SCORE*100: %d", s1, s2, list[a][0], list[a][1]);
		else if ( list[a][1]>0 && a>t)fprintf ( stdout, "\nRAPDB QUANTILE KEEP   S1: %3d S2: %3d COL: %3d SCORE*100: %d", s1, s2, list[a][0], list[a][1]);
	      }
	  }

	fprintf ( stdout, "\n# MINIMISATION FILTER ON: NiRMSD minimsation resulted in the removal of %d [out of %d] Columns On the alignment %s Vs %s\n", t, n, A->name[s1], A->name[s2]);
	for ( a=0; a<=t; a++)
	  {

	    residues[s1][s2][list[a][0]][0]=0;
	    residues[s1][s2][list[a][0]][1]=0;
	    residues[s1][s2][list[a][0]][2]=0;
	    residues[s1][s2][list[a][0]][3]=0;
	    residues[s1][s2][list[a][0]][4]=-1;

	  }

	free_int (list, -1);
      }
  free_float (normalized_len, -1);
  return residues;
}
float **** quantile_apdb_filtration ( Alignment *A, float ****residues, Constraint_list *CL, Pdb_param *PP,FILE *fp)
{
  int s1, s2, a,col1, n, t;
  int **pos, **list;

  pos=aln2pos_simple (A, A->nseq);
  for ( s1=0; s1< A->nseq; s1++)
    for ( s2=0; s2<A->nseq; s2++)
      {
	if ( s1==s2) continue;
	else if (!(CL->T[A->order[s1][0]]) || !(CL->T[A->order[s2][0]]))continue;

	list=declare_int (A->len_aln, 2);

	for ( n=0,col1=0; col1< A->len_aln; col1++)
	  {
	    if ( islower (A->seq_al[s1][col1]) || islower ( A->seq_al[s2][col1]))continue;
	    else if ( pos[s1][col1]<=0 || pos[s2][col1]<=0 ) continue;

	    list[n][0]=pos[s1][col1]-1;
	    list[n][1]=(int)100*residues[s1][s2][pos[s1][col1]-1][4];
	    n++;

	  }

	sort_int_inv (list, 2, 1,0, n-1);

	t=quantile_rank ( list,1, n,PP->filter);

	if ( PP->print_rapdb)
	  {
	    for ( a=0; a<n; a++)
	      {
		if      ( list[a][1]>0 && a<t)fprintf ( stdout, "\nRAPDB QUANTILE REMOVE S1: %3d S2: %3d COL: %3d SCORE*100: %d", s1, s2, list[a][0], list[a][1]);
		else if ( list[a][1]>0 && a>t)fprintf ( stdout, "\nRAPDB QUANTILE KEEP   S1: %3d S2: %3d COL: %3d SCORE*100: %d", s1, s2, list[a][0], list[a][1]);
	      }
	  }

	for ( a=0; a<t; a++)
	  {

	    residues[s1][s2][list[a][0]][0]=0;
	    residues[s1][s2][list[a][0]][1]=0;
	    residues[s1][s2][list[a][0]][2]=0;
	    residues[s1][s2][list[a][0]][3]=0;
	    residues[s1][s2][list[a][0]][4]=-1;

	  }

	free_int (list, -1);
      }

  return residues;
}
Alignment * analyse_pdb ( Alignment *A, Alignment *ST, char *results)
      {
	int s1,s2,r1, r2,b, p;
	int **pos;
	float **normalize_len;
	float m2, m4;
	float pair_tot=0, pair_m1, pair_m2, pair_m3, pair_m4, pair_m5, pair_len=0;
	float seq_tot, seq_m1, seq_m2, seq_m3, seq_m4, seq_m5,seq_len;
	float msa_tot, msa_m1, msa_m2, msa_m3, msa_m4, msa_m5, msa_len;
	float iRMSD_unit, iRMSD_max, iRMSD_min;
	float ****residues;
	Pdb_param *PP=NULL;
	Constraint_list *CL;
	char *average_file, *pairwise_file, *total_file, *irmsd_file=0;
	FILE *fp, *average,*pairwise, *total, *irmsd_graph=0;


	fp      =vfopen ( results, "w");
	pairwise=vfopen ((pairwise_file=vtmpnam (NULL)),"w");
	average =vfopen ((average_file =vtmpnam (NULL)),"w");
	total   =vfopen ((total_file   =vtmpnam (NULL)),"w");


	CL=A->CL;

	for ( s1=0; s1< (A->S)->nseq; s1++)
	  if ( CL->T[s1]){PP=(CL->T[s1])->pdb_param;break;}

	if (PP->irmsd_graph)irmsd_graph   =vfopen ((irmsd_file =vtmpnam (NULL)),"w");

	fprintf ( fp, "\nAPDB_RESULT_FORMAT_02\n");
	residues=analyse_pdb_residues ( A, A->CL,PP);
	if ( PP->filter>=0)residues=quantile_apdb_filtration (A, residues, A->CL,PP, fp);
	else if ( PP->filter<0)residues=irmsdmin_apdb_filtration (A, residues, A->CL,PP, fp);

	pos=aln2pos_simple (A, A->nseq);





	/*Compute the alignment length for normalization*/
	  normalize_len=declare_float (A->nseq+1, A->nseq+1);
	  for (s1=0; s1<A->nseq; s1++)
	    {
	      int l1, l2, r1, r2;
	      for (s2=0; s2<A->nseq; s2++)
		{
		  for ( l1=l2=p=0; p< A->len_aln; p++)
		    {
		      r1=A->seq_al[s1][p];
		      r2=A->seq_al[s2][p];
		      if (!is_gap(r1) && isupper(r1))l1++;
		      if (!is_gap(r2) && isupper(r2))l2++;
		    }
		  normalize_len[s1][s2]=MIN(l1,l2);
		}
	    }

	  msa_len=msa_tot=msa_m1=msa_m2=msa_m3=msa_m4=msa_m5=0;

	  for ( s1=0; s1< A->nseq; s1++)
	    {
	      if ( !(CL->T[A->order[s1][0]]))continue;
	      seq_len=seq_tot=seq_m1=seq_m2=seq_m3=seq_m4=seq_m5=0;
	      for ( s2=0; s2< A->nseq; s2++)
		{
		  if ( s1==s2)continue;
		  if ( !(CL->T[A->order[s2][0]]))continue;
		  pair_tot=pair_m1=pair_m2=pair_m3=pair_m4=pair_m5=0;
		  for ( p=0; p< A->len_aln; p++)
		    {
		      r1=A->seq_al[s1][p];
		      r2=A->seq_al[s2][p];
		      b=pos[s1][p]-1;


		      if (PP->filter_aln)
			{
			   if (is_gap(r1) || is_gap(r2) || residues[s1][s2][b][0]==0)
			     {
			       A->seq_al[s1][p]=tolower(r1);
			       A->seq_al[s2][p]=tolower(r2);
			     }
			   else
			     {
			        A->seq_al[s1][p]=toupper(r1);
				A->seq_al[s2][p]=toupper(r2);
			     }

			}

		      if ( PP->irmsd_graph && ( is_gap(r1) || is_gap(r2) || residues[s1][s2][b][0]==0))
			{

			  fprintf ( irmsd_graph, "\n_IRMSD_GRAPH %10s %10s ALN: %c%c iRMSD: -1.00", A->name[s1], A->name[s2],A->seq_al[s1][p], A->seq_al[s2][p]);
			}

		      if (is_gap(r1) || is_gap(r2) || residues[s1][s2][b][0]==0)continue;
		      pair_tot++;

		      /*APDB*/
		      m2=(residues[s1][s2][b][2]*100)/residues[s1][s2][b][0];
		      if (m2>PP->similarity_threshold){pair_m3++;}

		      /*iRMSD*/

		      m4=residues[s1][s2][b][4];

		      if ( PP->irmsd_graph )
			{
			  fprintf ( irmsd_graph, "\nIRMSD_GRAPH %10s %10s ALN: %c%c iRMSD: %.2f", A->name[s1], A->name[s2],A->seq_al[s1][p], A->seq_al[s2][p], m4);
			  }
		      pair_m4+=m4;
		    }
		  pair_len=normalize_len[s1][s2];
		  if ( s1>s2)
		    {

		      fprintf ( pairwise, "\n\n#PAIRWISE: %s Vs %s",A->name[s1], A->name[s2]);
		      fprintf ( pairwise, "\n\tPAIRWISE EVALUATED: %6.2f %%    [%s Vs %s] ",  (pair_len==0)?-1:(pair_tot*100)/pair_len,A->name[s1], A->name[s2]);
		      fprintf ( pairwise, "\n\tPAIRWISE APDB:      %6.2f %%    [%s Vs %s] ",  (pair_tot==0)?-1:(pair_m3*100)/pair_tot,A->name[s1], A->name[s2]);
		      fprintf ( pairwise, "\n\tPAIRWISE iRMSD:     %6.2f Angs [%s Vs %s]",  (pair_tot==0)?-1:pair_m4/pair_tot,A->name[s1], A->name[s2]);
		      fprintf ( pairwise, "\n\tPAIRWISE NiRMSD:    %6.2f Angs [%s Vs %s] [%d pos]", (pair_tot==0)?-1:(pair_m4*pair_len)/(pair_tot*pair_tot), A->name[s1], A->name[s2], (int)pair_tot);
		      fprintf ( pairwise, "\n\tRAPDB PAIRS PAIRWISE N_NONEMPTY_PAIRS %d N_MAXIMUM_PAIRS %d",(int) pair_tot, (int)pair_len);
		      msa_m3+=pair_m3;
		       msa_m4+=pair_m4;
		       msa_tot+=pair_tot;
		       msa_len+=pair_len;
		    }
		  seq_m3+=pair_m3;
		  seq_m4+=pair_m4;
		  seq_tot+=pair_tot;
		  seq_len+=pair_len;

		}

	      fprintf ( average, "\n\n#AVERAGE For Sequence %s", A->name[s1]);
	      fprintf ( average, "\n\tAVERAGE  EVALUATED: %6.2f %%    [%s]", (seq_len==0)?-1:(seq_tot*100)/seq_len, A->name[s1]);
	      fprintf ( average, "\n\tAVERAGE  APDB:      %6.2f %%    [%s]", (seq_tot==0)?-1:(seq_m3*100)/seq_tot, A->name[s1]);
	      fprintf ( average, "\n\tAVERAGE  iRMSD:     %6.2f Angs [%s]", (seq_tot==0)?-1:seq_m4/seq_tot, A->name[s1]);
	      fprintf ( average, "\n\tAVERAGE  NiRMS:     %6.2f Angs [%s]", (seq_tot==0)?-1:(seq_m4*seq_len)/(seq_tot*seq_tot), A->name[s1]);
	      if ( strm (PP->color_mode, "apdb"))ST->score_seq[s1]=(seq_tot==0)?-1:(seq_m3*100)/pair_tot;
	      if (PP->print_rapdb)fprintf (average, "\n\tRAPDB PAIRS AVERAGE N_NONEMPTY_PAIRS %d N_MAXIMUM_PAIRS %d", (int)pair_tot, (int)pair_len);

	      if ( strm (PP->color_mode, "irmsd"))ST->score_seq[s1]=(seq_tot==0)?-1:10*((seq_m4*pair_len)/(seq_tot*seq_tot));

	    }
	  fprintf ( total, "\n\n#TOTAL for the Full MSA");
	  fprintf ( total, "\n\tTOTAL     EVALUATED: %6.2f %%  ", (msa_len==0)?-1:(msa_tot*100)/msa_len);
	  fprintf ( total, "\n\tTOTAL     APDB:      %6.2f %%  ", (msa_tot==0)?-1:(msa_m3*100)/msa_tot);
	  fprintf ( total, "\n\tTOTAL     iRMSD:     %6.2f Angs", (msa_tot==0)?-1:msa_m4/msa_tot);
	  fprintf ( total, "\n\tTOTAL     NiRMSD:    %6.2f Angs", (msa_tot==0)?-1:(msa_m4*msa_len)/(msa_tot*msa_tot));
	  if (PP->print_rapdb)fprintf (total, "\n\tRAPDB PAIRS TOTAL N_NONEMPTY_PAIRS: %d N_MAXIMUM_PAIRS %d", (int)msa_tot, (int)msa_len);

	  if ( strm (PP->color_mode, "apdb")) ST->score_aln=ST->score=A->score_aln=A->score=(msa_tot==0)?-1:(msa_m3*100)/msa_tot;
	  if ( strm (PP->color_mode, "irmsd"))ST->score_aln=ST->score=A->score_aln=A->score=(msa_tot==0)?-1:10*((msa_m4*msa_len)/(msa_tot*msa_tot));

	  vfclose (average);vfclose (total); vfclose (pairwise);if (PP->irmsd_graph)vfclose (irmsd_graph);
	  fp=display_file_content (fp, pairwise_file);
	  fp=display_file_content (fp, average_file);
	  fp=display_file_content (fp, total_file);
	  if ( PP->irmsd_graph)fp=display_file_content (fp, irmsd_file);

	  fprintf ( fp, "\n\n# EVALUATED: Fraction of Pairwise Columns Evaluated\n");
	  fprintf ( fp, "# APDB:      Fraction of Correct Columns according to APDB\n");
	  fprintf ( fp, "# iRMDS:     Average iRMSD over all evaluated columns\n");
	  fprintf ( fp, "# NiRMDS:    iRMSD*MIN(L1,L2)/Number Evaluated Columns\n");
	  fprintf ( fp, "# Main Parameter: -maximum_distance %.2f Angstrom\n", PP->maximum_distance);

	  fprintf ( fp, "# Undefined values are set to -1 and indicate LOW Alignment Quality\n");
	  fp=print_program_information (fp, NULL);




	  /*Color Output*/
	  for (iRMSD_max=0,iRMSD_min=10000,s1=0; s1<A->nseq; s1++)
	    for ( s2=0; s2< A->nseq; s2++)
	      for (p=0; p<A->len_aln; p++)
		{
		  if ( residues[s1][s2][p][4]>0)
		    {
		    iRMSD_max=MAX(iRMSD_max, residues[s1][s2][p][4]);
		    iRMSD_min=MAX(iRMSD_min, residues[s1][s2][p][4]);
		    }

		}
	  iRMSD_unit=iRMSD_max/8;

	  for (p=0; p< A->len_aln; p++)
	    for ( s1=0; s1< A->nseq; s1++)
	      {

		for ( p=0; p< A->len_aln; p++)
		  {
		    r1=A->seq_al[s1][p];
		    b=pos[s1][p]-1;
		    if ( is_gap(r1) ||  !(CL->T[A->order[s1][0]]))
		      ST->seq_al[s1][p]=NO_COLOR_RESIDUE;
		    else
		      {
			float tot_m2=0, tot_m4=0, v=0;
			seq_m2=seq_m4=0;

			for (s2=0; s2< A->nseq; s2++)
			  {
			    r2=A->seq_al[s1][p];
			    if ( s1==s2) continue;
			    if (is_gap(r2) || !(CL->T[A->order[s1][0]]) || residues[s1][s2][b][0]==0)continue;

			    seq_m2+=m2=(residues[s1][s2][b][2]*100)/residues[s1][s2][b][0];
			    tot_m2++;

			    m4=residues[s1][s2][b][4];
			    if (m4>=0)
			      {
				seq_m4+=m4;
				tot_m4++;
			      }
			  }

			if (strm ( PP->color_mode, "apdb"))
			  {
			    if (tot_m2==0)v=NO_COLOR_RESIDUE;
			    else v=MIN((seq_m2/(10*tot_m2)),9);
			  }
			else if ( strm (PP->color_mode, "irmsd"))
			  {
			    if ( tot_m4==0)v=NO_COLOR_RESIDUE;
			    else v=(8-(int)((seq_m4/(iRMSD_unit*tot_m4))))+1;
			  }
			ST->seq_al[s1][p]=v;

		      }
		  }
	      }
	  for ( p=0; p<A->len_aln; p++) ST->seq_al[A->nseq][p]=NO_COLOR_RESIDUE;


	  ST->generic_comment=(char*)vcalloc ( 100, sizeof (int));
	  if ( strm (PP->color_mode, "apdb"))
	    {
	      sprintf ( ST->generic_comment, "# APDB Evaluation: Color Range Blue-[0 %% -- 100 %%]-Red\n# Sequence Score: APDB\n# Local Score: APDB\n\n");
	    }
	  else if ( strm (PP->color_mode, "irmsd"))
	    {
	      sprintf ( ST->generic_comment, "\n# iRMSD Evaluation:\n# Sequence score: NiRMSD (Angstrom*10)\n# Local Score: iRMSD, Blue-[%.2f Ang. -- 0.00 Ang.]-Red \n", iRMSD_max);
	    }

	  fprintf ( fp, "\n");
	  vfclose (fp);
	  free_int (pos, -1);
	  return ST;
      }
float **** analyse_pdb_residues ( Alignment *A, Constraint_list *CL, Pdb_param *pdb_param)
     {

	 int **pos;
	 int s1, s2, rs1, rs2;
	 int col1, col2;
	 float ****distances;

	      /*Distances[Nseq][len_aln][4]
                distances...[0]: Number of residues within the bubble
                distances...[1]: Absolute difference of distance of residues within Bubble
	        distances...[2]: Number of residues within the bubble with Delta dist < md_threshold
		distances ..[3]: Sum of squared difference of distances
		distances ..[4]: iRMSD
	      */
	 float d1, d2,delta;
	 int wd1, wd2;
	 int in_bubble=0;
	 int real_res1_col1=0;
	 int real_res1_col2;
	 int real_res2_col1;
	 int real_res2_col2;
	 Pdb_param *PP;
	 int print_rapdb;
	 float nrapdb, rapdb;
	 Alignment *BA=NULL;

	 PP=pdb_param;
	 print_rapdb=PP->print_rapdb;

	 distances=(float****)declare_arrayN(4, sizeof (float), A->nseq, A->nseq, 0, 0);

	 /*Pre-computation of the internal distances----> T[seq]->ca_dist[len][len]*/
	 /*Can be avoided if distance_on_request set to 1 */

	 for ( s1=0; s1< A->nseq; s1++)
	   {
	     rs1=A->order[s1][0];
	     if (CL->T[rs1] &&  !(CL->T[rs1])->ca_dist)(CL->T[rs1])->ca_dist=trace2ca_distances(CL->T[rs1],0);
	     for ( s2=0; s2< A->nseq; s2++)
	       {
		 distances[s1][s2]=declare_float ( A->len_aln, 6);
	       }
	   }
	 pos=aln2pos_simple (A, A->nseq);

	 for ( s1=0; s1< A->nseq; s1++)
	   for ( col1=0; col1< A->len_aln; col1++)
	     for ( s2=0; s2<A->nseq; s2++)
	       {
		 rs1=A->order[s1][0];
		 rs2=A->order[s2][0];
		 rapdb=0;
		 nrapdb=0;
		 if ( s1==s2) continue;
		 else if (!(CL->T[rs1]) || !(CL->T[rs2]))continue;
		 else if ( islower (A->seq_al[s1][col1]) || islower ( A->seq_al[s2][col1]))continue;
		 else if ( pos[s1][col1]<=0 || pos[s2][col1]<=0 ) continue;

		 if ( print_rapdb && s2>s1)
		   {

		     fprintf ( stdout, "RAPDB S1: %s S2: %s POS %d %d %c %d %c ", A->name[s1], A->name[s2], col1+1, pos[s1][col1],A->seq_al[s1][col1], pos[s2][col1],A->seq_al[s2][col1]);
		     BA=copy_aln (A, BA);
		     lower_string (BA->seq_al[s1]);
		     lower_string (BA->seq_al[s2]);
		     BA->seq_al[s1][col1]=toupper (BA->seq_al[s1][col1]);
		     BA->seq_al[s2][col1]=toupper (BA->seq_al[s2][col1]);
		   }

		 for ( col2=0; col2<A->len_aln; col2++)
		   {

		     if (pos[s1][col2]<=0 || pos[s2][col2]<=0 )continue;
		     else if ( FABS((pos[s1][col2]-pos[s1][col1]))<=PP->n_excluded_nb)continue;
		     else if ( FABS((pos[s2][col2]-pos[s2][col1]))<=PP->n_excluded_nb)continue;
		     else if ( islower (A->seq_al[s1][col2]) || islower ( A->seq_al[s2][col2]))continue;

		     real_res1_col1=pos[s1][col1]-1;
		     real_res1_col2=pos[s1][col2]-1;

		     real_res2_col1=pos[s2][col1]-1;
		     real_res2_col2=pos[s2][col2]-1;

		     d1=(CL->T[rs1])->ca_dist[real_res1_col1][real_res1_col2];
		     d2=(CL->T[rs2])->ca_dist[real_res2_col1][real_res2_col2];

		     if ( d1==UNDEFINED || d2 == UNDEFINED) continue;



		     if ( strm ( PP->local_mode, "sphere"))
		       {
			 in_bubble= (d1<PP->maximum_distance && d2<PP->maximum_distance)?1:0;		   ;
		       }
		     else if ( strm ( PP->local_mode, "window"))
		       {
			 wd1=FABS((pos[s1][col2]-pos[s1][col1]));
			 wd2=FABS((pos[s2][col2]-pos[s2][col1]));
			 in_bubble= (wd1<PP->maximum_distance && wd2<PP->maximum_distance)?1:0;		   ;
		       }

		     if (in_bubble)
		       {
			 if ( print_rapdb && s2 >s1)
			   {
			     fprintf ( stdout, "NB %d %d %c %d %c ", col2, pos[s1][col2], A->seq_al[s1][col2], pos[s2][col2], A->seq_al[s2][col2]);
			     BA->seq_al[s1][col2]=toupper (BA->seq_al[s1][col2]);
			     BA->seq_al[s2][col2]=toupper (BA->seq_al[s2][col2]);
			   }
			 delta=FABS((d1-d2));
			 if (delta<PP->md_threshold)
			   distances[s1][s2][real_res1_col1][2]++;
			 distances[s1][s2][real_res1_col1][1]+=delta;
			 distances[s1][s2][real_res1_col1][0]++;
			 distances[s1][s2][real_res1_col1][3]+=delta*delta;
			 nrapdb++;
			 rapdb+=delta*delta;
		       }
		   }

		 if ( nrapdb==0)distances[s1][s2][real_res1_col1][4]=-1;
		 else distances[s1][s2][real_res1_col1][4]=(float)sqrt((double)(rapdb/nrapdb));

		 if ( print_rapdb && s2>s1)
		   {
		     if (nrapdb==0)
		       {
			 fprintf ( stdout, "APDB: UNDEFINED\n");
		       }
		     else
		       {

			 fprintf ( stdout, " APDB: %.2f ",(float)sqrt((double)(rapdb/nrapdb)));
			 BA->residue_case=KEEP_CASE;unalign_residues (BA, s1, s2);
			 fprintf ( stdout,"SEQ1: %s %s SEQ2: %s %s\n", BA->name[s1], BA->seq_al[s1], BA->name[s2], BA->seq_al[s2]);
		       }
		   }

	       }

	 free_aln (BA);
	 free_int (pos, -1);
	 return distances;
     }

float drmsd ( Alignment *A, float max_distance, float delta)
     {

       int **pos;
       int s1, s2,col1,col2;
       float score=0;
       float tot=0;
       float max=0;
       Constraint_list *CL;
       if (A)CL=A->CL;
       else return -1;
       
       pos=aln2pos_simple (A, A->nseq);
       
       for (s1=0; s1<A->nseq; s1++)
	 for (col1=0; col1<A->len_aln; col1++)
	   for (col2=0; col2<A->len_aln; col2++)
	     for (s2=0; s2<A->nseq; s2++)
	       {
		 int rs1=A->order[s1][0];
		 int rs2=A->order[s2][0];
		 int rr1_s1=pos[s1][col1]-1;
		 int rr2_s1=pos[s1][col2]-1;
		 		 
		 if (s1==s2)continue;
		 if (col1==col2)continue;
		 
		 //gap in master pair
		 if (rr1_s1<0 || rr2_s1<0)continue;
		 //Missing Structure
		 else if (!(CL->T[rs1]) || !(CL->T[rs2]))continue;
		 //Master pair to be ignored
		 else if ( islower (A->seq_al[s1][col1]) || islower ( A->seq_al[s1][col2]))continue;//master pair to be ignore
		 else
		   {
		     double d1;
		     
		     if (CL->T[rs1] &&  !(CL->T[rs1])->ca_dist)(CL->T[rs1])->ca_dist=trace2ca_distances(CL->T[rs1],0);
		     if (CL->T[rs2] &&  !(CL->T[rs2])->ca_dist)(CL->T[rs2])->ca_dist=trace2ca_distances(CL->T[rs2],0);
		     d1=(CL->T[rs1])->ca_dist[rr1_s1][rr2_s1];;
		     
		     max++;
		     
		     
		     if (d1<= UNDEFINED  || d1>max_distance)continue;
		     else
			 {
			   int rr1_s2=pos[s2][col1]-1;
			   int rr2_s2=pos[s2][col2]-1;
			   
			   tot++;
			   if (rr1_s2<0 || rr2_s2<0);
			   else
			     {
			       
			       double d2=(CL->T[rs2])->ca_dist[rr1_s2][rr2_s2];
			       
			       
			       if ( d2<=UNDEFINED);
			       else if (d2>(d1-(d1*delta)) && d2< (d1+(d1*delta)))
				 {
				   score++;
				 }
			     }
			 }
		   }
	       }

       HERE ("MAX=%.0f Tot=%.0f Ratio: %.2f", max, tot, (max==0)?0:tot/max);
       return score/tot;
     }




int pair_res_suitable4trmsd    (int s1,int col1, int col2, Alignment *A, int **pos, Pdb_param *PP, Constraint_list *CL,int *s);
int aln_column_contains_gap (Alignment *A, int c);
float aln2ncol4trmsd(Alignment *A, int **pos, Constraint_list *CL, int **lc);
int pair_columns_suitable4trmsd(int col1, int col2, Alignment *A, int **pos, Pdb_param *PP, Constraint_list *CL,int *s);
int column_is_suitable4trmsd(int col1,Alignment *A, int **pos, Pdb_param *PP, Constraint_list *CL,int *s);



NT_node trmsdmat2tree (float **dm, int **count,Alignment *A, int colN);
Alignment * msa2struc_dist ( Alignment *A, Alignment *ST, char *results, char *output, int gapped, int min_ncol4trmsd)
     {

       int **pos, c;
       FILE *tl;
	 int s1, s2, rs1, rs2;
	 int col1, col2;
	 float ****distances;
	 float **dm,**tdm;
	 int **count,**tcount;
	 int print_subtrees=0;
	 float min, max;

	      /*Distances[Nseq][len_aln][4]
                distances...[0]: Number of residues within the bubble
                distances...[1]: Absolute difference of distance of residues within Bubble
	        distances...[2]: Number of residues within the bubble with Delta dist < md_threshold
		distances ..[3]: Sum of squared difference of distances
		distances ..[4]: iRMSD
	      */
	 Pdb_param *pdb_param;
	 Constraint_list *CL;
	 int a, b, ncol, npos,npos2,n;
	 float d1, d2,delta;
	 int wd1, wd2;
	 int in_bubble=0;
	 int real_res1_col1=0;
	 int real_res1_col2;
	 int real_res2_col1;
	 int real_res2_col2;
	 Pdb_param *PP;
	 int print_rapdb;
	 float nrapdb, rapdb;
	 Alignment *BA=NULL;
	 NT_node *T0,*T1,*T2,*PT, *POS;
	 NT_node BT0, BT10,BT50, BT100=NULL,RBT;
	 char **pair_pos_list;

	 int ntree=0, ntree2;

	 Alignment *B;
	 char *pos_list;
	 char *tot_pos_list;

	 char *struc_tree10;
	 char *struc_tree100;
	 char *struc_tree50;
	 char *struc_tree0;
	 char *consense_file;
	 
	 char *color_struc_tree;
	 int **score;
	 int proceed=1;
	 int **lc;
	 int used;

	 int *keep;
	 
	 if (min_ncol4trmsd<0)
	   {
	     min_ncol4trmsd*=-1;
	     min_ncol4trmsd=(min_ncol4trmsd*A->len_aln)/100;
	   }
	 else if ( min_ncol4trmsd>=A->len_aln)
	   {
	     min_ncol4trmsd=A->len_aln-1;
	   }

	 //prepare the list of positions to keep
	 
	
	
	 

	 
	 lc=declare_int (A->nseq, 2);
	 for (a=0; a<A->nseq; a++)lc[a][0]=a;

	 declare_name(tot_pos_list);
	 sprintf ( tot_pos_list, "%s.struc_tree.list", results);

	 declare_name(consense_file);
	 sprintf (consense_file, "%s.struc_tree.consense_output", results);

	 declare_name(pos_list);
	 sprintf ( pos_list, "%s.pos_list", results);

	 declare_name(struc_tree0);
	 sprintf ( struc_tree0, "%s.struc_tree.consensus",results);

	 declare_name(struc_tree10);
	 sprintf ( struc_tree10, "%s.struc_tree10",results);

	 declare_name(struc_tree100);
	 sprintf ( struc_tree100, "%s.struc_tree100",results);

	 declare_name(struc_tree50);
	 sprintf ( struc_tree50, "%s.struc_tree50",results);

	 declare_name(color_struc_tree);
	 sprintf ( color_struc_tree, "%s.struc_tree.html", results);

	 pair_pos_list=declare_char (A->len_aln*A->len_aln+1, 100);
	 T1=(tnode**)vcalloc (A->len_aln*A->len_aln+1, sizeof (NT_node));
	 T2=(tnode**)vcalloc (A->len_aln+1, sizeof (NT_node));

	 PT=(tnode**)vcalloc (A->len_aln*A->len_aln+1, sizeof (NT_node));
	 POS=(tnode**)vcalloc (A->len_aln+1, sizeof (NT_node));

	 CL=A->CL;

	 //Check all sequences have a PDB structure

	 for (used=0,a=0; a<A->nseq; a++)
	   {
	     if ( ! seq2P_template_file(A->S,a))
	       {
		 add_warning (stderr, "Sequence %s removed from the dataset [No Usable Structure]", A->name[a]);
	       }
	     else
	       {
		 if (used!=a)
		   {
		     sprintf (A->name[used], "%s", A->name[a]);
		     sprintf (A->seq_al[used], "%s", A->seq_al[a]);
		     for (b=0; b<4; b++)A->order[used][b]=A->order[a][b];
		   }
		 used++;
	       }
	   }

	 A->nseq=used;

	 if (A->nseq<2)myexit (fprintf_error(stderr, "Two sequences at least must have a known structure"));

	 for ( s1=0; s1< (A->S)->nseq; s1++)
	   if ( CL->T[s1]){PP=(CL->T[s1])->pdb_param;break;}

	 for ( s1=0; s1< A->nseq; s1++)
	   {
	     rs1=A->order[s1][0];
	     if (CL->T[rs1] &&  !(CL->T[rs1])->ca_dist)(CL->T[rs1])->ca_dist=trace2ca_distances(CL->T[rs1],0);
	   }
	 pos=aln2pos_simple (A, A->nseq);

	 dm=declare_float (A->nseq, A->nseq);
	 count=declare_int (A->nseq, A->nseq);
	 tdm=declare_float (A->nseq, A->nseq);
	 tcount=declare_int (A->nseq, A->nseq);

	 PP->maximum_distance=1000;
	 sprintf ( PP->local_mode, "sphere");

	 while ((npos=aln2ncol4trmsd(A,pos,CL,lc))<min_ncol4trmsd && A->nseq>1)
	   {

	     sort_int_inv (lc,2, 1, 0,A->nseq-1);
	     add_information (stderr, "Remove Sequence [%s] that contains %d un-suitable positions", A->name[lc[0][0]], lc[0][1]);
	     A=remove_seq_from_aln (A, A->name[lc[0][0]]);
	     ungap_aln (A);
	     pos=aln2pos_simple (A, A->nseq);
	   }
	 if (!A->nseq)
	   {
	     myexit (fprintf_error(stderr,"No suitable pair of column supporting a tree"));
	   }
	 else
	   fprintf ( stderr, "\n---- Number of usable positions: %d [%.2f %%]\n", npos, ((float)npos*100)/(float)A->len_aln);

	 tl=vfopen (tot_pos_list, "w");

	 //compile list of positions to keep
	 keep=(int*)vcalloc (A->len_aln, sizeof (int));
	 for (npos2=0,col1=0; col1<A->len_aln; col1++)
	   {
	     if (!gapped && aln_column_contains_gap (A, col1))keep[col1]=0;
	     else {keep[col1]=1;npos2++;}
	     
	   }
	 
	 if ((c=atoigetenv ("TRMSD_RANDOM_COLUMNS")))
	   {
	     int ** keep1=declare_int (A->len_aln, 2);
	     for (a=0; a<A->len_aln; a++)keep1[a][0]=a;
	     
	     for (col1=0; col1<A->len_aln; col1++)keep1[col1][1]=rand()%10000;
	     sort_int (keep1, 2, 1,0, A->len_aln-1);
	     c=(int)(((float)c/(float)100)*(float)(npos2));
	     for (col1=0; col1<A->len_aln; col1++)if (keep[col1])keep[col1]=2;
	     
	     for (npos2=0,col1=0; col1<A->len_aln && npos2<c; col1++)
	       {
		 col2=keep1[col1][0];
		 if (keep[col2]){keep[col2]=1; npos2++;}
	       }
	     for (col1=0; col1<A->len_aln; col1++)if (keep[col1]==2)keep[col1]=0;
	     
	     free_int (keep1, -1);
	   }

	 if (getenv ("TRMSD_COLUMNS2FILE"))
	   {
	     FILE *fp=vfopen (getenv ("TRMSD_COLUMNS2FILE"), "w");
	     for (a=0; a< A->nseq; a++)
	       {
		 fprintf ( fp, ">%s\n", A->name[a]);
		 for (b=0; b<A->len_aln; b++)if (keep[b])fprintf (fp, "%c", A->seq_al[a][b]);
		 fprintf (fp, "\n");
	       }
	     vfclose (fp);
	   }
	 
	 if (c) fprintf ( stderr, "\n---- Number of randomly selected position: %d [%.2f %%]\n", npos2, ((float)npos2*100)/(float)A->len_aln);
	 
	 for (ncol=0,ntree=0, col1=0; col1< A->len_aln; col1++)
	   {
	     int w,tree, cont;
	     //output_completion (stderr, col1, A->len_aln,1, "Sample Columns");
	     if (!keep[col1])continue;
	     for ( cont=1,ntree2=0,col2=0; col2<A->len_aln; col2++)
	       {
		 if (atoigetenv ("TRMSD_STRICT")==1)if (!keep[col2])continue;
		 for (s1=0; s1< A->nseq-1; s1++)
		   {
		     rs1=A->order[s1][0];
		     if (!pair_res_suitable4trmsd (s1,col1, col2, A, pos, PP, CL, &w))
		       {
			 continue;
		       }
		     for ( s2=s1+1; s2<A->nseq; s2++)
		       {
			 if (!pair_res_suitable4trmsd (s2,col1, col2, A, pos, PP, CL, &w))continue;
			 rs2=A->order[s2][0];
			 real_res1_col1=pos[s1][col1]-1;
			 real_res1_col2=pos[s1][col2]-1;
			 real_res2_col1=pos[s2][col1]-1;
			 real_res2_col2=pos[s2][col2]-1;

			 d1=(CL->T[rs1])->ca_dist[real_res1_col1][real_res1_col2];
			 d2=(CL->T[rs2])->ca_dist[real_res2_col1][real_res2_col2];

			 delta=FABS((d1-d2));
			 dm[s1][s2]=dm[s2][s1]+=delta;
			 tdm[s1][s2]=tdm[s2][s1]+=delta;
			 tcount[s1][s2]++;
			 tcount[s2][s1]++;

			 count[s1][s2]++;
			 count[s2][s1]++;
		       }
		   }
	       }



	     if ((POS[col1]=trmsdmat2tree (dm, count, A, col1+1)))
	       {
		 T1[ntree]=POS[col1];
		 fprintf (tl, "\n>Tree_%d Column\n", col1+1);
		 print_tree (T1[ntree], "newick", tl);
		 ntree++;
	       }
	   }
	 vfclose (tl);
	 vfree (keep);
	 if (!ntree){fprintf ( stderr, "\nERROR: No suitable pair of column supporting a tree. Consider removing the most distantly related sequences [FATAL]\n"); exit (EXIT_SUCCESS);}

	 //consensus tree: hijack all the output formats and print only the consensus tree
	
	 if (strm (output, "no"))
	   {
	     score=treelist2avg_treecmp (T1, NULL);
	     display_output_filename( stderr,"TreeList","newick",tot_pos_list, CHECK);
	     
	     if (treelist_file2consense (tot_pos_list, NULL, consense_file))
	       {
		 display_output_filename( stderr,"ConsenseTree","phylip",consense_file, CHECK);
	       }
	     else
	       {
		 fprintf ( stderr, "\nPhylip is not installed: the program could not produce the consense output. This is not mandatory but useful");
	       }
	   }
	 else
	   {
	     if ((BT100=treelist2filtered_bootstrap (T1, NULL,score, 1.0)))
	       {
		 vfclose (print_tree (BT100,"newick", vfopen (output, "w")));
		 display_output_filename( stderr,"Tree","newick",output, CHECK);
	       }
	     fprintf ( stderr, "\n");
	     free_int (pos, -1);
	     exit (EXIT_SUCCESS);
	   }
	 
	 //Default outpout
	 //consensus tree
	 if ((BT100=treelist2filtered_bootstrap (T1, NULL,score, 1.0)))
	       {
		 vfclose (print_tree (BT100,"newick", vfopen (struc_tree0, "w")));
		 display_output_filename( stderr,"Tree","newick",struc_tree0, CHECK);
	       }
	 if (print_subtrees)
	   {

	     if ( (BT0=trmsdmat2tree (tdm, tcount, A, 0)))
	       {
		 vfclose (print_tree (BT0,"newick", vfopen (struc_tree0, "w")));
		 display_output_filename( stderr,"Tree","newick",struc_tree0, CHECK);
	       }
	     if ((BT10=treelist2filtered_bootstrap (T1, NULL,score, 0.1)))
	       {
		 vfclose (print_tree (BT10,"newick", vfopen (struc_tree10, "w")));
		 display_output_filename( stderr,"Tree","newick",struc_tree10, CHECK);
	       }

	     if ((BT50=treelist2filtered_bootstrap (T1, NULL, score,0.5)))
	       {
		 vfclose (print_tree (BT50,"newick", vfopen (struc_tree50, "w")));
		 display_output_filename( stderr,"Tree","newick",struc_tree50, CHECK);
	       }
	   }


	 if (!BT100)BT100=treelist2filtered_bootstrap (T1, NULL,score, 1.0);

	 RBT=BT100;
	 if (RBT)
	   {
	     B=copy_aln (A, NULL);
	     for (a=0; a<A->len_aln; a++)
	       {
		 int score;
		 Tree_sim *S=NULL;

		 if (POS[a])
		   {
		     S=tree_cmp (POS[a], RBT);
		     score=S->uw/10;
		   }
		 else
		   {
		     score=NO_COLOR_RESIDUE;
		   }

		 for (b=0; b<B->nseq; b++)
		   {
		     if ( is_gap (B->seq_al[b][a]) || score == NO_COLOR_RESIDUE)
		       {
			 B->seq_al[b][a]=NO_COLOR_RESIDUE;
		       }
		     else
		       {
			 B->seq_al[b][a]=S->uw/10;
		       }
		   }
		 if (S)vfree (S);
	       }

	     output_format_aln ("score_html", A,B,color_struc_tree);
	     display_output_filename( stderr,"Colored MSA","score_html",color_struc_tree, CHECK);
	     free_aln (BA);
	     fprintf ( stderr, "\n");
	   }
	 fprintf ( stderr, "\n");
	 free_int (pos, -1);
	 exit (EXIT_SUCCESS);
	 return NULL;
     }
NT_node trmsdmat2tree (float **dm, int **count,Alignment *A, int colN)
{
  float min, max;
  int s1, s2;
  NT_node T;
  int ns;
  int **dm_int;
  char *dataset;
  ns=A->nseq;
  for (s1=0; s1<ns-1; s1++)
    for (s2=s1+1; s2<ns; s2++)
      {
	if ( count [s1][s2])dm[s1][s2]=dm[s2][s1]=dm[s1][s2]/(float)count[s1][s2];
	else
	  {
	    return NULL;
	  }
	if (s1==0 && s2==1)min=max=dm[s1][s2];
	min=MIN(dm[s1][s2], min);
	max=MAX(dm[s1][s2], max);
      }
  dm_int=declare_int (ns, ns);
  
  
  for (s1=0; s1<A->nseq-1; s1++)
    for (s2=s1+1; s2<A->nseq; s2++)
      {
	dm_int[s1][s2]=dm_int[s2][s1]=((dm[s1][s2])/(max))*100;
      }

  if ( colN && (dataset=get_string_variable ("dm"))!=NULL)
    {
      FILE *dmf;
      char dmfn[1000];
      sprintf (dmfn,"%s.column_%d.dm", dataset,colN);
      dmf=vfopen (dmfn, "w");
      fprintf (dmf, "%d\n",A->nseq);
      for (s1=0; s1<A->nseq; s1++)
	{
	  fprintf (dmf, "%-20s ", A->name[s1]);
	  for ( s2=0; s2<A->nseq; s2++)
	    fprintf (dmf, "%3d ", dm_int[s1][s2]);
	  fprintf (dmf, "\n"); 
	}
      vfclose (dmf);
    }
  
  T=compute_std_tree_2(A, dm_int, "_TMODE_upgma");
    
  free_int (dm_int, -1);
  for (s1=0; s1<ns; s1++)for ( s2=0; s2<ns; s2++){dm[s1][s2]=count[s1][s2]=0;}
  return T;
}

int pair_res_suitable4trmsd    (int s1,int col1, int col2, Alignment *A, int **pos, Pdb_param *PP, Constraint_list *CL,int *s)
{
   int rs;
   rs=A->order[s1][0];
   if ( !(CL->T[rs])){s[0]=s1; return 0;}
   else if (is_gap (A->seq_al[s1][col1])){s[0]=s1;return 0;}
   else if (is_gap (A->seq_al[s1][col2])){s[0]=s1;return 0;}

   else if (islower(A->seq_al[s1][col1])){s[0]=s1; return 0;}
   else if (islower(A->seq_al[s1][col2])){s[0]=s1; return 0;}

   else if ( FABS(((pos[s1][col2])-(pos[s1][col1])))<=PP->n_excluded_nb){s[0]=s1;return 0;}
   else if ((CL->T[rs])->ca_dist[pos[s1][col1]-1][pos[s1][col2]-1]==UNDEFINED){s[0]=s1;return 0;}
   return 1;
}
int pair_columns_suitable4trmsd(int col1, int col2, Alignment *A, int **pos, Pdb_param *PP, Constraint_list *CL,int *s)
{
  int s1;
  if (!column_is_suitable4trmsd (col1, A, pos, PP, CL,s))return 0;
  if (!column_is_suitable4trmsd (col2, A, pos, PP, CL,s))return 0;
  for (s1=0; s1<A->nseq; s1++)
    {
      int rs, rr1, rr2;

      rs=A->order[s1][0];
      if ( FABS(((pos[s1][col2])-(pos[s1][col1])))<=PP->n_excluded_nb){s[0]=s1;return 0;}
      if ((CL->T[rs])->ca_dist[pos[s1][col1]-1][pos[s1][col2]-1]==UNDEFINED){s[0]=s1;return 0;}
      rr1=pos[s1][col1]-1;
      rr2=pos[s1][col2]-1;
      if ((CL->T[rs])->ca_dist[rr1][rr2]>PP->maximum_distance){s[0]=s1;return 0;}
    }
  return 1;
}
int column_is_suitable4trmsd(int col1,Alignment *A, int **pos, Pdb_param *PP, Constraint_list *CL,int *s)
{
  int s1;
  for ( s1=0; s1<A->nseq; s1++)
    {
      int rs;
      rs=A->order[s1][0];
      if ( !(CL->T[rs])){s[0]=s1; return 0;}
      else if (is_gap (A->seq_al[s1][col1])){s[0]=s1;return 0;}
      else if (islower(A->seq_al[s1][col1])){s[0]=s1; return 0;}
    }
  return 1;
}
int aln_column_contains_gap (Alignment *A, int c)
{
  int a, b;
  if ( !A || c>=A->len_aln || c<0)
    {
      printf ( "\nERROR: values out of range in aln_column_contains_gap [FATL:%s]\n", PROGRAM);
      exit (EXIT_FAILURE);
    }
  for ( a=0; a<A->nseq; a++) if ( is_gap(A->seq_al[a][c]))return 1;
  return 0;
}


float aln2ncol4trmsd(Alignment *A, int **pos, Constraint_list *CL, int **lc)
{
  //This function estimates the number of columns suitable for constructing a trmsd
  int col1, s1, ncol, n, rs1, real_res1_col1;

  for (s1=0; s1<A->nseq; s1++){lc[s1][0]=s1; lc[s1][1]=0;}
  for (ncol=0,col1=0; col1< A->len_aln; col1++)
    {
      for (n=0,s1=0; s1<A->nseq; s1++)
	{
	  real_res1_col1=pos[s1][col1]-1;
	  rs1=A->order[s1][0];

	  if (real_res1_col1<0)lc[s1][1]++;
	  else if (!((CL->T[A->order[s1][0]])->ca[real_res1_col1]))lc[s1][1]++;
	  else n++;
	}
      if (n==A->nseq)
	{
	  ncol++;
	}
    }
  return ncol;
}

float square_atom ( Atom *X)
{

  return X->x*X->x + X->y*X->y + X->z*X->z;
}
Atom* reframe_atom ( Atom *X, Atom*Y, Atom *Z, Atom *IN, Atom *R)
     {
       float new_x, new_y, new_z;

       if ( R==NULL)R=(Atom*)vcalloc ( 1, sizeof (Atom));


        new_x= X->x*IN->x + Y->x*IN->y +Z->x*IN->z;
	new_y= X->y*IN->x + Y->y*IN->y +Z->y*IN->z;
	new_z= X->z*IN->x + Y->z*IN->y +Z->z*IN->z;

	R->x=new_x;
	R->y=new_y;
	R->z=new_z;
       return R;
     }

Atom* add_atom ( Atom *A, Atom*B, Atom *R)
{
  if ( R==NULL)R=(Atom*)vcalloc ( 1, sizeof (Atom));

  R->x=A->x+B->x;
  R->y=A->y+B->y;
  R->z=A->z+B->z;

  return R;
}
Atom* diff_atom ( Atom *A, Atom*B, Atom *R)
{
  if ( R==NULL)R=(Atom*)vcalloc ( 1, sizeof (Atom));

  R->x=A->x-B->x;
  R->y=A->y-B->y;
  R->z=A->z-B->z;

  return R;
}

Atom * copy_atom ( Atom *A, Atom*R)
{
  if ( R==NULL)R=(Atom*)vcalloc ( 1, sizeof (Atom));
  R->num=A->num;
  R->res_num=A->res_num;
  R->x=A->x;
  R->y=A->y;
  R->z=A->z;

  sprintf( R->type, "%s", A->type);
  return R;
}
 void print_atom (Atom *A)
{
  fprintf ( stdout, "%.2f %.2f %.2f", A->x, A->y, A->z);
}
/************************************************************************/
/*                                                                      */
/*                            NUSSINOV                                  */
/*                                                                      */
/************************************************************************/

/*---------prototypes ----------*/
static void computeBasePairMatrix(int**M,char*S,int l, int T);
static int backtrack(int a,int b,int**M,char*S,char*P, int T);



static int basePair(char x, char y)
{
  static short **mat;

  if (!mat)
    {
      char alp[20];
      int a, b, c1, c2, lc1, lc2;
      mat=declare_short (256, 256);
      sprintf ( alp, "AGCTUagctu");
      for (a=0; a<strlen (alp); a++)
	{
	  for (b=a; b<strlen (alp)-1; b++)
	    {
	      c1=alp[a];c2=alp[b];
	      lc1=tolower(c1); lc2=tolower(c2);
	      if ( lc1=='g' && lc2=='c')
		mat[c1][c2]=1;
	      else if ( lc1=='a' && lc2=='u')
		mat[c1][c2]=1;
	      else if ( lc1=='u' && lc2=='g')
		mat [c1][c2]=1;
	      mat[c2][c1]=mat[c1][c2];
	    }
	}
    }
  return (int)mat[(int)x][(int)y];
}



/* ------------------------------------------------------------ */

char *nussinov(char *S, int THRESHOLD)
{
  char *paren;
  int i;

  /*-------------------------------
    S is RNA sequence
    paren is parenthesis expression for
    optimal RNA secondary structure
    THRESHOLD: Min distance between two paired residues
    -------------------------------*/

  int **numBasePairs;
   int n;

   /*----- initialization  --*/
   n = strlen(S);
   paren=(char*)vcalloc (n+1, sizeof (char));
   numBasePairs=declare_int (n,n);

   for (i=0;i<n;i++) paren[i]='.';
   paren[n]='\0'; // paren is string of same length as S
   computeBasePairMatrix(numBasePairs,S,n, THRESHOLD);
   backtrack(0,n-1,numBasePairs,S,paren, THRESHOLD);
   free_int (numBasePairs, -1);
   return paren;
}

static void computeBasePairMatrix(int** numBasePairs,char *S,int n, int THRESHOLD)
{
   int i,j,d,k,max,val,index;

   for (d = THRESHOLD; d < n; d++){
     for(i=0; i < n; i++)
       {
        j=i+d;
        if (j < n){
          max=0;
          index=n;
           /*-------------------------------------
           if index<n at end of for-loop, then this
           means that index and j form a base pair,
           and this is noted by numBasePairs[j][i]=index.
           if index=n at end of for-loop, then this
           means that j is not base paired.
           -------------------------------------*/

          if ( numBasePairs[i][j-1]>max ){
             max = numBasePairs[i][j-1];
             index = n;
             // j not basepaired with some k such that i<k<j
          }

          val = basePair(S[i],S[j]) + numBasePairs[i+1][j-1];
          if ( j-i<= THRESHOLD && val > max ){
             max = val;
             index=i;
          }
          for(k=i; k<=j-THRESHOLD; k++){
             val = basePair(S[k],S[j]) + numBasePairs[i][k-1]
                   + numBasePairs[k+1][j-1];
             if (val > max) {
                max = val;
                index=k;
             }
          }
          numBasePairs[i][j]=max;
          if (index<n)
             numBasePairs[j][i]=index;
	  else
             numBasePairs[j][i]=-1;
       }
     }
   }

}




static int backtrack(int i, int j, int **numBasePairs,char *S, char *paren, int THRESHOLD)
{
  int k;

   k = numBasePairs[j][i];
   if(k != -1)
     {
     paren[k] = '(';
     paren[j] = ')';
     if( THRESHOLD <= (j-1)-(k+1) )
       backtrack(k+1,j-1,numBasePairs,S,paren, THRESHOLD);
     if (THRESHOLD <= k-1-i  )
       backtrack(i,k-1,numBasePairs,S,paren, THRESHOLD);
     }
   else{
     if( THRESHOLD <= j-1-i )
       {
	 backtrack(i,j-1,numBasePairs,S,paren, THRESHOLD);
       }
     else
       return 0;
   }
   return 0;}

int count;
char * rna_struc2rna_lib ( char *seq_name, char *seq, char *name)
{
  FILE *fp;
  char *st;


  st=nussinov (seq, 2);
  if ( name==NULL)name=vtmpnam(NULL);
  fp=vfopen ( name, "w");
  fprintf (fp, "! TC_LIB_FORMAT_01\n");
  fprintf (fp, "1\n%s %d %s\n", seq_name, (int)strlen (seq), seq);
  fprintf (fp, "#1 1\n");
  display_rna_ss (0, seq, st, fp);
  fprintf ( fp, "! SEQ_1_TO_N\n");
  vfclose (fp);
  vfree (st);
  //printf_system ( "cp %s test", name);
  return name;
}
int display_rna_ss ( int n, char *seq, char *st, FILE *fp)
{
  char p;
  char string[100];
  static int thread;

  while ((p=st[n])!='\0')
    {
      if ( p=='(')
	{
	  thread=count++;
	  sprintf (string, "%d",n+1);
	  n=display_rna_ss (n+1, seq, st, fp);
	  fprintf (fp, "%s %d 100\n", string, n+1);
	}
      else if (p=='.');
      else if (p==')')
	{
	  return n;
	}
      n++;
    }
  return n;
}

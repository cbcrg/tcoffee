#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <signal.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"
#include "t_coffee.h"
static void test();
static char * get_seq_type_from_cl (int argc, char **argv);
static char *get_defaults(char *buf, char *type);
static char *get_evaluate_defaults(char *buf, char *type);
static char *get_genome_defaults(char *buf, char *type);
static char *get_dali_defaults(char *buf, char *type);
static char *get_mcoffee_defaults(char *buf, char *type);
static char *get_fmcoffee_defaults(char *buf, char *type);
static char *get_t_coffee_defaults(char *buf, char *type);

static char *get_dmcoffee_defaults(char *buf, char *type);
static char *get_rcoffee_consan_defaults(char *buf, char *type);

static char *get_rmcoffee_defaults(char *buf, char *type);//Original R-Coffee Paper
static char *get_rcoffee_defaults(char *buf, char *type);//Original R-Coffee Paper
static char *get_rmcoffee_defaults_old(char *buf, char *type);//Original R-Coffee Paper
static char *get_rcoffee_defaults_old(char *buf, char *type);//Original R-Coffee Paper
static char *get_best4RNA_defaults(char *buf, char *type);

static char *get_very_fast_defaults(char *buf, char *type);
static char *get_precomputed_defaults(char *buf, char *type);
static char *get_3dcoffee_defaults(char *buf, char *type);
static char *get_expresso_defaults(char *buf, char *type);

static char *get_accurate_defaults(char *buf, char *type);
static char *get_accurate4PROTEIN_defaults(char *buf, char *type);
static char *get_accurate4DNA_defaults(char *buf, char *type);
static char *get_accurate4RNA_defaults(char *buf, char *type);

static char *get_psicoffee_defaults(char *buf, char *type);
static char *get_dna_defaults(char *buf, char *type);
static char *get_cdna_defaults(char *buf, char *type);
static char *get_repeat_defaults(char *buf, char *type);
static char *get_low_memory_defaults( char *buf, char *type);

static char *get_genepredx_defaults(char *buf, char *type);
static char *get_genepredpx_defaults(char *buf, char *type);

static int set_methods_limits (char **method_limits,int n_methods_limit,char **list_file, int n_list, int *maxnseq, int *maxlen);
static FILE *t_coffee_tip (FILE *fp,char *mode);

static int run_other_pg(int argc, char *argv[]);
static char* prepare_one2all (char *seq,Sequence *S, char *lib_file);
static char* prepare_subset2all (char *seq,Sequence *S, char *lib_file, Constraint_list *CL);

#define is_a_seq_file(file) (!is_matrix(file) && !is_matrix(file+1) && !is_method (file) && !is_method (file+1) &&(check_file_exists(file) || check_file_exists(file+1)))
static int NO_METHODS_IN_CL;
int batch_main ( int argc, char **argv);
int main (int argc, char *argv[])
{
  int r, a;

  if (argc>=2 && strcmp (argv[1], "-batch")==0)
    {
      char **list;
      list=file2lines (argv[2]);
      for (a=1; a<atoi (list[0]); a++)
	{
	  char **com;
	  com=string2list (list[a]);
	  r=batch_main (atoi (com[0])-1, com+1);
	  free_char (com, -1);
	}
    }
  else
    {
      r=batch_main (argc, argv);
    }
  myexit (r);
}
int batch_main ( int argc, char **argv)
	{

	int a, b, c, t;
	char *f;

	Sequence *S;
	char **new_order;
	char **initial_order;

	Fname *F=NULL;
	Constraint_list *CL;
	Alignment  *A=NULL, *EA=NULL;
	Alignment  **aln_list;

	FILE *OUT;


	NT_node **T=NULL;
	int tot_node;
	char *pc;
/*Parameters*/

	int check_configuration;
	int update;
	int garbage;
	int quiet;
	char *parameters;
	char *t_coffee_defaults;
	int t_coffee_defaults_flag;

	FILE *le=NULL;
	char *se_name;
	char *run_name;

	char *mem_mode;

	int do_extend;
	char *extend_mode;
	int max_n_pair;
	int  nseq_for_quadruplet;
	char **seq_name_for_quadruplet;

	int do_compact;
	int **do_list;
	int n_do;

	char *compact_mode;
	int do_clean;
	char *clean_mode;
	int do_self;
	int do_normalise;



	int n_list;
	char **list_file;
	
	char **setenv_list;
	int n_setenv;

	
	char **template_file_list;
	int n_template_file;

	char **template_mode_list;
	int n_template_mode;



	int remove_template_file;

	char **profile_template_file_list;
	int n_profile_template_file;

	int n_pdb;
	char **pdb_list;
	int pdb_start, pdb_end;
	char *pdb_name;

	char *out_lib;
	char *out_lib_mode;
	int shrink_lib;
	int relax_lib;
	int filter_lib;
	int processed_lib;

	int lib_only;
	char *outseqweight;



	char *seq_source;

	int cosmetic_penalty,gop,f_gop, gep,f_gep, nomatch;

	char *tree_file;
	char *ph_tree_file;

	char *use_tree;
	char *tree_mode;
	char *distance_matrix_mode;
	char *distance_matrix_sim_mode;

	int quicktree;
	char *out_aln;
	char **tot_out_aln;
	int maximise;
	char **out_aln_format;
	int  n_out_aln_format;
	char *infile;
	char *matrix;
	char *dp_mode;
	char *profile_mode;
	char *profile_comparison;


	int tg_mode;
	int   ktup;
	int   fasta_step;
	int   diag_threshold;
	int   diag_mode;
	char *sim_matrix;

	char *type;
	int check_type;
	int type_only;
	char *transform;
	char *outorder;
	char  *inorder;
	char *output_res_num;
	char *residue_case;
	int extra_cpu;


	char *weight;

	char *seq_weight;
	int do_align;
	char *evaluate_mode;
	char *method_evaluate_mode;
	int get_type;
	/*Post_processing*/
	int clean_aln;
	int clean_threshold;
	int clean_iteration;
	char *clean_evaluate_mode;
	/*Profile Alignment*/

	int n_seq_list;
	char **seq_list;

	char **method_list;
	int n_method_list;


	char **method_limits;
	int n_method_limits;

	char **aln_file_list;
	int n_aln_file_list;

	char **lib_file_list;
	int n_lib_file_list;

	char **profile_list;
	int n_profile_list;

	char *profile1;
	char *profile2;

	/*Domain Parameters*/
	int do_domain;
	int domain_start;
	int domain_len;
	int domain_scale;
	int domain_interactive;
	/* extended matrix analysis*/
	int do_extended_matrix;
	/*Action Parameters*/
	int do_evaluate;
	int do_genepred;
	int do_convert;
	int do_version;
	/*Genepred_score*/
	char *genepred_score;

	int maxnseq;
	int maxlen;
	/*Thread parameters*/
	int prot_min_sim;
	int prot_max_sim;
	int prot_min_cov;
	int pdb_min_sim;
	int pdb_max_sim;
	int pdb_min_cov;



	char *prot_blast_server;
	char *pdb_blast_server;


	char *pdb_db;
	char *prot_db;
	NT_node *SNL;

	/*
	char *dna_db;
	char *dna_blast_server;
	int dna_min_sim;
	int dna_max_sim;
	int dna_min_cov;
	*/


	/*Method log*/
	char * method_log;

	char **struc_to_use;
	int n_struc_to_use;

	/*Cache*/
	char * cache;

	/*align_pdb*/
	char *align_pdb_param_file;
	char *align_pdb_hasch_mode;

	/*msa_mode*/
	char *use_seqan;
	char *msa_mode;
	char *one2all;
	char *subset2all;

	int lalign_n_top;
	int iterate;
	/*split*/
	int split;
	int split_nseq_thres;
	int split_score_thres;
	char split_name[1000];
	char split_format[1000];
	Alignment *SPLIT_ALN;

	int check_pdb_status;
	/*trim*/
	int trim;
	Sequence *trim_subS=NULL;
	Sequence *trim_S=NULL;
	Sequence *SEQ_TO_KEEP=NULL;
	char *trimfile;
	char trim_format[1000];
	int clean_seq_name;
	int n_seq_to_keep;
	char **seq_to_keep;
	char **special_mode_list;
	int n_special_mode;
	char **special_mode_list1;
	int n_special_mode1;
	char **special_mode_list2;
	int n_special_mode2;
	/*dpa*/
	int dpa;
	char *dpa_master_aln;
	int dpa_min_score1;
	int dpa_min_score2;
	int dpa_maxnseq;
	int dpa_keep_tmpfile;
	int dpa_debug;
	/*error report*/
	char *full_log;

	/*Multithread*/
	char *multi_core;
	int   n_core;

	char *lib_list;
	char *prune_lib_mode;

	int no_error_report;
	int no_warning=0;
	char *tip;
	int run_local_script;
	char *plugins;
	char *email;
	char *proxy;

	char *rna_lib;
	/*over_aln*/
	char  **overaln_param;
	char    n_overaln_param;
	char  * exon_boundaries;
	char  * overaln_mode;
	char  *overaln_model;
	int    overaln_P1;
	int    overaln_P2;
	int    overaln_P3;
	int    overaln_P4;

	int     overaln_threshold;
	int     overaln_target;
	int     clean_overaln;

	
	argv=standard_initialisation (argv, &argc);
	set_string_variable ("t_coffee", argv[0]);
	
/*Running other programs via T-Coffee*/
	if (argc>=3 && strm (argv[1], "-other_pg"))
	  {
	    //standard_initialisation (NULL,NULL);
	    return run_other_pg (argc-2, argv+2);
	  }

/*PARAMETER PROTOTYPE:    READ PARAMETER FILE     */
	         get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-no_error_report"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Limit the maximum memory usage (in Megabytes). 0: no limit" ,\
			    /*Parameter*/ &no_error_report          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
					  );

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
			    /*DOC*/       "get bottom parameters" ,\
			    /*Parameter*/ &parameters          ,\
			    /*Def 1*/     "NULL"             ,\
			    /*Def 2*/     "stdin"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );


	       special_mode_list1=declare_char (100, STRING);

	       n_special_mode1=get_cl_param(			\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  100                ,\
			    /*DOC*/       "specifies a special mode: genome, quickaln, dali, 3dcoffee" ,\
			    /*Parameter*/ special_mode_list1          ,\
			    /*Def 1*/     "unspecified"             ,\
			    /*Def 2*/     "HARD_CODED"       ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );


	       special_mode_list2=declare_char (100, STRING);
	       n_special_mode2=get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-special_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  100                ,\
			    /*DOC*/       "[DEPRECATED ** -special_mode is deprected use -mode instead]" ,\
			    /*Parameter*/ special_mode_list2          ,\
			    /*Def 1*/     "unspecified"             ,\
			    /*Def 2*/     "HARD_CODED"       ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

	       special_mode_list=declare_char (n_special_mode1+n_special_mode2, STRING);
	       n_special_mode=0;
	       for (a=0; a<n_special_mode1; a++)
		 if (!strm (special_mode_list1[a], "unspecified"))sprintf ( special_mode_list[n_special_mode++], "%s", special_mode_list1[a]);
	       for (a=0; a<n_special_mode2; a++)
		 if (!strm (special_mode_list2[a], "unspecified"))sprintf ( special_mode_list[n_special_mode++], "%s", special_mode_list2[a]);
	       free_char (special_mode_list1, -1);
	       free_char (special_mode_list2, -1);

	       declare_name (t_coffee_defaults);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-t_coffee_defaults"        ,\
			    /*Flag*/      &t_coffee_defaults_flag     ,\
			    /*TYPE*/      "R_F"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "get top parameters" ,\
			    /*Parameter*/ &t_coffee_defaults          ,\
			    /*Def 1*/     "NULL"             ,\
			    /*Def 2*/     "NULL"       ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       /*PARAMETER PROTOTYPE:    -type_only: must stay here: needed by special_mode    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-type_only"        ,\
			    /*Flag*/      &type_only       ,\
			    /*TYPE*/      "FL"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "exit after checking the type and returning it to the stdout",\
			    /*Parameter*/ &type_only          ,\
			    /*Def 1*/    "0"              ,\
			    /*Def 2*/    "1"              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       /*PARAMETER PROTOTYPE:    CHECK_TYPE     */
	       get_cl_param(					\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,	\
 			    /*output*/    &le           ,	\
	 		    /*Name*/      "-check_type"  ,	\
		 	    /*Flag*/      &check_type    ,	\
			    /*TYPE*/      "FL"          ,	\
			    /*OPTIONAL?*/ OPTIONAL      ,	\
			    /*MAX Nval*/  0             ,		\
			    /*DOC*/       "Make sure that -type and the real type of the sequences agree"          , \
			    /*Parameter*/ &check_type    ,		\
			    /*Def 1*/    "0"            ,		\
			    /*Def 2*/    "1"            ,		\
			    /*Min_value*/ "any"         ,		\
			    /*Max Value*/ "any"				\
					  );
	       /*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (type);
	       get_cl_param(					\
			    /*argc*/      argc           ,	\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,	\
			    /*Name*/      "-type"        ,	\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "S"            ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,	\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "protein or dna. Automatically set, but can be forced with this flag"           , \
			    /*Parameter*/ &type          ,		\
			    /*Def 1*/    ""              ,		\
			    /*Def 2*/    ""              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
					  );
	

/* extra>prompt>special>parameters>defaults*/
 argv=break_list ( argv, &argc, "=;, \n");
 argv=merge_list ( argv, &argc);
 if (argc>1 && argv[1][0]!='-')argv=push_string ("-seq ", argv, &argc, 1);
 
 if ( name_is_in_list ("-method",argv, argc,100)==-1)
   {
     NO_METHODS_IN_CL=1;
   }
 
if (t_coffee_defaults_flag)
  {
    char *pname=NULL;

    pname=getenv ( "TCOFFEE_DEFAULTS");

    if (check_file_exists ( t_coffee_defaults))pname=t_coffee_defaults;
    else if ( getenv ( "TCOFFEE_DEFAULTS"))
      {
	pname=getenv ( "TCOFFEE_DEFAULTS");
	if (check_file_exists(pname));
	else pname=NULL;
      }
    else
      {
	declare_name(pname);sprintf (pname, "%s/.t_coffee_defaults",getenv ( "HOME") );
	if (!check_file_exists (pname)){vfree(pname);pname=NULL;}
      }

    if (pname)
      {
	argv=push_string (file2string(pname), argv, &argc, 1);
	t_coffee_defaults=pname;
      }
    else
      {
	t_coffee_defaults=NULL;
      }
  }
 
 if ( parameters && parameters[0])argv=push_string (file2string (parameters), argv, &argc, 1);


 if (n_special_mode && !type_only)
   {
     char *special_mode;
     char *lseq_type;
     declare_name(lseq_type);
     if (type && !strm (type, ""))
       sprintf (lseq_type,"%s",type);
     else
       sprintf (lseq_type,"%s",get_seq_type_from_cl (argc, argv));
     
     for ( a=0; a< n_special_mode; a++)
       {
	 char *new_arg=NULL;

	 special_mode=special_mode_list[a];

	 store_mode (special_mode);


	 if (special_mode && !special_mode[0]);
	 else if ( strm (special_mode, "genepredx"))new_arg=get_genepredx_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "genepredpx"))new_arg=get_genepredpx_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "regular") || strm (special_mode, "regular_fast")|| strm (special_mode, "default"))new_arg=get_defaults (NULL,lseq_type);
	 else if ( strm (special_mode, "genome"))new_arg=get_genome_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "quickaln"))new_arg=get_very_fast_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "dali"))new_arg=get_dali_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "evaluate"))new_arg=get_evaluate_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "precomputed"))new_arg=get_precomputed_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "3dcoffee"))new_arg=get_3dcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "expresso"))new_arg=get_expresso_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "repeats"))new_arg=get_repeat_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "psicoffee"))new_arg=get_psicoffee_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "accurate") || strm (special_mode, "accurate_slow") || strm (special_mode, "psicoffee_expresso"))new_arg=get_accurate_defaults(NULL, lseq_type);
	 else if ( strm (special_mode, "accurate4DNA"))new_arg=get_accurate4DNA_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "accurate4RNA"))new_arg=get_accurate4RNA_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "best4RNA"))new_arg=get_best4RNA_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "accurate4PROTEIN"))new_arg=get_accurate4PROTEIN_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "low_memory") || strm (special_mode, "memory"))new_arg=get_low_memory_defaults(NULL,lseq_type);


	 else if ( strm (special_mode, "dna"))new_arg=get_dna_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "cdna"))new_arg=get_dna_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "protein"))new_arg=get_low_memory_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "mcoffee"))new_arg=get_mcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "dmcoffee"))new_arg=get_dmcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "fmcoffee"))new_arg=get_fmcoffee_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "rcoffee_consan"))new_arg=get_rcoffee_consan_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rmcoffee") ||strm (special_mode, "mrcoffee") )new_arg=get_rmcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rcoffee"))new_arg=get_rcoffee_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "rcoffee_slow_accurate"))new_arg=get_rcoffee_consan_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rcoffee_fast_approximate"))new_arg=get_rmcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "t_coffee"))new_arg=get_t_coffee_defaults(NULL,lseq_type);


	 else if ( strm (special_mode, "unspecified"));
	 else
	   {
	     fprintf ( stderr, "\nERROR: special_mode %s is unknown [FATAL:%s]\n",special_mode, PROGRAM);
	     myexit (EXIT_FAILURE);
	   }

	 if (new_arg)argv=push_string (new_arg, argv, &argc, 1);
       }
   }

if ( getenv ("TCOFFEE_EXTRA_PARAM"))argv=push_string (getenv ("TCOFFEE_EXTRA_PARAM"), argv, &argc, argc);


argv=break_list ( argv, &argc, "=;, \n");
argv=merge_list ( argv, &argc);
/*check_cl4t_coffee ( argc, argv); */


/*PARAMETER PROTOTYPE:    VERSION      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-version"        ,\
			    /*Flag*/      &do_version        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "forces the program to output the version number and exit" ,\
			    /*Parameter*/ &do_version          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
	       


/*PARAMETER PROTOTYPE:    DO EVALUATE      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-score"        ,\
			    /*Flag*/      &do_evaluate        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "DEPRECATED: use -special_mode evaluate instead " ,\
			    /*Parameter*/ &do_evaluate          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
if ( !do_evaluate)
  {
/*PARAMETER PROTOTYPE:    DO EVALUATE      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-evaluate"        ,\
			    /*Flag*/      &do_evaluate        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Use -special_mode evaluate for a default behavior " ,\
			    /*Parameter*/ &do_evaluate          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
  }
/*PARAMETER PROTOTYPE:    DO EVALUATE      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-genepred"        ,\
			    /*Flag*/      &do_genepred        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Use -special_mode genepred for a default behavior " ,\
			    /*Parameter*/ &do_genepred          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
/*PARAMETER PROTOTYPE:    DO FORMAT      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-convert"        ,\
			    /*Flag*/      &do_convert        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "forces the program to make a conversion" ,\
			    /*Parameter*/ &do_convert          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );



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
			    /*DOC*/       "Defines the file in which the log output is written"          ,\
			    /*Parameter*/ &se_name      ,\
			    /*Def 1*/     "stderr"      ,\
			    /*Def 2*/     "/dev/null"   ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
     if (type_only==1)sprintf ( se_name, "/dev/null");

     /*PARAMETER PROTOTYPE:    DO FORMAT      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-check_configuration"        ,\
			    /*Flag*/      &check_configuration       ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "checks that the required programs are installed" ,\
			    /*Parameter*/ &check_configuration          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
			    );
	  /*PARAMETER PROTOTYPE:    UPDATE      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-update"        ,\
			    /*Flag*/      &update,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "checks the existence of an updated version" ,\
			    /*Parameter*/ &update          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
			    );



     if ( check_configuration)
       {

	 check_configuration4program();
	 return EXIT_SUCCESS;

       }
     if ( update)
       {
	 myexit (check_for_update(DISTRIBUTION_ADDRESS));
       }
     if ( do_version)
       {
	 fprintf ( stdout, "PROGRAM: %s (%s)\n",PROGRAM,VERSION);
	 return EXIT_SUCCESS;
       }


     le=vfopen ( se_name, "w");
     fprintf ( le, "\nPROGRAM: %s (%s)\n",PROGRAM,VERSION);

/*PARAMETER PROTOTYPE: RUN NAME*/
               declare_name (full_log);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-full_log"   ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Sets the prefix of all the output files"          ,\
			    /*Parameter*/ &full_log     ,\
			    /*Def 1*/     "NULL"        ,\
			    /*Def 2*/     "full_log"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       vremove(full_log);
/*PARAMETER PROTOTYPE: RUN NAME*/
               declare_name (genepred_score);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-genepred_score"   ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "nsd,tot, <seq_name>"          ,\
			    /*Parameter*/ &genepred_score     ,\
			    /*Def 1*/     "nsd"        ,\
			    /*Def 2*/     ""            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE: RUN NAME*/
               declare_name (run_name);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-run_name"   ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Sets the prefix of all the output files"          ,\
			    /*Parameter*/ &run_name     ,\
			    /*Def 1*/     "NULL"        ,\
			    /*Def 2*/     ""            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE: MEM MODE*/
	       declare_name(mem_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-mem_mode"   ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Deprecated" ,\
			    /*Parameter*/ &mem_mode     ,\
			    /*Def 1*/     "mem"         ,\
			    /*Def 2*/     ""            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE: EXTEND  */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-extend"     ,\
			    /*Flag*/      &do_extend    ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Do Library Extention On the Fly"          ,\
			    /*Parameter*/ &do_extend    ,\
			    /*Def 1*/     "1"           ,\
			    /*Def 2*/     "1"           ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE: EXTEND  */
	       declare_name (extend_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-extend_mode"     ,\
			    /*Flag*/      &garbage    ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Library extension mode"          ,\
			    /*Parameter*/ &extend_mode    ,\
			    /*Def 1*/     "very_fast_triplet"           ,\
			    /*Def 2*/     ""           ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       /*PARAMETER PROTOTYPE: EXTEND  */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-max_n_pair"     ,\
			    /*Flag*/      &garbage    ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Indicates the Number of Pairs to Compare when making prf Vs prf. 0<=>every pair "          ,\
			    /*Parameter*/ &max_n_pair    ,\
			    /*Def 1*/     "10"           ,\
			    /*Def 2*/     "3"           ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE: SEQUENCES TO EXTEND */
	seq_name_for_quadruplet=declare_char ( 200, STRING);
	nseq_for_quadruplet=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-seq_name_for_quadruplet"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "Indicates which sequence must be used to compute quadruplets"          ,\
			    /*Parameter*/ seq_name_for_quadruplet    ,\
			    /*Def 1*/     "all",\
			    /*Def 2*/     ""       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE: COMPACT */
	       declare_name (compact_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-compact"    ,\
			    /*Flag*/      &do_compact   ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Deprecated"          ,\
			    /*Parameter*/ &compact_mode ,\
			    /*Def 1*/     "default"      ,\
			    /*Def 2*/     "default"      ,\
			    /*Min_value*/ "0"           ,\
			    /*Max Value*/ "1"           \
		   );


/*PARAMETER PROTOTYPE:        CLEAN*/
	       declare_name ( clean_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-clean"      ,\
			    /*Flag*/      &do_clean     ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Deprecated"          ,\
			    /*Parameter*/ &clean_mode   ,\
			    /*Def 1*/     "no"          ,\
			    /*Def 2*/     "shadow"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:        DO SELF */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-do_self"    ,\
			    /*Flag*/      &do_self      ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "Make self extension. Used by Mocca"          ,\
			    /*Parameter*/ &do_self      ,\
			    /*Def 1*/     "0"           ,\
			    /*Def 2*/     "1"           ,\
			    /*Min_value*/ "0"           ,\
			    /*Max Value*/ "1"           \
		   );

/*PARAMETER PROTOTYPE:        DO NORMALISE */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-do_normalise"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Normalisation factor when computing scores"          ,\
			    /*Parameter*/ &do_normalise ,\
			    /*Def 1*/     "1000"           ,\
			    /*Def 2*/     "1000"       ,\
			    /*Min_value*/ "-10000"           ,\
			    /*Max Value*/ "10000"           \
		   );
/*PARAMETER PROTOTYPE:        IN */
	       template_file_list=declare_char (100, STRING);
	       n_template_file=get_cl_param(			\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-template_file"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1000           ,\
			    /*DOC*/       "List of templates file for the sequences",\
			    /*Parameter*/ template_file_list     ,	\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );

/*PARAMETER PROTOTYPE:    VERSION      */
	       setenv_list=declare_char (100, STRING);
	       n_setenv=get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-setenv"        ,\
			    /*Flag*/      &do_version        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  2                ,\
			    /*DOC*/       "Declares a parameter variable" ,\
			    /*Parameter*/ setenv_list          ,	\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
					  );
/*PARAMETER PROTOTYPE:        IN */
	       template_mode_list=declare_char (100, STRING);
	       n_template_mode=get_cl_param(			\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-template_mode"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1000           ,\
			    /*DOC*/       "List of template procedures",\
			    /*Parameter*/ template_mode_list     ,	\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );
	       for (a=0; a<n_template_mode; a++)
		 {
		   sprintf (template_file_list[n_template_file++], "%s", template_mode_list [a]);
		 }
/*PARAMETER PROTOTYPE:        remove_template_file*/

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-remove_template_file"      ,\
			    /*Flag*/      &remove_template_file     ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Remove all the template files: 0 keep all, 1: only remove the template files 2: remove template files AND template lists "          ,\
			    /*Parameter*/ &remove_template_file   ,\
			    /*Def 1*/     "0"          ,\
			    /*Def 2*/     "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
	       set_int_variable ("remove_template_file", remove_template_file);

	/*PARAMETER PROTOTYPE:        IN */
	       profile_template_file_list=declare_char (100, STRING);
	       n_profile_template_file=get_cl_param(			\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-profile_template_file"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1000           ,\
			    /*DOC*/       "List of templates files asscoaciated with profiles",\
			    /*Parameter*/ profile_template_file_list     , \
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );

	       /*PARAMETER PROTOTYPE:        IN */
	list_file=declare_char (2000, STRING);
	n_list=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-in"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  2000           ,\
			    /*DOC*/       "Reads the Ssequences, Mmethods, Llibraries,Xmatrices,Rprofiles,Pstructures,AAlignments"          ,\
			    /*Parameter*/ list_file     ,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
					  );


	/*PARAMETER PROTOTYPE:        IN */
	seq_list=declare_char (1000, STRING);
	n_seq_list=get_cl_param(				\
				/*argc*/      argc          ,	\
				/*argv*/      argv          ,	\
				/*output*/    &le           ,	\
				/*Name*/      "-seq"         ,	\
				/*Flag*/      &garbage      ,	\
				/*TYPE*/      "S"           ,	\
				/*OPTIONAL?*/ OPTIONAL      ,	\
				/*MAX Nval*/  1000           ,		\
				/*DOC*/       "List of sequences in any acceptable format", \
				/*Parameter*/ seq_list     ,		\
				/*Def 1*/    "",			\
				/*Def 2*/     "stdin"       ,		\
				/*Min_value*/ "any"         ,		\
				/*Max Value*/ "any"			\
					      );
	aln_file_list=declare_char (1000, STRING);
	n_aln_file_list=get_cl_param(				\
				/*argc*/      argc          ,	\
				/*argv*/      argv          ,	\
				/*output*/    &le           ,	\
				/*Name*/      "-aln"         ,	\
				/*Flag*/      &garbage      ,	\
				/*TYPE*/      "S"           ,	\
				/*OPTIONAL?*/ OPTIONAL      ,	\
				/*MAX Nval*/  1000           ,		\
				/*DOC*/       "List of sequences in any acceptable format", \
				/*Parameter*/ aln_file_list     ,		\
				/*Def 1*/    "",			\
				/*Def 2*/     "stdin"       ,		\
				/*Min_value*/ "any"         ,		\
				/*Max Value*/ "any"			\
					      );
	method_limits=declare_char (1000, STRING);
	n_method_limits=get_cl_param(				\
				/*argc*/      argc          ,	\
				/*argv*/      argv          ,	\
				/*output*/    &le           ,	\
				/*Name*/      "-method_limits"         ,	\
				/*Flag*/      &garbage      ,	\
				/*TYPE*/      "S"           ,	\
				/*OPTIONAL?*/ OPTIONAL      ,	\
				/*MAX Nval*/  1000           ,		\
				/*DOC*/       "List of limits for selected methods: method maxnseq maxlen (-1 = nolimit)", \
				/*Parameter*/ method_limits     ,		\
				/*Def 1*/    "",			\
				/*Def 2*/     ""       ,		\
				/*Min_value*/ "any"         ,		\
				/*Max Value*/ "any"			\
					      );
	method_list=declare_char (1000, STRING);
	n_method_list=get_cl_param(				\
				/*argc*/      argc          ,	\
				/*argv*/      argv          ,	\
				/*output*/    &le           ,	\
				/*Name*/      "-method"         ,	\
				/*Flag*/      &garbage      ,	\
				/*TYPE*/      "S"           ,	\
				/*OPTIONAL?*/ OPTIONAL      ,	\
				/*MAX Nval*/  1000           ,		\
				/*DOC*/       "List of sequences in any acceptable format", \
				/*Parameter*/ method_list     ,		\
				/*Def 1*/    "",			\
				/*Def 2*/     ""       ,		\
				/*Min_value*/ "any"         ,		\
				/*Max Value*/ "any"			\
					      );
	lib_file_list=declare_char (1000, STRING);
	n_lib_file_list=get_cl_param(				\
				/*argc*/      argc          ,	\
				/*argv*/      argv          ,	\
				/*output*/    &le           ,	\
				/*Name*/      "-lib"         ,	\
				/*Flag*/      &garbage      ,	\
				/*TYPE*/      "S"           ,	\
				/*OPTIONAL?*/ OPTIONAL      ,	\
				/*MAX Nval*/  1000           ,		\
				/*DOC*/       "List of sequences in any acceptable format", \
				/*Parameter*/ lib_file_list     ,		\
				/*Def 1*/    "",			\
				/*Def 2*/     "stdin"       ,		\
				/*Min_value*/ "any"         ,		\
				/*Max Value*/ "any"			\
					      );
	profile_list=declare_char ( 2000, STRING);
	n_profile_list=get_cl_param(					\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-profile"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  2000                ,\
			    /*DOC*/       "Input one or many MSA that will be treated as profiles" ,\
			    /*Parameter*/ profile_list          ,\
			    /*Def 1*/    ""             ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
					  );
	declare_name (profile1);
	get_cl_param(						\
		     /*argc*/      argc             ,		\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-profile1"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Input one profile (ClustalW option)" ,\
			    /*Parameter*/ &profile1          ,\
			    /*Def 1*/    ""             ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	declare_name (profile2);
	get_cl_param(						\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-profile2"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Input a profile (ClustalW option)" ,\
			    /*Parameter*/ &profile2          ,\
			    /*Def 1*/    ""             ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

	pdb_list=declare_char ( 200, STRING);
	n_pdb=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-pdb"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "Reads/fetch a pdb file: PDBID(PDB_CHAIN)[opt] (FIRST,LAST)[opt],"          ,\
			    /*Parameter*/ pdb_list     ,\
			    /*Def 1*/    "",\
			    /*Def 2*/     ""       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    OUT_LIB     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-relax_lib"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "self extend the library, without adding new positions", \
			    /*Parameter*/&relax_lib ,\
			    /*Def 1*/    "1"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-filter_lib"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Removes from the library every value below the threshold",\
			    /*Parameter*/&filter_lib ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "10"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    SHRINK_LIB     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-shrink_lib"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Runks linked_pairwise on the lib to remove every useless diagonal"          ,\
			    /*Parameter*/&shrink_lib       ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */
	declare_name (out_lib);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-out_lib"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Prompts the program to write the computed library file"          ,\
			    /*Parameter*/&out_lib       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       /*PARAMETER PROTOTYPE:    OUT_LIB_MODE     */
	       declare_name (out_lib_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-out_lib_mode"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Save the primary or the extended library:[primary|extende]extended_[pair|lib]_[raw|pc]"          ,\
			    /*Parameter*/&out_lib_mode       ,\
			    /*Def 1*/    "primary"      ,\
			    /*Def 2*/    "extended",\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    LIB_ONLY     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-lib_only"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Only Compute the library",\
			    /*Parameter*/&lib_only       ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    OUT_LIB     */
	declare_name (outseqweight);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-outseqweight"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Prompts the program to write the sequuence weight values"          ,\
			    /*Parameter*/&outseqweight       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       /*PARAMETER PROTOTYPE:    DPA    */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
 			    /*output*/    &le           ,\
	 		    /*Name*/      "-dpa"  ,\
		 	    /*Flag*/      &dpa    ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "Use DPA mode"          ,\
			    /*Parameter*/ &dpa    ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    SEQ TO ALIGN     */
	       declare_name (seq_source);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-seq_source",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Indicates the files that will be used as sequence sources, important for dpa. With the default mode alignments must be provided with the Sflag as well as tye Aflag if they contribute novel sequences",\
			    /*Parameter*/ &seq_source   ,\
			    /*Def 1*/    "ANY"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

/*PARAMETER PROTOTYPE:    COSMETIC PENALTY     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-cosmetic_penalty" ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "A very low Gap Opening Penalty.It only affects the non stable portions of the alignmnent.Negative values penalize gaps, positive values reward them"          ,\
			    /*Parameter*/ &cosmetic_penalty          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "0"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    GAPOPEN     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-gapopen"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Gap opening penalty. Must be negative, best matches get a score of 1000"          ,\
			    /*Parameter*/ &gop          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    GAPEXT     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-gapext"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Gap Extension Penalty. Positive values give rewards to gaps and prevent the alignment of unrelated segments"          ,\
			    /*Parameter*/ &gep          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:    F_GAPOPEN     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-fgapopen"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Deprecated"          ,\
			    /*Parameter*/ &f_gop          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    F_GAPEXT     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-fgapext"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Deprecated"          ,\
			    /*Parameter*/ &f_gep          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
/*PARAMETER PROTOTYPE:    NEW_TREE     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-nomatch"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Deprecated"          ,\
			    /*Parameter*/ &nomatch          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "0"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
	       declare_name ( tree_file);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-newtree"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Name of the output guide tree"          ,\
			    /*Parameter*/&tree_file     ,\
			    /*Def 1*/    "default"      ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
	       declare_name ( ph_tree_file);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-tree"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Name of the output guide tree"          ,\
			    /*Parameter*/&ph_tree_file     ,\
			    /*Def 1*/    "NO"      ,\
			    /*Def 2*/    "default"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:    USETREE     */
	       declare_name ( use_tree);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-usetree",\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "R_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Use an existing guide tree"          ,\
			    /*Parameter*/ &use_tree     ,\
			    /*Def 1*/    "NULL"         ,\
			    /*Def 2*/    "NULL"         ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
	       /*PARAMETER PROTOTYPE:    */
	       declare_name ( tree_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-tree_mode"  ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "nj, upgma, cwph",\
			    /*Parameter*/ &tree_mode    ,\
			    /*Def 1*/    "nj"         ,\
			    /*Def 2*/    "nj"         ,\
			    /*Min_value*/ "1"         ,\
			    /*Max Value*/ "1"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */
	       declare_name ( distance_matrix_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-distance_matrix_mode"  ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Computation of the distances for the tree: slow, fast, very_fast, ktup"          ,\
			    /*Parameter*/ &distance_matrix_mode    ,\
			    /*Def 1*/    "ktup"         ,\
			    /*Def 2*/    "idscore"         ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */
	       declare_name ( distance_matrix_sim_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-distance_matrix_sim_mode"  ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Choice of the distance measure: <mat>_sim1, _sim2, _sim3, _cov, _gap"          ,\
			    /*Parameter*/ &distance_matrix_sim_mode    ,\
			    /*Def 1*/    "idmat_sim1"         ,\
			    /*Def 2*/    "idmat_sim1"         ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
 			    /*output*/    &le           ,\
	 		    /*Name*/      "-quicktree"  ,\
		 	    /*Flag*/      &quicktree    ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "Use distance_matrix_mode=very_fast"          ,\
			    /*Parameter*/ &quicktree    ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       if ( quicktree)sprintf ( distance_matrix_mode, "very_fast");
/*PARAMETER PROTOTYPE:    OUTFILE     */
	       declare_name ( out_aln);
	       tot_out_aln=declare_char (200, STRING);
	        get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-outfile"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Name of the output alignment"          ,\
			    /*Parameter*/ &out_aln      ,\
			    /*Def 1*/    "default"      ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    MAXIMISE     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maximise"   ,\
			    /*Flag*/      &maximise     ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "Deprecated"          ,\
			    /*Parameter*/ &maximise     ,\
			    /*Def 1*/    "1"            ,\
			    /*Def 2*/    "1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
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
			    /*DOC*/       "Specifies one or many formats that must be output: clustalw_aln, msf_aln. The file extension is the output format"           ,\
			    /*Parameter*/ out_aln_format,\
			    /*Def 1*/    "aln,html"           ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (infile);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-infile"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "R_F"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "input a pre-computed alignment, or a file to reformat"           ,\
			    /*Parameter*/ &infile        ,\
			    /*Def 1*/    ""              ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (matrix);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-matrix"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Specifies the substitution matrix.",\
			    /*Parameter*/ &matrix        ,\
			    /*Def 1*/    "default"              ,\
			    /*Def 2*/    "default"    ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    TG_MODE    */

	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-tg_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "0: Penalise Term gap with gapopen and gapext\n1: gapopen only\n2: No penalty\n",\
			    /*Parameter*/ &tg_mode        ,\
			    /*Def 1*/    "1",\
			    /*Def 2*/    "0",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    DP_MODE    */
	       declare_name (profile_mode);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-profile_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Function used to compute profile2profile scores",\
			    /*Parameter*/ &profile_mode        ,\
			    /*Def 1*/    "cw_profile_profile",\
			    /*Def 2*/    "cw_profile_profile",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

	       declare_name (profile_comparison);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-profile_comparison"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Method used to compare two profiles: full<N>: compares <every | N best> pair of sequence and every pair of structure if a structure method is used,profile: compares only the profiles.  ",\
			    /*Parameter*/ &profile_comparison        ,\
			    /*Def 1*/    "profile",\
			    /*Def 2*/    "full50",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    DP_MODE    */
	       declare_name (dp_mode);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-dp_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Type of alignment algorithm used by T-Coffee: gotoh_pair_wise, myers_millers_pair_wise, "           ,\
			    /*Parameter*/ &dp_mode        ,\
			    /*Def 1*/    "linked_pair_wise",\
			    /*Def 2*/    "cfasta_pair_wise",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    KTUP    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-ktuple"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Word size when using the heursitic dynamic programming modes fasta_pair_wise and cfasta_pair_wise "           ,\
			    /*Parameter*/ &ktup          ,\
			    /*Def 1*/    "1",\
			    /*Def 2*/    "1",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    FASTA_STEP    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-ndiag"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Number of diagonals to consider when using the heursitic dynamic programming modes fasta_pair_wise and cfasta_pair_wise"           ,\
			    /*Parameter*/ &fasta_step          ,\
			    /*Def 1*/    "0",\
			    /*Def 2*/    "10",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    FASTA_STEP    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-diag_threshold"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &diag_threshold ,\
			    /*Def 1*/    "0",\
			    /*Def 2*/    "10",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    diag_mode    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-diag_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "0: Use the whole Diag\n1: Use the best match\n"           ,\
			    /*Parameter*/ &diag_mode          ,\
			    /*Def 1*/    "0",\
			    /*Def 2*/    "1",
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    SIM_MATRIX    */
	       declare_name (sim_matrix);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-sim_matrix"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Degenerated matrix used to compute a similarity"           ,\
			    /*Parameter*/ &sim_matrix        ,\
			    /*Def 1*/    "vasiliky",\
			    /*Def 2*/    "idmat",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (transform);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-transform"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "dna2rna, rna2dna, dna2prot",	\
			    /*Parameter*/ &transform          ,\
			    /*Def 1*/    ""              ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (outorder);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-outorder"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Specifies the order of the sequences in the msa: input or aligned"           ,\
			    /*Parameter*/ &outorder       ,\
			    /*Def 1*/    "input"          ,\
			    /*Def 2*/    "input"        ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (inorder);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-inorder"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "aligned: sort the sequences in alphabetic order before starting thus making the input order irrelevant but delivering a library in arbitratry order, keep: input order is used in the library but results become input order dependant"           ,\
			    /*Parameter*/ &inorder       ,\
			    /*Def 1*/    "aligned"          ,\
			    /*Def 2*/    "input"        ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
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
			    /*DOC*/       "Adds Residue Numbers to the MSA"           ,\
			    /*Parameter*/ &output_res_num ,\
			    /*Def 1*/    "off"            ,\
			    /*Def 2*/    "on"             ,\
			    /*Min_value*/ "any"           ,\
			    /*Max Value*/ "any"            \
		   );
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (residue_case);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-case"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Causes the case to be: kept:lower:upper."           ,\
			    /*Parameter*/ &residue_case         ,\
			    /*Def 1*/    "keep"            ,\
			    /*Def 2*/    "upper"             ,\
			    /*Min_value*/ "any"           ,\
			    /*Max Value*/ "any"            \
		   );

/*PARAMETER PROTOTYPE:    CPU     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-cpu"        ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Makes it possible to add a pre-specified amount of cpu time to the measured usage"          ,\
			    /*Parameter*/ &extra_cpu    ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "0"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    MAXNSEQ     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maxnseq"        ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Maximum number of sequences (-1=no max)"          ,\
			    /*Parameter*/ &maxnseq    ,\
			    /*Def 1*/    "1000"            ,\
			    /*Def 2*/    "0"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    MAXLEN     */

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maxlen"        ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Maximum length of a sequence (-1=no max)"          ,\
			    /*Parameter*/ &maxlen    ,\
			    /*Def 1*/    "-1"            ,\
			    /*Def 2*/    "-1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );


/*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name ( weight);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-weight"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Defines the library weight: sim OR  sim_(matrix) OR winsim" ,\
			    /*Parameter*/ &weight          ,\
			    /*Def 1*/    "default"             ,\
			    /*Def 2*/    "sim"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	  /*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name ( seq_weight);

	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-seq_weight"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Defines the sequences weighting scheme t_coffee" ,\
			    /*Parameter*/ &seq_weight          ,\
			    /*Def 1*/    "t_coffee"             ,\
			    /*Def 2*/    "t_coffee"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

/*PARAMETER PROTOTYPE:    DO ALIGN      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-align"        ,\
			    /*Flag*/      &do_align        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "forces the program to make the alignment" ,\
			    /*Parameter*/ &do_align          ,\
			    /*Def 1*/     "1"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
/*PARAMETER PROTOTYPE:    DO DOMAIN      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-mocca"        ,\
			    /*Flag*/      &do_domain        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "forces the program to extract domains" ,\
			    /*Parameter*/ &do_domain          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       if ( !do_domain)
		 {
/*PARAMETER PROTOTYPE:    DO DOMAIN      */
		   get_cl_param(				\
			    /*argc*/      argc             ,	\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-domain"        ,\
			    /*Flag*/      &do_domain        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "forces the program to extract domains" ,\
			    /*Parameter*/ &do_domain          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
					  );
		 }
/*PARAMETER PROTOTYPE:    Domain Param      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-start"        ,\
			    /*Flag*/      &domain_start        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "start of the master domain in the mocca mode" ,\
			    /*Parameter*/ &domain_start          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-len"        ,\
			    /*Flag*/      &domain_len        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "length of the master domain in the mocca mode" ,\
			    /*Parameter*/ &domain_len          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-scale"        ,\
			    /*Flag*/      &domain_scale        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Decreases the t_coffee score by Scale, so that non match get negative values" ,\
			    /*Parameter*/ &domain_scale          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-mocca_interactive"        ,\
			    /*Flag*/      &domain_interactive        ,\
			    /*TYPE*/      "FL"             ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "Runs Mocca in an interactive manneer" ,\
			    /*Parameter*/ &domain_interactive,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
/*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name (method_evaluate_mode);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-method_evaluate_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Specifies which method should be used to evaluate the score at the pairwise level" ,\
			    /*Parameter*/ &method_evaluate_mode          ,\
			    /*Def 1*/    "default"             ,\
			    /*Def 2*/    "default"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
					  );
 /*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name (evaluate_mode);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-evaluate_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Mode used to produce the color output:t_coffee_fast,t_coffee_slow  " ,\
			    /*Parameter*/ &evaluate_mode          ,\
			    /*Def 1*/    "t_coffee_fast"             ,\
			    /*Def 2*/    "dali"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-get_type"        ,\
			    /*Flag*/      &get_type        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "forces t_coffee top get the type of the sequences" ,\
			    /*Parameter*/ &get_type          ,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_aln"        ,\
			    /*Flag*/      &clean_aln        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Forces weak portion of aln to be realigned" ,\
			    /*Parameter*/ &clean_aln          ,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_threshold"        ,\
			    /*Flag*/      &clean_threshold        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Threshold for the portions of the MSA that will are realigned by '-clean_evaluate_mode'. The threshold refers to the CORE score set by '-evaluate_mode'" ,\
			    /*Parameter*/ &clean_threshold          ,\
			    /*Def 1*/     "1"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_iteration"        ,\
			    /*Flag*/      &clean_iteration        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Number of rounds for '-clean_aln'" ,\
			    /*Parameter*/ &clean_iteration          ,\
			    /*Def 1*/     "1"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
declare_name (clean_evaluate_mode);
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_evaluate_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Mode used to score residues (see evaluate_mode)" ,\
			    /*Parameter*/ &clean_evaluate_mode          ,\
			    /*Def 1*/    "t_coffee_fast"             ,\
			    /*Def 2*/    "t_coffee_fast"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

/*PARAMETER PROTOTYPE:    DO EXTENDED MATRIX      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-extend_matrix"        ,\
			    /*Flag*/      &do_extended_matrix        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "Deprecated" ,\
			    /*Parameter*/ &do_extended_matrix          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
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
 // HERE ("%s", blast_server);
 if ( strm (prot_blast_server, "env"))prot_blast_server=get_env_variable ("blast_server_4_TCOFFEE",IS_FATAL);
 set_string_variable ("blast_server", prot_blast_server);



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
 if ( strm (pdb_db, "env"))pdb_db=get_env_variable ("pdb_db_4_TCOFFEE", IS_FATAL);
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
 if ( strm (prot_db, "env"))prot_db=get_env_variable ("protein_db_4_TCOFFEE", IS_FATAL);
 set_string_variable ("prot_db", prot_db);

 declare_name (method_log);
 get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-method_log"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&method_log       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "default"      ,\
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
			    /*DOC*/       "Specifies the structures that must be used when combining sequences and structures. The default is to use all the structures."          ,\
			    /*Parameter*/ struc_to_use     ,\
			    /*Def 1*/    "",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
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
			    /*DOC*/       "Specifies that a cache must be used to save the structures and their comparison, as well as the blast searches.\navailable modes are: use,ignore,update,local, directory name"          ,\
			    /*Parameter*/ &cache       ,\
			    /*Def 1*/    "use"      ,\
			    /*Def 2*/    "update"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name (align_pdb_param_file);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-align_pdb_param_file"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "parameter_file"          ,\
			    /*Parameter*/ &align_pdb_param_file       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "no"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name (align_pdb_hasch_mode);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-align_pdb_hasch_mode"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "parameter_file"          ,\
			    /*Parameter*/ &align_pdb_hasch_mode       ,\
			    /*Def 1*/    "hasch_ca_trace_bubble"      ,\
			    /*Def 2*/    "hasch_ca_trace_bubble"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name (use_seqan);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-external_aligner"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Use seqan to compute the MSA",\
			    /*Parameter*/ &use_seqan      ,\
			    /*Def 1*/    "NO"      ,\
			    /*Def 2*/    "seqan_tcoffee"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name (msa_mode);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-msa_mode"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Algorithm used to compute the MSA: tree | graph"          ,\
			    /*Parameter*/ &msa_mode      ,\
			    /*Def 1*/    "tree"      ,\
			    /*Def 2*/    "tree"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name (one2all);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-one2all"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Align all the sequences to the master sequence"          ,\
			    /*Parameter*/ &one2all      ,\
			    /*Def 1*/    "NULL"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name (subset2all);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-subset2all"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Align all the sequences to the master sequence"          ,\
			    /*Parameter*/ &subset2all      ,\
			    /*Def 1*/    "NULL"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-lalign_n_top"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Number of local alignments reported by the local method (lalign) when building the library"          ,\
			    /*Parameter*/ &lalign_n_top      ,\
			    /*Def 1*/    "10"      ,\
			    /*Def 2*/    "10"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-iterate"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "NUmber of iteration on the progressive alignment [0: no iteration, -1: Nseq iterations]",\
			    /*Parameter*/ &iterate      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "100"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-trim"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "trim dataset",\
			    /*Parameter*/ &trim      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-split"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "split dataset",\
			    /*Parameter*/ &split      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
declare_name(trimfile);
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-trimfile"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "trim dataset filename",\
			    /*Parameter*/ &trimfile      ,\
			    /*Def 1*/    "default"      ,\
			    /*Def 2*/    ""      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-split"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "split dataset",\
			    /*Parameter*/ &split      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

if (trim && !split)split=trim;

get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-split_nseq_thres"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Maximum Number of sequences within a subgroup",\
			    /*Parameter*/ &split_nseq_thres      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-split_score_thres"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Minimum score within a split dataset",\
			    /*Parameter*/ &split_score_thres      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-check_pdb_status"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Reports the existance of a PDB file",\
			    /*Parameter*/ &check_pdb_status      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-clean_seq_name"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Remove Special Char from sequence names",\
			    /*Parameter*/ &clean_seq_name      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );


/*PARAMETER PROTOTYPE:    SEQ TO ALIGN     */
	       seq_to_keep=declare_char ( 2000, STRING);
	       n_seq_to_keep=get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-seq_to_keep",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  200              ,\
			    /*DOC*/       "File containing the name of the sequences to keep when triming OR a list of names)",\
			    /*Parameter*/ seq_to_keep   ,\
			    /*Def 1*/    "NULL"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*******************************************************************************************************/
/*                                                                                                     */
/*                           TCoffee_dpa Parameter:START                                               */
/*                                                                                                     */
/*******************************************************************************************************/
/*PARAMETER PROTOTYPE:    dpa_master_aln     */
	       declare_name (dpa_master_aln);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-dpa_master_aln",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Approximate Alignment: File|method",\
			    /*Parameter*/ &dpa_master_aln   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       /*PARAMETER PROTOTYPE:    dpa_maxnseq    */

	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-dpa_maxnseq",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Maximum number of sequences to be aligned with DPA",\
			    /*Parameter*/ &dpa_maxnseq   ,\
			    /*Def 1*/    "0"       ,\
			    /*Def 2*/    "50"              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    dpa_min_score1    */

	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-dpa_min_score1",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  200              ,\
			    /*DOC*/       "Minimum percent ID to merge sequences in the approximate alignment",\
			    /*Parameter*/ &dpa_min_score1   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );
/*PARAMETER PROTOTYPE:    dpa_min_score2    */

	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-dpa_min_score2",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  200              ,\
			    /*DOC*/       "Threshold for aligning a group in the slow double progressive alignment (automatically readjusted)",\
			    /*Parameter*/ &dpa_min_score2   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );
/*PARAMETER PROTOTYPE:    dpa_keep_tmp_file     */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-dpa_keep_tmpfile"        ,\
			    /*Flag*/      &dpa_keep_tmpfile        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Prevents deletion of the tmpfile generated by t_coffee_dpa",\
			    /*Parameter*/ &do_version          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \

			    );
/*PARAMETER PROTOTYPE:    dpa_debug     */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-dpa_debug"        ,\
			    /*Flag*/      &dpa_debug        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "DEbug mode for DPA ( causes dpa tmp files to be kept)",\
			    /*Parameter*/ &do_version          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \

			    );

/*PARAMETER PROTOTYPE:    multi_core    */
	       declare_name (multi_core);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-multi_core",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Multi core: template_jobs_relax_msa",\
			    /*Parameter*/ &multi_core   ,\
			    /*Def 1*/    "templates_jobs_relax_msa"       ,\
			    /*Def 2*/    "templates_jobs_relax_msa"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );
	       if (multi_core[0])set_string_variable ("multi_core",multi_core);
/*PARAMETER PROTOTYPE:    multi_core    */
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-n_core",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Number of cores to be used by machine [default=0 => all those defined in the environement]",\
			    /*Parameter*/ &n_core   ,\
			    /*Def 1*/    "0"       ,\
			    /*Def 2*/    "0"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );
	       if (n_core)set_int_variable ("n_core",n_core);


/*PARAMETER PROTOTYPE:    lib_list    */
	       declare_name (lib_list);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-lib_list",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "A File that contains every pair/group of sequence to process when computing the lib, Format:<nseq> <index1><index2>",\
			    /*Parameter*/ &lib_list   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    "default"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );

	       /*PARAMETER PROTOTYPE:    lib_list    */
	       declare_name (prune_lib_mode);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-prune_lib_mode",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "A File that contains every pair/group of sequence to process when computing the lib, Format:<nseq> <index1><index2>",\
			    /*Parameter*/ &prune_lib_mode   ,\
			    /*Def 1*/    "5"       ,\
			    /*Def 2*/    "5"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );
	       set_string_variable ("prune_lib_mode",prune_lib_mode);

	       /*PARAMETER PROTOTYPE:    multi_thread    */
	       declare_name (tip);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-tip",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Controls The Output of A TIP When Computation is over [one,all,none]",\
			    /*Parameter*/ &tip   ,\
			    /*Def 1*/    "one"       ,\
			    /*Def 2*/    "all"              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       /*PARAMETER PROTOTYPE:    RNA LIB    */
	       declare_name (rna_lib);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-rna_lib",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "",\
			    /*Parameter*/ &rna_lib   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-no_warning",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Suppresses all Warnings",\
			    /*Parameter*/ &no_warning   ,\
			    /*Def 1*/    "0"       ,\
			    /*Def 2*/    "1"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "1"           \
		   );
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,\
			    /*Name*/      "-run_local_script",\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "D"          ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,	\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "Run Local Script if in current directory",	\
			    /*Parameter*/ &run_local_script   ,		\
			    /*Def 1*/    "0"       ,			\
			    /*Def 2*/    "1"              ,		\
			    /*Min_value*/ "0"          ,		\
			    /*Max Value*/ "1"				\
		   );
	       set_int_variable ("run_local_script", run_local_script);
	       declare_name (plugins);
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,\
			    /*Name*/      "-plugins",\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "S"          ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,	\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "Set the directory containing the plugins [no if no plugin]",	\
			    /*Parameter*/ &plugins   ,		\
			    /*Def 1*/    "default"       ,			\
			    /*Def 2*/    ""              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
		   );
	       


	       declare_name (proxy);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-proxy",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "proxy used to access to webservices, when required",\
			    /*Parameter*/ &proxy   ,\
			    /*Def 1*/    "unset"       ,\
			    /*Def 2*/    " "              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if ( !strm (proxy, "unset"))set_string_variable ("cl_proxy",proxy);
	       declare_name (email);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-email",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "email provided to webservices, when required",\
			    /*Parameter*/ &email   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if ( strstr (email, "@"))
		 {
		   set_string_variable ("email", email);
		   set_string_variable ("cl_email", email);
		 }

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-clean_overaln",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &clean_overaln   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    "1"              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
					  );
	       overaln_param=declare_char ( 10, STRING);
	       n_overaln_param=get_cl_param(			\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_param",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  10              ,\
			    /*DOC*/       "Parameters for the overaln",\
			    /*Parameter*/ overaln_param   ,\
			    /*Def 1*/    "NULL"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       declare_name (overaln_mode);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_mode",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "lower || uanlaign",\
			    /*Parameter*/ &overaln_mode   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if (overaln_mode[0])set_string_variable ("overaln_mode", overaln_mode);
	       declare_name (overaln_model);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_model",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "fsa1 (no exon boundaries), fsa2 (exon boundaries)",\
			    /*Parameter*/ &overaln_model   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if (overaln_mode[0])set_string_variable ("overaln_model", overaln_model);

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_threshold",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &overaln_threshold   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
			    );
	       set_int_variable ("overaln_threshold", overaln_threshold);

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_target",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &overaln_target   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       set_int_variable ("overaln_target", overaln_threshold);

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_P1",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &overaln_P1   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if (overaln_P1)set_int_variable ("overaln_P1", overaln_P1);

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_P2",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &overaln_P2   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
		if (overaln_P2)set_int_variable ("overaln_P2", overaln_P2);

	        get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_P3",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &overaln_P3   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if (overaln_P3)set_int_variable ("overaln_P3", overaln_P3);

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-overaln_P4",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Ratio between overaligned exon id Vs legitimates *100",\
			    /*Parameter*/ &overaln_P4   ,\
			    /*Def 1*/    "0"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if (overaln_P4)set_int_variable ("overaln_P4", overaln_P4);


	       declare_name (exon_boundaries);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-exon_boundaries",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "exon_boundaries [EBI boj format]",\
			    /*Parameter*/ &exon_boundaries   ,\
			    /*Def 1*/    ""       ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
	       if ( exon_boundaries[0])set_string_variable ("exon_boundaries", exon_boundaries);





/*******************************************************************************************************/
/*                                                                                                     */
/*                           TCoffee_dpa Parameter:END                                                 */
/*                                                                                                     */
/*******************************************************************************************************/



	       if (argc==1  )
		 {
		   display_method_names ("display", stdout);
		   return EXIT_SUCCESS;
		 }
	       get_cl_param( argc, argv,&le, NULL,NULL,NULL,0,0,NULL);
	       prepare_cache (cache);
/*******************************************************************************************************/
/*                                                                                                     */
/*                           FILL list_file (contains seq, aln and meth)                               */
/*                                                                                                     */
/*******************************************************************************************************/



/*Re-introduce the sequences introduced with -infile*/
/*Standard*/

	       if ( infile[0] && !do_evaluate)
		   {
		   sprintf ( list_file[n_list++], "%s",infile);
		   }
/*DO EVALUATE: The aln to evaluate must be provided via -infile*/
	       else  if (do_evaluate)
	           {
		     if (!infile[0] ||  !(main_read_aln ( infile, NULL)))
		       {
			 fprintf ( stderr,"\nERROR: When using -evaluate, Provide a multiple sequence alignment via the -infile flag [FATAL:%s]\n", PROGRAM);
			 myexit (EXIT_FAILURE);
		       }
		     else if (! main_read_aln ( infile,NULL))
		       {
			 fprintf ( stderr,"\nERROR: FILE %s is NOT a valid alignment [FATAL:%s]\n", infile, PROGRAM);
			 myexit (EXIT_FAILURE);
		       }
		     else if ( infile[0]=='A' ||infile[0]=='S')
		       {
			 sprintf ( list_file[n_list++], "S%s",infile+1);
		       }
		     else sprintf ( list_file[n_list++], "S%s",infile);
		   }



/*Make Sure -infile is set*/
	       if (!infile[0]&& (do_evaluate || do_convert))
		 {

		   if ( do_evaluate || do_convert)sprintf ( infile, "%s",seq_list[0]);
		 }

/*EXPAND -in*/
	       /*Introduce the sequences from the -profile flag*/
	       if ( profile1 && profile1[0])
		 {
		   sprintf ( list_file[n_list++], "R%s",profile1);
		 }
	       if ( profile2 && profile2[0])
		 {
		   sprintf ( list_file[n_list++], "R%s",profile2);
		 }

	       for ( a=0; a< n_profile_list; a++)
		 {
		   FILE *fp;
		   if ( (fp=find_token_in_file (profile_list[a], NULL, "FILE_LIST"))!=NULL)
		     {
		       int z;
		       char rname[1000];
		       vfclose (fp);
		       fp=vfopen (profile_list[a], "r");

		       while ( (z=fgetc(fp))!=EOF)
			 {
			   ungetc(z, fp);
			   fscanf (fp, "%s\n", rname);
			   if ( check_file_exists(rname))sprintf ( list_file[n_list++], "R%s", rname);
			 }
		       vfclose (fp);
		     }
		   else if (format_is_conc_aln (profile_list[a]))
		     {
		       Alignment *P;
		       char *cname;

		       P=input_conc_aln (profile_list[a],NULL);
		       while (P)
			 {
			   cname=vtmpnam (NULL);
			   output_fasta_aln (cname, P);
			   P=P->A;
			   sprintf ( list_file[n_list++], "R%s",cname);
			 }
		       free_aln (P);
		     }

		   else
		     {
		       sprintf ( list_file[n_list++], "R%s",profile_list[a]);
		     }
		 }
	       /*Introduce the sequences from the -seq flag*/
	       for (a=0; a<n_seq_list; a++)
		 {
		   if (check_file_exists(seq_list[a]))
		     sprintf (list_file[n_list++], "S%s",seq_list[a]);
		   else if ( check_file_exists (seq_list[a]+1))
		     sprintf (list_file[n_list++], "%s",seq_list[a]);
		 }
	       /*introduce the alignments from the -aln flag*/
	       //Importnat: Must be introduced AFTER the profiles
	       for (a=0; a<n_aln_file_list; a++)
		 {
		   sprintf (list_file[n_list++], "A%s",aln_file_list[a]);
		 }
	       /*introduce the alignments from the -method flag*/
	       for (a=0; a<n_method_list; a++)
		 {
		   sprintf (list_file[n_list++], "M%s",method_list[a]);
		 }
	       /*introduce the alignments from the -library flag*/
	       for (a=0; a<n_lib_file_list; a++)
		 {
		   sprintf (list_file[n_list++], "L%s",lib_file_list[a]);
		 }
	       /*introduce sequences from the exon_boundaries flag flag*/
	       if ( exon_boundaries && exon_boundaries[0] && check_file_exists (exon_boundaries))
		 {
		   Sequence *ExS;
		   Alignment *ExA;
		   char *tmpf;
		   //make sure boundaries do not get into the sequences*/

		   ExS=main_read_seq (exon_boundaries);
		   ExS=seq2clean_seq (ExS, "BOJboj");
		   main_output_fasta_seq (tmpf=vtmpnam (NULL),ExA=seq2aln (ExS,NULL,RM_GAP), NO_HEADER);
		   sprintf (list_file[n_list++], "S%s",tmpf);
		   free_sequence (ExS, ExS->nseq);
		   free_aln (ExA);
		 }
	       /*FETCH THE STRUCTURES INTRODUCED WITH -pdb and add them to -in*/
	       if ( n_pdb)
		 {
		   for ( a=0; a< n_pdb; a++)
		     {
		       if ( is_number (pdb_list[a]));
		       else
			 {
			 pdb_start=pdb_end=0;
			 if ( a+1< n_pdb && is_number (pdb_list[a+1]))pdb_start=atoi (pdb_list[a+1]);
			 if ( a+2< n_pdb && is_number (pdb_list[a+2]))pdb_end=atoi (pdb_list[a+2]);

			 pdb_name=get_pdb_struc ( pdb_list[a],pdb_start, pdb_end);
			 if (pdb_name){sprintf (list_file[n_list++], "P%s", pdb_name);}
			 /*Warning: do not free pdb_name: it is statically allocated by get_pdb_struc*/
			 }
		     }
		 }

	       /*Check That Enough Methods/Libraries/Alignments Have been Chiped in*/

	       if (list_file)
		 {
		   int *nn;
		   nn=vcalloc ( 256, sizeof (int));
		   for (a=0; a<n_list; a++)
		     {
		       if ( !check_file_exists(list_file[a]))nn[(int)list_file[a][0]]++;
		       else
			 {
			   if (is_seq (list_file[a]))nn['S']++;
			   else if ( is_aln (list_file[a]))nn['A']++;
			   else if ( is_lib (list_file[a]))nn['L']++;
			   else if ( is_method (list_file[a]))nn['M']++;
			   else
			     add_warning (stderr, "\nWARNING: File %s was not properly tagged. Potential ambiguity\n",list_file[a]);
			 }
		     }


		   if ( (nn['A']+nn['L']+nn['M'])==0)
		     {
		       sprintf ( list_file[n_list++], "Mproba_pair"); //new default
		       //sprintf ( list_file[n_list++], "Mlalign_id_pair");
		       //sprintf ( list_file[n_list++], "Mslow_pair");
		     }
		   vfree (nn);
		 }

/*FILL THE F STRUCTURE (Contains Information for Output names For the defaults)*/
	       if (n_list==0 || argc<=1)
		 {
		   fprintf ( stderr, "\nERROR: You have NOT provided enough arguments [FATAL:%s]", PROGRAM);
		   myexit (EXIT_FAILURE);
		 }


	       else if ( argv[1][0]!='-' && (check_file_exists( argv[1]) || check_file_exists(argv[1]+1)))
		 {
		   if (check_file_exists(argv[1]))F=parse_fname(argv[1]);
		   else if ( check_file_exists(argv[1]+1))F=parse_fname(argv[1]+1);

		 }
	       else if (infile[0])
		 {

		   if ( check_file_exists (infile))F=parse_fname(infile);
		   else if (check_file_exists (infile+1))F =parse_fname(infile+1);
		 }
	       else if ( exon_boundaries && exon_boundaries[0])
		 {
		   if ( check_file_exists (exon_boundaries))F=parse_fname(exon_boundaries);
		   else if (check_file_exists (exon_boundaries+1))F =parse_fname(exon_boundaries+1);
		 }
	       else
	          {

		  for ( a=0; a< n_list; a++)
		      {
			if (!is_method(list_file[a]))
			  {


			     if ( check_file_exists( list_file[a])){F=parse_fname(list_file[a]);break;}
			     else if ( is_in_set ( list_file[a][0], "ASLX") && check_file_exists( list_file[a]+1)){F=parse_fname(list_file[a]+1);break;}
			     else if ( is_in_set ( list_file[a][0], "R") && check_file_exists( list_file[a]+1))
			       {
				 char lname[100];
				 F=parse_fname(list_file[a]+1);
				 sprintf ( lname, "%s_1", F->name);
				 sprintf ( F->name, "%s", lname);
				 break;
			       }

			     else if ( is_in_set ( list_file[a][0], "P") && is_pdb_struc (list_file[a]+1))
			       {
				 F=parse_fname(is_pdb_struc (list_file[a]+1));break;

			       }
			  }
		      }

		  }


	       /*Get Structures*/
	       for ( a=0; a< n_list; a++)
		 {
		   if ( list_file[a][0]=='P' && !check_file_exists(list_file[a]))
		     {
		       char buf[1000];
		       sprintf(buf, "%s", list_file[a]+1);
		       sprintf(list_file[a], "P%s",is_pdb_struc (buf));
		     }
		 }

	       /*FATAL: NO SEQUENCES*/
	       if (!F)
		 {
		   myexit (fprintf_error(stderr,"You have not provided any sequence"));
		 }
	       if (run_name)F=parse_fname(run_name);
	       else F->path[0]='\0';


	       identify_list_format      (list_file, n_list);


	       fprintf (le, "\nINPUT FILES\n");
	       for ( a=0; a< n_list; a++)
		   {
		     fprintf (le, "\tInput File (%c) %s ",list_file[a][0],list_file[a]+1);
		     if ( list_file[a][0]=='A' || list_file[a][0]=='S' || list_file[a][0]=='P'|| list_file[a][0]=='R' )
		       {
			 fprintf (le, " Format %s\n", f=identify_seq_format ( list_file[a]+1));

			 if (!f || f[0]=='\0')
			   {
			     myexit (fprintf_error(stderr,"The format of %s is not supported", list_file[a]+1));
				     
			   }
			 vfree (f);
		       }
		     else fprintf (le, "\n");
		   }


/*CONVERT, ALIGN OR EVALUATE: CHOSE THE RIGHT VERB*/
	       /*Set the Hierarchy of the verbs*/
	       /*The first one decides...*/


	       do_list=vcalloc ( 100, sizeof (int*));
	       n_do=0;
	       do_list[n_do++]=&do_genepred;
	       do_list[n_do++]=&do_extended_matrix;
	       do_list[n_do++]=&do_convert;
	       do_list[n_do++]=&do_evaluate;
	       do_list[n_do++]=&do_domain;
	       do_list[n_do++]=&do_align;


	       for ( a=0; a< n_do; a++)
		 {
		 if ( do_list[a][0])
		   {
		   for ( b=0; b< n_do; b++)if ( b!=a)do_list[b][0]=0;
		   break;
		   }
		 }



/*SET THE DEFAULT NAMES*/
	       if ( do_convert)
	           {
		     if ( strm (tree_file, "default"))sprintf ( tree_file, "no");
		   }



	       if (  do_evaluate)
		   {
		   sprintf ( out_lib, "no");
		   sprintf ( tree_file, "no");
		   clean_aln=0;
		   }
	       if (do_genepred)
		 {
		   sprintf ( tree_file, "no");
		   clean_aln=0;
		 }

	       if ( F && strm ( tree_file, "default"))sprintf ( tree_file ,"%s%s.dnd",F->path     ,F->name);
	       if ( F && strm ( ph_tree_file, "default"))sprintf ( ph_tree_file ,"%s%s.ph",F->path     ,F->name);

	       for (a=0; a< n_out_aln_format; a++)
		 {
		   if (is_out_format_list (out_aln_format[a]));
		   else
		     {
		       fprintf (stderr, "\n%s is not a valid format [FATAL:%s]\n", out_aln_format[a], PROGRAM);
		       myexit (EXIT_FAILURE);
		     }
		 }

	       for (a=0; a<n_out_aln_format; a++)
		 {
		   out_aln_format[a]=format_name2aln_format_name(out_aln_format[a]);
		 }

	       if ( F && strm ( out_aln  , "default"))
	          {
		  for (a=0; a< n_out_aln_format; a++)
		    {

		      sprintf ( tot_out_aln[a]   ,"%s%s.%s"      ,F->path,F->name,out_aln_format[a]);
		    }
		  }
	       else
	          {
		    sprintf ( tot_out_aln[0], "%s", out_aln);
		    for (a=1; a< n_out_aln_format; a++)
		      sprintf ( tot_out_aln[a]   ,"%s%s.%s", F->path  ,out_aln, out_aln_format[a]);
		  }



	       if ( F && strm ( out_lib  , "default"))sprintf ( out_lib   ,"%s%s.tc_lib",F->path     , F->name);

	       if ( type && type[0])
	          {
		      if (strm2 (type,"Protein", "protein"))sprintf ( type, "PROTEIN");
		      if (strm2 (type,"DNA", "dna"))sprintf ( type, "DNA");
		      if (strm2 (type,"RNA", "rna"))sprintf ( type, "RNA");

		  }


	       if (   !use_tree && check_file_exists (tree_file))vremove (tree_file);
	       else if ( !use_tree || (use_tree && strm (use_tree, "default")));
	       else sprintf ( tree_file, "%s", use_tree);
	       

/*******************************************************************************************************/
/*                                                                                                     */
/*                           Input Sequences and Library                                               */
/*                                                                                                     */
/*******************************************************************************************************/

	       set_methods_limits (method_limits,n_method_limits,list_file, n_list, &maxnseq, &maxlen);
	       /*Set Global Values*/
	     


/*START*/

	       /*1 READ THE SEQUENCES*/
	       
	       S=read_seq_in_n_list   (list_file, n_list, type,seq_source);
	       
	       if ( check_type)
		 {
		   if (!strm (S->type, get_array_type (S->nseq, S->seq)))
		     {
		       fprintf ( stderr, "\nINCORRECT SEQUENCE TYPE (USE %s ONLY) [FATAL:%s]", S->type, PROGRAM);
		       myexit (EXIT_FAILURE);
		     }
		 }
	       
	       if (S->nseq<=1 && !do_domain)
		 {
		   printf_exit (EXIT_FAILURE,stderr, "\nERROR: Your Dataset Contains %d Sequence. For multiple alignments you need at least 2 sequences[FATAL:%s]", S->nseq,PROGRAM);
		 }

	       store_seq_type (S->type);

	       if ( type_only==1)
		 {
		   fprintf ( stdout, "%s", S->type);
		   return EXIT_SUCCESS;
		 }
	       /*Translate Sequences*/
	       if ( transform && transform[0])
		 {
		   S=transform_sequence (S, transform);
		 }

	       /*Abort if the sequences are too long */
	       if (maxlen!=-1 && S->max_len>maxlen)
		  {
		      fprintf ( stderr, "\nSEQUENCES TOO LONG [Longuest=%d][MAX=%d][FATAL:%s]\n", S->max_len,maxlen, PROGRAM);
		      myexit (EXIT_FAILURE);

		  }

	      


	    
	       S=seq2template_seq(S, "SELF_S_",F);
	       /* Get the Templates*/
	       if ( n_template_file)
		 {
		   fprintf ( le, "\nLooking For Sequence Templates:\n");
		   for ( a=0; a< n_template_file; a++)
		     {
		       //correct for missing extension modes
		       if (strm (template_file_list[a],"RNA") && !strstr (extend_mode, "rna"))sprintf ( extend_mode, "rna2");


		       fprintf ( le, "\n\tTemplate Type: [%s] Mode Or File: [%s] [Start", template_type2type_name(template_file_list[a]), template_file_list[a]);
		       S=seq2template_seq(S, template_file_list[a], F);
		       fprintf ( le, "]");

		       if (S==NULL)
			 {
			   add_warning (stderr, "\nImpossible to find %s Templates\nCheck that your blast server is properly installed [See documentation][FATAL:%s]\n", template_file_list[a],PROGRAM);
			   myexit (EXIT_FAILURE);
			 }
		     }
		   
		   if (seq2n_X_template ( S, "_*_"))
		     {
		       
		       sprintf (S->template_file, "%s",seq2template_file (S, NULL));
		     }
		 }
	       else
		 {
		   int ptf=0;
		   for ( a=0; a<S->nseq; a++)
		     {
		       if ( seq_has_template ( S, a, "_P_"))ptf=1;
		     }
		   if (ptf)
		     {
		       int j;
		       sprintf ( S->template_file   ,"%s%s.template_file",F->path     , F->name);
		       seq2template_file (S,S->template_file);
		       display_output_filename ( stdout, "Template_List","fasta_seq", S->template_file, STORE);
		     }
		 }
	     

	       if (n_profile_template_file)
		 {
		   fprintf ( le, "\nLooking For Profile  Templates");
		   for ( a=0; a< n_profile_template_file; a++)
		     {
		       fprintf ( le, "\n\tTemplate Type: [%s] Mode Or File: [%s] [Start", template_type2type_name(profile_template_file_list[a]), profile_template_file_list[a]);
		       S=profile_seq2template_seq(S, profile_template_file_list[a], F);
		       fprintf ( le, "]");
		       if (S==NULL)
			 {
			   add_warning(stderr, "Impossible to find %s Templates\nCheck that your blast server is properly installed [See documentation][FATAL:%s]\n",profile_template_file_list[a], PROGRAM);
			   myexit (EXIT_FAILURE);
			 }
		     }
		 }
	       
	       S=seq2template_type (S);
	      
	       le=display_sequences_names   ( S, le, check_pdb_status, TEMPLATES);
	        
		 





	       if ( get_type)
		 {
		   S=get_sequence_type (S);
		   fprintf ( stdout , "%s\n", S->type);
		   free_sequence(S, S->nseq);
		   return 1;
		 }
	      

	       /*Reorder the sequences*/
	       new_order=duplicate_char (S->name, -1, -1);
	       if ( strm (inorder, "aligned"))new_order=sort_string_array   (new_order, S->nseq);

	       initial_order=duplicate_char (S->name, -1, -1);
	       S=reorder_seq(S,new_order,S->nseq);
	       free_char (new_order, -1);

	       

	       /*3 PREPARE THE CONSTRAINT LIST*/

	       CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?tmpfile():NULL, NULL);
	      
	       sprintf ( CL->method_evaluate_mode, "%s", method_evaluate_mode);
	     
	       (CL->TC)->use_seqan=use_seqan;
	       CL->local_stderr=le;
	       /*set the genepred parameters*/
	       sprintf ( CL->genepred_score, "%s", genepred_score);
	       /*Estimate the distance Matrix*/
	       CL->DM=cl2distance_matrix ( CL,NOALN,distance_matrix_mode, distance_matrix_sim_mode,1);

	       /*one to all alignment*/
	       if (one2all && one2all[0])prepare_one2all (one2all,S, lib_list);
	       else if ( subset2all)
		 {
		   prepare_subset2all (subset2all,S, lib_list,CL);
		 }

	       if (  matrix && matrix[0])
		 {
		   sprintf ( CL->method_matrix,"%s", matrix);

		 }
	       /*Set the filtering*/
	       CL->filter_lib=filter_lib;
	       /*Set the evaluation Functions*/
	       CL->profile_mode=get_profile_mode_function (profile_mode, NULL);
	       sprintf ( CL->profile_comparison, "%s", profile_comparison);
	       if ( n_struc_to_use)
		 {
		   CL->STRUC_LIST=declare_sequence (1,1,n_struc_to_use);
		   CL->STRUC_LIST->nseq=0;
		   for ( a=0; a< n_struc_to_use; a++)
		     {

		       sprintf ( (CL->STRUC_LIST)->name[(CL->STRUC_LIST)->nseq++],"%s",struc_to_use[a]);
		     }
		 }
	       sprintf (CL->align_pdb_param_file, "%s", align_pdb_param_file);
	       sprintf (CL->align_pdb_hasch_mode, "%s", align_pdb_hasch_mode);



	       /*Blast Parameters*/
	       (CL->Prot_Blast)->min_id=prot_min_sim;
	       (CL->Prot_Blast)->max_id=prot_max_sim;
	       (CL->Prot_Blast)->min_cov=prot_min_cov;
	       sprintf ( (CL->Prot_Blast)->blast_server, "%s", prot_blast_server);
	       sprintf ( (CL->Prot_Blast)->db, "%s", prot_db);

	       (CL->Pdb_Blast)->min_id=pdb_min_sim;

	       (CL->Pdb_Blast)->max_id=pdb_max_sim;
	       (CL->Pdb_Blast)->min_cov=pdb_min_cov;
	       sprintf ( (CL->Pdb_Blast)->blast_server, "%s", pdb_blast_server);
	       sprintf ( (CL->Pdb_Blast)->db, "%s", pdb_db);
	       CL->check_pdb_status=check_pdb_status;
	       /*split parameters */
	       CL->split=split;
	       CL->split_nseq_thres=split_nseq_thres;
	       CL->split_score_thres=split_score_thres;
	       /*Blast Parameters
	       (CL->DNA_Blast)->min_id=dna_min_sim;
	       (CL->DNA_Blast)->max_id=dna_max_sim;
	       (CL->DNA_Blast)->min_cov=dna_min_cov;
	       sprintf ( (CL->DNA_Blast)->blast_server, "%s", dna_blast_server);
	       sprintf ( (CL->DNA_Blast)->db, "%s", dna_db);
	       */

	       if ( method_log)
		 {
		   if ( strm (method_log, "default"))
		     {
		       sprintf ( CL->method_log, "%s%s.method_log",F->path, F->name);
		     }
		   else if ( !strm (method_log, "no"))
		     {
		      sprintf ( CL->method_log, "%s", method_log);
		     }
		   set_string_variable ("method_log", method_log);
		 }


	       CL->lalign_n_top=lalign_n_top;
	       sprintf ( CL->multi_thread, "%s", multi_core);
	       
	       sprintf ( CL->lib_list, "%s", lib_list);
	       sprintf (CL->rna_lib, "%s", rna_lib);
/* Important: This is where the library is compiled!!!!*/

	       if ((CL->S)->nseq>1 && !do_convert)
		 {
		   CL=read_n_constraint_list (list_file,n_list,NULL, mem_mode,weight,type, le, CL, seq_source);
		   //CL=post_process_constraint_list (CL); //needed when constraints are added, for instance the RNA modes
		 }
	       else if ( do_convert && out_lib[0])
		 {
		   if ( infile[0])
		     {sprintf (list_file[0], "%s", name2type_name(infile));
		     CL=read_n_constraint_list (list_file,1,NULL, mem_mode,weight,type, le, CL, seq_source);
		     }
		   else
		     {
		     CL=read_n_constraint_list (list_file,n_list,NULL, mem_mode,weight,type, le, CL, seq_source);
		     }
		 }
	       if ( CL->M)clean_aln=0;

	       if ( is_number (weight))set_weight4constraint_list (CL, atoi(weight));
	       
	       free_pair_wise ();//Free ststic memory allocated in some of the pairwise functions


	       //Shrink: re-run slow_pair using the library, remove everything


	       /*If the List is empty*/
	       if ( (CL->S)->nseq>1 && CL->ne==0 && !CL->M &&!(do_convert && n_list>0))
		 {
		   fprintf ( stderr, "\n******************ERROR*****************************************\n");

		   fprintf ( stderr, "\nYou have not provided any method or enough Sequences[FATAL]");
		   fprintf ( stderr, "\nIf you have used the '-in' Flag, ADD the methods you wish to use:");
		   fprintf ( stderr, "\n\t-in <your sequences> Mlalign_id_pair Mfast_pair\n");
		   fprintf ( stderr, "\nAnd make sure you provide at least TWO sequences\n");
		   for ( a=0; a< argc; a++)fprintf ( stderr, "%s ", argv[a]);
		   fprintf ( stderr, "\n*****************************************************************\n");
		   myexit(EXIT_FAILURE);
		 }


	       CL->normalise=do_normalise;

	       if ( type && type[0])sprintf ( (CL->S)->type, "%s", type);
	       CL->extend_jit=(do_extend>0)?1:0;

	       CL->extend_threshold=(do_extend==1)?0:do_extend;
	       CL->do_self=do_self;
	       sprintf (CL->extend_clean_mode,   "%s", clean_mode);
	       sprintf (CL->extend_compact_mode, "%s", compact_mode);
	       if ( CL->extend_jit && CL->extend_threshold !=0)filter_constraint_list (CL,WE,CL->extend_threshold);
	       CL->pw_parameters_set=1;
	       
	       

	       CL->nomatch=nomatch;
	       set_int_variable ("nomatch", nomatch);
	       /*Gep and Gop*/
	       if ( !gep &&  !gop && CL->M)
		   {
		    CL->gop=get_avg_matrix_mm ( CL->M, (strm3((CL->S)->type,"PROTEIN", "Protein", "protein")?AA_ALPHABET:"gcuta"))*10;
		    CL->gep=CL->gop/10;
		    fprintf ( CL->local_stderr, "\nAUTOMATIC PENALTIES: gapopen=%d gapext=%d", CL->gop, CL->gep);
		   }
	       else if ( !CL->M && cosmetic_penalty && !gep && !gop)
	           {
		    CL->gep=0;
		    CL->gop=cosmetic_penalty;
		   }
	       else
	           {
		     CL->gep=gep;
		     CL->gop=gop;
		     fprintf ( CL->local_stderr, "\nMANUAL PENALTIES: gapopen=%d gapext=%d", CL->gop, CL->gep);
		   }

	       /*Frame Penalties*/
	       CL->f_gep=f_gep;
	       CL->f_gop=f_gop;


	       CL->maximise=maximise;

	       if (strm(retrieve_seq_type(),"DNA")|| strm(retrieve_seq_type(),"RNA") )
		   CL->ktup=MAX(2,ktup);
	       else
		   CL->ktup=ktup;

	       CL->use_fragments=diag_mode;
	       CL->fasta_step=fasta_step;
	       CL->diagonal_threshold=diag_threshold;

	       sprintf ( CL->matrix_for_aa_group, "%s", sim_matrix);
	       sprintf ( CL->dp_mode, "%s", dp_mode);
	       CL->TG_MODE=tg_mode;

	       sprintf ( CL->evaluate_mode, "%s", evaluate_mode);
	       fprintf (le, "\n\n\tLibrary Total Size: [%d]\n", CL->ne);


	       CL=choose_extension_mode (extend_mode, CL);
	       CL->max_n_pair=max_n_pair;

	       processed_lib=0;
	       //use the L, vener touch it again
	       

	      if (CL->ne>0 && out_lib[0]!='\0' && !strm (out_lib, "no"))
	         {

		   if (strstr (out_lib_mode, "extended"))
		     {
		       char emode[1000];

		       //Do the processing before saving the extended lib*/
		       processed_lib=1;
		       if ( filter_lib) CL=filter_constraint_list (CL,CL->weight_field, filter_lib);
		       for (a=0; a<relax_lib; a++)CL=relax_constraint_list (CL);
		       for (a=0; a<shrink_lib; a++)CL=shrink_constraint_list (CL);
		       sprintf ( emode, "lib_%s", out_lib_mode);

		       OUT=vfopen (out_lib, "w");
		       OUT=save_extended_constraint_list(CL,emode,OUT);
		     }
		   else
		     {
		       OUT=save_constraint_list ( CL, 0, CL->ne, out_lib, NULL, "ascii",CL->S);
		     }
		   vfclose (OUT);
		   CL->local_stderr=display_output_filename (le, "TCLIB","tc_lib_format_01",out_lib, CHECK);
		 }

	      
	      if ( lib_only)return EXIT_SUCCESS;


	      if (!processed_lib)
		{
		   if ( filter_lib) CL=filter_constraint_list (CL,CL->weight_field, filter_lib);
		   if (atoigetenv ("EXTEND4TC")==1)CL=extend_constraint_list(CL);
		   for (a=0; a<relax_lib; a++)CL=relax_constraint_list (CL);
		   for (a=0; a<shrink_lib; a++)CL=shrink_constraint_list (CL);
		}

	      CL=evaluate_constraint_list_reference (CL);
	      sprintf ( CL->distance_matrix_mode, "%s", distance_matrix_mode);
	      sprintf ( CL->distance_matrix_sim_mode, "%s", distance_matrix_sim_mode);
	      
	      sprintf ( CL->tree_mode, "%s", tree_mode);
	      //Re-estimate the distance matrix with consistency//
	       if ( strm ("cscore", distance_matrix_mode))
		 {
		   CL->DM=cl2distance_matrix ( CL,NOALN,distance_matrix_mode, distance_matrix_sim_mode,1);
		 }
	      /*WEIGHT CONSTRAINT LIST*/

	      if ( !do_convert)
		{

		  CL->DM=cl2distance_matrix (CL, NOALN, NULL, NULL,0);

		  CL=weight_constraint_list(CL, seq_weight);
		  
		  if (output_seq_weights (CL->W, outseqweight))
		    CL->local_stderr=display_output_filename( CL->local_stderr,"WEIGHT","tc_weight",outseqweight, CHECK);
		  le=display_weights(CL->W, le);
		}



	      /*Prepare quadruplets*/
	      if ( nseq_for_quadruplet && !strm(seq_name_for_quadruplet[0], "all"))
		{
		  CL->nseq_for_quadruplet=nseq_for_quadruplet;
		  CL->seq_for_quadruplet=vcalloc ((CL->S)->nseq, sizeof (int));
		  for (a=0; a< CL->nseq_for_quadruplet; a++)
		    {
		      printf ( "\nquad: %s", seq_name_for_quadruplet[a]);
		      if ( (b=name_is_in_list (seq_name_for_quadruplet[a],(CL->S)->name,(CL->S)->nseq, 100))!=-1)CL->seq_for_quadruplet[b]=1;
		      else add_warning ( stderr, "\nWARNING: Sequence %s is not in the set and cannot be used for quadruplet extension\n",seq_name_for_quadruplet[a]);
		    }
		}
	      else if ( nseq_for_quadruplet && strm(seq_name_for_quadruplet[0], "all"))
		{

		  CL->nseq_for_quadruplet=(CL->S)->nseq;
		  CL->seq_for_quadruplet=vcalloc ((CL->S)->nseq, sizeof (int));
		  for (a=0; a< CL->nseq_for_quadruplet; a++)
		    {
		      CL->seq_for_quadruplet[a]=1;
		    }
		}

/*******************************************************************************************************/
/*                                                                                                     */
/*                           Prepare The Alignment                                                     */
/*                                                                                                     */
/*******************************************************************************************************/

	      if ( do_align )
		   {


		   A=seq2aln  ((CL->S),NULL,1);
		   ungap_array(A->seq_al,A->nseq);

		   /*Chose the right Mode for evaluating Columns*/

		   if ( A->nseq==1);
		   else if ( strm ( msa_mode, "seq_aln"))
		     {
		       A=seq_aln (A,(CL->S)->nseq, CL);
		     }
		   else if ( strm ( msa_mode, "sorted_aln"))
		     {
		       A=sorted_aln (A, CL);
		     }
		   else if ( strm ( msa_mode, "full_sorted_aln"))
		     {
		       full_sorted_aln (A, CL);
		       output_constraints (out_lib, "sim", A);
		       CL->local_stderr=display_output_filename (le, "TCLIB","tc_lib_format_01",out_lib, CHECK);
		       return EXIT_SUCCESS;
		     }

		   else if ( strm ( msa_mode, "profile_aln"))
		     {
		       A=iterative_tree_aln (A, 0, CL);
		       A=profile_aln (A, CL);
		     }
		   else if ( strm ( msa_mode, "iterative_aln"))
		     {
		       A=iterative_tree_aln (A, 0, CL);
		       A=iterative_aln (A,10, CL);
		     }
		    else if ( strm ( msa_mode, "iterative_tree_aln"))
		     {
		       A=iterative_tree_aln (A,1, CL);
		     }
		    else if ( strm ( msa_mode, "dpa_aln"))
		      {
			A=dpa_aln (A, CL);
		      }
		    else if ( strm ( msa_mode, "new_dpa_aln"))
		      {
			A=new_dpa_aln (A, CL);
		      }
		    else if ( strm ( msa_mode, "delayed_tree_aln"))
		     {
		       A=make_delayed_tree_aln (A,2, CL);
		     }
		   else if ( strm ( msa_mode, "groups"))
		     {
		       A=seq2aln_group (A,dpa_maxnseq, CL);
		       out_aln_format[0]="conc_aln";
		       n_out_aln_format=1;
		     }
		   else if ( strm ( msa_mode, "upgma"))
		     {
		       A=upgma_tree_aln (A, A->nseq, CL);
		     }
		   else if ( strm ( msa_mode, "graph"))
			{
			  fprintf ( stderr, "\nDO GRAPH ALIGNMENT");
			  A=graph_aln ( A, CL, (CL->S));
			}
		   else if ( strm ( msa_mode, "tsp"))
			{
			  fprintf ( stderr, "\nDO TSP ALIGNMENT");
			  A=tsp_aln ( A, CL, (CL->S));
			}
		   else if ( strm ( msa_mode, "precomputed"))
		     {
		       if (infile[0]) {free_aln (A);A=main_read_aln ( infile, declare_aln(CL->S));}
		       else{fprintf ( stderr, "\nERROR: distance_matrix_mode=aln requires an aln passed via the -infile flag [FATAL:%s]", PROGRAM);crash ("");}

		       sprintf ( CL->dp_mode, "precomputed_pair_wise");
		       sprintf ( CL->distance_matrix_mode, "aln");
		       CL->tree_aln=A=reorder_aln ( A, (CL->S)->name,(CL->S)->nseq);

		       pc=tree_file;
		       if ( strm (tree_file, "default") || !check_file_exists (tree_file))
			 T=make_tree ( A,CL,gop, gep,(CL->S),pc, maximise);
		       else if ( strm (tree_file, "no"))
			 T=make_tree ( A,CL,gop, gep,(CL->S),NULL, maximise);
		       else
			 {
			   T=read_tree (pc,&tot_node,(CL->S)->nseq,  (CL->S)->name);
			 }

		       SNL=tree_aln ((T[3][0])->left,(T[3][0])->right,A,(CL->S)->nseq, CL);
		     }
		   else if ( strm ( msa_mode, "tree"))
		     {
		       if ( strm (CL->distance_matrix_mode, "aln"))
			 {
			   if (infile[0]) {free_aln (A);A=main_read_aln ( infile, declare_aln(CL->S));}
			   else{fprintf ( stderr, "\nERROR: distance_matrix_mode=aln requires an aln passed via the -infile flag [FATAL:%s]", PROGRAM);crash ("");}
			   CL->tree_aln=A;
			 }
		       pc=tree_file;
		       if ( strm (tree_file, "default") || !check_file_exists (tree_file))
			 T=make_tree ( A,CL,gop, gep,(CL->S),pc,maximise);
		       else if ( strm (tree_file, "no"))
			 T=make_tree ( A,CL,gop, gep,(CL->S),NULL, maximise);
		       else
			 {
			   fprintf ( le, "\nREAD PRECOMPUTED TREE: %s\n", pc);
			   T=read_tree (pc,&tot_node,(CL->S)->nseq,  (CL->S)->name);
			 }
		       SNL=tree_aln ((T[3][0])->left,(T[3][0])->right,A,(CL->S)->nseq, CL);
		       A->nseq=(CL->S)->nseq;
		     }

		   else
		     {
		       fprintf ( stderr, "\nERROR: msa_mode %s is unknown [%s:FATAL]\n", msa_mode, PROGRAM);
		       crash ("");
		     }

		   }
	       else if ( (do_evaluate || do_convert))
	           {


		   A=(infile[0])?main_read_aln ( infile, declare_aln(CL->S)):NULL;

		   if (!A)A=seq2aln((CL->S), NULL,0);


		   A->S=CL->S;
		   A->nseq=(CL->S)->nseq;


		   }
	       else if ( do_genepred)
		 {
		   Alignment *A1, *A2;
		   A1=seq2aln(CL->S, NULL, 1);
		   A1->S=CL->S;
		   A2=coffee_seq_evaluate_output (A1, CL);
		   if (!A2->score_res)myexit(0);
		   for ( b=0; b< n_out_aln_format; b++)
		     {
		       Alignment *OUT;

		       OUT=copy_aln (A2,NULL);
		       output_format_aln (out_aln_format[b],OUT,NULL, tot_out_aln[b]);
		       le=display_output_filename( le,"MSA",out_aln_format[b], tot_out_aln[b], CHECK);
		       free_aln (OUT);
		     }
		   return EXIT_SUCCESS;
		 }

	      else if (do_domain)
		   {
		     CL->moca=vcalloc ( 1, sizeof ( Moca));
		     if (strm ( "cfasta_pair_wise", dp_mode))sprintf (CL->dp_mode, "%s","domain_pair_wise");
		     (CL->moca)->moca_start=domain_start;
		     (CL->moca)->moca_len  =domain_len;
		     (CL->moca)->moca_scale=(domain_scale==0)?-(CL->normalise/20):domain_scale;
		     (CL->moca)->moca_interactive=domain_interactive;



		     if (!cosmetic_penalty && !gep && !gop)
		       {
			 CL->gop=-200;
			 CL->gep=-100;
		       }

		     CL=prepare_cl_for_moca (CL);
		     aln_list=moca_aln (CL);
		     free_int ( CL->packed_seq_lu, -1);
		     CL->packed_seq_lu=NULL;

		     a=0;
		     while ( aln_list[a])
		       {
		       for ( b=0; b< n_out_aln_format; b++)
			  {

			    output_format_aln (out_aln_format[b],aln_list[a],EA=fast_coffee_evaluate_output(aln_list[a], CL), tot_out_aln[b]);
			    le=display_output_filename( le,"MSA",out_aln_format[b], tot_out_aln[b], CHECK);
			  }
		       a++;
		       }
		   return EXIT_SUCCESS;
		   }
	      else if ( do_extended_matrix)
		{
		  A=seq2aln(CL->S, NULL, 1);
		  A->CL=CL;
		  for ( a=0; a< n_out_aln_format; a++)
		    {
		      output_format_aln (out_aln_format[a],A,EA, tot_out_aln[a]);
		      le=display_output_filename( le,"MSA",out_aln_format[a], tot_out_aln[a], CHECK);
		    }

		  return EXIT_SUCCESS;
		}


/*******************************************************************************************************/
/*                                                                                                     */
/*                           PREPARE THE ALIGNMENT FOR OUTPUT                                          */
/*                                                                                                     */
/*******************************************************************************************************/

	      if (A)
	          {
		    /*
		    for ( a=0; a< A->nseq; a++)
		      {
			   for ( b=0; b< A->len_aln ; b++)
			     if ( A->seq_al[a][b]=='O' || A->seq_al[a][b]=='o')A->seq_al[a][b]='-';
			 }
		    */




		    if ( check_file_exists(outorder))
			{
			  Sequence *OS;
			  OS=get_fasta_sequence (outorder, NULL);
			  if ( prf_in_seq (CL->S))A->expanded_order=OS->name;
			  else A=reorder_aln ( A,OS->name,A->nseq);
			}
		      else if ( strm(outorder, "aligned") && T)
			{
			  A=reorder_aln ( A,A->tree_order,A->nseq);

			}
		      else
			{

			  A=reorder_aln ( A, (CL->S)->name,(CL->S)->nseq);
			  A=reorder_aln ( A, initial_order,(CL->S)->nseq);

			}

		      A->output_res_num=strm3 ( output_res_num, "on", "On", "ON");

		      if ( strm2 (residue_case, "keep", "retain"))A->residue_case=KEEP_CASE;
		      else if (strm3 (residue_case, "upper", "Upper", "UPPER"))A->residue_case=UPPER_CASE;
		      else if (strm3 (residue_case, "lower", "Lower", "LOWER"))A->residue_case=LOWER_CASE;
		      else A->residue_case=1;






		      if ( iterate)
			{
			  A=iterate_aln (A, iterate, CL);
			  A=ungap_aln(A);
			}

		      if ( clean_aln)
			{
			  EA=main_coffee_evaluate_output(A, CL,clean_evaluate_mode);
			  A=clean_maln(A, EA,clean_threshold,clean_iteration);
			  free_aln (EA);
			  A=ungap_aln(A);
			}

		      //overaln
		      if (clean_overaln)
			{
			  char *over_aln_tmp;
			  over_aln_tmp=vtmpnam(NULL);
			  output_format_aln ("overaln", A, NULL, over_aln_tmp);
			  A=main_read_aln (over_aln_tmp,A);
			}

		      EA=main_coffee_evaluate_output(A, CL, evaluate_mode);

		      //correct ascii file
		      if (clean_overaln)
			{
			  EA=overlay_alignment_evaluation (A,EA);
			}


		      if (A->A)A=A->A;
		      if (!strm2(out_aln, "stdout", "stderr") && le==stderr && !do_convert)output_format_aln ("aln",A,NULL,"stdout");


		      A->CL=CL;
		      for ( a=0; a< n_out_aln_format; a++)
			if ( !strstr ( out_aln_format[a], "expand"))output_format_aln (out_aln_format[a],A,EA, tot_out_aln[a]);
		      for ( a=0; a< n_out_aln_format; a++)
			if ( strstr (out_aln_format[a], "expand"))output_format_aln (out_aln_format[a],A,EA, tot_out_aln[a]);



		      fprintf (le, "\n\nOUTPUT RESULTS");
		      le=display_output_filename (le, "GUIDE_TREE","newick", tree_file, CHECK);

		      for ( a=0; a< n_out_aln_format; a++)
			le=display_output_filename( le,"MSA",out_aln_format[a], tot_out_aln[a], CHECK);
		      if (CL->ne>0 && out_lib[0]!='\0' && !strm (out_lib, "no"))
			CL->local_stderr=display_output_filename (le, "TCLIB","tc_lib_format_01",out_lib, CHECK);

		      if (!strm (ph_tree_file, "NO") && A->nseq>2)
			{
			  NT_node T;
			  FILE *tfp;
			  char **tmode;
			  tmode=declare_char (2, 100);
			  sprintf (tmode[0], "nj");
			  T=tree_compute (A, 1, tmode);
			  tfp=vfopen (ph_tree_file, "w");
			  tfp=print_tree (T, "newick", tfp);
			  vfclose (tfp);
			  le=display_output_filename (le, "PHYLOGENIC_TREE","newick", ph_tree_file, CHECK);
			}

		  }

	      if (split)
		{

		  if (trim && n_seq_to_keep)
		    {
		      if (n_seq_to_keep==1 && check_file_exists (seq_to_keep[0]))
			{

			  SEQ_TO_KEEP=read_sequences (seq_to_keep[0]);
			}
		      else
			{

			  SEQ_TO_KEEP=declare_sequence ( 1, 1,n_seq_to_keep);
			  for ( a=0; a< n_seq_to_keep; a++)sprintf ( SEQ_TO_KEEP->name[a], "%s", seq_to_keep[a]);
			}
		    }

		  sprintf ( CL->dp_mode, "precomputed_pair_wise");
		  sprintf ( CL->distance_matrix_mode, "aln");



		  CL->tree_aln=A=reorder_aln ( A, (CL->S)->name,(CL->S)->nseq);
		  CL->S=aln2seq ( A);

		  if (!T)
		    {

		      pc=tree_file;
		      if ( strm (tree_file, "default") || !check_file_exists (tree_file))
			T=make_tree ( A,CL,gop, gep,(CL->S),pc, maximise);
		      else if ( strm (tree_file, "no"))
			T=make_tree ( A,CL,gop, gep,(CL->S),NULL, maximise);
		      else
			{
			  T=read_tree (pc,&tot_node,(CL->S)->nseq,  (CL->S)->name);
			}
		    }

		  SNL=tree_aln ((T[3][0])->left,(T[3][0])->right,A,(CL->S)->nseq, CL);


		  for ( a=0, b=0; a<A->nseq; a++)b+=(SNL[a])?1:0;
		  fprintf ( le, "\n\nSPLIT DATASET: %d Groups\n", b);
		  /*Display Group Names*/

		  if ( trim && SEQ_TO_KEEP)
		    {
		      for ( a=0; a< SEQ_TO_KEEP->nseq; a++)
			{

			  trim_subS=extract_one_seq(SEQ_TO_KEEP->name[a],0,0,A,KEEP_NAME);
			  trim_S=add_sequence (trim_subS,trim_S,0);
			}
		    }
		  for ( a=0, b=0; a<A->nseq; a++)
		    {

		      if ( SNL[a])
			{
			  b++;
			  fprintf ( le, "\n\tSPLIT_GROUP %d ; Nseq %d ; Score %d ; List ",b, (SNL[a])->nseq, (int)(SNL[a])->score);
			  for ( c=0; c< (SNL[a])->nseq; c++)
			    {
			      fprintf ( le, "%s ",(CL->S)->name[(SNL[a])->lseq[c]]);
			    }

			  SPLIT_ALN=extract_sub_aln (A, (SNL[a])->nseq,(SNL[a])->lseq);
			  SPLIT_ALN->S=A->S;
			  ungap_aln (SPLIT_ALN);

			  if (!trim)
			    {
			      sprintf ( split_format, "%s", "clustalw");
			      sprintf ( split_name, "%s.split.%d.%s", F->name, b,split_format);
			      fprintf ( le, " ; File %s",  split_name);
			      output_format_aln (split_format,SPLIT_ALN,NULL,split_name);
			      le=display_output_filename( le,"SPLIT_SEQ",split_format,split_name, CHECK);
			    }
			  else if (trim)
			    {
			      t=aln2most_similar_sequence(SPLIT_ALN, "idmat");
			      trim_subS=extract_one_seq(SPLIT_ALN->name[t],0,0,SPLIT_ALN,KEEP_NAME);
			      trim_S=add_sequence (trim_subS,trim_S,0);
			      fprintf ( le, "\n\tTRIM_SEQ: Kept sequence %s",SPLIT_ALN->name[t]);
			    }
			  free_aln (SPLIT_ALN);
			  fprintf (le, "\n");
			}
		    }

		  if (trim)
		      {


		       SPLIT_ALN=seq2aln (trim_S,NULL, KEEP_GAP);
		       ungap_aln (SPLIT_ALN);
		       sprintf ( trim_format, "%s", "fasta_aln");
		       if ( strm (trimfile, "default"))sprintf ( trimfile, "%s.trim.%s", F->name,trim_format);

		       output_format_aln (trim_format,SPLIT_ALN,NULL,trimfile);
		       le=display_output_filename( le,"TRIM_SEQ",trim_format,trimfile, CHECK);
		      }
		}

	      if (remove_template_file){S=vremove_seq_template_files(S);}
	      else
		{
		  S=display_seq_template_files (S);
		}

	      //fLUSH OUT THE NAME OF ALL THE FILES THAT HAVE BEEN PRODUCED
	      le=display_output_filename (le, NULL, NULL, NULL, FLUSH);


	      fprintf (le, "\n\n");

	      free_char (list_file, -1);
	      free_Alignment (A);
	      free_Alignment (EA);


	      S=free_constraint_list (CL);
	      free_sequence (S, S->nseq);



	      vremove ( "core");

	      vfree_all();

	      le=t_coffee_tip (le, tip);
	      le=print_command_line ( le);
	      le=print_mem_usage (le, PROGRAM);
	      le=print_cpu_usage(le, PROGRAM);
	      le=print_program_information (le, NULL);


	      if (full_log && full_log[0])log_function(full_log);

	      return EXIT_SUCCESS;
	}

/*Specialized set of Parameters*/
char *get_defaults (char *buf, char *type)
{
  return NULL;
}
char *get_precomputed_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf," -msa_mode=precomputed ");
     buf=strcat (buf," -seq_weight=no ");
     buf=strcat (buf," -evaluate_mode no ");
     buf=strcat (buf," -in Xpam250mt ");
     return buf;
   }
char *get_evaluate_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf," -quiet=stdout ");
     /*buf=strcat (buf," -seq_weight=no ");*/
     buf=strcat (buf," -output score_ascii html ");
     buf=strcat (buf," -iterate 0 ");

     buf=strcat (buf," -evaluate ");



     return buf;
   }
char *get_genome_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf," -seq_weight=no ");
     buf=strcat (buf," -dp_mode sim_pair_wise_lalign ");
     buf=strcat (buf," -output glalign ");
     buf=strcat (buf," -iterate 0 ");
     buf=strcat (buf," -distance_matrix_mode ktup ");
     buf=strcat (buf," -evaluate_mode t_coffee_slow ");
     buf=strcat (buf," -gapopen 100 -gapext 20 -nomatch 30 ");
     buf=strcat (buf," -clean_aln 0 ");
     buf=strcat (buf,"-output clustalw,score_ascii ");


     return buf;
   }
char *get_dali_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-cosmetic_penalty=-50 ");
     buf=strcat (buf,"-distance_matrix_mode=slow ");
     buf=strcat (buf,"-output clustalw,score_ascii ");
     buf=strcat (buf,"-evaluate_mode=non_extended_t_coffee ");
     buf=strcat (buf,"-clean_aln 0 ");

     return buf;
   }

char *get_very_fast_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));
     
     buf=strcat (buf,"-in Xblosum62mt ");
     buf=strcat (buf,"-distance_matrix_mode ktup ");
     buf=strcat (buf,"-maxnseq 10000 ");
     buf=strcat (buf,"-dpa_maxnseq 0 ");
     buf=strcat (buf,"-dp_mode fasta_pair_wise ");
     buf=strcat (buf,"-extend_mode matrix ");
     buf=strcat (buf,"-gapopen -10 ");
     buf=strcat (buf,"-gapext -1 ");
     buf=strcat (buf,"-iterate 0 ");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }

char *get_low_memory_defaults(char *buf, char *type)
   {
     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     if (NO_METHODS_IN_CL)buf=strcat (buf,"-distance_matrix_mode=idscore -method lalign_id_pair slow_pair -dp_mode=linked_pair_wise ");
     else buf=strcat (buf,"-distance_matrix_mode=idscore -dp_mode=linked_pair_wise ");
     return buf;
   }
char *get_dna_defaults(char *buf, char *type)
{

  return buf;
}
char *get_cdna_defaults(char *buf, char *type)
{
  buf=strcat (buf,"-distance_matrix_mode=idscore -dp_mode=fasta_cdna_pair_wise ");
  return buf;
}
char *get_3dcoffee_defaults(char *buf, char *type)
   {
     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Msap_pair  -template_file SELF_P_ -profile_template_file SELF_P_");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_expresso_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Msap_pair -template_file EXPRESSO -profile_template_file EXPRESSO");

     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_psicoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Mproba_pair -template_file BLAST   ");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_accurate_defaults ( char *buf, char *type)
{

  if ( strm (type, "PROTEIN")) return get_accurate4PROTEIN_defaults(buf, type);
  else if ( strm (type, "DNA")) return get_accurate4DNA_defaults(buf, type);
  else if ( strm (type, "RNA")) return get_accurate4RNA_defaults(buf, type);
  else return get_defaults(buf, type);
}
char *get_accurate4PROTEIN_defaults(char *buf, char *type)
   {
     if (buf==NULL)buf=vcalloc (1000, sizeof (char));
     if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mbest_pair4prot -template_file BLAST -template_file EXPRESSO  ");
     else buf=strcat (buf,"-template_file BLAST -template_file EXPRESSO  ");
     buf=strcat (buf,"-output aln, expanded_fasta_aln ");

     return buf;
   }




char *get_accurate4DNA_defaults(char *buf, char *type)
{
  return get_low_memory_defaults (buf,type);
}
char *get_accurate4RNA_defaults(char *buf, char *type)
{
  return get_rcoffee_defaults (buf,type);
}
char *get_t_coffee_defaults(char *buf, char *type)
{
  return buf;
}
char *get_fmcoffee_defaults(char *buf, char *type)
   {
     //Fast Mcoffee
     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     if (NO_METHODS_IN_CL) buf=strcat (buf,"-in Mclustalw2_msa  Mmuscle_msa Mmafft_msa -multi_core methods_relax_msa");

     /*buf=strcat (buf,"-in ");*/

     return buf;
   }

char *get_mcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));


     if (NO_METHODS_IN_CL)   buf=strcat (buf,"-in Mclustalw2_msa Mt_coffee_msa Mpoa_msa Mmuscle_msa Mmafft_msa Mdialignt_msa Mpcma_msa Mprobcons_msa -multi_core methods_relax_msa  ");
     /*buf=strcat (buf,"-in ");*/
     return buf;
   }
char *get_dmcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mkalign_msa Mt_coffee_msa Mpoa_msa Mmuscle_msa Mmafft_msa Mdialignt_msa Mprobcons_msa Mamap_msa -multi_core methods_relax_msa");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_rcoffee_consan_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
	 if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mconsan_pair -multi_core templates_relax_msa -dp_mode myers_miller_pair_wise -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -type DNA -relax_lib 0");
     else buf=strcat (buf,"-dp_mode myers_miller_pair_wise -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -type DNA -relax_lib 0");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_rmcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
     if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mprobcons_msa Mmafft_msa Mmuscle_msa -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -check_type -type DNA -relax_lib 0");
     else buf=strcat (buf,"-extend_mode rna2 -template_file RCOFFEE -transform dna2rna -check_type -type DNA -relax_lib 0");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }

//    if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mbest_pair4prot -template_file BLAST -template_file EXPRESSO  ");
   char *get_best4RNA_defaults(char *buf, char *type)
   {

	   if (buf==NULL)buf=vcalloc (1000, sizeof (char));

	   check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
	   buf=strcat (buf," -extend_mode rna2 -template_file PDB,RNA -in Mbest_pair4rna -transform dna2rna -relax_lib 0");
	   /*buf=strcat (buf,"-in ");*/

	   return buf;
   }

char *get_rcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
     buf=strcat (buf," -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -check_type -type DNA -relax_lib 0");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_genepredx_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf, "-method tblastx_msa -evaluate_mode sequences  -genepred  -relax_lib 0 -output fasta_seq,exons,texons,wexons -seq_weight no -check_type -type DNA -out_lib");
     return buf;
   }
char *get_genepredpx_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf, "-method tblastpx_msa -evaluate_mode sequences  -genepred  -relax_lib 0 -output fasta_seq,exons,texons,wexons -seq_weight no -check_type -type DNA -out_lib");
     return buf;
   }

char *get_repeat_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in slow_pair -matrix idmat -out_lib -profile_comparison profile -profile_mode channel_profile_profile -dp_mode myers_miller_pair_wise ");
     /*buf=strcat (buf,"-in ");*/
     return buf;
   }


int check_configuration4program()
   {
     return 1;

   }

/*Chose the right Mode for comparing residues*/

void test ()
{
  char command[1000];
  char *c2;

  c2=vcalloc ( 100, sizeof (char));

  sprintf (command, "cat hmgt_mouseVsnrl3d.blast_result |blast_aln2fasta_aln.pl | fasta_aln2fasta_aln_unique_name.pl > my_test");

  fprintf ( stderr, "C1: %d, C2:%d", is_dynamic_memory (c2), is_dynamic_memory (c2));


  myexit (0);
}

int run_other_pg ( int argc, char *argv[])
{
  //make minimum initialization
  

  if ( strm (argv[0], "seq_reformat") || strm (argv[0], "saltt"))
    {
      return seq_reformat (argc, argv);
    }
  else if ( strm (argv[0], "aln_compare"))
    {
      return aln_compare (argc, argv);
    }
  else if ( strm (argv[0], "analyse_pdb") || strm (argv[0], "apdb") || strm (argv[0], "irmsd") || strm (argv[0], "trmsd"))
    {
      return apdb ( argc, argv);
    }
  else if ( strm (argv[0], "quantile"))
    {
      return quantile ( argc, argv);
    }
  else if ( strstr ( argv[0], "unpack_"))
    {
      unpack_all_perl_script (argv[0]+strlen ("unpack_"));
    }
	else if ( strstr ( argv[0], "fastal"))
	{
		return fastal_main(argc, argv);
	}
  else
    {
      return my_system_cl (argc, argv);
    }
  return EXIT_FAILURE;
}

FILE * t_coffee_tip (FILE *fp,char *mode)
{
  static char **tip;
  static int n;
  int a;
  if ( !tip)
    {
      tip=declare_char ( 100, 300);
      sprintf ( tip[n++],"Get the most accurate protein alignments with: t_coffee <yourseq> -special_mode accurate [Slow]\n");
      sprintf ( tip[n++],"Change the Width of your MSA with the environement variable ALN_LINE_LENGTH (all formats)");
      sprintf ( tip[n++],"Align 2 or more profiles with -profiles= aln1, aln2");
      sprintf ( tip[n++],"-special_mode=expresso to fetch your structures automatically");
      sprintf ( tip[n++],"-special_mode=psicoffee to expand your sequences");
      sprintf ( tip[n++],"-special_mode=accurate The best we can do [slow]");

      sprintf ( tip[n++],"-special_mode=3dcoffee to combine sequences and structures");
      sprintf ( tip[n++],"-special_mode=mcoffee to combine alternative msa methods");
      sprintf ( tip[n++],"-special_mode=dmcoffee to combine alternative msa methods on debian");

      sprintf ( tip[n++],"-usetree=<file> to use your own guide tree");
      sprintf ( tip[n++],"-infile=<aln> -special_mode=evaluate to evaluate your own alignment");
      sprintf ( tip[n++],"-other_pg seq_reformat to access seq_reformat");
      sprintf ( tip[n++],"-other_pg extract_from_pdb to use our pdb retriever");
      sprintf ( tip[n++],"All the latest versions on www.tcoffee.org");
      sprintf ( tip[n++],"-version to check for updates");
      sprintf ( tip[n++],"-output=html will produce a colored output");
      sprintf ( tip[n++],"-outorder=aligned will order the sequences according to the guide tree in newick");
      sprintf ( tip[n++],"-special_mode=quickaln will produce fast/low accuracy alignments");
      sprintf ( tip[n++],"-other_pg seq_reformat -in <aln> -action +trim %%50 Will reduce the redundancy of your MSA");
      sprintf ( tip[n++],"-tip=all to see all the tips, tip=no will prevent them all");
      sprintf ( tip[n++],"-other_pg unpack_all will unpack all the perl scripts within t_coffee");
    }

  if ( strm (mode, "none"))return fp;

  fprintf ( fp, "\n# TIP :See The Full Documentation on www.tcoffee.org\n");

  if (strm ( mode, "all"))
    {
      for ( a=0; a< n; a++)
	{
	  fprintf (fp, "# TIP %d: %s\n", a+1,tip[a]);
	}
    }
  else if ( strm ( mode, "one"))
    {
      int b;
      vsrand(0);
      b=(rand()%(n-1))+1;
      fprintf (fp, "# TIP %2s:  %s","1", tip[0]);
      fprintf (fp, "# TIP %2d:  %s\n", b+1, tip[b]);

    }

  fprintf ( fp, "\n");
  return fp;
}

char* prepare_one2all (char *seq,Sequence *S, char *lib_file)
 {
   int a, n, i;
   FILE *fp;
   char **name, *use_tree;


   if ( S->nseq==2) return NULL;


   if ((i=name_is_in_list (seq,S->name,S->nseq, 100))!=-1);
   else if ( is_number (seq))
     i=atoi(seq)-1;
   else
     return NULL;


   declare_name (use_tree);


   name=declare_char (S->nseq+1, 100);
   for (a=0; a<S->nseq; a++)
     sprintf (name[a], "%s", S->name[a]);
   n=S->nseq;

   if (i!=0)
     {
       sprintf (name[n], "%s", name[i]);
       sprintf (name[i], "%s", name[0]);
       sprintf (name[0], "%s", name[n]);
     }
   sprintf (lib_file, "%s", vtmpnam (NULL));
   fp=vfopen (lib_file, "w");
   for ( a=1; a<n; a++)
     {
       fprintf ( fp, "2 %s %s\n", name[0],name[a]);
     }
   vfclose (fp);


   sprintf ( use_tree, "%s", vtmpnam (NULL));
   fp=vfopen (use_tree, "w");
   vfclose (create_linear_tree (name, n, fp));
   free_char (name, -1);

   return use_tree;
 }
char* prepare_subset2all (char *mode, Sequence *S, char *lib_file, Constraint_list *CL)
 {
   int a,b,s1, n, i;
   FILE *fp;
   char **name;
   int **score;
   int **done;
   int nseq=0;

   if ( S->nseq==2) return NULL;
   name=declare_char (S->nseq+1, 100);
   done=declare_int (S->nseq, S->nseq);
   for (a=0; a<S->nseq; a++)done[a][a]=1;

   if ( check_file_exists (mode))
     {
       Sequence *L;
       L=main_read_seq (mode);
       for (a=0; a< L->nseq; a++)
	 if ( (b=name_is_in_list (L->name[a], S->name,S->nseq, 100))!=-1)
	   {
	     sprintf ( name[nseq++], "%s", L->name[a]);
	   }
     }
   else if ( strm (mode, "_P_"))
     {
       for (a=0; a<S->nseq; a++)
	 {
	   if (seq_has_template (S, a, "_P_"))
	     {
	       sprintf (name[nseq++], "%s", name [a]);
	     }
	 }
     }
   else if ( is_number (mode))
     {
       Sequence *LS;

       nseq=atoi (mode);
       if ( nseq<0)
	 nseq=((float)S->nseq*((float)nseq/(float)100.0)*(float)-1);

       nseq=MIN(nseq,S->nseq);
       if ( nseq>=S->nseq)LS=S;
       else
	 {
	   Alignment *A, *SA;
	   char tmode[1000];
	   A=very_fast_aln (seq2aln (S, NULL, RM_GAP), 0, NULL);
	   sprintf (tmode, "_aln_n%d", nseq);
	   SA=simple_trimseq (A, NULL, tmode, NULL);
	   LS=aln2seq(SA);
	   free_aln (A);
	   free_aln (SA);
	 }
       for (a=0; a<LS->nseq; a++)
	 {
	   sprintf (name[a], "%s", LS->name[a]);
	   fprintf ( stderr, "\n\tMaster Sequence: %s", name[a]);
	 }
       if (LS!=S)free_sequence (LS, LS->nseq);
     }
   else
     {
       printf_exit (EXIT_FAILURE, stderr, "ERROR: %s is neither a file nor a method nor a number for subset2all [FATAL:%s]\n",mode,PROGRAM);
     }

   sprintf (lib_file, "%s", vtmpnam (NULL));
   fp=vfopen (lib_file, "w");
   for (a=0; a<nseq; a++)
     {
       for ( b=0; b<S->nseq; b++)
	 {
	   if(!done[a][b] && !strm (name[a],S->name[b]))fprintf ( fp, "2 %s %s\n", name[a],S->name[b]);
	   done[a][b]=done[b][a]=1;
	 }
     }
   vfclose (fp);

   return NULL;
 }
int set_methods_limits (char ** method,int nl,char **list, int n, int *maxnseq, int *maxlen)
{
  int a,ns, ml, nm=0;
  char string[1000];

  nl/=3;
  for (a=0; a<nl; a+=3)
    {

      sprintf ( string, "M%s", method[a]);
      if ( name_is_in_list (string,list, n, 100)!=-1)
	{

	  ns=atoi(method[a+1]);
	  ml=atoi(method[a+2]);

	  if (ns!=-1 && (maxnseq[0]==-1 || maxnseq[0]>ns))maxnseq[0]=ns;
	  if (ml!=-1 && (maxlen[0]==-1 || maxlen[0]>ml))maxlen[0]=ml;
	  nm++;
	}
    }
  return nm;
}



char * get_seq_type_from_cl (int argc, char **argv)
{
  char *buf, *r;
  char file[100];
  int a;
  int seq=0;
 
  sprintf (file, "%d.tmp", rand()%10000);
  buf=vcalloc ( 1000, sizeof (char));
  sprintf ( buf, "%s ", get_string_variable ("t_coffee"));
  for (a=1, seq=0; a<argc; a++)
    {
      if ( check_file_exists (argv[a]))seq=1;
    }
  if (!seq) return "";
  
  for (a=1; a< argc; a++)
    {
      buf=strcat (buf, argv[a]);
      buf=strcat (buf, " ");
    }

  buf=strcat (buf, " -type_only >");
  buf=strcat (buf, file);
  
  my_system ( buf);
  
  r=file2string (file);
  vremove (file);
  return r;
}

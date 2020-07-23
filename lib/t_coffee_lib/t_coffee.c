#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"
#include "t_coffee.h"
#include "km_coffee_header.h"

/**
 * \file
 * Main file of T-Coffee
 * \author Cedric Notredame
 */


/** \mainpage T-Coffee Index Page
 *
 * \section premise Premise
 * 
 * This documentation is far from comprehensive... in fact only a few functions have been documented yet. 
 * But maybe it will be extended during the years by people working with the T-Coffe code ;)
 * 
 * \section start A Place to Start With
 *
 * In order to get an idea of the main procedure of T-Coffee when calculating a multiple alignment,
 * you should have a look at the ::batch_main function.
 *
 *
 *  \section documented_files Files Already Containing Documentation:
 *
 *	- t_coffee.c
 *	- util.c
 *	- util_constraints_list.c
 *	- aln_convertion_util.c
 *	- util.c
 *	- reformatp.c
 *	- evaluate.c
 *
 *	For datastructures, these files are interesting:
 *	- io_structures.h contains ::Sequence, ::Alignment and ::Template.
 *	- util_constraints_list.h contains the central structure ::Constraint_list, but also ::Distance_matrix and ::TC_method.
 *
 *
 */







/**
 * \file t_coffee.c
 *
*/

static void test();
static char * get_seq_type_from_cl (int argc, char **argv);
static char *get_defaults(char *buf, char *type);
static char *get_evaluate_defaults(char *buf, char *type);
static char *get_genome_defaults(char *buf, char *type);
static char *get_dali_defaults(char *buf, char *type);
static char *get_mcoffee_defaults(char *buf, char *type);
static char *get_xcoffee_defaults(char *buf, char *type);
static char *get_fmcoffee_defaults(char *buf, char *type);
static char *get_t_coffee_defaults(char *buf, char *type);

static char *get_dmcoffee_defaults(char *buf, char *type);
static char *get_rcoffee_consan_defaults(char *buf, char *type);
static char *get_rmcoffee_defaults(char *buf, char *type);//Original R-Coffee Paper
static char *get_rcoffee_defaults(char *buf, char *type);//Original R-Coffee Paper
static char *get_rmcoffee_defaults_old(char *buf, char *type);//Original R-Coffee Paper
static char *get_rcoffee_defaults_old(char *buf, char *type);//Original R-Coffee Paper
static char *get_best4RNA_defaults(char *buf, char *type);
static char *get_saracoffee_defaults(char *buf, char *type);
static char *get_rsapcoffee_defaults(char *buf, char *type);

static char *get_very_fast_defaults(char *buf, char *type);
static char *get_precomputed_defaults(char *buf, char *type);
static char *get_3dcoffee_defaults(char *buf, char *type);
static char *get_expresso_defaults(char *buf, char *type);

static char *get_accurate_defaults(char *buf, char *type);
static char *get_accurate4PROTEIN_defaults(char *buf, char *type);
static char *get_accurate4DNA_defaults(char *buf, char *type);
static char *get_accurate4RNA_defaults(char *buf, char *type);

static char *get_procoffee_defaults(char *buf, char *type);
static char *get_blastr_defaults(char *buf, char *type);
static char *get_psicoffee_defaults(char *buf, char *type);
static char *get_dna_defaults(char *buf, char *type);
static char *get_cdna_defaults(char *buf, char *type);
static char *get_repeat_defaults(char *buf, char *type);
static char *get_sample_defaults(char *buf, char *type);
static char *get_highlow_defaults(char *buf, char *type);

static char *get_low_memory_defaults( char *buf, char *type);

static char *get_genepredx_defaults(char *buf, char *type);
static char *get_genepredpx_defaults(char *buf, char *type);

static int set_methods_limits (char **method_limits,int n_methods_limit,char **list_file, int n_list, int *maxnseq, int *maxlen);
static FILE *t_coffee_tip (FILE *fp,char *mode);

static int run_other_pg (int argc, char *argv[]);

static Sequence* prepare_master (char *seq,Sequence *S,Constraint_list *CL, char *dmode);
char ** km_coffee (int argc, char **argv);
Alignment * t_coffee_dpa (int argc, char **argv);
Alignment *kmir(Alignment *A);
int fastal_main(int argc, char **argv); //this function was declared in fastal.c

#define is_a_seq_file(file) (!is_matrix(file) && !is_matrix(file+1) && !is_method (file) && !is_method (file+1) &&(check_file_exists(file) || check_file_exists(file+1)))
static int NO_METHODS_IN_CL;

int batch_main ( int argc, char **argv);

/**
 * Pass arguments to ::batch_main.
 *
 * \callgraph
 */
int main (int argc, char *argv[])
{
// printf("RUNNING DEBUG\n");
  int r, a;
  

  if (argv[0][0]=='\0'){argv[0]="t_coffee";}
      
  
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


/**
 * Brings everything together, but does not contain computations.
 *
 * \callgraph
 * 
 * \section remark Remark
 * 
 * The follwowing text is just a summation of short lines of documentation appearing throughout
 * the code of this function. Often they make more sense in their specific context,
 * that means while reading the source code. 
 *
 *
 * \section prep Preparation
 *
 * 
 */

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

	int check_configuration;
	int update;
	int garbage;
	int quiet;
	char *parameters;
	char *t_coffee_defaults;
	int t_coffee_defaults_flag;

	FILE *le=NULL;
	char *se_name;
	char *clean_list;
	int debug=0;
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

	char **export_list;
	int n_export;


	char **template_file_list;
	int n_template_file;

	char **template_mode_list;
	int n_template_mode;


	int flip;
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
	int len;
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
	int   extend_seq;
	char *outorder;
	char  *inorder;
	char *output_res_num;
	char *residue_case;
	int extra_cpu;


	char *weight;
	int sample_dp=0;

	char *seq_weight;
	int do_align;
	char *evaluate_mode;
	char *color_mode;
	int aln_line_length=0;
	char *method_evaluate_mode;
	int get_type;
	/*Post_processing*/
	int clean_aln;
	int clean_threshold;
	int clean_iteration;
	char *clean_evaluate_mode;
	/*Profile Alignment*/

	int n_seq_list=0;
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
	int ulimit;
	/*Thread parameters*/
	int psiJ=0;
	int psitrim=0;
	char *psitrim_mode;
	char *psitrim_tree;
	int prot_min_sim;
	int prot_max_sim;
	int prot_min_cov;
	int pdb_min_sim;
	int pdb_max_sim;
	int pdb_min_cov;
	char *pdb_type;


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
	int print_cache;
	
	/*align_pdb*/
	char *align_pdb_param_file;
	char *align_pdb_hasch_mode;

	/*msa_mode*/
	char *use_seqan;
	char *msa_mode;
	char *master_mode;
	Sequence *MASTER_SEQ=NULL;
	Sequence *TEMP_SEQ=NULL;
	char *et_mode;

	int blast_maxnseq;

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
	int trim=0;
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
	int dpa=0;
	int reg=0;
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
	int max_n_proc; //legacy for Nature protocol
	int n_thread;
	
	char *lib_list;
	char *prune_lib_mode;

	int no_error_report;
	int no_warning=0;
	char *tip;
	int run_local_script;
	char *plugins;
	char *plugins_order;
	char *email;
	char *proxy;
	char *tmp_4_tcoffee;

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
	int     display;
	int sand_box=0;



	if (sand_box==1)
	  {
	    sb();
	    exit (0);
	  }
	
	/**
	 * Before anything else, check if we want to use the \b other_pg option of T-Coffee.
	 *        If so, redirect to ::run_other_pg.
	 *        (In case of kmcoffee, load the kmcoffee arguments from ::km_coffee)
	 */

	argv=standard_initialisation (argv, &argc);
	
	set_string_variable ("t_coffee", argv[0]);

	if (argc>=3 && strm (argv[1], "-other_pg"))
	  {
	    set_string_variable ("-other_pg", "set");
	    myexit(run_other_pg (argc-2, argv+2));
	  }
	else if ( name_is_in_list ("kmcoffee", argv, argc, 100)!=-1)
	  {
	    argv=km_coffee(argc, argv); 
	  }
	
	
	 /**
	  * Read all the parameters of T_Coffee using ::get_cl_param
	  *
	 */

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
			    /*Def 1*/     "1"             ,\
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
			    /*DOC*/       "PROTEIN, DNA or RNA. Automatically set, but can be forced with this flag"           , \
			    /*Parameter*/ &type          ,		\
			    /*Def 1*/    ""              ,		\
			    /*Def 2*/    ""              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
					  );
	       /*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (type);
	       get_cl_param(					\
			    /*argc*/      argc           ,	\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,	\
			    /*Name*/      "-dpa"        ,	\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "D"            ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,	\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "Run DPA mode"           ,	\
			    /*Parameter*/ &dpa          ,		\
			    /*Def 1*/    "0"              ,		\
			    /*Def 2*/    "1"              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
					  );
	        /*PARAMETER PROTOTYPE:    INFILE    */
	       
	       declare_name (type);
	       get_cl_param(					\
			    /*argc*/      argc           ,	\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,	\
			    /*Name*/      "-reg"        ,	\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "D"            ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,	\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "Run DPA mode"           ,	\
			    /*Parameter*/ &reg          ,		\
			    /*Def 1*/    "0"              ,		\
			    /*Def 2*/    "1"              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
					  );



	       declare_name (plugins);
	       get_cl_param(					\
			    /*argc*/      argc           ,	\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,	\
			    /*Name*/      "-plugins",		\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "S"          ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,		\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "Set the directory containing the plugins [no if no plugin]",	\
			    /*Parameter*/ &plugins   ,			\
			    /*Def 1*/    "default"       ,		\
			    /*Def 2*/    ""              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
					  );
	       if (!strm (plugins,"default"))
		 {
		   cputenv ("LOCAL_PLUGINS_4_TCOFFEE=%s", plugins);// single use environement variable
		   get_plugins_4_tcoffee();
		 }
	      
	       declare_name (plugins_order);
	       get_cl_param(					\
			    /*argc*/      argc           ,	\
			    /*argv*/      argv           ,	\
			    /*output*/    &le            ,	\
			    /*Name*/      "-plugins_order",	\
			    /*Flag*/      &garbage       ,	\
			    /*TYPE*/      "S"          ,	\
			    /*OPTIONAL?*/ OPTIONAL       ,		\
			    /*MAX Nval*/  1              ,		\
			    /*DOC*/       "first or last. Set the order of the plugins for T-Coffee. First means the plugins are used first. Last means that local installations are used first and the plugins only used if the local installation cannot", \
			    /*Parameter*/ &plugins_order   ,		\
			    /*Def 1*/    "last"       ,		\
			    /*Def 2*/    ""              ,		\
			    /*Min_value*/ "any"          ,		\
			    /*Max Value*/ "any"				\
					  );
	       if (strm (plugins_order,"first"))cputenv4pathFirst(get_plugins_4_tcoffee());
	       else cputenv4pathLast(get_plugins_4_tcoffee());
	       
	       if (reg || dpa)reg=dpa=1;
	       

	       /**
		* Load T-Coffee default parameters.
		* Next, decide whether a special mode of T_Coffee has been called and load the
		* corresponing parameters from a function. Here's an example:
		* \code
		* else if ( strm (special_mode, "psicoffee")) new_arg = get_psicoffee_defaults( NULL,lseq_type );
		* \endcode
		*
		*/
	       

if (dpa)
  {
    dump_io_start (NULL);
    t_coffee_dpa (argc, argv);
  }




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
	register_file4dump(parameters, "r");
	argv=push_string (file2string(pname), argv, &argc, 1);
	t_coffee_defaults=pname;
      }
    else
      {
	t_coffee_defaults=NULL;
      }
  }

 if ( parameters && parameters[0])
   {
     register_file4dump(parameters, "r");
     argv=push_string (file2string (parameters), argv, &argc, 1);
   }
 
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
	 char *new_arg=NULL;//to be pushed after t_coffee
	 char *new_arg2=NULL;//to be pushed last
	 
	 
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
	 else if ( strm (special_mode, "sample"))new_arg=get_sample_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "highlow"))new_arg=get_highlow_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "psicoffee"))new_arg=get_psicoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "procoffee"))new_arg=get_procoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "blastr"))new_arg=get_blastr_defaults(NULL,lseq_type);
	 
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
	 else if ( strm (special_mode, "xcoffee"))new_arg=get_xcoffee_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "dmcoffee"))new_arg=get_dmcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "fmcoffee"))new_arg=get_fmcoffee_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "rcoffee_consan"))new_arg=get_rcoffee_consan_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rmcoffee") ||strm (special_mode, "mrcoffee") )new_arg=get_rmcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rcoffee"))new_arg=get_rcoffee_defaults(NULL,lseq_type);

	 else if ( strm (special_mode, "rcoffee_slow_accurate"))new_arg=get_rcoffee_consan_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rcoffee_fast_approximate"))new_arg=get_rmcoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "t_coffee"))new_arg=get_t_coffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "saracoffee"))new_arg2=get_saracoffee_defaults(NULL,lseq_type);
	 else if ( strm (special_mode, "rsapcoffee"))new_arg2=get_rsapcoffee_defaults(NULL,lseq_type);
	 

	 else if ( strm (special_mode, "unspecified"));
	 else
	   {
	     fprintf ( stderr, "\nERROR: special_mode %s is unknown [FATAL:%s]\n",special_mode, PROGRAM);
	     myexit (EXIT_FAILURE);
	   }
	 
	 if (new_arg)argv=push_string (new_arg, argv, &argc, 1);
	 else if (new_arg2)argv=push_string (new_arg2, argv, &argc,argc);
       }
   }

if ( getenv ("TCOFFEE_EXTRA_PARAM"))argv=push_string (getenv ("TCOFFEE_EXTRA_PARAM"), argv, &argc, argc);

	       
dump_io_start (NULL);
argv=break_list ( argv, &argc, "=;, \n");
argv=merge_list ( argv, &argc);

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
     
     get_cl_param(						\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-debug"    ,\
			    /*Flag*/      &quiet        ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "0 [default]: no dump; 1: dump the input, 2: dump input and keep tmp files", \
			    /*Parameter*/ &debug      ,\
			    /*Def 1*/     "0"      ,\
			    /*Def 2*/     "1"   ,\
			    /*Min_value*/ "0"         ,\
			    /*Max Value*/ "2"          \
		   );
     if (debug==0)
       {
	 cputenv ("ERRORFILE_4_TCOFFEE=NO");
       }
     else if (debug==2)
       {
	 cputenv ("DEBUG_TMP_FILE=1");
       }
     /*PARAMETER PROTOTYPE: MEM MODE*/
     declare_name(clean_list);
     get_cl_param(						\
		  /*argc*/      argc          ,			\
		  /*argv*/      argv          ,			\
		  /*output*/    &le           ,			\
		  /*Name*/      "-clean"   ,		\
		  /*Flag*/      &garbage      ,			\
		  /*TYPE*/      "S"          ,			\
		  /*OPTIONAL?*/ OPTIONAL      ,			\
		  /*MAX Nval*/  1             ,			\
		  /*DOC*/       "Will delete cached data and exit: all, cache, lock, tmp. It is possible to specify a list: cache_lock_tmp" ,			\
		  /*Parameter*/ &clean_list     ,			\
		  /*Def 1*/     "no"         ,			\
		  /*Def 2*/     "all"            ,			\
		  /*Min_value*/ "any"         ,			\
		  /*Max Value*/ "any"				\
				);
     


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


	 /**
	  * Special commands, that usually end T-Coffee are handled next:
	  * check_configuation, update, version etc.
	  */
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
	 fprintf ( stdout, "PROGRAM: %s %s (%s)\n",PROGRAM,VERSION,BUILD_INFO);
	 return EXIT_SUCCESS;
       }

     le=get_stdout1(se_name);
     fprintf ( le, "\nPROGRAM: %s %s (%s)\n",PROGRAM,VERSION,BUILD_INFO);
     

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
	       if (n_template_file)cputenv ("template_file_4_TCOFFEE=%s",template_file_list[0]);
	       
/*PARAMETER PROTOTYPE:    VERSION      */
	       setenv_list=declare_char (100, STRING);
	       n_setenv=get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-setenv"        ,\
			    /*Flag*/      &garbage        ,\
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

/*PARAMETER PROTOTYPE:    VERSION      */
	       export_list=declare_char (100, STRING);
	       n_export=get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-export"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  2                ,\
			    /*DOC*/       "Declares a parameter variable" ,\
			    /*Parameter*/ export_list          ,	\
			    /*Def 1*/     "0"             ,		\
			    /*Def 2*/     "1"             ,		\
			    /*Min_value*/ "0"            ,		\
			    /*Max Value*/ "1"				\
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
/*PARAMETER PROTOTYPE:        flip*/

	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-flip"      ,\
			    /*Flag*/      &flip     ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "flip sequences"          ,\
			    /*Parameter*/ &flip   ,\
			    /*Def 1*/     "0"          ,\
			    /*Def 2*/     "50"      ,\
			    /*Min_value*/ "0"         ,\
			    /*Max Value*/ "100"         \
		   );
	       set_int_variable ("flip", flip);

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
	profile_list=declare_char ( 2001, STRING);
	n_profile_list=get_cl_param(				\
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
	
	if (n_profile_list)
	  {
	    profile_list[n_profile_list]=NULL;
	    profile_list=list2expanded_flist(profile_list,&n_profile_list, "FILE::");
	  }
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
			    /*DOC*/       "nj, upgma, cwph,kmeans",\
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
			    /*DOC*/       "Specifies one or many formats that must be output: clustalw_aln, msf_aln, tcs_[residue,column]_[filter,lower][0-9], tcs_[weighted,replicate][Nreplicates],sp_ascii, score_ascii . The file extension is the output format "           ,\
			    /*Parameter*/ out_aln_format,\
			    /*Def 1*/    "aln,html"           ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );

	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-len"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Line Length\n",\
			    /*Parameter*/ &len        ,\
			    /*Def 1*/    "0",\
			    /*Def 2*/    "100",\
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
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-extend_seq"        ,\
			    /*Flag*/      &extend_seq       ,\
			    /*TYPE*/      "FL"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "extend the sequences",	\
			    /*Parameter*/ &extend_seq          ,\
			    /*Def 1*/    "0"              ,\
			    /*Def 2*/    "1"              ,\
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

/*PARAMETER PROTOTYPE:    -ulimit     */
	       ulimit=-1;
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-ulimit"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Maximum amount of memory to be used. Kill job otherwise"          ,\
			    /*Parameter*/ &ulimit       ,\
			    /*Def 1*/    "-1"            ,\
			    /*Def 2*/    "0"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       if (ulimit!=-1)set_max_mem (ulimit);

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
			    /*Def 1*/    "-1"            ,\
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

	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-sample_dp"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "defines the tie breaking strategy (only with gotoh_pair_wise)" ,\
			    /*Parameter*/ &sample_dp          ,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "2"             \
		   );
	       if ( sample_dp)cputenv ("SAMPLE_DP_4_TCOFFEE=%d", sample_dp);


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
			    /*Def 2*/    "default"             ,\
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
			    /*Def 1*/    "no"             ,\
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
	       declare_name (color_mode);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-color_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Mode used to produce the color output:new (default) or old  " ,\
			    /*Parameter*/ &color_mode          ,\
			    /*Def 1*/    "new"             ,\
			    /*Def 2*/    "old"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       if (!strm (color_mode, "default"))cputenv ("COLOR_4_TCOFFEE=%s", color_mode);
	     
/*PARAMETER PROTOTYPE:    WEIGHT      */
	       
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-aln_line_length"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Mode used to produce the color output:t_coffee_fast,t_coffee_slow  " ,\
			    /*Parameter*/ &aln_line_length          ,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "0"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       if (aln_line_length)cputenv ("ALN_LINE_LENGTH=%d", aln_line_length);

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
			    /*Def 1*/    "triplet"             ,\
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
			    /*DOC*/       "Minimum similarity between a sequence and its BLAST relatives" ,\
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
			    /*Def 1*/     "100"             ,\
			    /*Def 2*/     "50"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
set_int_variable ("prot_max_sim", prot_max_sim);

get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-psiJ"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Defines the number of iteration of psiblast (-j)",\
			    /*Parameter*/&psiJ       ,\
			    /*Def 1*/    "3"      ,\
			    /*Def 2*/    "3",\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 
cputenv ("psiJ_4_TCOFFEE=%d", psiJ);
cputenv ("num_iterations_4_TCOFFEE=%d", psiJ);

declare_name (psitrim_mode);
get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-psitrim_mode"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Mode used to trim profiles, regtrim or trim (def)",\
			    /*Parameter*/&psitrim_mode       ,\
			    /*Def 1*/    "regtrim"      ,\
			    /*Def 2*/    "regtrim",\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

cputenv ("psitrim_mode_4_TCOFFEE=%s", psitrim_mode);

declare_name (psitrim_tree);
get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-psitrim_tree"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Mode used to compute the tree when using regtree to trim profiles (codnd def)",\
			    /*Parameter*/&psitrim_tree       ,\
			    /*Def 1*/    "codnd"      ,\
			    /*Def 2*/    "codnd",\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
cputenv ("psitrim_tree_4_TCOFFEE=%s", psitrim_tree);


get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-psitrim"        ,\
			    /*Flag*/      &psitrim        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Maximum number of sequences to keep when building a profile [0 to keep everything, negative value to keep X%%, positive value to keep X Sequences]" ,\
			    /*Parameter*/ &psitrim          ,\
			    /*Def 1*/     "100"             ,\
			    /*Def 2*/     "100"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
set_int_variable ("psitrim", psitrim);
cputenv ("psitrim_4_TCOFFEE=%d",psitrim);

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
			    /*Def 1*/     "90"             ,\
			    /*Def 2*/     "0"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
set_int_variable ("prot_min_cov", prot_min_cov);
declare_name(pdb_type);
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-pdb_type"        ,\
			    /*Flag*/      &pdb_min_sim        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "d: diffraction, n: nmr, m:model" ,\
			    /*Parameter*/ &pdb_type          ,\
			    /*Def 1*/     "d"             ,\
			    /*Def 2*/     "d"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
set_string_variable ("pdb_type", pdb_type);
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
 if ( strm (prot_blast_server, "env"))prot_blast_server=get_env_variable ("blast_server_4_TCOFFEE",IS_FATAL);
 set_string_variable ("blast_server", prot_blast_server);
 cputenv ("blast_server_4_TCOFFEE=%s",prot_blast_server);


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
 cputenv ("pdb_db_4_TCOFFEE=%s",pdb_db);

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
	 		    /*Def 1*/    "uniref50"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 if ( strm (prot_db, "env"))prot_db=get_env_variable ("protein_db_4_TCOFFEE", IS_FATAL);
 set_string_variable ("prot_db", prot_db);
 cputenv ("protein_db_4_TCOFFEE=%s",prot_db);
 
 
 
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
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-print_cache"    ,\
			    /*Flag*/      &print_cache      ,\
			    /*TYPE*/      "FL"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "print the cache dir to stdout and exit\n",\
			    /*Parameter*/ &print_cache       ,		\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "1"      ,\
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
declare_name (et_mode);
 get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-et_mode"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Algorithm used to the et score: id, et, sankoff, sp"          ,\
			    /*Parameter*/ &et_mode      ,\
			    /*Def 1*/    "et"      ,\
			    /*Def 2*/    "et"      ,		\
			    /*Min_value*/ "any"         ,	\
			    /*Max Value*/ "any"			\
					  );
 
declare_name (master_mode);
get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-master"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Align all the sequences to the master sequences: file or number"          ,\
			    /*Parameter*/ &master_mode      ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "_LONG_n_100_kmeans_"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
 get_cl_param(							\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-blast_nseq"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Maximum number of querries for BLAST (0: all)"          ,\
			    /*Parameter*/ &blast_maxnseq      ,\
			    /*Def 1*/    "0"      ,\
			    /*Def 2*/    "0"      ,\
			    /*Min_value*/ "0"         ,\
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
split=0;

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
			    /*DOC*/       "Multi core: template_jobs_relax_[msa|pairwise]_evaluate",\
			    /*Parameter*/ &multi_core   ,\
			    /*Def 1*/    "templates_jobs_relax_msa_evaluate"       ,\
			    /*Def 2*/    "templates_jobs_relax_msa_evaluate"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "100"           \
		   );
	       if (multi_core[0])set_string_variable ("multi_core",multi_core);
/*PARAMETER PROTOTYPE:    multi_core    */
	       n_core=0;
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-n_core",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Number of cores to be used by machine [default=1, 0=> all those defined in the environement]",\
			    /*Parameter*/ &n_core   ,\
			    /*Def 1*/    "1"       ,\
			    /*Def 2*/    "1"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "10000"           \
		   );
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-thread",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Number of cores to be used by machine [default=1, 0=> all those defined in the environement]",\
			    /*Parameter*/ &n_thread   ,\
			    /*Def 1*/    "1"       ,\
			    /*Def 2*/    "1"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "10000"           \
		   );
	      
	       

/*PARAMETER PROTOTYPE:    max_n_proc:: legacy for the nature protocol    */
	       max_n_proc=0;
	       get_cl_param(					\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-max_n_proc",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Number of cores to be used by machine [default=1, 0=> all those defined in the environement]",\
			    /*Parameter*/ &max_n_proc   ,\
			    /*Def 1*/    "1"       ,\
			    /*Def 2*/    "1"              ,\
			    /*Min_value*/ "0"          ,\
			    /*Max Value*/ "10000"           \
		   );

	       //-n_core and -max_n_core are deprecated. everything should arrive via -thread
	       //if any of the three is set 0, ALL cores are used
	       //if any of the three are set top >1 , the highest value is used
	       
	       if (!n_core || !max_n_proc || !n_thread){n_core=get_nproc();}
	       else n_core=MAX3(n_core, max_n_proc,n_thread);
	       set_int_variable ("n_core",n_core);
	       set_nproc (n_core);
	       cputenv ("thread_4_TCOFFEE=%d", get_nproc());
	       
	       

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
	       register_file4dump(lib_list, "r");
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
			    /*Def 1*/    "none"       ,\
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
	       
	       overaln_P1=0;
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
	       
	       overaln_P2=0;
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

		overaln_P3=0;
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
	       
	       overaln_P4=0;
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


	       

	       display=100;
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-display",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "Sets the threshold (nseq) for the full display of the groups. -1 results in a full display, and 0 in no display at all",\
			    /*Parameter*/ &display   ,\
			    /*Def 1*/    "100"       ,\
			    /*Def 2*/    "-1"              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
			   );
	       if (display)set_int_variable ("display",display);
	       
/*******************************************************************************************************/
/*                                                                                                     */
/*                           TCoffee_dpa Parameter:END                                                 */
/*                                                                                                     */
/*******************************************************************************************************/
	      

	       
	       if (argc==1 )
		 {
		   display_method_names ("display", stdout);
		   return EXIT_SUCCESS;
		 }
	       get_cl_param( argc, argv,&le, NULL,NULL,NULL,0,0,NULL);
	       prepare_cache (cache);

	       if (print_cache)
		 {
		   fprintf (stdout, "%s\n", get_cache_dir());
		   myexit (EXIT_SUCCESS);
		 }
	       

/*******************************************************************************************************/
/*                                                                                                     */
/*                           TCoffee clean mode                                                        */
/*                                                                                                     */
/*******************************************************************************************************/
	       
	       if (!strstr (clean_list, "no"))
		 {
		   char command[10000];
		   if (strstr (clean_list, "all") || strstr (clean_list, "cache"))
		     {
		       sprintf (command, "rm -rf %s", getenv ("CACHE_4_TCOFFEE")); 
		       if ( !strstr (command, "cache"))
			 {
			   fprintf ( stderr, "For security reasons the cache dir must contain the string cache.\nYour cached data seems to be stored in [%s]\nYou must delete it manually.",getenv ("CACHE_4_TCOFFEE")); 
			 }
		       else
			 {
			   
			   printf_system_direct ("%s",command);
			 }
		     }
		   if (strstr (clean_list, "all") || strstr (clean_list, "lock"))
		     {
		       
		       sprintf (command, "rm -rf %s", getenv ("LOCKDIR_4_TCOFFEE")); 
		       if ( !strstr (command, "lck") &&!strstr (command, "lock")  )
			 {
			   fprintf ( stderr, "For security reasons the lock dir must contain the string lck.\nYour lock data seems to be stored in [%s]\nYou must delete it manually.",getenv ("LOCKDIR_4_TCOFFEE")); 
			 }
		       else
			 {
			   printf_system_direct ("%S",command);
			 }
		       
		     }
		   if (strstr (clean_list, "all") || strstr (clean_list, "tmp"))
		     {
		       sprintf (command, "rm -rf %s", getenv ("ROOT_TMP_4_TCOFFEE"));
		       if ( !strstr (command, "tmp"))
			 {
			   fprintf ( stderr, "For security reasons the tmp dir must contain the string tmp.\nYour tmp data seems to be stored in [%s]\nYou must delete it manually.",getenv ("ROOT_TMP_4_TCOFFEE")); 
			 }
		       else
			 {
			   printf_system_direct ("%S",command);
			 }
		     }
		   exit (EXIT_SUCCESS);
		 }
	       
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
	       
	       /**
	        *  Array \c list_file[] contains the names of all input files like sequences, profiles and templates,
	        *  but also methods.
	        *  There is a convention in T-Coffee to start the filename of with a
	        *  capital letter specifying its type:
	        *      - \b M for methods
	        *      - \b R for profiles
	        *      - \b A for alignments
	        *      - \b L for libraries (via the -lib flag)
	        *      - \b P for PDB structures
	        *   
	        *   See ::read_constraint_list for more information on how these files will be processed.
	        */

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
		   if (format_is_conc_aln (profile_list[a]))
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
		   //
		   
		   if (check_file_exists(seq_list[a]))
		     sprintf (list_file[n_list++], "S%s",seq_list[a]);
		   else if ( check_file_exists (seq_list[a]+1) && seq_list[a][0]=='S')
		     sprintf (list_file[n_list++], "%s",seq_list[a]);
		   else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: %s does not exist [FATAL]",seq_list[a]);

		 }
	      
	       /*introduce the alignments from the -aln flag*/
	       //Important: Must be introduced AFTER the profiles
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
		   nn=(int*)vcalloc ( 256, sizeof (int));
		   for (a=0; a<n_list; a++)
		     {
		       
		       if ( !check_file_exists(list_file[a]))
			 {
			   nn[(int)list_file[a][0]]++;
			   
			 }
		       else
			 {

			   if (is_seq (list_file[a]))nn['S']++;
			   else if ( is_aln (list_file[a]))nn['A']++;
			   else if ( is_lib (list_file[a]))nn['L']++;
			   else if ( is_method (list_file[a]))nn['M']++;
			   else
			     {
			       add_warning (stderr, "File %s was not properly tagged. Potential ambiguity\n",list_file[a]);
			     }
			 }
		     }
		   
		   
		   /**
		    * If no alignment A, library L or method M is in the list of infiles (that means in \c list_files[]),
		    * T-Coffee will add a default method here. The default method is currently \b Mproba_pair.
		    */
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
	       
	       
	       /**
	     * Check the content of each input file and report possible errors.
	     */
	       else if ( argv[1][0]!='-' && (check_file_exists( argv[1]) || check_file_exists(argv[1]+1)))
		 {

		   if (check_file_exists(argv[1]))F=parse_fname(argv[1]);
		   else if ( check_file_exists(argv[1]+1))F=parse_fname(argv[1]+1);
		 
		 }
	       else if (infile[0])
		 {
		   
		   
		   if (isvtmpnam (infile))
		     {
		       F=parse_fname("output");
		     }
		   else if ( strstr (infile,"stdin"))F=parse_fname(infile);
		   else if ( check_file_exists (infile))F=parse_fname(infile);
		   else if (check_file_exists (infile+1))F =parse_fname(infile+1);
		   
		 }
	       else if ( exon_boundaries && exon_boundaries[0])
		 {
		   if ( check_file_exists (exon_boundaries))F=parse_fname(exon_boundaries);
		   else if (check_file_exists (exon_boundaries+1))F =parse_fname(exon_boundaries+1);
		 }
	       else if (n_seq_list)//use sequence names for the default name, in priority
		 {
		   F=parse_fname(seq_list[0]);
		 }
	       else if (n_aln_file_list)//use sequence names for the default name, in priority
		 {
		   F=parse_fname(aln_file_list[0]);
		 }

	       else
	          {

		  for ( a=0; a< n_list; a++)
		      {
			if (!is_method(list_file[a]))
			  {
			    if ( check_file_exists( list_file[a]))
			       {
				 F=parse_fname(list_file[a]);break;
			       }
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
	       

	       /**
	        * First non-error output: List the Input files, which are usually a sequence
	        * file (S) and one or several methods (M), as shown in the example:
	        * \code
	        * INPUT FILES
        	*         Input File (S) sh3.fasta  Format fasta_seq
			*         Input File (M) proba_pair
	        * \endcode
	        */
	       
	       fprintf (le, "\nINPUT FILES\n");
	       for ( a=0; a< n_list; a++)
		   {
		     fprintf (le, "\tInput File (%c) %s ",list_file[a][0],list_file[a]+1);
		     if ( list_file[a][0]=='A' || list_file[a][0]=='S' || list_file[a][0]=='P'|| list_file[a][0]=='R' )
		       {
			 fprintf (le, " Format %s\n", f=identify_seq_format ( list_file[a]+1));

			 if (!f || f[0]=='\0')
			   {
			     printf_system_direct ("cp %s wrong.file", list_file[a]+1);
			     myexit (fprintf_error(stderr,"The format of %s is not supported", list_file[a]+1));

			   }
			 vfree (f);
		       }
		     else fprintf (le, "\n");
		   }

	       
/*CONVERT, ALIGN OR EVALUATE: CHOSE THE RIGHT VERB*/
	       /*Set the Hierarchy of the verbs*/
	       /*The first one decides...*/


	       do_list=(int**)vcalloc ( 100, sizeof (int*));
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
	      




	       /**
	        * \section input Input Sequences and Libraries
	        *
	        * Read the Sequences (::read_seq_in_n_list). Moreover, for these methods that need Blast scans on the sequences,
	        * ::precompute_blast_db is called.
	        * \note The ::precompute_blast_db method is currently invoked only for methods that contain the substring "blast".
	        *       As far as I know, it's main purpose is the method \b blastr_pair which is used in the special mode blastr.
	        *       See the source code of ::precompute_blast_db for mmore information.
	        * */
	       
	       S=read_seq_in_n_list   (list_file, n_list, type,seq_source);
	       
	       if (maxnseq!=-1 && S->nseq>maxnseq)
		 {
		   printf_exit ( EXIT_FAILURE,stderr, "\nNumber of sequences (%d) exceeds the allowed maximum (%d) [FATAL:%s]", S->nseq,maxnseq, PROGRAM);
		 }
	       
	       S=precompute_blast_db(S,method_list, n_method_list);

	       if ( check_type)
		 {
		   if (!strm (S->type, get_array_type (S->nseq, S->seq)))
		     {
		       fprintf ( stderr, "\nINCORRECT SEQUENCE TYPE (USE %s ONLY) [FATAL:%s]", S->type, PROGRAM);
		       myexit (EXIT_FAILURE);
		     }
		 }


	       if (S->nseq<1 && !do_domain)
		 {
		   printf_exit (EXIT_FAILURE,stderr, "\nERROR: Your Dataset Contains %d Sequence. For multiple alignments you need at least 2 sequences[FATAL:%s]", S->nseq,PROGRAM);
		 }

	       store_seq_type (S->type);

	       if ( type_only==1)
		 {
		   fprintf ( stdout, "%s", S->type);
		   return EXIT_SUCCESS;
		 }
	       /**
	        * Translate Sequences from RNA to DNA or vice versa,
	        * if wanted (::transform_sequence).
	        *
	        */
	       
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







	       if ( get_type)
		 {
		   S=get_sequence_type (S);
		   fprintf ( stdout , "%s\n", S->type);
		   free_sequence(S, S->nseq);
		   return 1;
		 }

	       /**
	        * Reorder the sequences, if wanted, by using ::reorder_seq.
	        * You can make T-Coffee order the sequences lexicographically by their names by
	        * using the \c -inorder=\c aligned flag. If it is set, the order of the sequences
	        * is specified by ::sort_string_array(new_oder), otherwise the order from the
	        * input file (i.e. the order in the Sequence::name array) is kept.
	        */
	       new_order=duplicate_char (S->name, -1, -1);
	       if ( strm (inorder, "aligned"))new_order=sort_string_array   (new_order, S->nseq);

	       initial_order=duplicate_char (S->name, -1, -1);
	       S=reorder_seq(S,new_order,S->nseq);
	       free_char (new_order, -1);



	       /**
	        * Prepare the global ::Constraint_list object using ::declare_constraint_list.
	        * This function allocates memory for a Constraint_list and sets the Constraint_list::residue_index
	        * to default values.
	        */
	      
	       CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?tmpfile():NULL, NULL);

	       sprintf ( CL->method_evaluate_mode, "%s", method_evaluate_mode);

	       (CL->TC)->use_seqan=use_seqan;
	       CL->local_stderr=le;
	       /*set the genepred parameters*/
	       sprintf ( CL->genepred_score, "%s", genepred_score);
	       /*Estimate the distance Matrix*/


	       if (extend_seq)extend_seqaln(CL->S,NULL);
	       //removed to prevent systematic pw distance computation
	       //CL->DM=cl2distance_matrix ( CL,NOALN,distance_matrix_mode, distance_matrix_sim_mode,1);
	       if (extend_seq)unextend_seqaln(CL->S,NULL);



	       /**
	        * Identify Master Sequences for the master mode (option \c -master).
	        * Per default, master mode is not used (\c -master=\c no), but there are
	        * several ways to specify master sequences. The master mode is experimental only,
	        * the idea is not to run all pairwise alignments, but only some of them.
	        * See ::prepare_master (best in the soruce code) for more details.
		* Note that when running master, the library must be relaxed in order to compensate for missing pairs
	        */
	       fprintf ( le, "\nIdentify Master Sequences [%s]:\n", master_mode);
	       (CL->S)->MasterS=MASTER_SEQ=prepare_master(master_mode,S,CL, "_kmeans_");
	       fprintf ( le, "\nMaster Sequences Identified");
	       if (!blast_maxnseq)CL->o2a_byte=(CL->S)->nseq;
	       else CL->o2a_byte=blast_maxnseq;


	       if (MASTER_SEQ)
		 {
		   TEMP_SEQ=S;
		   S=MASTER_SEQ;
		 }







	       /**
	        * \section templates Get the Templates
	        *
	        * Call ::seq2template_seq for each template file in \c template_file_list. This will read/fetch template
	        * files and incorporate them into the ::Seqeunce object. If no template file is given but the Sequence
	        * object already contains templates, these are written into a \c .template_file file.
	        * Next,the same procedure is repeated for profile_template_files (\c -profile_template_file), using
	        * ::profile_seq2template_seq
	        *
	        * \sa ::Template for more information on how templates are stored within T-Coffee.
	        */

	       if (name_is_in_list("BLAST",template_file_list, n_template_file,-1)!=-1  && strstr ("LOCAL", prot_blast_server))
		 {
		   cputenv ("db_4_BLAST=%s", prot_db);
		   cputenv ("num_iterations_4_BLAST=%d", psiJ);
		   cputenv ("outdir_4_BLAST=%d", get_cache_dir());
		   cputenv ("thread_4_BLAST=%d",get_nproc());
		   
		   fprintf ( le, "\nPrecompute the Blasts -- Use Cache if available\n");
		   seq2blast (S);
		 }
		     
		     
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
			   add_warning (stderr, "Impossible to find %s Templates\nCheck that your blast server is properly installed [See documentation][FATAL:%s]\n", template_file_list[a],PROGRAM);
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

	       //Release Master Sequences
	       if (MASTER_SEQ )
		 {
		   int i;
		   S=TEMP_SEQ;
		   for (a=0; a< MASTER_SEQ->nseq; a++)
		     if ((i=name_is_in_list (MASTER_SEQ->name[a], S->name, S->nseq, 100))!=-1)
		       {
			 S->T[i]=MASTER_SEQ->T[a];
		       }
		 }

	       /**
	        * Add the sequences as their own template via the generic ::seq2template_seq and the option \c "SELF_S".
	        * Finally, specify the types of used templates in the Sequence object by ::seq2template_type.
	        */
	       S=seq2template_seq(S, "SELF_S_",F);
	       S=seq2template_type (S);

	       le=display_sequences_names   ( S, le, check_pdb_status, TEMPLATES);







		  /**
		   * Set further options like the potential substitution matrix, the \c -filter_lib flag, the profile_mode
		   * (see ::get_profile_mode_function), some BLAST parameters, some PDB-BLAST parameters and the multicore mode,
		   * before we can go on to the real stuff.
		   */
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








	       /**
	        * \section libcompilation Compile the Library
	        */


	       /**
	        * \subsection read_CL Build Up the Initial Constraint_list
	        *
	        * Finally we can start with the actual algorithm to build an alignment. In this step
	        * we build the initial ::Constraint_list from various sources. Usually, the initial
	        * Constraint_list will contain the information from pairwise alignment methods (per
	        * default \c proba_pair), but of course there are different modes of T-Coffee that use
	        * structural pairwise alignments or the output of multiple aligners. However, this is
	        * the place where all this information is produced and collected. The function
	        * that summarizes all this is ::read_n_constraint_list and its documentation is recommended
	        * reading.
	        *
	        * Just to give a short overview: It will distibute the process of generating pairwise alignments
	        * (i.e. calling external aligners or built-ins like proba_pair) over several CPUs and incorporate
	        * their outputs in one single ::Constraint_list object called \b CL. There are some maybe
	        * non-intuitive decisions done in the process of reading the outputs, which should be documented
	        * in the function.
	        * \sa Constraint_list::residue_index to understand how \b CL stores edges/constraints.
	        *
	        */
	       
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
	       //This very important step insures the CL is symetrical
	       //This is essential for the remaining computation
	       //It should not be done before because some methods may not be symetrical for specific reasons
	       //And the CL may be used in different contexts
	       CL=CL2simCL (CL);

	       
	       if ( CL->M)clean_aln=0;

	       if ( is_number (weight))set_weight4constraint_list (CL, atoi(weight));

	       /**
	        * Some methods, like \c proba_pair, need a cleanup afterwards to free static memory allocated during their process.
	        * That's why ::free_pair_wise is called here.
	        */
	       free_pair_wise ();
	       
	       
	       
	
	       /**
	        * If the Constraint_list should, for some reason, be empty afterwards, report an Error.
	        */
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

	       /**
	        * If an \c -extend has been set to a value greater than 0 an initial filtering
	        * of the Constraint_list is performed, see ::filter_constraint_list.
	        */
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
		    CL->gop=get_avg_matrix_mm ( CL->M, const_cast<char*>( (strm3((CL->S)->type,"PROTEIN", "Protein", "protein")?AA_ALPHABET:"gcuta")) )*10;
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

	       /**
	        * Choose the extension mode given by the \c -extend flag, see ::choose_extension_mode.
	        * Default values seems to be \c very_fast_triplet (??)
	        * \sa Constraint_list::evaluate_residue_pair for more about the Extension mode.
	        */
	       CL=choose_extension_mode (extend_mode, CL);
	       CL->max_n_pair=max_n_pair;

	       processed_lib=0;
	       //use the L, vener touch it again



	       /**
			* \subsection extend_CL Consistency: Relaxation, Filtering, Extension!
			*
			* This is where the magic happens!
			*
			* First, if \c -out_lib is used, the library will be written into the specified file
			* via the ::save_constraint_list function. If you demand an extended library (by
			* \c -out_lib_mode = \c extended) then first the extension steps are performed and
			* then the library is written to a file via ::save_extended_constraint_list.
			* If this was all you wanted to do (\c lib_only), exit the program afterwards.
			*
			*/
	       
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
	      

	      /**
	       * Extension step:
	       * -# If \c filter_lib > 0 then reduce the Constraint_list::residue_index by
	       *    removing edges with weight smaller or equal to \c filter_lib.
	       *    Usually not used right here. See ::filter_constraint_list.
	       * -# Only if the environment variable \c EXTEND4TC is set to 1, call the ::extend_constraint_list.
	       * -# For a=0 to \c relax_lib (which is 1 per default) call ::relax_constraint_list.
	       * -# For a=0 to \c shrink_lib call ::shrink_constraint_list, but \c schrink_lib is usually set to 0.
	       * -# A call to ::evaluate_constraint_list_reference determines the maximum edge weights.
	       *
	       * We recommend reading the documentation of each of these functions to get to grips with the extension.
	       */

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


	       /**
	        * If \c weight is not set to no (but to \c sim or \c matrix for example),
	        * then compute a new distance matrix (::cl2distance_matrix) and ::weight_constraint_list
	        * with the \c seq_weight flag.
	        */
	       if ( !do_convert && !strm (weight, "no"))
	       	   {
	    	   CL->DM=cl2distance_matrix (CL, NOALN, NULL, NULL,0);

	    	   CL=weight_constraint_list(CL, seq_weight);

	    	   if (output_seq_weights (CL->W, outseqweight))
	    		   CL->local_stderr=display_output_filename( CL->local_stderr,"WEIGHT","tc_weight",outseqweight, CHECK);

	    	   le=display_weights(CL->W, le);
	       	   }



	      /**
	       * Quadruplet mode (\b experimental). Only if \c nseq_for_quadruplet > 0.
	       *
	       */
	      if ( nseq_for_quadruplet && !strm(seq_name_for_quadruplet[0], "all"))
			{
			  CL->nseq_for_quadruplet=nseq_for_quadruplet;
			  CL->seq_for_quadruplet=(int*)vcalloc ((CL->S)->nseq, sizeof (int));
			  for (a=0; a< CL->nseq_for_quadruplet; a++)
				{
				  printf ( "\nquad: %s", seq_name_for_quadruplet[a]);
				  if ( (b=name_is_in_list (seq_name_for_quadruplet[a],(CL->S)->name,(CL->S)->nseq, 100))!=-1)CL->seq_for_quadruplet[b]=1;
				  else add_warning ( stderr, "Sequence %s is not in the set and cannot be used for quadruplet extension\n",seq_name_for_quadruplet[a]);
				}
			}
	      else if ( nseq_for_quadruplet && strm(seq_name_for_quadruplet[0], "all"))
			{

			  CL->nseq_for_quadruplet=(CL->S)->nseq;
			  CL->seq_for_quadruplet=(int*)vcalloc ((CL->S)->nseq, sizeof (int));
			  for (a=0; a< CL->nseq_for_quadruplet; a++)
				{
				  CL->seq_for_quadruplet[a]=1;
				}
			}



	      /**
	       * \section prog_aln Progressive Alignment
	       *
	       * At first, transform the ::Seqeunce into an ::Alignment via the small function ::seq2aln and
	       * use ::ungap_array to remove unwanted gaps.
	       * Then you have to choose the right \c msa_mode for evaluating columns. Right now,
	       * there are several options in the code, but only a few of them are used in general.
	       *    - \c sorted_aln
	       *    - \c seq_aln
	       *
	       * \todo The documentation of this function is not yet finished
	       */

	      if ( do_align )
		   {

		   A=seq2aln  ((CL->S),NULL,1);
		   ungap_array(A->seq_al,A->nseq);

		   /*Chose the right Mode for evaluating Columns*/

		   if ( A->nseq==1);
		   else if ( strm ( msa_mode, "etcoffee"))
		     {
		       set_string_variable ("et_mode", et_mode);
		       //CL->get_dp_cost=slow_get_dp_cost;
		       if (!CL->M)
			 {
			   CL->M=read_matrice ("blosum62mt");
			 }
		       if (!CL->gop)CL->gop=-20;
		       if (!CL->gep)CL->gep=-2;
		       
		       CL->get_dp_cost=get_dp_cost_sankoff_tree;
		       CL->pair_wise=gotoh_pair_wise;
		       
		       pc=tree_file;
		       if ( strm (tree_file, "default") || !check_file_exists (tree_file))
			 {
			   T=make_tree ( A,CL,gop, gep,(CL->S),pc,maximise);
			 }
		       else if ( strm (tree_file, "no"))
			 T=make_tree ( A,CL,gop, gep,(CL->S),NULL, maximise);
		       else
			 {
			   fprintf ( le, "\nREAD PRECOMPUTED TREE: %s\n", pc);
			   T=read_tree (pc,&tot_node,(CL->S)->nseq,  (CL->S)->name);
			 }
		       vfclose (tree2file ((T[3][0]),CL->S, "newick",vfopen (pc, "w")));
		       A->tname=(char*)vcalloc ( strlen (pc)+1, sizeof(char));
		       sprintf (A->tname, "%s", pc);
		       SNL=tree_aln ((T[3][0])->left,(T[3][0])->right,A,(CL->S)->nseq, CL);
		       A->nseq=(CL->S)->nseq;
		     }
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
		       strcpy(out_aln_format[0],"conc_aln");
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
			 {

			   T=make_tree ( A,CL,gop, gep,(CL->S),pc,maximise);
			 }
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
		     CL->moca=( Moca*)vcalloc ( 1, sizeof ( Moca));
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
		    int ii;
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
			  kmir(A);
			  //A=iterate_aln (A, iterate, CL);
			  //A=ungap_aln(A);
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
		      if (!strm2(out_aln, "stdout", "stderr") && le==stderr && !do_convert && A->nseq < MAX_NSEQ_4_DISPLAY)output_format_aln ("aln",A,NULL,"stdout");


		      A->CL=CL;
		      for ( a=0; a< n_out_aln_format; a++)
			if ( !strstr ( out_aln_format[a], "expand"))output_format_aln (out_aln_format[a],A,EA, tot_out_aln[a]);
		      for ( a=0; a< n_out_aln_format; a++)
			if ( strstr (out_aln_format[a], "expand"))output_format_aln (out_aln_format[a],A,EA, tot_out_aln[a]);



		      fprintf (le, "\n\nOUTPUT RESULTS");
		      if ((CL->S)->nseq>2)
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

	      vfree_all(NULL);

	      le=t_coffee_tip (le, tip);
	      le=print_command_line ( le);
		  le=print_mem_usage (le, PROGRAM);
	      //le=print_cpu_usage(le, PROGRAM);
	      le=print_program_information (le, NULL);


	      if (full_log && full_log[0])log_function(full_log);

	      return EXIT_SUCCESS;
	      /**
	       * Done.
	       * \tableofcontents
	       */
	}

/*Specialized set of Parameters*/
char *get_defaults (char *buf, char *type)
{
  return NULL;
}
char *get_precomputed_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf," -msa_mode=precomputed ");
     buf=strcat (buf," -seq_weight=no ");
     buf=strcat (buf," -evaluate_mode no ");
     buf=strcat (buf," -in Xpam250mt ");
     return buf;
   }
char *get_evaluate_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf," -quiet=stdout ");
     /*buf=strcat (buf," -seq_weight=no ");*/
     buf=strcat (buf," -output score_ascii html ");
     buf=strcat (buf," -iterate 0 ");

     buf=strcat (buf," -evaluate ");



     return buf;
   }
char *get_genome_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

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

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-cosmetic_penalty=-50 ");
     buf=strcat (buf,"-distance_matrix_mode=slow ");
     buf=strcat (buf,"-output clustalw,score_ascii ");
     buf=strcat (buf,"-evaluate_mode=non_extended_t_coffee ");
     buf=strcat (buf,"-clean_aln 0 ");

     return buf;
   }

char *get_very_fast_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Xblosum62mt ");
     buf=strcat (buf,"-distance_matrix_mode ktup ");
     buf=strcat (buf,"-maxnseq 10000 ");
     buf=strcat (buf,"-dpa_maxnseq 0 ");
     buf=strcat (buf,"-dp_mode gotoh_pair_wise_lgp ");

     buf=strcat (buf,"-extend_mode matrix ");
     buf=strcat (buf,"-gapopen -25 ");
     buf=strcat (buf,"-gapext -1 ");
     buf=strcat (buf,"-iterate 0 ");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }

char *get_low_memory_defaults(char *buf, char *type)
   {
     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

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
     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Msap_pair  -template_file SELF_P_ -profile_template_file SELF_P_");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_expresso_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-method sap_pair  -template_file EXPRESSO -profile_template_file EXPRESSO");

     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_procoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Mpromo_pair -extend_seq  ");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_blastr_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Mblastr_pair -extend_seq   ");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_psicoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Mproba_pair -template_file BLAST  ");
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
     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));
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
     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     if (NO_METHODS_IN_CL) buf=strcat (buf,"-in Mkalign_msa  Mmuscle_msa Mmafft_msa -multi_core methods_relax_msa");

     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_xcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf, " -master _kmeans_N_25_LONG_ ");
     buf=strcat (buf, " -dp_mode linked_pair_wise");
     buf=strcat (buf, " -tree_mode kmeans ");
     //buf=strcat (buf, " -method blastp_o2a proba_pair");

     //if (NO_METHODS_IN_CL)   buf=strcat (buf,"-method blastp_o2a ");

     return buf;
   }
char *get_mcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));


     if (NO_METHODS_IN_CL)   buf=strcat (buf,"-in Mclustalw2_msa Mt_coffee_msa Mpoa_msa Mmuscle_msa Mmafft_msa Mdialignt_msa Mpcma_msa Mprobcons_msa -multi_core methods_relax_msa  ");
     /*buf=strcat (buf,"-in ");*/
     return buf;
   }
char *get_dmcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mkalign_msa Mt_coffee_msa Mpoa_msa Mmuscle_msa Mmafft_msa Mdialignt_msa Mprobcons_msa Mamap_msa -multi_core methods_relax_msa");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_rcoffee_consan_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
	 if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mconsan_pair -multi_core templates_relax_msa -dp_mode myers_miller_pair_wise -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -type DNA -relax_lib 0");
     else buf=strcat (buf,"-dp_mode myers_miller_pair_wise -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -type DNA -relax_lib 0");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_rmcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
     if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mprobcons_msa Mmafft_msa Mmuscle_msa -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -check_type -type DNA -relax_lib 0");
     else buf=strcat (buf,"-extend_mode rna2 -template_file RCOFFEE -transform dna2rna -check_type -type DNA -relax_lib 0");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }

//    if (NO_METHODS_IN_CL)buf=strcat (buf,"-in Mbest_pair4prot -template_file BLAST -template_file EXPRESSO  ");
   char *get_best4RNA_defaults(char *buf, char *type)
   {

	   if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

	   check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);
	   buf=strcat (buf," -extend_mode rna2 -template_file PDB,RNA -in Mbest_pair4rna -transform dna2rna -relax_lib 0");
	   /*buf=strcat (buf,"-in ");*/

	   return buf;
   }

char *get_rcoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     check_program_is_installed (RNAPLFOLD_4_TCOFFEE,NULL, NULL,RNAPLFOLD_ADDRESS, INSTALL_OR_DIE);

     buf=strcat (buf," -extend_mode rna2 -template_file RCOFFEE -transform dna2rna -check_type -type DNA -relax_lib 0");
     /*buf=strcat (buf,"-in ");*/

     return buf;
   }
char *get_rsapcoffee_defaults(char *buf, char *type)
{ 
  if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));
  check_program_is_installed (SAP_4_TCOFFEE,NULL, NULL,SAP_ADDRESS, INSTALL_OR_DIE);
  check_program_is_installed (MUSTANG_4_TCOFFEE,NULL, NULL,MUSTANG_ADDRESS, INSTALL_OR_DIE);
  
  buf=strcat (buf,"-method sap_pair -template_file=RNA -extend_mode rna2 -output clustalw,html -transform dna2rna");
  
  return buf;
   }
char *get_saracoffee_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));
     
     check_program_is_installed ("sara.py",NULL, NULL,"http://structure.biofold.org/sara",INSTALL_OR_DIE);
     buf=strcat (buf,"-method sara_pair -template_file=RNA -extend_mode rna2 -relax_lib 0 -output clustalw,html -transform dna2rna");
     
     return buf;
   }

char *get_genepredx_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf, "-method tblastx_msa -evaluate_mode sequences  -genepred  -relax_lib 0 -output fasta_seq,exons,texons,wexons -seq_weight no -check_type -type DNA -out_lib");
     return buf;
   }
char *get_genepredpx_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf, "-method tblastpx_msa -evaluate_mode sequences  -genepred  -relax_lib 0 -output fasta_seq,exons,texons,wexons -seq_weight no -check_type -type DNA -out_lib");
     return buf;
   }

char *get_repeat_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in slow_pair -matrix idmat -out_lib -profile_comparison profile -profile_mode channel_profile_profile -dp_mode myers_miller_pair_wise ");
     /*buf=strcat (buf,"-in ");*/
     return buf;
   }

char *get_sample_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-extend_mode matrix -dp_mode gotoh_pair_wise -sample_dp -in Xblosum62mt");
     
     return buf;
   }
char *get_highlow_defaults(char *buf, char *type)
   {

     if (buf==NULL)buf=(char*)vcalloc (1000, sizeof (char));

     buf=strcat (buf,"-in Xblosum62mt -dp_mode gotoh_pair_wise_test");
     
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

  c2=(char*)vcalloc ( 100, sizeof (char));

  sprintf (command, "cat hmgt_mouseVsnrl3d.blast_result |blast_aln2fasta_aln.pl | fasta_aln2fasta_aln_unique_name.pl > my_test");

  fprintf ( stderr, "C1: %d, C2:%d", is_dynamic_memory (c2), is_dynamic_memory (c2));


  myexit (0);
}


/**
 * Use one of the other_pg options.
 *
 * If T-Coffee is called with the \b other_pg option, this function determines what to do.
 *
 * \callgraph
 * \see ::seq_reformat, ::aln_compare
 */
int run_other_pg ( int argc, char *argv[])
{
  //make minimum initialization


  if ( strm (argv[0], "seq_reformat") || strm (argv[0], "saltt"))
    {
      return seq_reformat (argc, argv);
    }
  else if ( strm (argv[0], "unistrap"))
    {
      return unistrap (argc, argv);
    }
  else if ( strm (argv[0], "aln_compare"))
    {
      return aln_compare (argc, argv);
    }
  else if ( strm (argv[0], "analyse_pdb") || strm (argv[0], "apdb") || strm (argv[0], "irmsd") || strm (argv[0], "trmsd") || strm (argv[0], "strike"))
    {
      return apdb ( argc, argv);
    }
  else if ( strm (argv[0], "quantile"))
    {
      return quantile ( argc, argv);
    }
  else if ( strstr ( argv[0], "unpack_"))
    {
      char buf;
      unpack_all_perl_script (argv[0]+strlen ("unpack_"));
    }
  else if ( strstr ( argv[0], "mat2process"))
    {
      return mat2process (argc-1, argv+1);
    }
  else
    {
      return my_system_cl (argc, argv);
    }
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



Sequence* prepare_master (char *seq, Sequence *S, Constraint_list *CL, char *dmode)
 {
   int a,b,s1, n, i;
   FILE *fp;
   int trim_mode=2;
   int **lu=NULL;
   char ttag;
   
   char tmode[100];
   int nseq=0;
   CL->master=(int*)vcalloc (S->nseq+1, sizeof(int));

   if ( check_file_exists (seq))
     {
       Sequence *L;
       L=main_read_seq (seq);
       for (a=0; a< L->nseq; a++)
	 if ( (b=name_is_in_list (L->name[a], S->name,S->nseq, 100))!=-1)CL->master[b]=1;
     }
   else if ( S->nseq==2 || strm (seq, "no") || strm (seq, "default"))
     {
       for (a=0; a<S->nseq; a++)CL->master[a]=1;
       return NULL;
     }
   
   //figure out the number of sequences to keep as a master
   nseq=0;
   if ( is_number (seq) || strstr (seq, "_N") || strstr (seq, "_n"))
     {
       char ttag;
       
       if      ( strstr (seq, "_N_")){nseq=atoi (strstr(seq, "_N_")+strlen ("_N_"));ttag='N';}
       else if ( strstr (seq, "_n_")){nseq=atoi (strstr(seq, "_n_")+strlen ("_n_"));ttag='n';}
       else if ( strstr (seq, "_N")){nseq=atoi (strstr(seq, "_N")+strlen ("_N"));ttag='N';}
       else if ( strstr (seq, "_n")){nseq=atoi (strstr(seq, "_n")+strlen ("_n_"));ttag='n';}
       else {nseq=atoi (seq);ttag='n';}
       
       if ( ttag=='N')nseq=((float)S->nseq*((float)nseq/(float)100.0));  
     }
  
   //if no seqiuences or all, keep everything and return
   if ( nseq>=S->nseq || nseq==0)
     {
       for (a=0; a<(CL->S)->nseq; a++)CL->master[a]=1;
       return NULL;
     }
  
   //If keep tghe nlonguest make a sorted list of sequence indexes
   if ( strstr (seq, "_NLONG_"))
     {
       lu=declare_int ((CL->S)->nseq, 2);
       for (a=0; a< (CL->S)->nseq; a++)
	 {
	   lu[a][0]=a;
	   lu[a][1]=strlen ((CL->S)->seq[a]);
	   sort_int_inv(lu,2,1,0, (CL->S)->nseq-1);
	 }
       for (a=0; a<nseq; a++)
	 {
	   CL->master[lu[a][0]]=1;
	 }
       lu=free_int (lu, -1);
     }
   
   //prepare a special library mode using the most commected sequences only
   if ( strstr (seq, "_PLIB_"))
     {
       set_int_variable ("N_4_PLIB", atoi (strstr(seq, "_PLIB_")+strlen ("_PLIB_")));
     }
   
   //Keep all sequences with a known structure
   if ( strstr (seq, "_P_"))
     {
       for (a=0; a<S->nseq; a++)
	 {
	   if (seq_has_template (S, a, "_P_"))CL->master[a]=1;
	 }
     }
 
   //Keep only the most informative sequences according to simple trim
   //Criteria depends on a pairwise distance estimate provided by dmode
   if (strstr (seq, "_TRIM_"))
     {
       Alignment *A,*SA;
       char **name;
       
       A=(strm (dmode, "msa"))?(very_fast_aln (seq2aln (S, NULL, RM_GAP), 0, NULL)):(seq2aln (S, NULL, RM_GAP));
       
       if (strstr (dmode,"_kmeans_"))sprintf (tmode, "_kmeans_n%d_", nseq);
       else sprintf (tmode, "_aln_%c%d_", ttag,nseq);
       
       fprintf ( CL->local_stderr, "------- Master Mode For Trimming: %s\n", tmode);
       SA=simple_trimseq (A, NULL, tmode, NULL, NULL);
       nseq=SA->nseq;
       name=SA->name;
       for (a=0; a<nseq;a++)
	 {
	   if (nseq==(CL->S)->nseq)CL->master[a]=1;
	   else if ((b=name_is_in_list (name[a], S->name,S->nseq, 100))!=-1)CL->master[b]=1;
	 }
       free_aln (A);
       free_aln (SA);
     }


   //keep the longuest sequence
   if ( strstr (seq, "_LONG_"))
     {
       int ml=0;
       int ls=0;
       for (a=0; a< (CL->S)->nseq; a++)
	 {
	   int l=strlen ((CL->S)->seq[a]);
	   if (l>ml){ml=l;ls=a;}
	 }
       CL->master[ls]=1; //keep the longest seqquence
     }
   
   fprintf ( CL->local_stderr, "\n");
   
   for (b=0,a=0; a<S->nseq; a++)
     {
       if ( CL->master[a])
	 {
	   fprintf (CL->local_stderr, "\tMaster_sequence: %s\n", S->name[a]);
	   b++;
	 }
     }
   fprintf (CL->local_stderr, "------- Selected a total of %d Master Sequences\n",b);
   if ( b==0)
     {
       printf_exit (EXIT_FAILURE, stderr, "ERROR: %s defines an invalid mode for -master [FATAL:%s]\n",seq,PROGRAM);
     }

   if (b!=(CL->S)->nseq)
     {
       Sequence *T, *MS;
       T=duplicate_sequence (CL->S);
       for(a=0; a<T->nseq; a++)
	 {
	   if (!CL->master[a]){vfree (T->seq[a]); T->seq[a]=NULL;}
	 }
       MS=duplicate_sequence (T);
       free_sequence (T, -1);
       MS->blastdbS=CL->S;
       return MS;
     }
   else
     {

       (CL->S)->blastdbS=CL->S;
       return CL->S;
     }
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
  buf=(char*)vcalloc ( 10000, sizeof (char));
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

  printf_system ( "%s",buf);

  r=file2string (file);
  vremove (file);
  return r;
}




/////////////////////////
double max_kmcoffee;
static int aligned;
static int toalign;
static int kmk;
Alignment* km_coffee_align3 (Sequence *S, char *km_tree, int k, char *outfile, int argc, char **argv);
Alignment* km_coffee_align2 (Sequence *S, char *km_tree, int k, int argc, char **argv);
Alignment* km_coffee_align1 (char *method,Alignment *A,int mn,int argc, char **argv, int nit,int round);
Alignment *km_align_profile (Alignment **AL,int n, int argc, char **argv,int nit,int round);
Alignment *km_align_kprofile (Alignment **AL,int k,int n, int argc, char **argv,int nit,int round);
Alignment *km_align_seq_slow (char *method,Alignment *A, int argc, char **argv,int nit,int round);
Alignment *km_align_seq_fast (Alignment *A, int argc, char **argv,int nit,int round);
Alignment *km_refine_msa (Alignment *A, int argc, char **argv, int k);
int km_coffee_count (Alignment *A,int k, int t);
static Fname *F;



char** km_coffee (int argc, char **argv)
{
	char **new_argv;
	int new_argc=0;
	Sequence *S;
	int a;
	int k=0;
	int k_leaf=0;
	int nit=0;
	Alignment *A;
	char *tm=NULL;
	int seqset=0;
	char *method=NULL;
	char *km_mode=NULL;
	char *km_tree=NULL;

	new_argv=(char**)vcalloc (argc+100, sizeof (char*));
	char *seq_f = NULL;
	char *out_f = NULL;
	char *km_init = NULL;
	int n_core=1;
	int gapopen=0;
	int gapext=0;
	for (a=0; a<argc; a++)
	{
		if ( strm (argv[a], "-seq"))
		  {
		    seq_f=argv[++a];
		  }
		
		else if ( strm (argv[a], "-mode"))
		  {++a;}
		else if ( strm (argv[a], "-km_mode"))
		  {km_mode=argv[++a];}
		else if ( strm (argv[a], "-km_k"))
		  {k=atoi (argv[++a]);}
		else if ( strm (argv[a], "-km_tree"))
		  {km_tree=argv[++a];}
		else if ( strm (argv[a], "-km_nit"))
		  {nit=atoi (argv[++a]);}
		else if ( strm (argv[a], "-km_kl"))
		  {k_leaf=atoi (argv[++a]);}
		else if ( strm (argv[a], "-gapopen"))
		  {gapopen=atoi (argv[++a]);
		    new_argv[new_argc++]=argv[a-1];
		    new_argv[new_argc++]=argv[a];
		  }
		else if ( strm (argv[a], "-gapext"))
		  {gapext=atoi (argv[++a]);
		    new_argv[new_argc++]=argv[a-1];
		    new_argv[new_argc++]=argv[a];
		  }
		else if ( strm (argv[a], "-tree_mode"))
		  {tm=argv[++a];}
		else if ( strm (argv[a], "-km_method"))
		  {method=argv[++a];}
		else if ( strm (argv[a], "-km_init"))
		  {km_init=argv[++a];}
		else if ( strm (argv[a], "-n_core"))
		  {
		    n_core =atoi (argv[++a]);
		    new_argv[new_argc++]=argv[a-1];
		    new_argv[new_argc++]=argv[a];
		  }
		else if ( strm (argv[a], "-outfile"))
		  {
		    out_f=argv[++a];
		    new_argv[new_argc++]=argv[a-1];
		    new_argv[new_argc++]=argv[a];
		  }
		else
		  {
		    new_argv[new_argc++]=argv[a];
		  }
		
	}
	
	
	//Set the default values
	if (!out_f){declare_name(out_f);sprintf (out_f, "%s.aln", seq_f);}
	if (!k)k=10;
	
	

	if (seq_f==NULL)
	  {
	    myexit(fprintf_error (stderr, "Error: with -mode=kmcoffee sequences MUST be provided with -seq"));
	  }


	//if (!tm){new_argv[new_argc]=vcalloc(100, sizeof (char)); sprintf ( new_argv[new_argc++], "kmeans");}


	if (strm (km_mode, "km_fast"))
	  {
	    if (method == NULL)
	      strcpy(method,"proba_pair");
	    if (km_init == NULL)
	      strcpy(km_init, "distributed");
	    if (k_leaf==0)
	      k_leaf=k;
	    km_coffee_align3(seq_f, k, k_leaf, method, out_f, n_core, gapopen, gapext, km_init );
	  }
	else
	  {
	    S=main_read_seq (seq_f);
	    F=parse_fname (seq_f);
	    if (!k)k=100;
	    if (S->nseq<=k)k=S->nseq;
	    
	    
	    if (!km_mode || strm (km_mode, "topdown"))
	      {
		A=seq2aln(S,NULL, RM_GAP);
		toalign=km_coffee_count (A,k, toalign);
		km_coffee_align1 (method, A, k, new_argc, new_argv,nit,0);
	      }
	    else if (strm (km_mode, "bottomup"))
	      {
		km_coffee_align2 (S,km_tree,k, new_argc,new_argv);
	      }
	    else if (strm (km_mode, "updown"))
	      {
		km_coffee_align3 (S,km_tree,k, out_f, new_argc,new_argv);
	      }
	    
	    else
	      myexit(fprintf_error (stderr,"Please specify km_mode (topdown/bottomup/km_fast/updown)!\n"));
	  }
	myexit (EXIT_SUCCESS);
}

Alignment* kmir (Alignment *A)
{
	int **prf;
	char *seq;
	int a;

	seq=(char*)vcalloc ( A->len_aln+1, sizeof (char));

	prf=aln2prf(A, "blosum62mt");

	for (a=0; a<A->nseq; a++)
	  {
	    ungap (A->seq_al[a]);

	    add_sequence2prf(A->seq_al[a],seq, prf, A->len_aln, -10, -1);
	    sprintf ( A->seq_al[a], "%s", seq);
	  }
	return A;
}

int km_coffee_count (Alignment *A,int k, int t)
{


  if (A->nseq<=k)t++;
  else
    {
      Alignment**AL;
      int n=0;
      int a;

      AL=seq2kmeans_subset(A,k,&n, "triaa");
      if (n==1)return t;
      else
	{
	  for (a=0; a<n; a++)
	    {
	      t=km_coffee_count (AL[a],k,t);
	    }
	  t++;
	  for (a=0; a<n; a++)free_aln (AL[a]); vfree(AL);
	}
    }
  return t;
}
Alignment* km_coffee_align1 (char *method, Alignment *A,int k,int argc, char **argv, int nit,int round)
{
  int tot=0;



  if (A->nseq<=k)
    {
      A=km_align_seq_slow (method,A, argc, argv,nit,round);
      aligned++;
    }
  else
    {
      Alignment**AL;
      int n,a;

      AL=seq2kmeans_subset(A,k,&n, "triaa");
      free_aln (A);

      if (n==1)
	{

	  AL=seq2id_subset (AL[0],k,&n, "100");
	  A=km_align_kprofile(AL,n,k,argc, argv,nit,round+1);
	  aligned++;
	}
      else
	{
	  for (a=0; a<n; a++)
	    {
	      AL[a]=km_coffee_align1 (method,AL[a],k,argc, argv,nit,round+1);
	    }

	  A=km_align_profile(AL,n,argc, argv, nit,round);
	  aligned++;
	  for (a=0; a<n; a++)free_aln (AL[a]); vfree(AL);
	}
    }


  fprintf ( stderr, "\r\t\tkmcoffee: Aligned: %4.2f %%", (float)aligned*100/(float)toalign);
  return A;
}


Alignment *km_align_kprofile (Alignment **AL,int n, int k, int argc, char **argv,int nit, int round)
{
	Alignment **NAL;
	Alignment *A;
	int nn=0;
	int a;
	NAL=(Alignment **)vcalloc (n, sizeof (Alignment *));

	if (n==1)return km_align_profile (AL,n,argc, argv, nit, round);
	else
	{
		while (n>1)
		{
			nn=0;
			for (a=0; a<n; a+=k)
			{
				int n2;
				n2=((a+k)>n)?(n-a):k;
				NAL[nn++]=km_align_profile (AL+a,n2,argc, argv, nit, round);
			}
			for (a=0; a<n; a++)
			{free_aln (AL[a]);AL[a]=NULL;}
			for (a=0; a<nn; a++)
				AL[a]=NAL[a];
			n=nn;
		}
		A=AL[0];
	}

	vfree (NAL);
	vfree (AL);
	return A;
}


Alignment *km_align_profile (Alignment **AL,int n, int argc, char **argv,int nit, int round)
{
	char hostname[150];
	gethostname(&hostname[0], 149);
	char *cl=NULL;
	static char *prf;
	char *aln;
	char *prff;
	int a;
	char buf [10000];
	FILE *fp;
	Alignment *A;
	static char *treefile;
	int tot_nseq=0;


	for (a=0; a<n; a++)
		tot_nseq+=(AL[a])->nseq;


	if (!treefile)treefile=vtmpnam(NULL);
	if (!prf)prf=(char*)vcalloc (500, sizeof (char));
	if (n==1 && round!=-1)
	{
		A=copy_aln (AL[0], NULL);
		return A;
	}

	buf[0]='\0';
	aln=vtmpnam(NULL);
	prff=vtmpnam(NULL);
	cl=list2string (argv, argc);
	strcat (buf, cl);
	strcat (buf, " ");



	strcat (buf, " -profile FILE::");
	strcat (buf, prff);

	if (round!=0)
	{
		strcat (buf, " -outfile ");
		strcat (buf, aln);
	}
	strcat (buf, " -output fasta_aln -quiet  ");
	strcat (buf, " -newtree ");
	strcat (buf, treefile);
	strcat (buf, " >/dev/null 2>/dev/null");

	int pid = getpid();
	fp=vfopen (prff, "w");
	for ( a=0; a<n; a++)
	{
		sprintf (prf, "%s/prf_%d_%d", get_tmp_4_tcoffee(), a+1,pid);
		fprintf (fp, "%s\n", prf);
		output_fasta_aln (prf, AL[a]);
	}

	vfclose (fp);

	printf_system_direct_check ("%s",buf);

	for ( a=0; a<n; a++)
	{
		sprintf (prf, "%s/prf_%d_%d", get_tmp_4_tcoffee(), a+1,pid);
		remove (prf);
	}

	A=main_read_aln (aln, NULL);

	if (round==-1)
	  {
	    sprintf (prf, "%s.%d",F->name,getpid());
	    output_fasta_aln (prf,A);
	    printf_system_direct_check ("%s -profile %s", cl, prf);

	    remove (prf);
	    myexit(EXIT_SUCCESS);
	  }
	vfree (cl);
	return A;
}
Alignment *km_align_seq_slow_withTC (Alignment *A, int argc, char **argv, int nit,int round);
Alignment *km_align_seq_slow (char *method, Alignment *A, int argc, char **argv, int nit,int round)
{
  if (A->nseq==1)return A;
  else if (!method || strm (method, "tcoffee"))return km_align_seq_slow_withTC (A, argc, argv, nit, round);
  else
    {
      static char *seq;
      static char *aln;
      if (!seq)seq=vtmpnam (NULL);
      if (!aln)aln=vtmpnam (NULL);
      output_fasta_seq (seq,A);
      free_aln (A);

      printf_system ("tc_generic_method.pl -mode=seq_msa -infile=%s -method=%s -outfile=%s", seq, method, aln);

      if (round==-1)
	{
	  char *cl;
	  char buf[10000];
	  cl=list2string (argv, argc);
	  sprintf ( buf, "%s -profile %s ", cl, aln);
	  myexit (EXIT_SUCCESS);
	}
      else
	{
	  return main_read_aln (aln, NULL);
	}
    }
}
Alignment *km_align_seq_slow_withTC (Alignment *A, int argc, char **argv, int nit,int round)
{

	static char *seq;
	static char *aln;
	char *cl=NULL;
	char buf[100000];
	static char *treefile;
	if (!treefile)treefile=vtmpnam(NULL);

	if (!seq)
	{
		seq=vtmpnam(NULL);
		aln=vtmpnam(NULL);
	}


	if ( A->nseq <=1)
		return A;
	else
	{
		cl=list2string (argv, argc);
		output_fasta_seq (seq,A);
		sprintf ( buf, "%s -seq %s ", cl, seq);

		if (round!=-1)
		  {
		    strcat (buf, " -output fasta_aln -quiet -outfile ");
		    strcat (buf, aln);
		    strcat (buf, " -newtree ");
		    strcat (buf, treefile);
		    strcat (buf, " >/dev/null 2>/dev/null");
		  }

		printf_system_direct_check ("%s",buf);

		if (round==-1)
		  myexit (EXIT_SUCCESS);

		A=main_read_aln (aln,A);
		vfree (cl);
	}
	return A;
}


Alignment *km_align_seq_fast (Alignment *A, int argc, char **argv, int nit, int round)
{
	static Constraint_list *CL;
	int *ns;
	int **ls;
	int n;
	int a;

	if (!CL)
	  {
	    CL=(Constraint_list*)vcalloc ( 1, sizeof (Constraint_list));
	    CL->pw_parameters_set=1;
	    CL->matrices_list=declare_char (10, 10);
	    CL->M=read_matrice ("blosum62mt");
	    CL->evaluate_residue_pair=evaluate_matrix_score;
	    CL->get_dp_cost=get_dp_cost;
	    CL->extend_jit=0;
	    CL->maximise=1;
	    CL->gop=-10;
	    CL->gep=-2;
	    CL->TG_MODE=2;
	    sprintf (CL->matrix_for_aa_group, "vasiliky");
	    CL->use_fragments=0;
	    CL->ktup=3;
	    if ( !CL->use_fragments)CL->diagonal_threshold=0;
	    else CL->diagonal_threshold=6;
	    sprintf (CL->dp_mode, "myers_miller_pair_wise");
	    CL->S=declare_sequence (1,1,1);
	}

	ungap (A->seq_al[0]);
	n=A->nseq;

	if (n<=1)
	  return A;

	CL->S=fill_sequence_struc(n,A->seq_al, A->name, NULL);
	ns=(int*)vcalloc (2, sizeof (int));
	ls=declare_int (2, n);
	ns[0]=ns[1]=1;
	ls[0][0]=0;

	//Alignment
	A->len_aln=strlen (A->seq_al[0]);


	for (a=1; a< n; a++)
	{
		ungap(A->seq_al[a]);
		ls[1][0]=a;
		pair_wise (A,ns,ls,CL);
		A->len_aln=strlen (A->seq_al[a]);
		ls[0][ns[0]++]=a;
		A->nseq=ns[0];
	}
	free_int (ls,-1);
	vfree (ns);

	return A;
}




Alignment *km_refine_msa (Alignment *A,int argc, char **argv, int k)
{
  char *seq;
  char *prf;
  char *tree;
  char *aln;
  char buf [10000];

  prf=vtmpnam(NULL);
  seq=vtmpnam(NULL);
  aln=vtmpnam(NULL);
  tree=vtmpnam (NULL);

  if (A->nseq<=k)k=A->nseq/2;

  aln2gap_trimmed (A,k, prf, seq);
  printf_system_direct_check ("t_coffee -seq %s -profile %s -outfile %s -newtree %s ",seq, prf, aln, tree);
  free_aln (A);
  A=main_read_aln (aln, NULL);
  remove (prf);
  remove (seq);
  remove (tree);

  return A;
}

Alignment * km_coffee_align2 (Sequence *S, char *km_tree, int k, int argc, char **argv)
{
	char *km_tree2=vtmpnam (NULL);
	NT_node T;
	if (strm (km_tree, "kmeans")){
		seq2km_tree_old (S, km_tree2);}
	else if (!km_tree)
	{
		seq2km_tree_old (S, km_tree2);
	}
	T=main_read_tree (km_tree2);
	tree_aln_N(T,S, k, argc, argv);
	
	myexit (EXIT_SUCCESS);
	return NULL;
}

Alignment * km_coffee_align3 (Sequence *S, char *km_tree, int k, char *out_f, int argc, char **argv)
{

  
  char *km_tree2=vtmpnam (NULL);
  NT_node T;
  int n;
  
  //This insures that the function aln2cons_cov is used to generate the consensus
  //cputenv ("KM_COFFEE_CONS_COV=1");
  
  //This insures that  get_tot_prob estimates prf/prf alignments using the cons and not the voectorized MSA
  //cputenv ("KM_COFFEE_PRF_CONS=1");
  
  //Split profiles can be triguered by adding to any method:-method <method>@EP@PRFMODE@prf1-2 or3
  //@EP@PRFMODE@prf1: pick the N sequences covering entirely the prf
  //@EP@PRFMODE@prf2: make a consensus sequence patching the most covering sequences
  //@EP@PRFMODE@prf3: make a consensus sequence using blosum62mt
  

  if (strm (km_tree, "kmeans")){
    seq2km_tree_old (S, km_tree2);}
  else if (!km_tree)
    {
      seq2km_tree_old (S, km_tree2);
    }
  else if (check_file_exists (km_tree))
    {
      printf_system ("cp %s %s", km_tree, km_tree2);
    }
  	
  T=main_read_tree (km_tree2);
  updown_tree_aln (T,S, k,&n, argc, argv);
  
  
  if (!check_file_exists (T->alfile))
    {
      printf_exit ( EXIT_FAILURE, stderr, "kmcoffee did not manage to align the provided dataset\n");
    }
  
  printf_system ("mv %s %s", T->alfile, out_f);
  display_output_filename(stdout,"MSA","ALN",out_f, CHECK);
  fprintf (stdout, "\n\n");
  myexit (EXIT_SUCCESS);
  return NULL;
}

Alignment * t_coffee_dpa (int argc, char **argv)
{
  NT_node T;
  FILE *fp;
  Sequence *S=NULL;
  char *seqfile=NULL;
  char *dpa_tree=NULL;
  char *usetree=NULL;
  char *newtree=NULL;
  char *outfile=NULL;
  char *alnfile=NULL;
  char *dpa_weight=NULL;
  int dpa_nseq=0;
  char *dpa_aligner=NULL;
  char *command;
  char *run_name=NULL;
  Fname *F=NULL;
  char *cache=NULL;
  float *w;
  int a;
  int output_dpa_tree=1;
  int seqflag=0;
  FILE *le;
  char *se_name;
  char *homoplasy=NULL;
  int reg_dynamic=1;
  int reg_pool=0;
  int n_core=1;
  
  /* This is used for the dump function see -dump option*/
  declare_name (se_name);
  sprintf (se_name, "stderr");
  le=get_stdout1(se_name);
  
  
  command=(char*)vcalloc (10000, sizeof (char));
  sprintf ( command, "#");
  
  //default values
  set_int_variable ("swlN",50);
  set_string_variable ("output", "fasta_aln");
  set_int_variable    ("reg_dnd_nseq",100);
  set_int_variable    ("reg_dnd_depth",3);
  set_string_variable ("reg_dnd_mode", "codnd");
  //*_4_CLTCOFFEE gets added to the T-Coffee command line of any slave;

  //Set Some Defaults that can be over-written by the CL parsing below
  //Note that by default ALL flag get added to the env variable <name>_4_CLTCOFFEE. These variables are then used to construct the slave CL in dynamic.pl
  cputenv ("blast_server_4_CLTCOFFEE=LOCAL");
 
  for (a=1; a<argc; a++)
    {
      
      if ( argv[a][0]!='-')
	{
	  myexit (fprintf_error (stderr, "%s is an unknown flag of the -reg mode [FATAL:%s]", argv[a],PROGRAM));
	}
      
      if (strm (argv[a], "-seq" ))
	{
	  
	  seqfile=argv[++a];
	  S=quick_read_seq (seqfile);
	  seqflag=1;
	}
      else if (strm ( argv[a], "-output"))
	{
	  char *dpa_output=argv[++a];
	  if (strm (dpa_output, "fastaz_aln"))set_string_variable ("output", "fastaz_aln");
	  else if (strm (dpa_output, "fasta_aln"))set_string_variable ("output", "fasta_aln");
	  else
	    {
	       myexit (fprintf_error (stderr, "regressive mode only supports fasta_aln and fastaz_aln with the -output flag [FATAL:%s]", argv[a],PROGRAM));
	    }
	}
      else if (strm (argv[a],"-tree") || strm (argv[a],"-dpa_tree") ||strm (argv[a],"-reg_tree") )
	{
	  dpa_tree=argv[++a];
	}
      else if (strm (argv[a],"-child_tree"))
	{
	  cputenv ("child_tree_4_TCOFFEE=%s", argv[++a]);	  
	}
      else if (strm (argv[a],"-child_thread"))
	{
	  cputenv ("thread_4_TCOFFEE=%d", atoi(argv[++a]));
	}
      else if (strm (argv[a],"-dynamic_config"))
	{
	  cputenv ("dynamic_config_4_TCOFFEE=%s", (fname2abs(argv[a+1])));
	  a++;
	}
      
      else if ( strm(argv[a], "-reg_chaindnd_mode"))
	{
	  set_string_variable ("reg_chaindnd_mode",argv[++a]);
	}
      else if ( strm(argv[a], "-reg_dnd_nseq"))
	{
	  set_int_variable ("reg_dnd_nseq",atoi(argv[++a]));
	}
      else if ( strm (argv[a],"-reg_dnd_depth"))
	{
	  set_int_variable ("reg_dnd_depth",atoi(argv[++a]));
	}
      else if ( strm (argv[a],"-reg_dnd_mode"))
	{
	  set_string_variable ("reg_dnd_mode",argv[++a]);
	}
	
      else if (strm (argv[a],"-dpa_swlN")|| strm (argv[a],"-reg_swlN") )
	{
	  set_int_variable ("swlN",atoi (argv[++a]));
	}
      else if (strm (argv[a],"-dpa_weight") || strm (argv[a],"-reg_weight"))
	{
	  dpa_weight=argv[++a];
	}

      else if (strm (argv[a],"-thread") || strm (argv[a],"-dpa_thread") || strm (argv[a],"-dpa_n_core") || strm (argv[a],"-reg_thread") || strm (argv[a],"-reg_n_core"))
	{
	  n_core=atoi(argv[++a]);
	}
       
      else if (strm (argv[a],"-dpa")|| strm (argv[a],"-reg"));
      
      else if (strm (argv[a],"-usetree"))
	{
	  usetree=argv[++a];
	}
      else if (strm (argv[a],"-newtree") ||strm (argv[a],"-outtree") )
	{
	  newtree=argv[++a];
	}
      else if (strm (argv[a],"-run_name"))
	{
	  run_name=argv[++a];
	}
      else if (strm (argv[a],"-outfile"))
	{
	  outfile=argv[++a];
	}
      else if (strm (argv[a], "-nseq") || strm (argv[a],"-dpa_nseq") || strm (argv[a], "-reg_nseq")|| strm (argv[a], "-N") )
	{
	  dpa_nseq=atoi(argv[++a]);
	}
      else if (strm (argv[a],"-dynamic") || strm (argv[a],"-reg_dynamic") )
	{
	  reg_dynamic=atoi(argv[++a]);
	  
	}
      else if (strm (argv[a],"-pool")  )
	{
	  reg_pool=1;
	  
	}
      else if (strm (argv[a], "-method") || strm (argv[a], "-dpa_method") || strm (argv[a], "-reg_method"))
	{
	  dpa_aligner=argv[++a];
	}
      else if (strm (argv[a], "-cache"))
	{
	  cache=argv[++a];
	}
      else if (strm (argv[a],"-in") || strm (argv[a],"-infile"))
	{
	  myexit (fprintf_error (stderr, "%s is not supported when using -dpa [FATAL:%s]", argv[a],PROGRAM));
	}
     
      else if ( strstr (argv[a], "reg_homoplasy"))
	{
	  homoplasy=(char*)vcalloc ( 1000, sizeof (char));
	 
	}
      else if (argv[a][0]=='-' &&(a==argc-1 || argv[a+1][0]=='-'))
	{
	  cputenv ("%s_4_CLTCOFFEE=FLAGSET", argv[a]+1);
	}
      else if (strm (argv[a], "-child_method"))
	{
	  char *ml=NULL;
	  a++;
	  while (a<argc && argv[a][0]!='-')
	    {
	      ml=vcat(ml," ");
	      ml=vcat(ml, argv[a]);
	      a++;
	    }
	  if (ml)cputenv ("method_4_CLTCOFFEE=%s", ml);
	}
      else if (argv[a][0]=='-')       
	{
	  if (getenv("DEBUG_4_TCOFFEE"))
	    {
	      HERE ("Parse Arg: %s -> %s -> %s", argv[a], argv[a]+1,fname2abs(argv[a+1]));
	    }
	  cputenv ("%s_4_CLTCOFFEE=%s", argv[a]+1,fname2abs(argv[a+1]));
	  a++;
	  
	}
      else 
	{
	   myexit (fprintf_error (stderr, "%s is an unknown flag of the -reg mode [FATAL:%s]", argv[a],PROGRAM));
	}
    }

  if (getenv("DEBUG_4_TCOFFEE"))
    {
      HERE ("******* ARGUMENTS*****: start");
      for (a=1; a<argc; a++)fprintf ( stderr, "%s\n", argv[a]);
      HERE ("******* ARGUMENTS*****: end");
      HERE ("******* Environement: start*****");
      system ("printenv");
      HERE ("******* Environement: end*****");
    }
      

  //Core management
 
  if (n_core==0)n_core=get_nproc();
  set_nproc (n_core);
  set_int_variable ("n_core",n_core);
  

  //Dynamic regression
  set_int_variable ("reg_dynamic",reg_dynamic);
  set_int_variable ("reg_pool",reg_pool);
  
  //Set the cache
  if (!cache)cache=csprintf (NULL,"use");
  prepare_cache (cache);
  cputenv ("cache_4_TCOFFEE=%s", get_cache_4_tcoffee());
  cputenv ("cache_4_CLTCOFFEE=%s", get_cache_4_tcoffee());
  
  //prepare the aligner CL
  if (dpa_aligner)
    {
      
      command=(char*)vcalloc (10000, sizeof (char));
      sprintf ( command, "%s", dpa_aligner);
    }
  else
    {
      command=(char*)vcalloc (10000, sizeof (char));
      sprintf (command, "clustalo_msa");
    }
  

  //prepare output names
  //output the MSA
  F=parse_fname(seqfile);
  if (run_name){vfree(F->name); F->name=run_name;F->path[0]='\0';}
  
  if (!outfile)
    {
      outfile=(char*)vcalloc ( 1000, sizeof (char));
      sprintf (outfile, "%s.%s", F->name, get_string_variable("output"));
    }
  //output the MSA
  if (homoplasy)
    {
      sprintf (homoplasy, "%s.homoplasy", F->name);
      set_string_variable("homoplasy", homoplasy);
    }
  
  
  //check Sequences are here
  
  if (!seqflag)
    myexit (fprintf_error ( stderr, "\nERROR: When using -reg, sequences must be provided via -seq [FATAL:%s]", PROGRAM));
  else if (!S)
    myexit (fprintf_error ( stderr, "\nERROR: Could not read %s [FATAL:%s]", seqfile,PROGRAM));
  //decide on bucket sizes
  if (dpa_nseq==0)
    {
      if (S->nseq>10000)dpa_nseq=1000;
      else dpa_nseq=MAX((S->nseq/10),2);
    }
  fprintf ( stdout, "PROGRAM: %s %s (%s) -- regressive mode\n",PROGRAM,VERSION,BUILD_INFO);
  fprintf ( le, "!Maximum N Threads --- %d\n",get_nproc());
  //prepare the guide tree
  fprintf ( le, "!Compute Guide Tree --- ");
  
  if (usetree){dpa_tree=usetree;}
  if (!dpa_tree)dpa_tree="codnd";
  if (strm (dpa_tree, "dpa"))
    {
      
      T=seq2dnd (S, "blength");
      w=seq2dpa_weight (S, "longuest");
      T=node2master (T, S, w);
      T=tree2dnd4dpa(T, S, dpa_nseq, command);
      
    }
  else
    T=seq2dnd (S, dpa_tree);
  fprintf ( le, " reg_tree %s\n", dpa_tree);
  
  if (dpa_tree && check_file_exists (dpa_tree))output_dpa_tree=0;
  else output_dpa_tree=1;
  
  //Save Guide Tree for Children Aligners
  
  fprintf ( le, "!Compute Guide Tree ---  done\n");
  
 
  
  
  //get the weight
  fprintf (le, "!Compute Weights --- ");
  if (dpa_weight){w=seq2dpa_weight (S, dpa_weight); fprintf ( le, "%s\n", dpa_weight);}
  else { w=seq2dpa_weight (S, "longuest");fprintf ( le, "default longuest\n", dpa_weight);}
  fprintf ( le, "!Compute Weights --- done\n");
  
  //run the alignment
  fprintf (le, "!Compute MSA --- reg_method %s -- reg_nseq %d -- start\n", command, dpa_nseq);
  T=node2master (T, S, w);
  
  alnfile=tree2msa4dpa(T, S, dpa_nseq, command);
  fprintf ( le, "!Compute MSA --- done\n");
 
  printf_system ("mv %s %s", alnfile, outfile);
  display_output_filename (le, "MSA",get_string_variable ("output"),outfile, CHECK);
  if (homoplasy)display_output_filename (le, "HOMOPLASY","homoplasy",homoplasy, CHECK);
  
  //output The tree
  if (output_dpa_tree)
    {
      if (!newtree)
	{
	  newtree=(char*)vcalloc ( 1000, sizeof (char));
	  sprintf (newtree, "%s.%s", F->name, dpa_tree);
	}
      print_newick_tree (T,newtree);
      display_output_filename (le, "TREE","newick",newtree, CHECK);
    }
  
  //terminate
  fprintf (le,"\n\n");
  print_command_line (le);
  print_mem_usage (le, "REG memory Usage");
  print_program_information (le, NULL);
  fprintf ( le, "\n");
  myexit (EXIT_SUCCESS);
  return NULL;
}

//

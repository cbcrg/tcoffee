#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"

//-----------MARIA ADDED CHANGES------------//

FILE * compare_al_to_lib ( Constraint_list *CL, int start, char *fname, Sequence *S)
{
    
	int a, b, s1, s2, r1, r2;
	int len;
	static int* translation;
	FILE *fp;


	if ( translation!=NULL)vfree(translation);
	  translation=(int*)vcalloc ( (CL->S)->nseq+1, sizeof (int));
	    for ( b=0,a=0; a< (CL->S)->nseq; a++)
		{
			if ( name_is_in_list((CL->S)->name[a],S->name,S->nseq, 100)==-1)
			{
				(CL->S)->len[a]=-1;
				translation [a]=-1;
			}
			else
			{
				translation[a]=b++;
			}
		}

	
	fp=vfopen ( fname, "w");
	fp=save_list_header (fp,CL);
		
	len=CL->ne;
	
	if (len==start && CL->cpu!=-1)
	{
		fprintf (fp, "! CPU %d\n",get_time());
		return fp;
	}
	else
	{
		Sequence *S=CL->S;
		int ***cacheI=(int***)vcalloc (S->nseq, sizeof (int**));
		int **cacheR=(int**)vcalloc (S->nseq, sizeof (int*));
		int *tot=(int*)vcalloc ( S->nseq, sizeof (int));
		int *max=(int*)vcalloc ( S->nseq, sizeof (int));
		for (a=0;a<S->nseq; a++)
		  {
		    cacheI[a]=(int**)vcalloc (S->len[a], sizeof (int*));
		    cacheR[a]=(int*)vcalloc (S->len[a], sizeof (int));
		  }
		for (a=0; a<S->nseq; a++)max[a]=S->len[a];

		for (s1=0; s1<S->nseq; s1++)
		  {
		  for (r1=1; r1<=S->len[s1]; r1++)
		    {
		      for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{

			  s2=CL->residue_index[s1][r1][b+SEQ2];
			  if (tot[s2]>=max[s2])
			    {
			      max[s2]+=100;
			      cacheI[s2]=(int**)vrealloc (cacheI[s2], max[s2]*sizeof (int*));
			      cacheR[s2]=(int*)vrealloc (cacheR[s2], max[s2]*sizeof (int*));
			    }

			  cacheI[s2][tot[s2]]=CL->residue_index[s1][r1]+b;
			  cacheR[s2][tot[s2]]=r1;
			  tot[s2]++;
			}
		    }

			for (s2=0;s2<S->nseq; s2++)
			{
				int x1=translation[s1];
				int x2=translation[s2];
				if (tot[s2] && x1!=-1 && x2!=-1 && x1<x2)
				{
					fprintf ( fp, "#%d %d\n", x1+1, x2+1);
					for ( r2=0; r2<tot[s2];r2++)
					{
						int *v=cacheI[s2][r2];
						fprintf (fp, "%5d %5d %5d %5d %5d\n",cacheR[s2][r2], v[R2], v[WE], v[CONS], v[MISC]);
					}
				}
				tot[s2]=0;
			}
		}
		for (a=0; a<S->nseq; a++)
		{
			vfree (cacheI[a]);
			vfree (cacheR[a]);
		}
		vfree (cacheI); vfree (cacheR);
		vfree (tot);
		vfree (max);
	}

	return save_list_footer (fp, CL);

}
//-----------END OF CHANGES------------//
 
 
 
int aln_compare ( int argc, char *argv[])
    {
    int a, b, c, f;
	
    Alignment *A, *B;
    Sequence *SA, *SB, *TOT_SEQ=NULL;
    Sequence *defined_residueA;
    Sequence *defined_residueB;
    char **seq_list;
    int n_seq_file;
    
    
    Sequence  *S=NULL;
    Structure *ST=NULL;
    Result *R=NULL;
    char *buf1;
    char *buf2;

/*PARAMETERS*/
    char ***grep_list;
    int n_greps;
    char compare_mode[STRING];
    char sim_aln[STRING];
    char *alignment1_file;
    char *alignment2_file;
    //---Maria added this---
    char *lib_file; //---
    
    char *pep1_file;
    char *pep2_file;
    
    char io_format[STRING];

    int n_structure;
    char **struct_file;
    char **struct_format;
    int *n_symbol;
    char ***symbol_list;
    int pep_compare;
    int aln_compare;
    int count;
    int print_pair;
    int output_aln;
    int output_aln_threshold;
    char output_aln_file[LONG_STRING];
    char output_aln_format[LONG_STRING];
    char output_aln_modif[LONG_STRING];
/*LIST VARIABLES*/
    Constraint_list *CL_A;
    Constraint_list *CL_B;
    //---Maria added this---//
    Constraint_list *CL;

    int pos_in_clist;
    Sequence *CLS;


/*Column Comparison Variables*/
 int **posA;
 int **posB;
 int **seq_cache;
 int is_same;
 int *correct_column;
/*RESULTS_VARIABLES*/
    int **tot_count;
    int **pos_count;
    int ***pw_tot_count;
    int ***pw_pos_count;
    int *glob;
    int **pw_glob;
    int s1, r1, s2, r2;
/*IO VARIABLES*/
    int n_categories;
    char ***category;
    char *category_list;
    int *n_sub_categories;
    char sep_l;
    char sep_r;
    char io_file[STRING];
    FILE *fp;
    int **aln_output_count;
    int **aln_output_tot;
    
/*Sims VARIABLES*/
    float **sim;
    float **sim_param;
    char sim_matrix[STRING];
    
    int sim_n_categories;
    char ***sim_category;
    char *sim_category_list;
    int *sim_n_sub_categories;
    
    
    if ( argc==1|| strm6 ( argv[1], "h", "-h", "help", "-help", "-man", "?"))
      {
	output_informations();
      }
    
    argv=standard_initialisation (argv, &argc);
    /*Declarations and Initializations*/
    alignment1_file=(char*)vcalloc ( LONG_STRING, sizeof (char));
    alignment2_file=(char*)vcalloc ( LONG_STRING, sizeof (char));
    //---Maria added this---//
    lib_file=(char*)vcalloc ( LONG_STRING, sizeof (char));
    
    pep1_file=(char*)vcalloc ( LONG_STRING, sizeof (char));
    pep2_file=(char*)vcalloc ( LONG_STRING, sizeof (char));
    

    sprintf (compare_mode, "sp");
    count=0;
    pep_compare=0;
    aln_compare=1;

    
    grep_list=(char***)vcalloc ( STRING, sizeof (char**));
    for ( a=0; a< STRING; a++)grep_list[a]=declare_char (3, STRING);
    n_greps=0;

    n_structure=0;
    struct_file=declare_char ( MAX_N_STRUC, LONG_STRING);
    struct_format=declare_char (MAX_N_STRUC, STRING);
    
    n_symbol=(int*)vcalloc ( MAX_N_STRUC, sizeof (int));
    symbol_list=(char***)vcalloc (MAX_N_STRUC, sizeof (char**));
    for ( a=0; a< MAX_N_STRUC; a++)symbol_list[a]=declare_char ( 100, 100);
    

   
    
    n_categories=1;
    category=(char***)vcalloc ( MAX_N_CATEGORIES, sizeof (char**));
    for ( a=0; a< MAX_N_CATEGORIES; a++)category[a]=declare_char(100, STRING);
    n_sub_categories=(int*)vcalloc ( 100, sizeof (int));
    category_list=(char*)vcalloc ( LONG_STRING, sizeof (char));
    sprintf ( category_list, "[*][*]=[ALL]");
    

    sim_n_categories=1;
    sim_category=(char***)vcalloc ( MAX_N_CATEGORIES, sizeof (char**));
    for ( a=0; a< MAX_N_CATEGORIES; a++)sim_category[a]=declare_char(100, STRING);
    sim_n_sub_categories=(int*)vcalloc ( 100, sizeof (int));
    sim_category_list=(char*)vcalloc ( LONG_STRING, sizeof (char));
    sprintf ( sim_category_list, "[*][*]=[ALL]");
    sprintf ( sim_aln, "al1");
    sim_matrix[0]='\0';
    sprintf ( sim_category_list, "[*][*]=[ALL]");
    
    sprintf ( io_format, "ht");
    sprintf ( io_file, "stdout");
    sep_l='[';
    sep_r=']';

    output_aln=0;
    print_pair=0;
    output_aln_threshold=100;
    sprintf ( output_aln_file, "stdout");
    sprintf ( output_aln_format, "clustalw");
    sprintf ( output_aln_modif, "lower");
/*END OF INITIALIZATION*/
    
   

    
/*PARAMETERS INPUT*/

    for ( a=1; a< argc; a++)
    	{
	if      (strcmp ( argv[a], "-f")==0)
	        {
		sprintf  (io_file,"%s", argv[++a]);
		}
        else if ( strcmp ( argv[a], "-sim_aln")==0)
    		{
    		sprintf  (sim_aln,"%s", argv[++a]);
    		}
    	else if ( strcmp ( argv[a], "-sim_matrix")==0)
    		{
    		sprintf  (sim_matrix,"%s", argv[++a]);
    		}
	else if ( strm ( argv[a], "-compare_mode"))
	        {
		sprintf ( compare_mode, "%s", argv[++a]);
		}
	else if ( strcmp ( argv[a], "-sim_cat")==0)
		{	
		if ( argv[++a][0]!='[')
		   {
		   if ( strcmp ( argv[a], "3d_ali")==0)sprintf ( sim_category_list, "[b][b]+[h][h]=[struc]");
		   else
		       {
		       fprintf ( stderr, "\n%s: Unknown category for distance measure", argv[a]);
		       } 
		   }
		else
		   {
		   sprintf ( sim_category_list, "%s", argv[a]);
		   }
    		}
    	else if ( strcmp ( argv[a], "-grep_value")==0)
    		{
    		sprintf ( grep_list[n_greps][0], "%s", argv[++a]);
		sprintf ( grep_list[n_greps][1], "%s", argv[++a]);
    		n_greps++;
	     
    		}
	//----- Maria added this----
	else if ( strcmp ( argv[a], "-lib")==0)
    		{
    		sprintf ( lib_file, "%s", argv[++a]);
		
		}	
	else if ( strcmp ( argv[a], "-al")==0)
		{
    		
    		sprintf ( alignment2_file, "%s", argv[++a]);
		
		}	
	//-----END of changes-----
	else if ( strcmp ( argv[a], "-al1")==0)
    		{
    		
    		sprintf ( alignment1_file, "%s", argv[++a]);
		
		}	
	else if ( strcmp ( argv[a], "-al2")==0)
		{
    		
    		sprintf ( alignment2_file, "%s", argv[++a]);
		
		}
	else if (  strcmp ( argv[a], "-pep1")==0)
	        {
		pep_compare=1;
		sprintf ( pep1_file, "%s", argv[++a]);
		}
	else if (  strcmp ( argv[a], "-pep2")==0)
	        {
		pep_compare=1;
		sprintf ( pep2_file, "%s", argv[++a]);
		}
	else if (  strcmp ( argv[a], "-pep")==0)
	        {
		pep_compare=1;
		}
	else if (  strcmp ( argv[a], "-count")==0)
	        {
		count=1;
		}
	else if (  strcmp ( argv[a], "-print_pair")==0)
	        {
		print_pair=1;
		}		
	else if (  strcmp ( argv[a], "-output_aln")==0)
	        {
		output_aln=1;
		}		
	else if (  strcmp ( argv[a], "-output_aln_threshold")==0)
	        {
		output_aln_threshold=atoi(argv[++a]);
		}
	else if (  strcmp ( argv[a], "-output_aln_file")==0)
	        {
		sprintf (output_aln_file,"%s",argv[++a]);
		}
	else if (  strcmp ( argv[a], "-output_aln_format")==0)
	        {
		sprintf (output_aln_format,"%s",argv[++a]);
		}
	else if (  strcmp ( argv[a], "-output_aln_modif")==0)
	        {
		sprintf (output_aln_modif,"%s",argv[++a]);
		}
	else if ( strcmp ( argv[a], "-st")==0)
		{
    		sprintf ( struct_file  [n_structure], "%s", argv[++a]);
		if (!NEXT_ARG_IS_FLAG && is_a_struc_format (argv[a+1]))
		    sprintf ( struct_format[n_structure], "%s", argv[++a]);
		else
		    sprintf ( struct_format[n_structure], "%s", "pep");
		
		if ( !NEXT_ARG_IS_FLAG && strcmp ( argv[a+1], "conv")==0)
		   {
		   a++;
		   while(!NEXT_ARG_IS_FLAG)
		         { 
			 sprintf ( symbol_list[n_structure][n_symbol[n_structure]], "%s", argv[++a]);
			 n_symbol[n_structure]++; 
			 }
		   }
		else if (!NEXT_ARG_IS_FLAG)
		         {
			 symbol_list[n_structure]=make_symbols ( argv[++a], &n_symbol[n_structure]);
			 }

		else 
		         {
                         symbol_list[n_structure]=make_symbols ( "any", &n_symbol[n_structure]);
			 }
		
	  	n_structure++;
		
		}
	else if ( strcmp  (argv[a], "-sep")==0)
	        {
		if ( !NEXT_ARG_IS_FLAG)
		    get_separating_char ( argv[++a][0], &sep_l, &sep_r);
		else
		    sep_l=sep_r=' ';
		}
	else if ( strncmp ( argv[a], "-io_format",5)==0)
	        {
		sprintf ( io_format, "%s", argv[++a]);
		}
	else if ( strcmp ( argv[a], "-io_cat")==0)
		{
		if ( argv[++a][0]!='[')
		   {
		   if ( strcmp ( argv[a], "3d_ali")==0)sprintf ( category_list, "[b][b]+[h][h]=[struc];[*][*]=[tot]");
		   }
		else
		   {
		   sprintf ( category_list, "%s", argv[a]);
		   }
		}
	else
		{
		fprintf ( stdout, "\nOPTION %s UNKNOWN[FATAL]\n", argv[a]);
		myexit (EXIT_FAILURE);
		}
	}
 
/*PARAMETER PROCESSING*/
 //---Maria added this---//
       CL=NULL;

if ( pep_compare==1 || count==1)aln_compare=0;
if ( aln_compare==1)pep_compare=0; 

/*READ THE TOTAL SEQUENCES*/
       seq_list=declare_char ( 100,STRING);
       n_seq_file=0;

        if (  alignment1_file[0] && !check_file_exists ( alignment1_file))
	 {
	   fprintf (stderr, "\nERROR: %s DOES NOT EXIST[FATAL:%s]\n",  alignment1_file, PROGRAM);
	   myexit(EXIT_FAILURE);
	 } 
	if (  alignment2_file[0] && !check_file_exists ( alignment2_file))
	 {
	   fprintf (stderr, "\nERROR: %s DOES NOT EXIST[FATAL:%s]\n",  alignment2_file, PROGRAM);
	   myexit(EXIT_FAILURE);
	 }
	if (  pep1_file[0] && !check_file_exists ( pep1_file))
	 {
	   fprintf (stderr, "\nERROR: %s DOES NOT EXIST[FATAL:%s]\n",  pep1_file, PROGRAM);
	   myexit(EXIT_FAILURE);
	 }
	if (  pep2_file[0] && !check_file_exists ( pep2_file))
	 {
	   fprintf (stderr, "\nERROR: %s DOES NOT EXIST[FATAL:%s]\n",  pep2_file, PROGRAM);
	   myexit(EXIT_FAILURE);
	 }
	
       if ( alignment1_file[0])sprintf ( seq_list[n_seq_file++], "A%s", alignment1_file);
       if ( alignment2_file[0])sprintf ( seq_list[n_seq_file++], "A%s", alignment2_file);
       if ( pep1_file[0])sprintf ( seq_list[n_seq_file++], "S%s", pep1_file);
       if ( pep2_file[0])sprintf ( seq_list[n_seq_file++], "S%s", pep2_file);
       if ( lib_file[0]){		 //----------------Maria added this-----------------// 
	   
	  int *entry, *aln_pos;
	  char **tmp_pos;
	  int tot_we=0, al_we=0, len, i, j;
	  double comp_Score=0;
	
	  CL=read_constraint_list_file( CL, lib_file);
	  printf("%s  %s\n",lib_file, alignment2_file);
	  
	  TOT_SEQ=read_seq_in_n_list ( seq_list, n_seq_file, NULL, NULL);
	  B=declare_aln (TOT_SEQ);
	  main_read_aln ( alignment2_file, B);
	  
	  tmp_pos= (char **) malloc ( CL->S->nseq * sizeof(char *) );
	  for (i=0; i< B->nseq; i++){
	    for (j=0; j< CL->S->nseq; j++){
		if( strcmp(CL->S->name[j], B->name[i])==0 ){
		      tmp_pos[j]=B->seq_al[i]; 						//printf(" %s   \n", tmp_pos[j]);
		      break;
		}
	    }
	  }
	  B->seq_al = tmp_pos;								//printf("\n");   for (i=0; i< B->nseq; i++){ printf(" %s   \n", B->seq_al[i]);  }
	    
	  SB=aln2seq ( B);
	  CL_B=declare_constraint_list ( SB, NULL, NULL, 0, NULL, NULL);	  
	  CL_B=aln2constraint_list_full (B,CL_B, "1");
	  
	  while ( ( entry=extract_entry (CL) ) ){ tot_we+=entry[WE]; }
	 
	  while ( ( entry=extract_entry (CL_B) ) ){
	   
	      if ((main_search_in_list_constraint (entry,&pos_in_clist,4,CL))!=NULL){   //printf ("-%d- S1:%d  S2:%d    R1:%d R2:%d", ++c, entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]); 
		
		len=CL->residue_index[ entry[SEQ1] ][ entry[R1] ][0];
		i=0;
		while ((i<len) && (CL->residue_index[ entry[SEQ1] ][ entry[R1] ][1+(i*5)])!=entry[SEQ2])
		  ++i;
		
		if (i<len)
		{									 //printf("     %d\n", CL->residue_index[ entry[SEQ1] ][ entry[R1] ][(i*5)+3]);
		    al_we+= (CL->residue_index[ entry[SEQ1] ][ entry[R1] ][(i*5)+3]);
		} 
	      } 
	  }
	  
	  comp_Score=(double)al_we/tot_we; 
	  printf("SCORE:  %d   %d    %.3lf\n", al_we, tot_we, comp_Score );
       }
       //--------------------------------------------------------------------------------//
       
       else{
       
       TOT_SEQ=read_seq_in_n_list ( seq_list, n_seq_file, NULL, NULL);
             
       A=declare_aln (TOT_SEQ);
       B=declare_aln (TOT_SEQ);

/*1 COMPARISON OF THE SEQUENCES*/
    if ( pep_compare==1 || count==1)
       {
       f=0; 

       if ( pep1_file[0]!='\0')SA=main_read_seq (pep1_file);
       else if (alignment1_file[0]!='\0')  
           {	       
	   main_read_aln ( alignment1_file, A);
	   SA=aln2seq ( A);
	   }
       else 
	   {
	   main_read_aln ("stdin", A);
	   sprintf ( alignment1_file, "stdin");
	   SA=aln2seq ( A);
	   }
       if ( pep2_file[0]!='\0')SB=main_read_seq (pep2_file);
       else if (alignment2_file[0]!='\0') 
           {
	   main_read_aln ( alignment2_file, B);
	   
	   SB=aln2seq (B);
	   }
       else 
           {
	   SB=SA;
	   sprintf ( alignment2_file, "%s", alignment1_file );
	   }
       buf1=(pep1_file[0]!='\0')?pep1_file: alignment1_file;
       buf2=(pep2_file[0]!='\0')?pep2_file: alignment2_file;
       /*Output of the Results*/    

       fp=vfopen ( io_file, "w");      
       

       if ( count==1)
          {
	  fprintf (fp, "Number of seq: %d %d\n", SA->nseq,SA->nseq);
	  for ( a=0; a< SA->nseq; a++) fprintf (fp, "%-15s %d %d\n", SA->name[a], (int)strlen (SA->seq[a]), (int)strlen (SA->seq[a]));
	  }
	      
       if (SA->nseq!=SB->nseq)
	   {
	   
	   fprintf ( fp, "DIFFERENCE TYPE 1: Different number of sequences %3d/%3d",SA->nseq,SB->nseq);
	   f=1;
	   }

       trim_seq ( SA, SB);
       for ( a=0; a< SA->nseq; a++)
           {
	   lower_string (SA->seq[a]);
	   lower_string (SB->seq[a]);
	   ungap ( SA->seq[a]);
	   ungap ( SB->seq[a]);

	   if ( strcmp ( SA->seq[a], SB->seq[a])!=0) 
	       {
	       fprintf ( fp, "DIFFERENCE TYPE 2: %s is different in the 2 files\n", SA->name[a]);
	       f=1;
	       }
	   }
       for ( a=0; a< SA->nseq; a++)
           {
	   lower_string (SA->seq[a]);
	   lower_string (SB->seq[a]);
	   ungap ( SA->seq[a]);
	   ungap ( SB->seq[a]);

	   if ( strlen ( SA->seq[a])!= strlen (SB->seq[a]))
	       {
		 fprintf ( fp, "DIFFERENCE TYPE 3: %s has != length in the 2 files (%d-%d)\n", SA->name[a],(int)strlen ( SA->seq[a]), (int)strlen (SB->seq[a]));
	       f=1;
	       }
	   }
       if ( f==1)
           {
	   fprintf ( fp, "\nDIFFERENCES found between:\n\t%s\n\t%s\n**********\n\n",buf1, buf2);
	   }
       fclose (fp);
       }

/*2 COMPARISON OF THE ALIGNMENTS*/			 
    else if ( aln_compare==1)
       {

       n_categories=parse_category_list ( category_list, category, n_sub_categories);
       sim_n_categories=parse_category_list ( sim_category_list, sim_category, sim_n_sub_categories);
       
       main_read_aln ( alignment1_file, A);
       main_read_aln ( alignment2_file, B);
       CLS=trim_aln_seq ( A, B);
       
       
       defined_residueA=get_defined_residues (A);
       defined_residueB=get_defined_residues (B);
       

       A=thread_defined_residues_on_aln(A, defined_residueA);
       A=thread_defined_residues_on_aln(A, defined_residueB);
       B=thread_defined_residues_on_aln(B, defined_residueA);
       B=thread_defined_residues_on_aln(B, defined_residueB);
      
     
       CL_A=declare_constraint_list ( CLS, NULL, NULL, 0, NULL, NULL);
       CL_B=declare_constraint_list ( CLS, NULL, NULL, 0, NULL, NULL);
             
      
       CL_A=aln2constraint_list_full (A,CL_A, "sim");
       CL_B=aln2constraint_list_full (B,CL_B, "sim");
       
       
       
       glob=(int*)vcalloc ( A->nseq+1, sizeof (int));
       pw_glob=declare_int ( A->nseq+1, A->nseq+1);


       if ( strm( compare_mode, "sp"))
	    {
	      int *entry;
	      while ( (entry=extract_entry (CL_A)) )
		{
// 		  HERE ("-%d-S1:%d S2:%d R1:%d R2:%d", ++c, entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]);
		  s1=entry[SEQ1];
		  s2=entry[SEQ2];
		  glob[0]++;
		  glob[s1+1]++;
		  glob[s2+1]++;
		  pw_glob[s1][s2]++;
		  pw_glob[s2][s1]++;
		  entry[MISC]=1;
		  
		  if ((main_search_in_list_constraint (entry,&pos_in_clist,4,CL_B))!=NULL)
		    add_entry2list (entry,CL_A);
		  
		}
	    }
       else if ( strm( compare_mode, "column") )
	    {
	      int *entry;
	      posA=aln2pos_simple_2(A);
	      posB=aln2pos_simple_2(B);
	      seq_cache=declare_int ( A->nseq, A->len_aln+1);
	      for ( a=0; a< A->len_aln; a++)
		for ( b=0; b<B->len_aln; b++)
		  {
		    is_same=compare_pos_column(posA, a, posB, b, A->nseq);    
		    if(is_same)
		    {
		      for (c=0; c< A->nseq;c++)
			  {
			    if(posA[c][a]>0)
			    {
				seq_cache[c][posA[c][a]]=1;
			    }
			  }
		    }
		  }
	      
	      while ( (entry=extract_entry (CL_A)) )
		{
		  s1=entry[SEQ1];
		  s2=entry[SEQ2];
		  r1=entry[R1];
		  entry[MISC]=0;
		  glob[0]++;
		  glob[s1+1]++;
		  glob[s2+1]++;
		  pw_glob[s1][s2]++;
		  pw_glob[s2][s1]++;
		  
		  if (seq_cache[s1][r1])
		  {
		     entry[MISC]=1;
		     add_entry2list (entry, CL_A);
		  }
		}
	      free_int (posA, -1);
	      free_int (posB, -1);
	      free_int (seq_cache, -1);      
	    }	
	else if ( strm( compare_mode, "tc") )
	{
	      correct_column = (int*)vcalloc(A->len_aln+1, sizeof (int));
	      posA=aln2pos_simple_2(A);
	      posB=aln2pos_simple_2(B);
	      for ( a=0; a< A->len_aln; a++)
		for ( b=0; b<B->len_aln; b++)
		  {
		    is_same = compare_pos_column(posA, a, posB, b, A->nseq);
		    if(is_same)
		      correct_column[a] = is_same;
		  }
	      free_int (posB, -1);
	}
	
       for ( a=0; a< n_structure; a++)
           {
	       ST=read_structure (struct_file[a],struct_format[a], A,B,ST,n_symbol[a], symbol_list[a]); 
	   }
       
      /*RESULT ARRAY DECLARATION*/    
    
       tot_count=declare_int (n_categories+1, A->nseq+1);
       pos_count=declare_int (n_categories+1, A->nseq+1);
       pw_tot_count=(int***)vcalloc ( A->nseq, sizeof (int**));
       for ( a=0; a< A->nseq; a++)pw_tot_count[a]=declare_int ( A->nseq, n_categories);

    
       pw_pos_count=(int***)vcalloc ( A->nseq, sizeof (int**));
       for ( a=0; a< A->nseq; a++)pw_pos_count[a]=declare_int ( A->nseq, n_categories);
 
     /*COMPARISON MODULE*/
       if (strm( compare_mode, "tc") )
	 {
	   int column_is_structure;
	   int pair_is_structure;
	   
	   for ( a=0; a< n_categories; a++)
	     {     
	       for(b = 0; b < A->len_aln; b++)
		 {
		   column_is_structure = 1;
		   pair_is_structure = 0;
		   for (s1=0; s1< A->nseq-1; s1++)
		     {
		       for(s2=s1+1; s2 < A->nseq; s2++)
			 {
			   if ((posA[s1][b]>0) && (posA[s2][b]>0))
			     {  
			       pair_is_structure=is_in_struct_category(s1,s2,posA[s1][b],posA[s2][b],ST,category[a], n_sub_categories[a]);     
			       column_is_structure = (column_is_structure && pair_is_structure);
			       if(pair_is_structure)
				 {
				   tot_count[a][s1+1]++;
				   tot_count[a][s2+1]++;
				   pw_tot_count[s1][s2][a]++;
				   pw_tot_count[s2][s1][a]++;
				   if (correct_column[b])
				     {
				       pw_pos_count[s1][s2][a]++;
				       pw_pos_count[s2][s1][a]++;
				       pos_count[a][s1+1]++;
				       pos_count[a][s2+1]++;
				     }
				 }
			     }
			 }
		     }	       
		   if(column_is_structure && pair_is_structure)
		     {
		       tot_count[a][0]++;
		       if(correct_column[b]) 
			 {
			   pos_count[a][0]++;
			 }
		     }
		 }
	     }
	   free_int (posA, -1);
	 }
       else
	 {
	   /*COMPARISON MODULE*/
	   int *entry;
	   for ( a=0; a< n_categories; a++)
	     {
	       while ( (entry=extract_entry (CL_A)) )
		 {
		   s1=entry[SEQ1];
		   s2=entry[SEQ2];
		   r1=entry[R1];
		   r2=entry[R2];
		   c= entry[MISC];
		   if ( is_in_struct_category ( s1, s2, r1, r2, ST, category[a], n_sub_categories[a]))
		     {
		       tot_count[a][0]++;
		       tot_count[a][s1+1]++;
		       tot_count[a][s2+1]++;
		       pw_tot_count[s1][s2][a]++;
		       pw_tot_count[s2][s1][a]++;
		       
		       if(print_pair)
		       {
			  HERE("%s,%s,%s",get_structure_residue(s1,r1,ST), get_structure_residue(s2,r2,ST),(c==1)?"Y":"N");
		       }
		       
		       if ( c==1)
			 {
			   pw_pos_count[s1][s2][a]++;
			   pw_pos_count[s2][s1][a]++;
			   pos_count[a][0]++;
			   pos_count[a][s1+1]++;
			   pos_count[a][s2+1]++;
			 }
		     }
		 }
	     }
	 }
       // 	printf("pos=%d,tot=%d\n", pos_count[0][0], tot_count[0][0]);
       /*Measure of Aligned Sequences Similarity*/
       
       sim=get_aln_compare_sim ((strcmp (sim_aln, "al1")==0)?A:B, ST,sim_category[0], sim_n_sub_categories[0], sim_matrix);
       sim_param=analyse_sim ((strcmp (sim_aln, "al1")==0)?A:B, sim);
       
   /*Fill the Result_structure*/
       R=(Result*)vcalloc ( 1, sizeof (Result));
    
       R->grep_list=grep_list;
       R->n_greps=n_greps;
       R->A=A;
       R->B=B;

       R->S=S;
       R->ST=ST;
       R->sim_aln=sim_aln;
       R->alignment1_file=alignment1_file;
       R->alignment2_file=alignment2_file;
       R->io_format=io_format;
       R->n_structure=n_structure;
       R->struct_file=struct_file;
       R->struct_format=struct_format;
       R->n_symbol=n_symbol;
       R->symbol_list=symbol_list;

       
       R->tot_count=tot_count;
       R->pos_count=pos_count;
       R->pw_tot_count=pw_tot_count;
       R->pw_pos_count=pw_pos_count;
       R->glob=glob;
       R->pw_glob=pw_glob;
       R->n_categories=n_categories;
       R->category=category;
       R->category_list=category_list;
       R->n_sub_categories=n_sub_categories;
       R->sim=sim;
       R->sim_param=sim_param;
       R->sim_matrix=sim_matrix;
       R->sim_n_categories=sim_n_categories;
       R->sim_category=sim_category;
       R->sim_category_list=sim_category_list;
       R->sim_n_sub_categories=sim_n_sub_categories; 
       R->sep_r=sep_r;
       R->sep_l=sep_l;

      /*Output of the Results*/    

       fp=vfopen ( io_file, "w");      
       fp=output_format (io_format, fp, R);
       vfclose ( fp); 


     /*Rewriting of Alignment A*/
       if ( output_aln)
	 {
	   int *entry;
	   A->residue_case=2;
	   aln_output_tot  =declare_int ( A->nseq, A->len_aln+1);
	   aln_output_count=declare_int ( A->nseq, A->len_aln+1);
	   
	   while ( (entry=extract_entry (CL_A)) )
	     {
	       aln_output_tot[entry[SEQ1]][entry[R1]]++;
	       aln_output_tot[entry[SEQ2]][entry[R2]]++;
	       
	       aln_output_count[entry[SEQ1]][entry[R1]]+=entry[MISC];
	       aln_output_count[entry[SEQ2]][entry[R2]]+=entry[MISC];
	     }
	  
	   for ( a=0; a< A->nseq; a++)
	     {
	       
	       for (c=0, b=0; b< A->len_aln; b++)
		 {
		   if ( !is_gap(A->seq_al[a][b]))
		     {
		       c++;
		       if ( aln_output_tot[a][c] && ((aln_output_count[a][c]*100)/aln_output_tot[a][c])<output_aln_threshold)
			 { 
			   if ( strm (output_aln_modif, "lower"))
				A->seq_al[a][b]=tolower(A->seq_al[a][b]);
			   else
			     A->seq_al[a][b]=output_aln_modif[0];
			 }
		       else A->seq_al[a][b]=toupper(A->seq_al[a][b]);
		     }
		     
		 }
	     }
	   A->score_aln=(int)(R->tot_count[0][0]==0)?0:((R->pos_count[0][0]*100)/(float)R->tot_count[0][0]);

	   output_format_aln (output_aln_format,A,NULL,output_aln_file);
	   
	   free_int ( aln_output_tot, -1);
	   free_int ( aln_output_count, -1  );
	 }
       }
       }
    return EXIT_SUCCESS;
    }   
/************************************************************************************/
/*                                                                                  */
/*                            OUTPUT                                                */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/
FILE *output_format (char *iof,FILE *fp, Result *R)
      {
      int a;
      int l;
      
      /*
      H: files Header
      h: basic header;
      s: sequence results
      t: total results (global);
      p: pairwise_results;
      */
      l=strlen ( iof);

      for ( a=0; a< l; a++)
          {
	  if ( iof[a]=='H')fp=output_large_header (fp,R);
	  else if ( iof[a]=='h')fp=output_header (fp,R);
	  else if ( iof[a]=='t')fp=output_total_results (fp, R);
	  else if ( iof[a]=='s')fp=output_sequence_results (fp,R);
	  else if ( iof[a]=='p')fp=output_pair_wise_sequence_results(fp,R);
	  }
      return fp;
      }

FILE *output_pair_wise_sequence_results (FILE *fp,  Result *R)
     {
     int a,c,d;
     
     
     for ( c=0; c<(R->A)->nseq-1; c++)
         {
	 for ( d=c+1; d< (R->A)->nseq; d++)    
	     {
	     fprintf (fp, "%-10s %-10s%s",(R->A)->name[c],(R->A)->name[d],SSPACE);
	     fprintf (fp, "%5.1f%s", R->sim[c][d], SSPACE);
	     
	     for (a=0; a< R->n_categories; a++)
		 {
		 fprintf ( fp, "%5.1f ",(R->pw_tot_count[c][d][a]==0)?0:((float)(R->pw_pos_count[c][d][a]*100)/(float)R->pw_tot_count[c][d][a]));
		 fprintf ( fp, "%c%5.1f%c%s",(R->sep_l),(R->pw_glob[c][d]==0)?0:((float)(R->pw_tot_count[c][d][a]*100)/(float)R->pw_glob[c][d]),(R->sep_r),SSPACE);
		 }
	     fprintf ( fp, "%c%5d%c\n",(R->sep_l), R->pw_glob[c][d],(R->sep_r));
	     }
	 }
	 
    return fp;
    }
FILE *output_sequence_results (FILE *fp,  Result *R)
     {
     int a,c;
    
     for ( c=1; c<=R->A->nseq; c++)
         {
	 fprintf (fp, "%-10s %-10s%s",(R->A)->name[c-1], "..",SSPACE);
	 fprintf (fp, "%5.1f%s", R->sim_param[c-1][0],SSPACE);
	 for (a=0; a< R->n_categories; a++)
		 {
		 fprintf ( fp, "%5.1f ",(R->tot_count[a][c]==0)?0:((float)(R->pos_count[a][c]*100)/(float)R->tot_count[a][c]));
		 fprintf ( fp, "%c%5.1f%c%s",(R->sep_l),(R->glob[c]==0)?0:((float)(R->tot_count[a][c]*100)/(float)R->glob[c]),(R->sep_r), SSPACE);
		 }
	 fprintf ( fp, "%c%5d%c\n",(R->sep_l), R->glob[c],(R->sep_r));
	 }
    return fp;
    }
	 

FILE *output_total_results (FILE *fp,  Result *R)
     {
     int a;
     
     
     fprintf ( fp, "%-13s %-7d%s",extract_suffixe (R->alignment1_file),(R->A)->nseq, SSPACE);
     fprintf (fp, "%5.1f%s", R->sim_param[(R->A)->nseq][0], SSPACE);
     for (a=0; a< R->n_categories; a++)
         {
	 fprintf ( fp, "%5.1f ",(R->tot_count[a][0]==0)?0:((float)(R->pos_count[a][0]*100)/(float)R->tot_count[a][0]));
	 fprintf ( fp, "%c%5.1f%c%s",(R->sep_l),(R->glob[0]==0)?0:((float)(R->tot_count[a][0]*100)/(float)R->glob[0]),(R->sep_r), SSPACE);
	 }
     fprintf ( fp, "%c%5d%c\n",(R->sep_l),  R->glob[0],(R->sep_r));
     return fp;
     }
FILE *output_header (FILE *fp, Result *R)
     {
     int a;

     
     fprintf ( fp, "%s\n",generate_string ( R->n_categories*(13+strlen(SSPACE))+31+2*strlen(SSPACE),'*'));
 
     fprintf ( fp, "%-10s %-10s %s%-3s%s", "seq1", "seq2",SSPACE,"Sim",SSPACE);
     
     for ( a=0; a< R->n_categories; a++)
	 fprintf ( fp, "%-12s%s ",R->category[a][0], SSPACE);
     fprintf (fp, "%-5s", "Tot");
     fprintf (fp, "\n");
     return fp;
     }
FILE *output_large_header ( FILE *fp, Result *R)
     {
     int a, b;

     fprintf ( fp, "AL1: %s\n", R->alignment1_file);
     fprintf ( fp, "AL2: %s\n", R->alignment2_file);
     for ( a=0; a< R->n_structure; a++)
          {
	   fprintf (fp, "ST %d: %s [%s]/[", a, R->struct_file[a], R->struct_format[a]);
	   for ( b=0; b< R->n_symbol[a]; b++)fprintf (fp, "%s ", R->symbol_list[a][b]);
	   fprintf ( fp, "]\n");
	  }
     return (fp);
     }

void get_separating_char ( char s, char *l, char *r)
    {
    if ( s=='{' || s=='}')
       {
       l[0]='{';
       r[0]='}';
       return;
       }
    else if ( s==']' || s=='[')
       {
       l[0]='[';
       r[0]=']';
       return;
       }
    else if ( s==')' || s=='(')
       {
       l[0]='(';
       r[0]=')';
       return;
       }
    else
       {
       l[0]=s;
       r[0]=s;
       return;
       }
    }

/************************************************************************************/
/*                                                                                  */
/*                            SIM MEASURE                                           */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/
float **get_aln_compare_sim ( Alignment *A, Structure *ST, char **cat, int n_cat, char *matrix)
      {
      int a, b, c;

      float **sim;
      char cr1, cr2;
      int r1, r2;
      int p1, p2;
      float pos,tot;

      sim=declare_float ( A->nseq, A->nseq);

      for ( a=0; a<A->nseq-1; a++)
          {
	  for (b=a+1; b< A->nseq; b++)
	      {
	      for ( r1=0, r2=0,tot=0, pos=0,c=0; c< A->len_aln; c++)
	          {
		  p1=is_gap(A->seq_al[a][c]);
		  p2=is_gap(A->seq_al[b][c]);
		  
		  r1+=1-p1;
		  r2+=1-p2;
		  cr1=A->seq_al[a][c];
		  cr2=A->seq_al[b][c];
		  if (!p1 && !p2)
		      {
		      if (is_in_struct_category (a, b, r1, r2, ST, cat, n_cat))
			  {
			  tot++;
		          pos+=is_in_same_group_aa ( cr1, cr2, 0, NULL, matrix);
			  }
		      }
		  }
	      sim[a][b]=sim[b][a]=(tot==0)?0:((pos*100)/tot);
	      }
	  }
      return sim;
      }
float **analyse_sim ( Alignment *A, float **sim)
	{
	int a,b,c;

	float **an, d;
	
	
	an=declare_float ( A->nseq+1,2);
	
	for (d=0, a=0; a< A->nseq;a++)
		{
		for ( b=0; b< A->nseq; b++)
			{
			if ( b!=a)
				{
				an[a][0]+=sim[a][b];
				an[A->nseq][0]+=sim[a][b];
				d++;
				}
			}
		an[a][0]=((float)an[a][0]/(float)(A->nseq-1));
		}
	
	
	an[A->nseq][0]=an[A->nseq][0]/d;
		
	for ( d=0,a=0; a< A->nseq; a++)
		{
		for ( b=0; b< A->nseq; b++)
			{
			if ( b!=a)
				{
				c=an[a][0]-sim[a][b];
				an[a][1]+=(c>0)?c:-c;
				an[A->nseq][1]+=(c>0)?c:-c;
				d++;
				}
			}
		an[a][1]=((float)an[a][1]/(float)(A->nseq-1));
		}
	
	an[A->nseq][1]=an[A->nseq][1]/d;
	
	return an;
	}
	    	







/************************************************************************************/
/*                                                                                  */
/*                            STRUC ANALYSE                                         */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/
int is_in_struct_category ( int s1, int s2, int r1, int r2, Structure *ST, char **cat, int n_cat)
    {
    int a;
    static char *struc_r1;
    static char *struc_r2;
    char first[STRING];
    char second[STRING];
    static int **r;
    
    

    
    if ( ST==NULL)return 1;
    if ( struc_r1!=NULL)
       {
       vfree ( struc_r1);
       vfree ( struc_r2);
       }

    if ( r==NULL)r=declare_int (2, 2);
    else r[0][0]=r[1][1]=r[1][0]=r[0][1]=0;
    
    
    struc_r1=get_structure_residue ( s1, r1, ST);
    struc_r2=get_structure_residue ( s2, r2, ST);

    
    for ( a=1; a< n_cat; a+=2)
	{
        sprintf ( first, "%s", cat[a]);
	sprintf ( second,"%s", cat[a+1]);
	r[0][0]=struc_matches_pattern ( struc_r1, first);
	r[0][1]=struc_matches_pattern ( struc_r2, first);
	r[1][0]=struc_matches_pattern ( struc_r1, second);
	r[1][1]=struc_matches_pattern ( struc_r2, second);
	
	if ( (r[0][0]&&r[1][1])||(r[1][0]&&r[0][1]))return 1;
	}	
   return 0;
   }

char * get_structure_residue (int seq, int res, Structure *ST)
    {
    int a;
    char *s;
    
    s=(char*)vcalloc ( ST->n_fields+1, sizeof (char));
    for (a=0; a< ST->n_fields; a++)
	s[a]=ST->struc[seq][res][a];
    s[a]='\0';
    return s;
    }

int struc_matches_pattern ( char *struc, char *pattern)
    {
    char p[STRING];
    char *y;
    
    int b,l;
    

    sprintf ( p, "%s", pattern);
    
    if ( strcmp (p, "*")==0)return 1;
    else
        {
	l=strlen ( struc);
	y=strtok ( p, ".");

	for ( b=0; b<l; b++)
	        {

		if ( strcmp ( y, "*")==0);		   
		else if ( strchr(y, struc[b])==NULL)return 0;
		y=strtok ( NULL, ".");
		
		if ( y==NULL && b<(l-1))
		   {
		   crash ( "\nA Structure is Missing[FATAL]");
		   }
		}
	}

    return 1;
    }
    

Structure * read_structure (char *fname, char *format, Alignment *A,Alignment *B, Structure *ST, int n_symbols, char **symbol_table)        
        {
               
        Alignment *StrucAln;
        Sequence *SEQ=NULL, *SA=NULL;
	
        


        if (ST==NULL)ST=declare_structure (A->nseq, A->seq_al);
        else
            ST=extend_structure ( ST);
        
	StrucAln=declare_Alignment(NULL); 
        if (strm ( format, "pep"))
                        {
                        SEQ=main_read_seq ( fname);
			StrucAln=seq2aln(SEQ,StrucAln,0);
			}
        
        else if ( strcmp ( format, "aln")==0)
                        {
			StrucAln=main_read_aln (fname, StrucAln);
			}

	reorder_aln(StrucAln, A->name, A->nseq);

	SA=aln2seq(StrucAln);
	string_array_convert (SA->seq, SA->nseq, n_symbols, symbol_table);
        seq2struc (SA, ST);
	
	free_aln(StrucAln);
	if (SEQ)free_sequence (SEQ, SEQ->nseq);
	
	return ST;
        }

int parse_category_list ( char *category_list_in, char ***category, int*n_sub_categories)
    {
    int n,a;
    char *category_list;
    char *y,*z;
    category_list=(char*)vcalloc ( strlen(category_list_in)+1, sizeof (char));
    sprintf (category_list, "%s", category_list_in);

    n=0;
    z=category_list;
    while ((y=strtok(z, ";"))!=NULL)
        {
	sprintf ( category[n++][2], "%s", y);
	z=NULL;
	}

    
    for ( a=0; a< n; a++)
	{
	 sprintf (category_list,"%s",strtok(category[a][2], "="));
	 sprintf (category[a][0],"%s",strtok (NULL, "="));
	 sprintf ( category[a][++n_sub_categories[a]],"%s", strtok(category_list, "[]+"));
	 while ( ((y=strtok(NULL, "[]+"))!=NULL))
	        {
		if ( strcmp (y, "#")==0)y=category[a][n_sub_categories[a]];
	        sprintf ( category[a][++n_sub_categories[a]],"%s",y);
		}
	}
    return n;
    }


int is_a_struc_format (char *format)
    {
    if ( strcmp ( format, "pep")==0)return 1;
    if ( strcmp ( format, "aln")==0)return 1;
    return 0;
    }
float aln_file2compare (char *aln1, char *aln2)
{
  Alignment *A, *B;
  float sim;
  A=main_read_aln ( aln1, NULL);
  B=main_read_aln ( aln2, NULL);
  sim=aln2compare (A, B);
  free_aln (A);
  free_aln (B);
  return sim;
}
float aln2compare (Alignment *A, Alignment *B)
{
  Sequence *S;
  Constraint_list *CL_A, *CL_B;
  int tot=0;
  int id=0;
  int r, s;
  
  S=trim_aln_seq ( A, B);
  CL_A=declare_constraint_list ( S, NULL, NULL, 0, NULL, NULL);
  CL_B=declare_constraint_list ( S, NULL, NULL, 0, NULL, NULL);
  
  
  CL_A=aln2constraint_list_full (A,CL_A, "sim");
  CL_B=aln2constraint_list_full (B,CL_B, "sim");

  S=CL_A->S;
  for (s=0; s<S->nseq;s++)
    for (r=1; r<=S->len[s]; r++)
      compare_residues_between_clists (s, r, CL_A, CL_B, &id, &tot);
  
  free_constraint_list (CL_A);
  free_constraint_list (CL_B);
  free_sequence (S, -1);
  
  return (id>0)?(100*(float)id/(float)tot):0;
}
/************************************************************************************/
/*                                                                                  */
/*                            Informations                                          */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/
void output_informations ()
{
fprintf ( stderr, "\nPROGRAM: %s (%s)\n",PROGRAM,VERSION);
fprintf ( stderr, "******INPUT***************************");
fprintf ( stderr, "\n-al1 al1_file");
fprintf ( stderr, "\n-al2 al2_file");
fprintf ( stderr, "\n-compare_mode [sp] or column or tc");
fprintf ( stderr, "\n-pep (compare only the sequences");
fprintf ( stderr, "\n-count");
fprintf ( stderr, "\n-pep1 pep1_file");
fprintf ( stderr, "\n-pep1 pep2_file");
fprintf ( stderr, "\n-st  str_file st_format conversion");
fprintf ( stderr, "\n   **st_format:  aln, pep");
fprintf ( stderr, "\n   **conversion: 3d_ali,  conv abcZ #X");
fprintf ( stderr, "\nNOTE: Several structures in a row are possible"); 

fprintf ( stderr, "\n\n****DISTANCE MEASURE*****************");
fprintf ( stderr, "\n-sim_cat category_format  or category_name");
fprintf ( stderr, "\n       **category_format: [*][*]=[tot]");
fprintf ( stderr, "\n       **category_name  : 3d_ali (<=>[h][h]+[e][e]=[Struc]");
fprintf ( stderr, "\n-sim_matrix matrix_name");
fprintf ( stderr, "\n          **matrix_name: idmat,pam250mt..");
fprintf ( stderr, "\n-sim_aln al1 or al2");
fprintf ( stderr, "\n\n****COMPARISON***********************");
fprintf ( stderr, "\n-io_cat category_format or category_name");
fprintf ( stderr, "\n      **category_format: [*][*]=[tot]");
fprintf ( stderr, "\n      **category_name  : 3d_ali(<=>[h][h]+[e][e]=[Struc];[*][*]=[Tot]");
fprintf ( stderr, "\n\nNOTE: if two structures:");
fprintf ( stderr, "\n[he.123][#]=[he123VShe123];[beh.*][he.2345]=[other]");
fprintf ( stderr, "\n\n****OUTPUT****************************");
fprintf ( stderr, "\n-f stdout");
fprintf ( stderr, "\n   stderr");
fprintf ( stderr, "\n   file_name");
fprintf ( stderr, "\n-io_format hts");
fprintf ( stderr, "\n           H ->large Header");
fprintf ( stderr, "\n           h ->small Header");
fprintf ( stderr, "\n           t ->global (average)results");
fprintf ( stderr, "\n           s ->average results for each sequence");
fprintf ( stderr, "\n           p ->results for each pair of sequences");
fprintf ( stderr, "\n-print_pair           Print structures of pair residues and correctness by ref");
fprintf ( stderr, "\n-output_aln           Outputs al1 with conserved bits in Upper");
fprintf ( stderr, "\n-output_aln_threshold [100]");
fprintf ( stderr, "\n-output_aln_file      [stdout]");
fprintf ( stderr, "\n-output_aln_format    [clustalw]");
fprintf ( stderr, "\n-output_aln_motif     [lower]");
fprintf ( stderr, "\n");
myexit (EXIT_SUCCESS);
}

/*********************************COPYRIGHT NOTICE**********************************/
/*� Centro de Regulacio Genomica */
/*and */
/*Cedric Notredame */
/*Tue Jul  1 10:00:54 WEST 2008. */
/*All rights reserved.*/
/*This file is part of T-COFFEE.*/
/**/
/*    T-COFFEE is free software; you can redistribute it and/or modify*/
/*    it under the terms of the GNU General Public License as published by*/
/*    the Free Software Foundation; either version 2 of the License, or*/
/*    (at your option) any later version.*/
/**/
/*    T-COFFEE is distributed in the hope that it will be useful,*/
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of*/
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*/
/*    GNU General Public License for more details.*/
/**/
/*    You should have received a copy of the GNU General Public License*/
/*    along with Foobar; if not, write to the Free Software*/
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA*/
/*...............................................                                                                                      |*/
/*  If you need some more information*/
/*  cedric.notredame@europe.com*/
/*...............................................                                                                                                                                     |*/
/**/
/**/
/*	*/
/*********************************COPYRIGHT NOTICE**********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

/******************************************************************/
/*                   MOCA DRIVER                                  */
/*                                                                */
/*                                                                */
/******************************************************************/
Constraint_list *prepare_cl_for_moca ( Constraint_list *CL)
   {
    int a, b, c;
    int tot_l, l;
    char **name, **seq;
    Sequence *NS=NULL;

   /*Prepare the constraint list*/
		     CL->do_self=1;
		     CL->get_dp_cost=moca_slow_get_dp_cost;
		     CL->evaluate_residue_pair=moca_residue_pair_extended_list;
		     
   /*Prepare the moca parameters*/
		     (CL->moca)->evaluate_domain=evaluate_moca_domain;
		     (CL->moca)->cache_cl_with_domain=cache_cl_with_moca_domain;
		     (CL->moca)->make_nol_aln=make_moca_nol_aln;
   
   /*Prepare the packing of the sequences*/		     
		     for ( a=0, b=1; a< (CL->S)->nseq; a++)b+=strlen ( (CL->S)->seq[a])+1;
			 
		     seq =declare_char ( 1,b+1);
		     name=declare_char(  1,30);
		     CL->packed_seq_lu  =declare_int ( b, 2);
		     
		     
		     for (tot_l=1,a=0; a< (CL->S)->nseq; a++)
		       {
			 strcat (seq[0], (CL->S)->seq[a]);
			 strcat (seq[0], "X");
			 l=strlen((CL->S)->seq[a]);
			 for ( c=1; c<= l; c++, tot_l++)
			   {
			     CL->packed_seq_lu[tot_l][0]=a;
			     CL->packed_seq_lu[tot_l][1]=c;
			   }
			 CL->packed_seq_lu[tot_l++][0]=UNDEFINED; 
		       }
		     sprintf ( name[0], "catseq");
		     NS=fill_sequence_struc(1, seq, name); 
		     CL->S=add_sequence (NS, CL->S, 0);
		     free_char( seq, -1);
		     free_char(name, -1);
		     free_sequence (NS, NS->nseq);
		     
		     
    return CL;
   }

Alignment ** moca_aln ( Constraint_list *CL)
   {
     /*
       function documentation: start
       
       Alignment ** moca_aln ( Constraint_list *CL)
       
       This function inputs CL and outputs a series of local multiple alignments
       contained in aln_list;
      
       The terminator of aln_list is set to NULL;
       
       function documentation: end
     */
       
     
   static int max_n_domains=1000;
   int n_domains=0;
   
   Alignment **aln_list;
   
   

   aln_list=vcalloc (max_n_domains, sizeof (Alignment *));
   if ((CL->moca)->moca_interactive)aln_list[n_domains++]=extract_domain ( CL);
   else
     {
       while ( (aln_list[n_domains++]=extract_domain ( CL))!=NULL)
	 {
	   if ((CL->moca)->moca_len)break;
	   if ( n_domains==max_n_domains)
	     {
	       n_domains+=1000;
	       aln_list=vrealloc (aln_list, max_n_domains*sizeof (Alignment*));
	     }
	 }
     }
   return aln_list;
   }
   
Alignment * extract_domain ( Constraint_list *CL)
   {
     /*
       function documentation: start
       Alignment * extract_domain ( Constraint_list *CL)
       
       given a CL, this function extracts the next best scoring local multiple alignment
       It returns a CL where the aligned residues have been indicated in (CL->moca)->forbiden_residues;
       
       the local alignment is extracted with the dp function indicated by
                           CL->dp_mode: (gotoh_sw_pair_wise)
       Evaluation:
                  CL->get_dp_cost=slow_get_dp_cost;
		  CL->evaluate_residue_pair=sw_residue_pair_extended_list;
       Continuation:
                 (CL->moca)->evaluate_domain=evaluate_moca_domain;       
       Cache of CL:		  
                 (CL->moca)->cache_cl_with_domain=cache_cl_with_moca_domain;
       Domain post processing:
                 (CL->moca)->make_nol_aln=make_moca_nol_aln;
       function documentation: end
     */
     int min_start, max_start, start,min_len, max_len, len, score;
     int step;
     Alignment *C=NULL;
     Alignment *RESULT=NULL;
     Alignment *EA=NULL;


     

     /*CASE 1: Non Automatic Domain Extraction*/
     if ((CL->moca)->moca_interactive)
	 {
	   return interactive_domain_extraction (CL);
	 }
     else if ((CL->moca)->moca_len)
       {
	 while ((C=extract_domain_with_coordinates (C,(CL->moca)->moca_start,(CL->moca)->moca_len,CL))->nseq==0)(CL->moca)->moca_scale=(CL->moca)->moca_scale*0.9;    
	 RESULT=copy_aln ( C, RESULT);
	 unpack_seq_aln (RESULT, CL);
	 output_format_aln ("mocca_aln",RESULT,EA=fast_coffee_evaluate_output(RESULT, CL),"stdout");
	 free_aln(EA);
	
	 return RESULT;
       }
     else if ( !(CL->moca)->moca_len)
       {
	 analyse_sequence (CL);
	 myexit (EXIT_FAILURE);
       }

     /*CASE 2: Automatic Domain Extraction: Find Coordinates*/
          
       
     start=500;

     step=10;
     min_start=0;
     max_start=strlen ((CL->S)->seq[0]);
     min_len=20;
     max_len=strlen ((CL->S)->seq[0]);
     
     C=extract_domain_with_coordinates (C,13,30,CL);     
     C->output_res_num=1;
     print_aln (C);

     (CL->moca)->moca_scale=-180;
     C=add_seq2aln (CL,C, CL->S);
     print_aln (C);

     (CL->moca)->moca_scale=-160;
     C=add_seq2aln (CL,C, CL->S);
     print_aln (C);
     
     myexit (EXIT_FAILURE);
     
     while ( step>0)
       {
	 C=approximate_domain (min_start,max_start,step,min_len,max_len, step,&start, &len, &score, CL);
	 min_start=start-step;
	 max_start=start+step;
	 min_len=len-step;
	 max_len=len+step;
	 step=step/2;
       }
   
     C=extract_domain_with_coordinates (C,start-10, len+20,CL);
     C->output_res_num=1;
     print_aln (C);

     myexit (EXIT_FAILURE);
     return C;


   }
Alignment * interactive_domain_extraction ( Constraint_list *CL)
   {
     int LEN=0;
     int START=1;
     int SCALE=2;
     int GOPP=3;

     int iteration=0;
     char *choice;
     int a,b, c;
     int index;
     char *s;
     char last_start[100];
     char out_format[100];
     Alignment *RESULT=NULL;
     Alignment *PREVIOUS=NULL;
     Alignment *C=NULL;
     Alignment *EA=NULL;
     
     int **parameters;
     
     
     choice=vcalloc ( 100, sizeof (char));
     parameters=declare_int (10000, 4);
     
     parameters[0][START]=(CL->moca)->moca_start;
     parameters[0][LEN]=  (CL->moca)->moca_len;
     parameters[0][SCALE]=(CL->moca)->moca_scale;
     parameters[0][GOPP]=CL->gop;
     iteration=0;
     sprintf ( last_start, "%d", (CL->moca)->moca_start);
     sprintf ( out_format, "mocca_aln");
     
     print_moca_interactive_choices ();
     while ( !strm4 (choice, "Q","X", "q", "x" ))
	     {
	       c=choice[0];
	       
	       if (c=='b' || c=='B')
		 {
		 iteration-=atoi(choice+1)+1;		   
		   
		 if (iteration<0)iteration=1;  
		 }
	       else
		 {
		 iteration++;
		 parameters[iteration][START]=parameters[iteration-1][START];
		 parameters[iteration][LEN]=parameters[iteration-1][LEN];
		 parameters[iteration][SCALE]=parameters[iteration-1][SCALE];
		 parameters[iteration][GOPP]=parameters[iteration-1][GOPP];

		 if ( c=='>')parameters[iteration][LEN]=atoi(choice+1);
		 else if ( c=='|')
		   {
		   sprintf ( last_start, "%s", choice);
		   parameters[iteration][START]=0;  
		   s=strrchr(choice, ':');
		   
		   if (s==NULL)
		     {
		       parameters[iteration][START]=atoi(choice+1);
		     }
		   else
		     {
		      
		       s[0]='\0';

		       if((index=name_is_in_list (choice+1,(CL->S)->name,(CL->S)->nseq,100))==-1)
			 {
			   fprintf ( stderr, "\n\tERROR: %s NOT in Sequence Set",choice+1);
			   continue;
			 }

		     for ( a=0; a< index; a++)
		       {
			 parameters[iteration][START]+=(CL->S)->len[a]+1;
		       }
		     parameters[iteration][START]+=atoi(s+1)-1;
		     }
		       
		   }
		 else if ( c=='C'||c=='c')parameters[iteration][SCALE]=atoi(choice+1);
		 else if ( c=='G'||c=='g')
		   {
		     parameters[iteration][GOPP]=atoi(choice+1);
		     CL->gop=parameters[iteration][GOPP];
		   }
		 else if (  c=='F'||c=='f')
		   {
		     sprintf ( out_format, "%s", choice+1);
		   }
		 else if ( c=='S'||c=='s')
		   {
		     if (choice[1]=='\0')sprintf ( choice, "default.domain_aln.%d", iteration); 
		     output_format_aln (out_format,RESULT,EA=fast_coffee_evaluate_output(RESULT, CL),choice+1);
		     fprintf (stderr, "\tOutput  file [%15s] in [%10s] format\n",choice+1,out_format);
		     free_aln (EA);
		   }
		 else if (c=='\0')
		   {
		     if ( parameters[iteration][SCALE]>0)
		       {
			 fprintf ( stderr, "\nWARNING: THRESHOLD RESET to 0");
			 parameters[iteration][SCALE]=0;
		       }
		     
		     (CL->moca)->moca_scale=parameters[iteration][SCALE];
		      CL->gop=parameters[iteration][GOPP];

		     C=extract_domain_with_coordinates (C,parameters[iteration][START],parameters[iteration][LEN],CL);
		     
		     if ( C==NULL)
		       {
			 fprintf ( stderr, "\nERROR: ILLEGAL COORDINATES! SEQUENCE BOUNDARY CROSSED\n");
			 for ( b=1,a=0; a< (CL->S)->nseq-1; a++)
			   {
			     
			     fprintf ( stderr, "\n\t%15s=> Abs:[%d %d] Rel:[0 %d]", (CL->S)->name[a],b, b+(CL->S)->len[a]-1,(CL->S)->len[a]);
			     b+=(CL->S)->len[a];
			   }
			 fprintf ( stderr, "\n");
		       }
		     else if (parameters[iteration][START]==0 && parameters[iteration][LEN]==0)
		       {
			 fprintf ( stderr, "\n\tEnter the following parameters:\n\n\t\tSTART  value: |x [Return]\n\t\tLENgth value: >y [Return]\n\t\ttype             [Return]\n\n");
			 fprintf ( stderr, "\n\n\tSTART is measured on the total length of the concatenated sequences\n\tx and y are positive integers\n\n");
		       }

		     else if ( C->nseq==0)
		       {
			 fprintf ( stderr, "\nNO MATCH FOUND: LOWER THE SCALE (C)\n");
		       }
		     else
		       {
			 RESULT=copy_aln ( C, RESULT);
			 unpack_seq_aln (RESULT, CL);
			 RESULT->output_res_num=1;
			 
			 output_format_aln (out_format,RESULT,EA=fast_coffee_evaluate_output(RESULT, CL),"stdout");
			 free_aln(EA);
			 PREVIOUS=copy_aln ( RESULT, PREVIOUS);
			 free_aln (C);
			 print_moca_interactive_choices ();
			
		       }
		   }
		
		 fprintf ( stderr, "\t[ITERATION %3d][START=%s][LEN=%3d][GOPP=%3d][SCALE=%4d]\t",iteration,last_start,parameters[iteration][LEN],parameters[iteration][GOPP],parameters[iteration][SCALE]);
		 a=0;
		 fprintf ( stderr, "Your Choice: ");
		 while ( (c=fgetc(stdin))!='\n')choice[a++]=c;
		 choice[a]=0;
		 }
	     }
     
     if (!RESULT)myexit(EXIT_SUCCESS);
     if ( RESULT)RESULT->output_res_num=0;
     return RESULT;
   }

int print_moca_interactive_choices ()
{
  fprintf ( stderr, "\n**************************************************************");
  fprintf ( stderr, "\n******************** MOCCA: %s           ***********",VERSION);
  fprintf ( stderr, "\n**************************************************************");


fprintf ( stderr, "\nMENU: Type Flag[number] and Return: ex |10");

fprintf ( stderr, "\n\t|x      -->Set     the  START to x");
fprintf ( stderr, "\n\t           100       start=100 on concatenated sequences");
fprintf ( stderr, "\n\t           human:100 start=100 on human sequence");
fprintf ( stderr, "\n\t>x      -->Set     the  LEN   to x");
fprintf ( stderr, "\n\tGx      -->Set     the  Gap Opening Penalty to x");
fprintf ( stderr, "\n\tCx      -->Set     the  sCale to x");
fprintf ( stderr, "\n\tSname   -->Save    the  Alignment ");
fprintf ( stderr, "\n\tFformat -->Save    the  Alignment Format");

fprintf ( stderr, "\n\treturn  -->Compute the  Alignment");

fprintf ( stderr, "\n\tX       -->eXit\n\n");

return 0;
}

Alignment * approximate_domain ( int min_start, int max_start, int step_start,int min_len, int max_len, int step_len, int *best_start, int *best_len, int *best_score, Constraint_list *CL)
   {    
     Alignment *C=NULL;
     int start;
     int len;
     int score;
         
     /*1 Extract the first*/
     best_score[0]=UNDEFINED;
     best_start[0]=min_start;
     best_len[0]=min_len;
     
     for (start=min_start; start< max_start; start+=step_start)
       {
	 for ( len=min_len; len<max_len; len+=step_len)
	   {
	     C=extract_domain_with_coordinates (C,start,len,CL);
	     if ( C==NULL)continue;
	     score=((CL->moca)->evaluate_domain)(C, CL); 
	     fprintf ( stderr, "\nSTART=%d LEN=%3d SCORE=%5d [%d]",start,len,score, C->nseq);
	     

	     if ( best_score[0]==UNDEFINED)best_score[0]=score;
	     if ( score>best_score[0])
	       {
		 best_score[0]=score;
		 best_start[0]=start;
		 best_len[0]=len;
	       }	     
	   }
       }
    
     C=extract_domain_with_coordinates (C,best_start[0], best_len[0],CL);
     C->output_res_num=1;
     return C;
   }
int measure_domain_length ( Constraint_list *CL,Alignment *IN, int start, int min_len, int max_len, int step)
   {
   Alignment *C=NULL;
   int score, best_score,best_len,a, b, l;
   int *score_matrix, *len_matrix;
   int n_val, best_val;

   score_matrix=vcalloc ( max_len, sizeof (int));
   len_matrix=vcalloc ( max_len, sizeof (int));
   
   
   l=strlen ( (CL->S)->seq[0]);
   
   min_len=MAX(0, min_len);
   min_len=MIN(l-start, min_len);
   
   if ( !IN)C=extract_domain_with_coordinates (C,start,min_len, CL);
   else 
     {
     C=copy_aln (IN, C);
     C->len_aln=min_len;
     for ( a=0; a< C->nseq; a++)C->seq_al[a][min_len]='\0';
     C=add_seq2aln (CL,C, CL->S);
     }
   
  best_score= score=((CL->moca)->evaluate_domain)(C, CL);
  
  
  min_len=MAX(0, min_len);
  for ( best_len=best_val=n_val=0,b=min_len; b<max_len && (start+b)<l; b+=step, n_val++)
       {
       if ( !IN)C=extract_domain_with_coordinates (C,start, b, CL);
       else
	   {
	   C=copy_aln (IN, C);
	   C->len_aln=min_len;
	   for ( a=0; a< C->nseq; a++)C->seq_al[a][b]='\0';
	   C=add_seq2aln (CL,C, CL->S);
	   }
       if ( C->len_aln>0 )score=((CL->moca)->evaluate_domain)(C, CL);
       else score=-1;
       
       if ( score< -3000)break;

       fprintf ( stderr, "\n\t%d %d=>%d (%d, %d)[%d]",start, b, score, C->nseq, C->len_aln, step);
       score_matrix[n_val]=score;
       len_matrix [n_val]=b;
       if ( score>best_score)
          {
          best_score=score;
          best_len=b;
	  best_val=n_val;
          }
      }
   free_aln(C);

   for ( a=best_val; a<n_val;a++)
     {
       if (score_matrix[a]>best_score/2)best_len=len_matrix[a];
       else break;
     }
   vfree ( score_matrix);
   vfree ( len_matrix);
   
   return best_len;
   }

Alignment *extract_domain_with_coordinates ( Alignment *RESULT,int start, int len, Constraint_list *CL)
{
  int a;
  char *buf;
  Alignment *SEQ_DOMAIN=NULL;
  

 
  

  /*ADJUST THE DIRECTION OF THE DOMAIN: len<0:left and len>0:right*/

  if (len>0);
  else if (len<0)
    {
      len=len*-1;
      start=start-len+1;
    }

  /*CHECK THAT THE BOUNDARY CONDITIONS*/
  
 
  if (start<0 || (!CL->packed_seq_lu && (start+len)>strlen((CL->S)->seq[0])) ||(CL->packed_seq_lu && (start+len)>strlen((CL->S)->seq[(CL->S)->nseq-1])) )return NULL;
  else
    {
      for ( a=start; a< start+len; a++)
	{
	  if ((CL->moca)->forbiden_residues && (CL->moca)->forbiden_residues[0][a+1]==UNDEFINED)
	    {
	      fprintf ( stderr, "*");
	      return NULL;
	    }
	}
    }
 
  /*EXTRACT THE DOMAIN*/
  
  SEQ_DOMAIN=add_seq2aln (CL,SEQ_DOMAIN, CL->S);
  buf=extract_char (SEQ_DOMAIN->seq_al[0], start, len);

  for (a=0; a<len; a++)if ( buf[a]=='X'){free_aln(SEQ_DOMAIN);return NULL;}
  sprintf ( SEQ_DOMAIN->seq_al[0], "%s", buf);
  SEQ_DOMAIN->order[0][1]=start;
  SEQ_DOMAIN=add_seq2aln (CL,SEQ_DOMAIN, CL->S);

 
  
  return SEQ_DOMAIN;
}




int get_starting_point ( Constraint_list *CL)
{
  int a;
  
  
  int l;
 
 
  
 
 
  int **seq;
  int start;
  int *entry=NULL;
  
  l=strlen ( (CL->S)->seq[0]);
  
  seq=declare_int ( l, 2);
  
  
  
   while (entry=extract_entry (CL))
     {
       seq[entry[R1]][1]=entry[R1];
       seq[entry[R2]][1]=entry[R2];
       if ((CL->moca) && (CL->moca)->forbiden_residues && ((CL->moca)->forbiden_residues[0][entry[R1]]==UNDEFINED||(CL->moca)->forbiden_residues[0][entry[R2]]==UNDEFINED ))continue; 
       else
	 {
	   seq[entry[R1]][0]+=entry[MISC];
	   seq[entry[R2]][0]+=entry[MISC];
	 }
     }
   
   sort_int_inv ( seq, 2, 0, 0, l-1);
   fprintf ( stderr, "\nStart=%d %d", seq[0][1], seq[0][0]);
   start=seq[0][1];
   
   
   free_int ( seq, -1);
   return start;
   
   
}


int * analyse_sequence ( Constraint_list *CL)
{
 int a, p;
 int len, start, n_dots;
 int left, right, tw, r, w;
 int best_tw, best_start=0, best_len=0;
 int l;
 int max_len=200;
 

 l=strlen (( CL->S)->seq[0]);

 for ( best_tw=UNDEFINED,start=0; start<l; start++)
    {
      for ( len=10; len< max_len && start+len<l; len++)
	{
	left=start;
	right=start+len;
	for (tw=0, p=0; p<len; p++)
	  {
	    n_dots=CL->residue_index[0][p+start+1][0];
	    
	    for ( a=1; a<n_dots; a+=ICHUNK)
	      {
		
		r=CL->residue_index[0][p+start+1][a+R2];
		w=CL->residue_index[0][p+start+1][a+WE];

		if (r<left || r>right)tw+=w;
	      }
	  }
	
	if ( tw> best_tw || best_tw==UNDEFINED)
	  {
	    best_tw=tw;
	    best_start=start;
	    best_len=len;
	  }
	}
    }
  fprintf ( stderr, "\nStart=%d Len=%d", best_start, best_len);
  return NULL;
}
			  
/*********************************COPYRIGHT NOTICE**********************************/
/*© Centre National de la Recherche Scientifique (CNRS) */
/*and */
/*Cedric Notredame */
/*Fri Aug  8 19:03:27 MDT 2003. */
/*All rights reserved.*/
/*NOTICE:                                                                                                                                     |*/
/*  This file is an integral part of the */
/*  T-COFFEE Software. */
/*  Its content is protected and all */
/*  the conditions mentioned in the licensing */
/*  agreement of the software apply to this file.*/
/*...............................................                                                                                      |*/
/*  If you need some more information, or if you */
/*  wish to obtain a full license, please contact: */
/*  cedric.notredame@europe.com*/
/*...............................................                                                                                                                                     |*/
/**/
/**/
/*	*/
/*********************************COPYRIGHT NOTICE**********************************/

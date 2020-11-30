#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"



/**
 * \file aln_convertion_util.c
 * Contains several auxiliary functions for alignments and templates.
 */
Constraint_list * seq2contacts (Sequence *S, Sequence *T,Constraint_list *CL, char *mode)
{
  //Turns an RNA into a 2D maps (CL) using either Zuker or bona fide PDB
  //Does nothing for proteins

  if (!S || !S->nseq) return CL;
  S=fast_get_sequence_type(S);
  if (!strm (S->type, "RNA"))return CL;//seq2contacts is only supported for RNA
  if (T)
    {
      char *template_file=vtmpnam(NULL);
      output_fasta_seqS (template_file,T);
      S=seq2template_seq(S,template_file,NULL);
    }
  
  cputenv ("SEQ2TEMPLATE4_F_=%s",(mode)?mode:"RNAplfold");
  cputenv ("PDB2TEMPLATE4_F_=no");
  
  if (seq2n_template (S, "_P_")==S->nseq)return CL;
  else if (strm (mode,"RNAplfold") || !mode)S=seq2template_seq (S,"RNA", NULL);
  else if (strm (mode, "join"));
  else printf_exit (EXIT_FAILURE,stderr, "+seq2contact %s: unknown mode [FATAL]", mode);
  
  return read_contact_lib (S,NULL,CL); 
}
Constraint_list * pdb2contacts (Sequence *S, Sequence *T,Constraint_list *CL, char *contact_mode, char *contact_type, float maxD)
{
  char *template_file=vtmpnam(NULL);
  char *lib_name=NULL;
 
  
  if (!T && strm (contact_type, "RNAplfold"));
  else if (!T)
    {
      if (!T)printf_exit (EXIT_FAILURE,stderr, "pdb2contact requires a PDB template  [FATAL]");
    }
  else
    {
      output_fasta_seqS (template_file,T);
      S=seq2template_seq(S,template_file,NULL);
    }
  
  S=fast_get_sequence_type(S);
  

  if (contact_type)
    {
      cputenv ("SEQ2TEMPLATE4_F_=no");
      cputenv ("PDB2TEMPLATE4_F_=%s", contact_type);
    }
  

  if (strm5 (contact_type, "contacts", "best", "count", "closest", "distances") || !contact_type)
    {
      
      lib_name=vtmpnam(NULL);
      pdb2contacts2lib(S,contact_type,maxD,lib_name, contact_mode);
    }
  else if (strm(S->type, "RNA") && strm4 (contact_type, "find_pair","RNAplfold","find_pair-p","x3dna-ssr"))
    {
      S=seq2template_seq (S,"RNA", NULL);
      lib_name=NULL;
    }
  else if (strm4 (contact_type, "find_pair","RNAplfold","find_pair-p","x3dna-ssr"))
    {
      printf_exit (EXIT_FAILURE,stderr, "+pdb2contact %s: is only compatible with RNA [FATAL]", contact_type);
    }
  else
    {
      printf_exit (EXIT_FAILURE,stderr, "+pdb2contact %s: unknown contact_type [FATAL]",contact_type );
    }
  return read_contact_lib (S,lib_name,CL); 
}
  
  
int vienna2template_file (char *outfile, Sequence *R, Sequence *ST)
{
  int a;
  char *seq, *str;
  FILE *fpout;
  char *lib;
  
  fpout=vfopen (outfile,"w");
  seq=vtmpnam(NULL);
  str=vtmpnam(NULL);
  lib=(char*)vcalloc (1000, sizeof (char));
  
  for (a=0; a<ST->nseq; a++)
    {
      int s=name_is_in_list (ST->name[a], R->name, R->nseq, 100);
      
      if (s!=-1)
	{
	  Sequence *STR1;
	  Sequence *RNA1;
	  FILE *fp;
	  
	  fp=vfopen (seq, "w");
	  fprintf (fp, ">%s\n%s\n", R->name[s], R->seq[s]);
	  vfclose (fp);
	  
	  fp=vfopen (str, "w");
	  fprintf (fp, ">%s\n%s\n", ST->name[a], ST->seq[a]);
	  vfclose (fp);

	  RNA1=get_fasta_sequence(seq, NULL);
	  STR1=get_fasta_sequence(str, NULL);
	  
	  sprintf (lib, "%s.fold_tc_lib", R->name[a]);
	  vienna2tc_lib (lib,RNA1, STR1);
	  fprintf (fpout, ">%s _F_ %s\n", R->name[a], lib);
	 
	  free_sequence (RNA1, -1);
	  free_sequence (STR1, -1);
	}
    }
  vfclose (fpout);
  
  
  return 1;
  
}
Constraint_list * vienna2tc_lib (char *out, Sequence *R, Sequence *ST)
{
  int a, nseq;
  int **lu;
  FILE *fp;
  char *outfile;

  outfile=vtmpnam (NULL);
  lu=declare_int (ST->nseq, 2);
  
  
  for (nseq=0,a=0; a<ST->nseq; a++)
    {
      int s=name_is_in_list (ST->name[a], R->name, R->nseq, 100);
      if (s!=-1)
	{
	  lu[nseq][0]=s;
	  lu[nseq][1]=a;
	  nseq++;
	}
    }
  fp=vfopen (outfile, "w");
  
  fprintf (fp, "! TC_LIB_FORMAT_01\n%d\n", nseq);
  for (a=0; a<nseq; a++)
    {
      fprintf (fp,"%s %d %s\n", R->name[lu[a][0]], (int)strlen ( R->seq[lu[a][0]]),  R->seq[lu[a][0]]);
    }
  for (a=0; a<nseq; a++)
    {
      int **list,i;
      fprintf ( fp, "#%d %d\n", a+1, a+1);
      list=vienna2list(ST->seq[lu[a][1]]);
      i=0;
      while (list[i])
	{
	  fprintf (fp, "%d %d %d\n", list[i][0]+1,list[i][1]+1,100);
	  i++;
	}
      free_int (list, -1);
    }
  
  free_int (lu, -1);
  fprintf (fp, "! SEQ_1_TO_N\n");
  vfclose (fp);
 
  if (out && out[0])
    {
      printf_system ("cp %s %s", outfile, out);
      return NULL;
    }
  else
    {
      Constraint_list *L;
      L=declare_constraint_list ( R,NULL, NULL, 0,NULL, NULL);
      return read_constraint_list_file (L, outfile);
    }
}

Alignment *trim_RNA (Alignment *RNA, Sequence *ST, int max)
{
  Alignment *SC;
  int **score;
  int a, ctot=0, tot=0, fraction;
  char *aln;
  FILE *fp;
  
  //evaluate
  
  SC=copy_aln (RNA, NULL);
  SC=sp3_evaluate (SC, ST);
  score=declare_int (RNA->nseq, 2);
  //sort and filter
  for (a=0; a<RNA->nseq; a++)
    {
      score[a][0]=a;
      score[a][1]=SC->score_seq[a];
      tot+=SC->score_seq[a];
    }
  sort_int_inv (score, 2, 1, 0, RNA->nseq-1);
  
  aln=vtmpnam (NULL);
  fp=vfopen (aln, "w");
  for (a=0; a<RNA->nseq; a++)
    {
      int s=score[a][0];
      ctot+=score[a][1];
      
      if (((ctot*100)/tot)<max)
	{
	  HERE ("%s %d", RNA->name[s],SC->score_seq[s]);
	  fprintf (fp, ">%s Score: %d\n%s\n", RNA->name[s],SC->score_seq[s],RNA->seq_al[s]);
	}
      else continue;
    }
 
  vfclose (fp);
  main_read_aln (aln,RNA);
  free_int (score, -1);
  free_aln (SC);
  return RNA;
}


  
Alignment *sp3_evaluate (Alignment *RNA, Sequence *ST)
{
  int a, b, c, d, i;
  int *npairs;

  float **max_res_sc;
  float **tot_res_sc;
  float *max_seq_sc;
  float *tot_seq_sc;
  
  float *max_col_sc;
  float *tot_col_sc;
  
  float tot_sc=0;
  float max_sc=0;
  Alignment *A, *OUT, *IN;

  if (!ST)
    {
      HERE ("sp3_evaluate requires an RNA secondary structure via -in2");
      exit (0);
    }
  
  
  thread_seq_struc2aln (RNA,ST);
  IN=copy_aln (RNA, NULL);
  
  A=copy_aln (IN, NULL);
  OUT=copy_aln (IN, NULL);
  
  max_res_sc=declare_float (A->nseq, A->len_aln);
  tot_res_sc=declare_float (A->nseq, A->len_aln);
  
  max_seq_sc=(float*)vcalloc (A->nseq, sizeof (float));
  tot_seq_sc=(float*)vcalloc (A->nseq, sizeof (float));
  

  max_col_sc=(float*)vcalloc (A->len_aln, sizeof (float));
  tot_col_sc=(float*)vcalloc (A->len_aln, sizeof (float));
  
  for (a=0; a<A->nseq; a++)
    {
      int **list=vienna2list (A->seq_al[a]);
      for (i=0; i<A->len_aln; i++)A->seq_al[a][i]=-1;
      
      i=0;
      while (list[i])
	{
	  A->seq_al[a][list[i][0]]=list[i][1];
	  A->seq_al[a][list[i][1]]=list[i][0];
	  i++;
	}
      free_int (list, -1);
    }
  
  
  
  
  for (a=0; a<A->nseq; a++)
    for (b=0; b<A->nseq; b++)
      {
	for (c=0; c<A->len_aln; c++)
	  {
	    int r1=A->seq_al[a][c];
	    int r2=A->seq_al[b][c];
	    
	    max_res_sc[b][c]+=(r1!=-1)?1:0;
	    tot_res_sc[b][c]+=(r1==r2 && r1!=-1)?1:0;
	    
	    max_col_sc[c]+=((r1+r2)!=-2)?1:0;
	    tot_col_sc[c]+=(r1==r2 && r1!=-1)?1:0;
	    
	    max_seq_sc[b]+=((r1+r2)!=-2)?1:0;
	    tot_seq_sc[b]+=(r1==r2 && r1!=-1)?1:0;
	    
	    max_sc+=((r1+r2)!=-2)?1:0;
	    tot_sc+=(r1==r2 && r1!=-1)?1:0;
	  }	
      }
  for (a=0; a<A->nseq; a++)
    {
      for (c=0; c<A->len_aln; c++)
	{
	  if (A->seq_al[a][c]!=-1)
	    {
	      int r1=(max_res_sc[a][c]==0)?0:(tot_res_sc[a][c]*(float)10/max_res_sc[a][c]);
	      r1=(r1>=10)?9:r1;
	      OUT->seq_al[a][c]=r1+'0';
	    }
	}
      OUT->score_seq[a]=(max_seq_sc[a]==0)?0:(tot_seq_sc[a]*(float)100)/max_seq_sc[a];
    }
  for (c=0; c<A->len_aln; c++)
    {
      int r1=(max_col_sc[c]==0)?0:(tot_col_sc[c]*(float)10)/max_col_sc[c];
      OUT->seq_al[A->nseq][c]=(r1>=10)?9:r1;
    }
 
  OUT->score=OUT->score_aln=(int)(max_sc==0)?0:(tot_sc*1000)/max_sc;



  free_aln (A);
  free_float (tot_res_sc, -1);
  free_float (max_res_sc, -1);
  vfree (tot_seq_sc);vfree (max_seq_sc);
  vfree (tot_col_sc);vfree (max_col_sc);
  copy_aln (OUT, RNA);
  free_aln (IN);
  free_aln (OUT);
  return RNA;
}
		

int aln_has_stockholm_structure (Alignment *A)
{
  return name_is_in_list ("#=GC SS_cons", A->name, A->nseq, 100);
}

int get_aln_stockholm_structure (Alignment *A)
{
  int i;
  if ((i=aln_has_stockholm_structure(A))==-1)
    A=add_alifold2aln (A, NULL);
  return aln_has_stockholm_structure(A);
}
int ** update_RNAfold_list (Alignment *A, int **pos, int s, int **l)
{
  int a=0;
  while (l[a])
    {
      if (!is_gap(A->seq_al[s][l[a][0]]) && !is_gap (A->seq_al[s][l[a][1]]))
	{
	  l[a][2]=pos[s][l[a][0]];
	  l[a][3]=pos[s][l[a][1]];
	}
      else
	{
	  l[a][2]=l[a][3]=-1;
	}
      a++;
    }
  return l;
}

Alignment *compare_RNA_fold ( Alignment *A, Alignment *B)
{
  int i1, i2, i;
  int **l1, **l2;
  int **pos1, **pos2;
  int a, b;
  int tot_ol=0, tot_l=0;

  i1=get_aln_stockholm_structure (A);
  i2=get_aln_stockholm_structure (B);

  l1=vienna2list (A->seq_al[i1]);
  l2=vienna2list (B->seq_al[i2]);

  pos1=aln2pos_simple(A, A->nseq);
  pos2=aln2pos_simple(B, B->nseq);



  for (a=0; a< A->nseq; a++)
    {
      char **lu;
      int ol=0, ll1=0, ll2=0;
      if ( A->name[a][0]=='#')continue;
      i=name_is_in_list (A->name[a], B->name, B->nseq, 100);
      if (i!=-1)
	{
	  l1=update_RNAfold_list (A,pos1,a, l1);
	  l2=update_RNAfold_list (B,pos2,i, l2);
	  lu=declare_char (A->len_aln, B->len_aln);

	  b=0;
	  while (l2[b])
	    {

	      if (l2[b][2]==-1 || l2[b][3]==-1);
	      else
		{
		  ll2++;
		  lu[l2[b][2]][l2[b][3]]=1;

		}
	      b++;
	    }
	  b=0;

	  while (l1[b])
	    {

	      if (l1[b][2]==-1 || l1[b][3]==-1);
	      else
		{
		  ll1++;
		  if (lu[l1[b][2]][l1[b][3]]==1)
		    {
		      A->seq_al[a][l1[b][0]]='6';
		      A->seq_al[a][l1[b][1]]='6';
		      ol++;
		    }
		  else
		    {
		      A->seq_al[a][l1[b][0]]='0';
		      A->seq_al[a][l1[b][1]]='0';
		    }
		}
	      b++;
	    }

	  free_char (lu, -1);
	}
      tot_ol+=ol;
      tot_l+=ll1;
      tot_l+=ll2;
      fprintf ( stdout, "@@ Seq: %s Overalp: %.2f Al1: %.2f Al2: %.2f \n", A->name[a], (float)(ol*200)/(ll1+ll2), (float)(ol*100)/ll1,(float)(ol*100)/ll2);
    }

  fprintf ( stdout, "@@ Seq: Tot Overalp: %.2f \n", (float)(tot_ol*200)/(tot_l));

  return A;
}
int is_neutral(char c1, char c2);
int is_watson (char c1, char c2);
int is_watson2 (char c1, char c2);
int is_watson (char c1, char c2)
{
  c1=tolower (c1);
  c2=tolower (c2);
  if ( is_watson2 (c1, c2)) return 1;
  else return is_watson2 (c2, c1);
}
int is_watson2 (char c1, char c2)
{

  if ( c1=='g' && c2=='c')return 1;
  else if (c1=='a' && (c2=='t' || c2=='u'))return 1;
  return 0;
}
int is_neutral (char c1, char c2)
{

  c1=tolower (c1);
  c2=tolower (c2);
  if (is_watson (c1, c2)) return 1;
  else if (c1=='g' && (c2=='t' || c2=='u'))return 1;
  else if ((c1=='t' || c1=='u') && c2=='g')return 1;
  return 0;
}

int ** vienna2list ( char *seq)
{
  int a, b, i, i2,l;
  int **list;
  l=strlen (seq);
  list=declare_int (l+1, 8);
  for (i=0,a=0; a<l; a++)
    {
      if ( seq[a]=='(')
	{
	  list[i][0]=a;
	  for (i2=0,b=a+1; b<l && i2>=0; b++)
	    {
	      if (seq[b]=='(')i2++;
	      else if (seq[b]==')')i2--;
	    }
	  list[i][1]=b-1;
	  i++;
	}
      else if (seq[a]=='<')
	{
	  list[i][0]=a;
	  for (i2=0,b=a+1; b<l && i2>=0; b++)
	    {
	      if (seq[b]=='<')i2++;
	      else if (seq[b]=='>')i2--;
	    }
	  list[i][1]=b-1;
	  i++;
	}
      else if (seq[a]=='[')
	{
	  list[i][0]=a;
	  for (i2=0,b=a+1; b<l && i2>=0; b++)
	    {
	      if (seq[b]=='[')i2++;
	      else if (seq[b]==']')i2--;
	    }
	  list[i][1]=b-1;
	  i++;
	}
       else if (seq[a]=='{')
	{
	  list[i][0]=a;
	  for (i2=0,b=a+1; b<l && i2>=0; b++)
	    {
	      if (seq[b]=='{')i2++;
	      else if (seq[b]=='}')i2--;
	    }
	  list[i][1]=b-1;
	  i++;
	}
    }

  list[i]=NULL;
  return list;
}
Alignment *aln2alifold(Alignment *A)
{
  char *tmp1;
  char *tmp2;

  print_aln (A);
  tmp1=vtmpnam (NULL);
  tmp2=vtmpnam (NULL);
  output_clustal_aln (tmp1,A);
  printf_system ("RNAalifold %s >%s 2>/dev/null", tmp1, tmp2);
  return alifold2aln (tmp2);
}

Alignment *add_alifold2aln  (Alignment *A, Alignment *ST)
{
  int a,b,c,d,p1,p2;
  int r1, rr1, r2, rr2;
  int watson, comp,tot;
  int **compmat;
  int max, p,k;
  int minseq=3;
  int **list;
  int ncomp=0, nwatson=0;
  int cons_l, fold_l;
  int i,l;

  if (!ST)
    {
      char *tmp1, *tmp2;
      int f;
      Alignment *T;
      T=copy_aln (A, NULL);
      tmp1=vtmpnam (NULL);
      tmp2=vtmpnam (NULL);
      cons_l=A->len_aln;
      for (a=0; a<A->len_aln; a++)
	{
	  for (f=0,b=0; b<A->nseq && f==0; b++)
	    {
	      if (is_gap (A->seq_al[b][a]))f=1;

	    }
	  if (f)
	    {
	      cons_l--;
	      for (b=0; b<A->nseq; b++)T->seq_al[b][a]='-';
	    }
	}
      ST=aln2alifold (T);
    }
  //add or Replace the structure
  l=strlen (ST->seq_al[1]);
  for (a=0; a< l; a++)if (ST->seq_al[1][a]==STOCKHOLM_CHAR)ST->seq_al[1][a]='.';
  if ((i=name_is_in_list ("#=GC SS_cons", A->name, A->nseq, 100))!=-1)
    {
       sprintf (A->seq_al[i], "%s", ST->seq_al[1]);
    }
  else
    {
      A=realloc_aln2 ( A, A->nseq+1, A->len_aln+1);
      sprintf (A->name[A->nseq], "#=GC SS_cons");
      sprintf (A->seq_al[A->nseq], "%s", ST->seq_al[1]);
      A->nseq++;
    }
  return A;
}
Alignment * alifold2analyze (Alignment *A, Alignment *ST, char *mode)
{
  int s;
  int **list;
  int usegap;

  s=name_is_in_list ("#=GC SS_cons", A->name,A->nseq, 100);

  if (s==-1)
    {
      A=add_alifold2aln (A,ST);
      s=name_is_in_list ("#=GC SS_cons", A->name,A->nseq, 100);
    }

  list=vienna2list (A->seq_al[s]);
  list=alifold_list2cov_list (A, list);

  usegap=0; //do not use gaped positions by default
  if (mode && strstr (mode, "usegap"))usegap=1;//count positions with gaps

  if (!mode)
    {
      A=alifold2cov_stat   (A, list,usegap);
    }
  else
    {
      if ( strstr (mode, "stat"))  A=alifold2cov_stat   (A, list, usegap);
      if ( strstr (mode, "list"))  A=alifold2cov_list   (A, list, usegap);
      if ( strstr (mode, "aln"))   A=alifold2cov_aln    (A, list, usegap);
      if ( strstr (mode, "color") )
	{
	  Alignment *C;
	  C=copy_aln (A, NULL);
	  C=alifold2cov_cache (C, list, usegap);
	  A=alifold2cov_aln (A, list, usegap);
	  if ( strstr ( mode, "ps"))
	    output_color_ps (A, C, "stdout");
	  else
	    output_color_html (A, C, "stdout");
	  myexit (EXIT_SUCCESS);
	}
    }
  return A;
}


int **    alifold_list2cov_list (Alignment *A, int **list)
{
   int a,b,c,d,p1,p2,s;
  int r1, rr1, r2, rr2;
  int neutral,watson, comp,tot, occupancy;
  int **compmat;
  int max, p,k;
  int minseq=3;

  int ncomp=0, nwatson=0, nneutral=0, ncomp_wc=0;
  int cons_l, fold_l;
  int nseq;


  for (nseq=0,a=0; a< A->nseq; a++)if ( A->name[a][0]!='#')nseq++;
  max=((nseq*(nseq-1))/2);
  a=0;
  while (list[a])
    {
      p1=list[a][0];
      p2=list[a][1];
      watson=0;
      comp=0;
      neutral=0;
      tot=0;
      occupancy=0;
      for (c=0; c<A->nseq-1; c++)
	{
	  if (A->name[c][0]=='#')continue;
	  r1=tolower(A->seq_al[c][p1]);
	  r2=tolower(A->seq_al[c][p2]);
	  if (is_gap(r1) || is_gap(r2))continue;
	  for (d=c+1; d<A->nseq; d++)
	    {
	      if (A->name[d][0]=='#')continue;
	      rr1=tolower(A->seq_al[d][p1]);
	      rr2=tolower(A->seq_al[d][p2]);
	      if (is_gap(rr1) || is_gap(rr2))continue;
	      if (is_watson (r1, r2))watson++;
	      if (is_watson (rr1, rr2))watson++;
	      if (is_neutral (r1, r2))neutral++;
	      if (is_neutral (rr1, rr2))neutral++;
	      if (r1!=rr1 && r2!=rr2)comp++;
	      occupancy++;
	    }
	}
      if (occupancy==0)
        {
	  a++;
          continue;
        }
      watson=(watson*100)/(occupancy*2);
      comp=(comp*100)/occupancy;
      neutral=(neutral*100)/(occupancy*2);
      occupancy=(occupancy*100)/max;
      list[a][3]=neutral;
      list[a][4]=watson;
      list[a][5]=comp;
      list[a][6]=occupancy;

      if (list[a][3]<100)list[a][7]='I';//incompatible pair
      else
	{
	  list[a][7]='N';//Neutral pair
	  if (list[a][4]==100)
	    {
	      list[a][7]='W';//Watson and Crick
	      if ( list[a][5]>0)list[a][7]='C'; //Watson and crick compensated
	    }
	  else if ( list[a][5]>0)
	    {
	      list[a][7]='c';//compensated
	    }
	}
      a++;
    }

  return list;
}
Alignment *alifold2cov_aln (Alignment *inA,int **list, int ug)
{
  int a=0;
  a=0;
  Alignment *A;

  A=copy_aln (inA, NULL);
  A=realloc_aln2 ( A, A->nseq+1, A->len_aln+1);
  sprintf (A->name[A->nseq], "#=GC SS_analyze");
  sprintf (A->seq_al[A->nseq], "%s", A->seq_al[A->nseq-1]);
  A->nseq++;
  while (list[a])
    {
      char s;
      if (list[a][6]<100 && !ug);
      else
	{
	  s=list[a][7];
	  A->seq_al[A->nseq-1][list[a][0]]=s;
	  A->seq_al[A->nseq-1][list[a][1]]=s;
	}
      a++;
    }
  return A;
}
Alignment *alifold2cov_stat (Alignment *A,int **list, int ug)
{
  int fold=0,watson=0, comp=0, compwc=0, incomp=0, neutral=0;
  int a;

  a=0;
  while (list[a])
    {
      int s;
      fold++;
      if (list[a][6]<100 && !ug);
      else
	{
	  s=list[a][7];
	  watson +=(s=='W')?1:0;
	  compwc +=(s=='C')?1:0;
	  comp   +=(s=='c')?1:0;
	  neutral+=(s=='N')?1:0;
	  incomp +=(s=='I')?1:0;
	}
      a++;
    }
  fprintf ( stdout, "@@ TOT Nseq:%d tot_len: %d  fold: %d neutral: %d watson: %d CorWC: %d cor: %d Incompatible: %d\n",A->nseq-1, A->len_aln,fold, neutral,watson, compwc,comp,incomp);
  return A;
}
Alignment *alifold2cov_cache (Alignment *inA, int **list, int ug)
{
  int a,b, c;
  Alignment *A;

  A=copy_aln (inA, NULL);
  a=0;
  while (list[a])
    {
      int v, s;
      if (list[a][6]<100 && !ug);
      else
	{
	  s=list[a][7];
	  if (s=='C')v=9; //red
	  else if ( s=='c')v=7; //orange
	  else if ( s=='W')v=5; //Yellow
	  else if ( s=='N')v=2; //green
	  else if ( s=='I')v=0; //blue;
	  for (b=0;b<A->nseq; b++)
	    {
	      if (A->name[b][0]=='#');
	      else
		{
		  for (c=0; c<2; c++)
		    {
		      A->seq_al[b][list[a][c]]='0'+v;
		    }
		}
	    }
	}
      a++;
    }
  return A;
}

Alignment *alifold2cov_list (Alignment *A,int **list, int ug)
{
  int a,b, s;

  a=0;
  while (list[a])
    {
      s=list[a][7];
      if (list[a][6]<100 && !ug);
      else if (s=='C')
	{
	  fprintf ( stdout, "@@ WC Compensated pair: %4d %4d =>", list[a][0]+1, list [a][1]+1);
	  for (b=0; b<A->nseq; b++)if (A->name[b][0]!='#')fprintf ( stdout, "[%c%c]", toupper (A->seq_al[b][list[a][0]]), toupper(A->seq_al[b][list[a][1]]));
	  fprintf (stdout,"\n");
	}
      else if (s=='c')
	{
	  fprintf ( stdout, "@@ Neural Compensated pair: %4d %4d =>", list[a][0]+1, list [a][1]+1);
	  for (b=0; b<A->nseq; b++)if (A->name[b][0]!='#')fprintf ( stdout, "[%c%c]", toupper (A->seq_al[b][list[a][0]]), toupper(A->seq_al[b][list[a][1]]));
	  fprintf (stdout,"\n");
	}
      else if (s=='W')
	{
	  fprintf ( stdout, "@@ WC pair: %4d %4d =>", list[a][0]+1, list [a][1]+1);
	  for (b=0; b<A->nseq; b++)if (A->name[b][0]!='#')fprintf ( stdout, "[%c%c]", toupper (A->seq_al[b][list[a][0]]), toupper(A->seq_al[b][list[a][1]]));
	  fprintf (stdout,"\n");
	}
      else if (s=='N')
	 {
	  fprintf ( stdout, "@@ Neutral pair: %4d %4d =>", list[a][0]+1, list [a][1]+1);
	  for (b=0; b<A->nseq; b++)if (A->name[b][0]!='#')fprintf ( stdout, "[%c%c]", toupper (A->seq_al[b][list[a][0]]), toupper(A->seq_al[b][list[a][1]]));
	  fprintf (stdout,"\n");
	}
      else if (s=='I')
	{
	  fprintf ( stdout, "@@ incompatible pair: %4d %4d =>", list[a][0]+1, list [a][1]+1);
	  for (b=0; b<A->nseq; b++)if (A->name[b][0]!='#')fprintf ( stdout, "[%c%c]", toupper (A->seq_al[b][list[a][0]]), toupper(A->seq_al[b][list[a][1]]));
	  fprintf (stdout,"\n");
	}
      a++;
    }

  return A;
}


Alignment *aln2sample (Alignment *A, int n)
{
  Alignment *B;
  int a, b, p;
  int **pos;

  B=copy_aln (A, NULL);

  vsrand(0);

  pos=declare_int (A->len_aln, 2);
  for (a=0; a<A->len_aln; a++){pos[a][0]=a;pos[a][1]=rand()%(1000*A->len_aln);}

  sort_int (pos, 2, 1, 0, A->len_aln-1);

  n=(n==0)?A->len_aln:(MIN (n, (A->len_aln)));
  for (a=0; a<n; a++)
    for (b=0; b<A->nseq; b++)
      A->seq_al[b][a]=B->seq_al[b][pos[a][0]];
  for (b=0; b<A->nseq; b++)
    A->seq_al[b][n]='\0';
  A->len_aln=n;

  free_aln (B);
  free_int (pos, -1);
  return A;
}
Alignment *aln2bootstrap (Alignment *A, int n)
{
  Alignment *B;
  int a, b, p;

  if (n==0)n=A->len_aln;
  else A=realloc_aln (A, n+1);
  vsrand(0);
  B=copy_aln (A, NULL);
  for (a=0; a<n; a++)
    {
      p=rand ()%A->len_aln;
      for (b=0; b<A->nseq; b++)
	A->seq_al[b][a]=B->seq_al[b][p];
    }
  for ( b=0; b<A->nseq; b++)A->seq_al[b][n]='\0';
  A->len_aln=n;

  free_aln (B);
  return A;

}


Alignment * aln2random_aln (Alignment *A, char *smode)

{
  int a, b, n, **res;
  int max;



  if ( smode==NULL)
    {
      smode=(char*)vcalloc (4, sizeof (char));
      sprintf ( smode, "SCR");//Sequences, Column Residues
    }
  else if ( strm (smode, "NO"))return A;

  vsrand(0);
  max=A->nseq*1000;

  if ( strstr ( smode, "S"))
    {
      A=aln2scramble_seq (A);
    }
  if ( strstr ( smode, "C"))
    {

      res=declare_int (A->nseq, 2);
      for (a=0; a< A->len_aln; a++)
	  {
	    for (n=0,b=0;b<A->nseq; b++)
	      {
		if ( !is_gap(A->seq_al[b][a]))
		  {
		    res[n][0]=A->seq_al[b][a];
		    res[n][1]=rand()%max;
		    n++;
		  }
		sort_int (res, 2, 1, 0, n-1);
	      }
	    for (n=0,b=0;b<A->nseq; b++)
	      {
		if ( !is_gap(A->seq_al[b][a]))A->seq_al[b][a]=res[n++][0];
	      }
	  }
      free_int (res, -a);
    }


  //Redistributes the residues randomly without changing the gap pattern
  if ( strstr ( smode, "R"))
    {
      max=A->len_aln*A->nseq;
      res=declare_int (max, 2);

      for (n=0,a=0; a< A->len_aln; a++)
	{
	  for (b=0;b<A->nseq; b++)
	    {
	      if ( !is_gap(A->seq_al[b][a]))
		{
		  res[n][0]=A->seq_al[b][a];
		  res[n][1]=rand()%max;
		  n++;
		}

	    }
	}
      sort_int (res, 2, 1, 0, n-1);
      for (n=0,a=0; a< A->len_aln; a++)
	{
	  for (b=0;b<A->nseq; b++)
	    {
	      if ( !is_gap(A->seq_al[b][a]))
		{
		  A->seq_al[b][a]=res[n++][0];
		}

	    }
	}

      free_int (res, -1);
    }

  return A;
}
Alignment *score_aln2score_ascii_aln (Alignment *A, Alignment *C)
{
  //Convert the output of T-Coffee evaluate into a printable score_ascii alignment*/
  //A and C must be sorted
  //sets to 0 lone residues
  int a, b;

  for (a=0; a<A->nseq; a++)
    for (b=0; b<A->len_aln; b++)
      {

	int rC=C->seq_al[a][b];
	int rA=A->seq_al[a][b];
	if ( !strm (A->name[a], C->name[a])){HERE ("Unsorted aln in score_aln2score_ascii"); myexit (EXIT_FAILURE);}

	if ( rA=='x' || rA=='X')C->seq_al[a][b]='9';
	else if ( rC >='0' && rC<='9');
	else if ( rC<10)C->seq_al[a][b]='0'+rC;
	else if ( rC==NO_COLOR_RESIDUE && !is_gap(rA)) C->seq_al[a][b]='0';
	else if ( rC==NO_COLOR_RESIDUE && is_gap(rA))C->seq_al[a][b]='-';
      }
  return C;
}
Alignment*aln2gap_cache (Alignment *A, int val)
{
  Alignment *B;
  int a, b, c, nr;

  B=copy_aln (A, NULL);
  for (b=0; b<A->len_aln; b++)
    {
      for (nr=0,a=0; a<A->nseq; a++)nr+=!is_gap (A->seq_al[a][b]);
      for (a=0; a<A->nseq; a++)if (!is_gap(A->seq_al[a][b]))B->seq_al[a][b]=(nr==1)?'0'+val:'1';
    }
  return B;
}

Alignment* aln2case_aln (Alignment *B, char *upper, char *lower)
{
  int a, b, c, up, lo;
  Alignment *A;

  A=copy_aln (B, NULL);

  up=(upper)?upper[0]:'u';
  lo=(lower)?lower[0]:'l';

  for (a=0; a<A->nseq; a++)
    for (b=0; b<A->len_aln; b++)
      {
	c=A->seq_al[a][b];

	if ( is_gap(c));
	else A->seq_al[a][b]=(isupper (c))?up:lo;
      }
  return A;
}
Alignment *aln2scale (Alignment *A, char *coffset)
{
  int a, b, t, v, n;
  char *s1, *s2;
  char s[1000];
  int offset;

  if (coffset)offset=atoi(coffset);
  else offset=0;

  sprintf (s, "%d", A->len_aln+offset);
  n=strlen (s);

  A=realloc_aln2 (A, A->nseq+n, A->len_aln+1);
  s1=(char*)vcalloc ( n+1, sizeof (char));
  s2=(char*)vcalloc ( n+1, sizeof (char));

  for (a=0; a<n; a++)
    {
      if (a==0)s2[a]='1';
      else strcat (s2, "0");
      sprintf (A->name[A->nseq+a], "%s", s2);
    }

  for (a=0; a<A->len_aln; a++)
    {
      sprintf (s1, "%d", a+1+offset);
      s2=invert_string (s1);
      t=strlen (s2);

      for (b=0; b<=n; b++)
	{
	  if (b>=t) v='0';
	  else v=s2[b];

	  A->seq_al[A->nseq+b][a]=v;
	}
    }

  A->nseq+=n;
  return A;
}




int * pos2list (int * pos, int len, int *nl)
{
  int *list;
  int a;
  nl[0]=0;
  list=(int*)vcalloc (len, sizeof (int));
  for (a=0; a<len; a++)if (pos[a])list[nl[0]++]=a;
  return list;
}
int *list2pos (int *list, int nl, int len)
{
  int *pos, a;
  pos=(int*)vcalloc (len, sizeof (int));
  for (a=0; a<nl; a++)pos[list[a]]=1;
  return pos;
}

int **aln2resindex ( Alignment *A, Alignment *B, FILE *fp)
{
  int *list, **pos;
  int a, b, n, s;


  list=(int*)vcalloc (A->nseq+((B)?B->nseq:0), sizeof (int));
  pos=aln2pos_simple_2 (A);
  if (B)
    {
      n=B->nseq;
      for ( a=0; a<B->nseq; a++)
	{
	  list[a]=name_is_in_list(B->name[a], A->name, A->nseq, 100);
	}
    }
  else
    {
      for ( a=0; a<A->nseq; a++)
	list[a]=a;
      n=A->nseq;
    }


  fprintf ( fp, "#");
  for ( b=0; b<n; b++)
    {
      s=list[b];
      if ( s!=-1)fprintf (fp, " %s",A->name[s]);
    }
  fprintf (fp, "\n");

  for ( a=0; a<A->len_aln; a++)
    {
      for ( b=0; b<n; b++)
	{
	  s=list[b];
	  if ( s==-1);
	  else if (pos[s][a]<0)
	    fprintf (fp, "%4d", -1);
	  else
	    fprintf (fp, "%4d", pos[s][a]);
	}
      fprintf (fp, "\n");
    }
  return pos;
}

int **index_seq_res      ( Sequence *S1, Sequence *S2, int **name_index)
{
  /*Index the residues of S1 according to S2
    index[seq1 of S1][z]->x, where x is the position of residue z of seq1/S1 in S2->seq[index[Seq1/S1]]
  */
  int a;
  int **index;
  char *seq1=NULL, *seq2=NULL;
  Alignment *Profile;

  index=(int**)vcalloc ( S1->nseq, sizeof (int*));

  for (a=0; a< S1->nseq; a++)
    {
      int  len1, len2, b, c;

      seq1=S1->seq[a];

      if (name_index[a][0]==-1)
	seq2=NULL;
      else if (name_index[a][1]==-1)
	{
	  seq2=S2->seq[name_index[a][0]];
	}
      else if ((Profile=seq2R_template_profile (S2, name_index[a][0])) !=NULL)
	{
	  seq2=Profile->seq_al[name_index[a][1]];
	}

      len1=(seq1)?strlen (seq1):0;
      len2=(seq2)?strlen (seq2):0;
      index[a]=(int*)vcalloc (len2, sizeof(int));


      for (c=0,b=0; b<len2; b++)if( !is_gap(seq2[b]))index[a][c++]=b;
      //index[a]=get_res_index ( seq1, seq2);
    }
  return index;
}

int **index_seq_name ( Sequence *S1, Sequence *S2)
{
  /*Index the names of S1 according to S2
    index[seq1 of S1][0]->x  if seq1 is the xth sequence of S2
                        ->-1 if seq1 is nowhere to be found
    index[seq1 of S1][1]->z if seq1 is the zth sequence within the xth profile of S2
  */
  int **index;
  int a, b, x, z;
  Alignment *Profile;
  index=declare_int (S1->nseq, 2);


  for ( a=0; a<S1->nseq; a++)
    {
      index[a][0]=index[a][1]=-1;
      x=name_is_in_list (S1->name[a],S2->name,S2->nseq,100);
      if ( x!=-1){index[a][0]=x;index[a][1]=-1;}
      for ( b=0; b<S2->nseq; b++)
	    {
	      if ((Profile=seq2R_template_profile (S2,b)))
		{
		  z=name_is_in_list (S1->name[a],Profile->name,Profile->nseq,100);
		  if ( z!=-1){index[a][0]=b;index[a][1]=z;b=S2->nseq;}
		}
	    }
    }
  return index;
}




int *get_name_index (char **l1, int n1, char **l2, int n2)
{
  int *r;
  int a;
  /*return Array[Index_L1]=Index_L2 */
  r=(int*)vcalloc ( n1, sizeof (int));
  for ( a=0; a< n1; a++)
    r[a]=name_is_in_list (l1[a],l2,n2,100);
  return r;
}

int* get_res_index (char *seq0, char *seq1)
{
  int *coor, a;

  if ( !seq0 || !seq1) return NULL;


  coor=(int*)vcalloc ( strlen (seq0)+1, sizeof (int));
  if (!strm (seq0, seq1))
    {
      int r0, r1 , isr0, isr1;
      int l0=0, l1=0;
      Alignment *A;
      A=align_two_sequences (seq0,seq1,"pam250mt",-5,-1, "myers_miller_pair_wise");

      for ( a=0; a< A->len_aln; a++)
	{
	  r0=A->seq_al[0][a];r1=A->seq_al[1][a];
	  isr0=!is_gap(r0);
	  isr1=!is_gap(r1);
	  l0+= isr0;
	  l1+= isr1;
	  if (isr0 && isr1)coor[l0-1]=l1-1;
	  else if (isr0)  coor[l0-1]=-1;
	}
      free_aln (A);
    }
  else
    {
      int l0;

      l0=strlen (seq0);
      for ( a=0;a< l0; a++)
	coor[a]=a;
    }

  return coor;
}

int change_residue_coordinate ( char *in_seq1, char *in_seq2, int v)
{
  /*Expresses the coordinate of a residue in seq1, in the coordinate system of seq2*/


  static char *seq1, *seq2;
  static int *coor;


  if ( seq1 !=in_seq1 || seq2 !=in_seq2)
    {
      int r0, r1 , isr0, isr1;
      int l0=0, l1=0;
      Alignment *A;
      int a;

      vfree (coor);
      seq1=in_seq1, seq2=in_seq2;
      A=align_two_sequences (seq1,seq2,"pam250mt", -14, -2, "myers_miller_pair_wise");

      coor=(int*)vcalloc ( A->len_aln, sizeof (int));
      for ( a=0; a< A->len_aln; a++)
	{
	  r0=A->seq_al[0][a];r1=A->seq_al[1][a];

	  isr0=!is_gap(r0);
	  isr1=!is_gap(r1);
	  l0+= isr0;
	  l1+= isr1;

	  if (isr0 && isr1)coor[l0-1]=l1-1;
	  else if (isr0)  coor[l0-1]=-1;
	}
      free_aln (A);
    }
  return coor[v];
}


int ** minimise_repeat_coor (int **coor, int nseq, Sequence *S)
    {
    int **new_coor;
    int a, min;
    new_coor=declare_int ( nseq, 3);
    min=return_min_int (coor, nseq, 2);
    for ( a=0; a< nseq; a++)
        {
	new_coor[a][0]=coor[a][0];
	new_coor[a][1]=coor[a][1];
	new_coor[a][2]=min;
	}
    return new_coor;
    }
int ** get_nol_seq ( Constraint_list *CL, int **coor, int nseq, Sequence *S)
    {
    int a, s, p, l, nl;
    int **buf;
    int **new_coor;

    new_coor=declare_int ( nseq+1, 3);


    buf=get_undefined_list ( CL);



    for ( a=0; a< nseq; a++)buf[coor[a][0]][coor[a][1]]=1;


    for ( a=0; a< nseq; a++)
        {
	s=coor[a][0];
	p=coor[a][1]+1;
	l=strlen(S->seq[s]);
	nl=0;
	while ( p<=l && !buf[s][p++])nl++;
	new_coor[a][0]=s;
	new_coor[a][1]=coor[a][1];
	new_coor[a][2]=nl;
	}
    free_int ( buf, -1);
    return new_coor;
    }



int compare_pos_column( int **pos1,int p1, int **pos2,int p2, int nseq)
    {
    int a,v1, v2;
    int identical=0;

    for ( a=0; a< nseq; a++)
	{

	v1=pos1[a][p1];
	v2=pos2[a][p2];

	if (v1>0 || v2>0)
	    {
	    if ( v1!=v2)return 0;
	    else identical=1;
	    }
	}

    return identical;
    }


char *seq2alphabet (Sequence *S)
{
  return array2alphabet (S->seq, S->nseq, "");
}

double *seq2tetraa (char *seq, double *v)
{
  static int ****diaa;
  int tot=0;
  int l=strlen (seq);
  int a,a1,a2,a3,a4,b,c,d,e;


  if (!diaa)
    {
      diaa=(int****)vcalloc ( 26, sizeof(int***));
      for (a=0; a<26; a++)
	{
	  diaa[a]=(int***)vcalloc (26, sizeof (int**));
	  for (b=0; b<26; b++)
	    {
	      diaa[a][b]=(int**)vcalloc (26, sizeof(int*));
	      for (c=0; c<26; c++)
		diaa[a][b][c]=(int*)vcalloc (26, sizeof(int));
	    }
	}
    }

  for (a=0; a<l-3;a++)
    {
      a1=tolower(seq[a])-'a';
      a2=tolower(seq[a+1])-'a';
      a3=tolower(seq[a+2])-'a';
      a4=tolower(seq[a+3])-'a';

      if (a1<0 || a2<0 || a3<0 || a4<0)continue;
      diaa[a1][a2][a3][a4]++;
      tot++;
    }

  if (!v)v=(double*)malloc (26*26*26*26*sizeof (double));
  for (e=0,a=0; a<26; a++)
    for (b=0; b<26; b++)
      for (c=0; c<26; c++)
	for (d=0;d<26; d++,e++)
	{
	  v[e]=(double)diaa[a][b][c][d]/(double)tot;
	  diaa[a][b][c][d]=0;
	}

  return v;
}
double *seq2swv   (char *seq, char **seql, int n)
{
  double *v=(double*)vcalloc (n+10, sizeof (double));
  int a;
  static int **matrix=read_matrice ("blosum62mt");
  for (a=0; a<n; a++)
    {
      int l=strlen (seql[a]);
      int s =seq2swl (seq, seql[a], -4, -1, matrix);
      v[a]=(double)s/(double)l;
    }
  return v;
}
double *seq2swr   (char *seq, char **seql, int n)
{
  double *v=(double*)vcalloc (n+10, sizeof (double));
  int l1=strlen (seq);
  int a,l;
  static int **matrix=read_matrice ("blosum62mt");
  for (a=0; a<n; a++)
    {
      v[a]=(double)seq2swl (seq, seql[a], -4, -1, matrix);
    }
  return v;
}
            
double *seq2triaa (char *seq, double *v)
{
  static int ***diaa;
  int tot=0;
  int l=strlen (seq);
  int a,a1,a2,a3,b,c,d;


  if (!diaa)
    {
      diaa=(int***)vcalloc ( 26, sizeof(int**));
      for (a=0; a<26; a++)
	{
	  diaa[a]=(int**)vcalloc (26, sizeof (int*));
	  for (b=0; b<26; b++)
	    diaa[a][b]=(int*)vcalloc (26, sizeof(int));
	}
    }
  for (a=0; a<l-2;a++)
    {
      a1=tolower(seq[a])-'a';
      a2=tolower(seq[a+1])-'a';
      a3=tolower(seq[a+2])-'a';
      if (a1<0 || a2<0 || a3<0)continue;
      diaa[a1][a2][a3]++;
      tot++;
    }

  if (!v)v=(double*)malloc (26*26*26*sizeof (double));
  for (d=0,a=0; a<26; a++)
    for (b=0; b<26; b++)
      for (c=0; c<26; c++,d++)
	{
	  v[d]=(double)diaa[a][b][c]/(double)tot;
	  diaa[a][b][c]=0;
	}
  return v;
}
double *seq2diaa (char *seq, double *v)
{
  static int **diaa;
  int tot=0;
  int l=strlen (seq);
  int a,a1,a2,b,c;

  if (!diaa)diaa=declare_int (256, 256);
  for (a=0; a<l-1;a++)
    {
      a1=tolower(seq[a])-'a';
      a2=tolower(seq[a+1])-'a';
      if (a1<0 || a2<0)continue;
      diaa[a1][a2]++;
      tot++;
    }

  if (!v)v=(double*)malloc (26*26*sizeof (double));
  for (c=0,a=0; a<26; a++)
    for (b=0; b<26; b++, c++)
      {
	v[c]=(double)diaa[a][b]/(double)tot;
	diaa[a][b]=0;
      }

  return v;
}
char *aln2alphabet (Alignment *A)
{
  return array2alphabet (A->seq_al, A->nseq, "");
}

char *array2alphabet (char **array, int n, char *forbiden)
{
  int a, b, l;
  int *hasch;
  char *alphabet;

  hasch=(int*)vcalloc (256, sizeof (int));
  alphabet=(char*)vcalloc ( 257, sizeof (char));


  for ( a=0; a<n; a++)
    {
      l=strlen (array[a]);
      for ( b=0; b<l; b++)
	hasch[tolower(array[a][b])]++;
    }

  for ( a=0, b=0; a< 256; a++)
    {
      if (hasch[a] && !strrchr(forbiden,a))alphabet[b++]=a;
    }

  alphabet[b]='\0';
  vfree (hasch);
  return alphabet;
}


//***************************************************************
//
  //                          TM PRED
//***************************************************************

char* alnpos2hmmtop_pred (Alignment *A,Alignment *Pred, int pos, int mode)
{
  static char *result;
  static Alignment *Cache;
  static int *score;
  int a, tot, cons;

  if (!score)
    {
      score=(int*)vcalloc (256, sizeof (int));
      result=(char*)vcalloc (100, sizeof (char));
    }

  if (!Pred && !Cache)
    {
      Cache=aln2hmmtop_pred (A);
    }
  if (!Pred) Pred=Cache;


  for (tot=0,a=0; a<A->nseq; a++)
    {
      char s;
      s=Pred->seq_al[a][pos];
      if (!is_gap(s))
	{
	  score[tolower(s)]++;
	  tot++;
	}
    }

  if ( score['h']>score['i'] && score['h']>score['o'])cons='h';

  else if ( score['i']>score['o'])cons='i';
  else cons='o';
  if (tot==0) return "";


  if (mode==VERBOSE)sprintf (result, " H: %3d I: %3d O: %3d P: %c", (score['h']*100)/tot, (score['i']*100)/tot, (score['o']*100)/tot, cons);
  else if (mode == SHORT)sprintf ( result, "%c", cons);
  score['h']=score['o']=score['i']=0;
  return result;
}


Alignment * aln2hmmtop_pred (Alignment *A)
  {
    int a, b, c;
    char *buf, *pred;
    Alignment *PA;

    PA=copy_aln (A, NULL);
    buf=(char*)vcalloc ( A->len_aln+1, sizeof (char));

    for ( a=0; a< A->nseq; a++)
      {
	sprintf (buf, "%s", A->seq_al[a]);
	pred=seq2tmstruc (buf);
	for (c=0,b=0; b<A->len_aln; b++)
	  {
	    if (!is_gap (PA->seq_al[a][b]))PA->seq_al[a][b]=pred[c++];
	  }
	vfree (pred);
      }
    vfree (buf);
    return PA;
  }

char * seq2tmstruc ( char *seq)
   {
     static Sequence *S;
     char *seqfile, *predfile, *buf;
     FILE *fp;

     seqfile=vtmpnam (NULL);
     predfile=vtmpnam (NULL);

     fp=vfopen (seqfile, "w");
     fprintf ( fp, ">seq1\n%s", seq);
     vfclose (fp);


     printf_system ( "fasta_seq2hmmtop_fasta.pl -in=%s -out=%s -arch=%s/%s -psv=%s/%s", seqfile, predfile, get_mcoffee_4_tcoffee(), "hmmtop.arch", get_mcoffee_4_tcoffee(), "hmmtop.psv");
     S=get_fasta_sequence (predfile, NULL);
     buf=(char*)vcalloc ( strlen (S->seq[0])+1, sizeof (char));
     sprintf ( buf, "%s", S->seq[0]);

     free_sequence (S, S->nseq);

     return buf;
   }

void set_blast_default_values()
{
  set_string_variable ("blast_server", const_cast<char*>( (getenv ("blast_server_4_TCOFFEE"))?getenv ("blast_server_4_TCOFFEE"):"EBI") );
  set_string_variable ("pdb_db", const_cast<char*>( (getenv ("pdb_db_4_TCOFFEE"))?getenv ("pdb_db_4_TCOFFEE"):"pdb") );
  set_string_variable ("prot_db", const_cast<char*>( (getenv ("protein_db_4_TCOFFEE"))?getenv ("protein_db_4_TCOFFEE"):"uniprot") );
				  //Maria added this to cast a const char* to char*
  set_int_variable ("prot_min_sim", 0);
  set_int_variable ("prot_max_sim", 100);

  set_int_variable ("prot_min_cov", 0);
  set_int_variable ("prot_max_cov", 100);

  set_int_variable ("pdb_min_sim", 0);
  set_int_variable ("pdb_max_sim", 100);
  set_int_variable ("pdb_min_cov", 0);
  set_int_variable ("pdb_max_cov", 100);

}

char * seq2pdb   (Sequence *S)
{
  set_blast_default_values();
  S->nseq=1;

  
  S=seq2template_seq (S, "PDB", NULL);
  return seq2P_pdb_id(S,0);
}

char* clean_sname     ( char *in)
{
  int l=0;
  while (in && in[l])
    {
      char c=in[l];
      if (!(isalnum(c) || c=='.' || c=='_'))in[l]='_';
      l++;
    } 
  return in;
}
Sequence * seq2blast ( Sequence *S)
{
  int thread=0;
  int   *pid_list;
  char **pid_tmpfile;
  char **blist;
  int pid, npid, njob, max_nproc,a;
  int max=max_n_pid();
  int nseq,nseq2, nnseq;
  int i, b;
  FILE *fp;
  char *db=NULL;
  char *dbn=NULL;
  char *outdir=NULL;
  int num_iterations, outfmt;
  int compress;
  char *in=NULL;
  char *inz=NULL;
  char *sname=NULL;
  
  if (!S->blastfile)S->blastfile=(char**)vcalloc (S->nseq, sizeof (char*));
  max_nproc=max*2;
  if (getenv("thread_4_TCOFFEE"))thread=atoi(getenv("thread_4_TCOFFEE"));
  else thread=1;
  
  nnseq=MAX((S->nseq/thread),1);
  pid_tmpfile=(char**)vcalloc (max_nproc+1, sizeof (char*));
  pid_list   =(int *)vcalloc (max_nproc, sizeof (int *));
  blist=(char**)vcalloc (max_nproc+1, sizeof (char*));;
  npid=0;
  nseq=0;
  
  
  if (getenv("protein_db_4_TCOFFEE"))db=csprintf (db, "%s", getenv("protein_db_4_TCOFFEE"));
  if (!getenv("protein_db_4_TCOFFEE"))myexit(fprintf_error (stderr,"Database has not been set for Blast"));
  if (!isfile(db))myexit(fprintf_error (stderr,"Database %s not available",db));
  dbn=path2filename(db);
  
  if (getenv("num_iterations_4_TCOFFEE"))num_iterations=atoi(getenv("num_iterations_4_TCOFFEE"));
  else num_iterations=1;
  
  if (getenv("outfmt_4_TCOFFEE"))outfmt=atoi(getenv("outfmt_4_TCOFFEE"));
  else outfmt=5;
  
  if (getenv("cache_4_TCOFFEE"))
    {
      outdir=csprintf (outdir, "%s", getenv("cache_4_TCOFFEE"));
      if (outdir[0]!='/')outdir=csprintf (outdir, "./%s", getenv("cache_4_TCOFFEE"));
    }
  else  outdir=csprintf (outdir, "./");
  my_mkdir (outdir);
      
  if (getenv("compress_4_TCOFFEE"))compress=1;
  else compress=0;

  fprintf ( stderr, "\n!Run Blast in batch -- db: %s number_iterations: %d outdir: %s outfmt: %d compression: %s\n", db, num_iterations, outdir, outfmt, (compress)?"yes":"no"); 
  for (b=0,a=0; a<S->nseq; a++)
    {
      sname=clean_sname(csprintf (sname, "%s", S->name[a]));
      in=csprintf (in,"%s/%s.blastp.%s.LOCAL.%d.tmp", outdir,sname,dbn,num_iterations);
      inz=csprintf (inz, "%s.gz", in);
      
      if ( isfile (in) || isfile (inz))
	{
	  fprintf ( stderr, "!\tLOOKUP BLAST %s vs %s\n", S->name[a], dbn); 
	  b++;
	}
      else 
	{
	  fprintf ( stderr, "!\tCOMPUTE BLAST %s vs %s\n", S->name[a], dbn); 
	}
    }
  
    

  npid=0;
  while (nseq<S->nseq)
    {
      int nseq2;
      char *tmp=vtmpnam (NULL);
      
      blist[npid]=vtmpnam (NULL);
      FILE *fp=vfopen  (tmp, "w");
      for (nseq2=0,a=0; a<nnseq && nseq<S->nseq; a++, nseq++, nseq2++)
	{
	  
	  fprintf ( fp, ">%s\n%s\n", S->name[nseq], S->seq[nseq]);
	}
     
      vfclose (fp);
      pid=vvfork(NULL);
      if ( pid==0)
	{
	  int b;
	  initiate_vtmpnam(NULL);
	  vfclose(vfopen (pid_tmpfile[a],"w"));
	  S=seq2blast_thread(get_fasta_sequence (tmp, NULL));
	  fp=vfopen (blist[npid], "w");
	  for (b=0; b<nseq2; b++)
	      fprintf (fp,">%s %s\n%s\n", S->name[b], S->blastfile[b], S->seq[b]);
	  vfclose (fp);

	  myexit (EXIT_SUCCESS);
	}
      else
	{
	  pid_list[pid]=npid;
	  npid++;
	}
      
    }
  for (a=0; a<npid; a++)
    {
      pid=vwait(NULL);
      vremove(pid_tmpfile[pid_list[pid]]);
    }
  for (a=0; a<npid; a++)
    {
      int b;
      Sequence *S2=get_fasta_sequence (blist[a], NULL);
      for (b=0; b<S2->nseq; b++)
	{
	  
	  if ((i=name_is_in_hlist (S2->name[b], S->name,-1))!=-1)
	    {
	      S->blastfile[i]=csprintf (S->blastfile[i], "%s", S2->seq_comment[b]);
	      S->seq_comment[i] =strcatf (S->seq_comment[i], " _BLAST_ %s", S2->seq_comment[b]);
	      S->aln_comment[i] =strcatf (S->aln_comment[i], " _BLAST_ %s", S2->seq_comment[b]);
	    }

	}
      free_sequence(S2, -1);
    }
  
  vfree (pid_list);
  free_char (pid_tmpfile, -1);
  return S;
}
  
Sequence * seq2blast_thread (Sequence *S)
{
  char *db=NULL;
  char *dbn=NULL;
  int num_iterations=1;
  int outfmt=5;
  char *FullBlastF=vtmpnam (NULL);
  char *FullSeqF=vtmpnam (NULL);
  char *sname=NULL;
  char *outdir=NULL;
  char *buf=NULL;
  FILE *fp;
  FILE *out;
  char *header=NULL;
  char *lheader=NULL;
  char *ilen=NULL;
  char *olen=NULL;
  char *in=NULL;
  char *inz=NULL;
  int  *dolist;
  char *footer=NULL;
  int a,b;
  int compress=0;
  
  free_char (S->file, -1);
  S->file=(char**)vcalloc (S->nseq, sizeof (char*));

  free_char (S->blastfile, -1);
  S->blastfile=(char**)vcalloc (S->nseq, sizeof (char*));
  
  if (getenv("protein_db_4_TCOFFEE"))db=csprintf (db, "%s", getenv("protein_db_4_TCOFFEE"));
  if (!getenv("protein_db_4_TCOFFEE"))myexit(fprintf_error (stderr,"Database has not been set for Blast"));
  if (!isfile(db))myexit(fprintf_error (stderr,"Database %s not available",db));
  dbn=path2filename(db);
  
  if (getenv("num_iterations_4_TCOFFEE"))num_iterations=atoi(getenv("num_iterations_4_TCOFFEE"));
  else num_iterations=1;
  
  if (getenv("outfmt_4_TCOFFEE"))outfmt=atoi(getenv("outfmt_4_TCOFFEE"));
  else outfmt=5;
  
  if (getenv("cache_4_TCOFFEE"))
    {
      outdir=csprintf (outdir, "%s", getenv("cache_4_TCOFFEE"));
      if (outdir[0]!='/')outdir=csprintf (outdir, "./%s", getenv("cache_4_TCOFFEE"));
    }
  else  outdir=csprintf (outdir, "./");
  my_mkdir (outdir);
      
  if (getenv("compress_4_TCOFFEE"))compress=1;
  else compress=0;
  
  dolist=(int*)vcalloc ( S->nseq, sizeof (int));
  
  for ( b=0,a=0; a<S->nseq; a++)
    {
      char *in;
      sname=clean_sname(csprintf (sname, "%s", S->name[a]));
      S->file[a]=csprintf (S->file[a],"%s/%s.blastp.%s.LOCAL.%d.infile.tmp", outdir,sname,dbn,num_iterations);
      in=S->blastfile[a]=csprintf (S->blastfile[a],"%s/%s.blastp.%s.LOCAL.%d.tmp", outdir,sname,dbn,num_iterations);
      inz=csprintf (inz, "%s.gz",in);

      if ( isfile (in) || isfile (inz))
	{
	  dolist[a]=0;
	}
      else 
	{
	  dolist[a]=1;
	  b++;
	}
    }
  
  if (!b)
    {
      vfree (db);vfree(dbn);vfree(sname);vfree(outdir);
      vfree(buf); 
      vfree(header);vfree(lheader); vfree(ilen);vfree (olen);vfree(in); vfree(inz);
      vfree(dolist );vfree(footer );
      return S;
    }

  fp=vfopen ( FullSeqF, "w");
  for ( a=0; a<S->nseq; a++)
    {
      if (dolist[a])
	{
	  fprintf ( fp, ">%s\n%s\n", S->name[a], S->seq[a]);
	  b++;
	}
    }
  vfclose (fp);
  
  printf_system ("psiblast -db %s -query %s -num_iterations %d -out %s -outfmt %d", db, FullSeqF, num_iterations, FullBlastF, outfmt);
  
  fp=vfopen (FullBlastF, "r");
  while ((buf=vfgets (buf, fp))!=NULL && !strstr (buf, "<BlastOutput_iterations>"))
    header=vcat (header, buf);
  footer=csprintf (NULL, "</BlastOutput_iterations>\n</BlastOutput>\n");
  

  for ( a=0; a<S->nseq; a++)
    {
      if (!dolist[a])
	{
	  HERE ("...Skip %s", S->name[a]);
	  continue;
	}
      lheader=csprintf (lheader, "%s", header);
      if (a>0)
	{
	  lheader=substitute(lheader, S->name[0], S->name[a]);
	  ilen=csprintf (ilen,"<BlastOutput_query-len>%d</BlastOutput_query-len>", strlen (S->seq[0]));
	  olen=csprintf (olen,"<BlastOutput_query-len>%d</BlastOutput_query-len>", strlen (S->seq[a]));
	}
      
      string2file(S->file[a],"w", ">A\n%s\n", S->seq[a]);//Note: name is removed
      out=vfopen (S->blastfile[a], "w");
      fprintf ( out, "%s", lheader);
      while ((buf=vfgets (buf, fp))!=NULL && !strstr (buf, "</Iteration>"))
	{
	  fprintf (out, "%s", buf);
	}
      fprintf (out, "%s", buf);
      fprintf ( out,"%s", footer);
      vfclose (out);
      if ( compress)
	{
	  ;
	  printf_system_direct ("gzip %s", (S->file[a]));
	  printf_system_direct ("gzip %s", (S->blastfile[a]));
	}
    }
  vfclose (fp);
  vfree (db);vfree(dbn);vfree(sname);vfree(outdir);
  vfree(buf); 
  vfree(header);vfree(lheader); vfree(ilen);vfree (olen);vfree(in); vfree(inz);
  vfree(dolist );vfree(footer );
    
  return S;
}
Sequence * seq2prf ( Sequence *S)
{
  int prot_min_cov, prot_min_sim, prot_max_sim, psitrim;
  char *outdir=NULL;
  char *psitrim_mode=NULL;
  char *psitrim_tree=NULL;
  char *sname=NULL;
  char *prfF=NULL;
  char *seqF=vtmpnam (NULL);
  int a;
  FILE *fp;

  seq2blast (S);
  
  if (getenv("prot_min_cov_4_TCOFFEE"))prot_min_cov=atoi(getenv("prot_min_cov_4_TCOFFEE"));
  else prot_min_cov=40;
  
  if (getenv("prot_min_sim_4_TCOFFEE"))prot_min_sim=atoi(getenv("prot_min_sim_4_TCOFFEE"));
  else prot_min_sim=50;
  
  if (getenv("prot_max_sim_4_TCOFFEE"))prot_max_sim=atoi(getenv("prot_max_sim_4_TCOFFEE"));
  else prot_max_sim=90;
  
  if (getenv("psitrim_mode_4_TCOFFEE"))psitrim_mode=getenv("psitrim_mode_4_TCOFFEE");
  else psitrim_mode=csprintf (psitrim_mode, "regtrim");
  
  if (getenv("psitrim_tree_4_TCOFFEE"))psitrim_mode=getenv("psitrim_tree_4_TCOFFEE");
  else psitrim_tree=csprintf (psitrim_tree, "codns");
  
  if (getenv("psitrim_4_TCOFFEE"))prot_max_sim=atoi(getenv("psitrim_4_TCOFFEE"));
  else psitrim=100;
  
  
  if (getenv("cache_4_TCOFFEE"))
    {
      outdir=csprintf (outdir, "%s", getenv("cache_4_TCOFFEE"));
      if (outdir[0]!='/')outdir=csprintf (outdir, "./%s", getenv("cache_4_TCOFFEE"));
    }
  else  outdir=csprintf (outdir, "./");
  my_mkdir (outdir);
  
  for ( a=0; a<S->nseq; a++)
    {
      sname=clean_sname(csprintf (sname, "%s", S->name[a]));
      prfF=csprintf (prfF, "%s/%s.prf", outdir,sname);
      S->seq_comment[a]=strcatf(S->seq_comment[a], " _R_ %s ", prfF);
      S->aln_comment[a]=strcatf(S->aln_comment[a], " _R_ %s ", prfF);
      string2file (seqF,"w", ">%s\n%s\n", S->name[a], S->seq[a]);
      printf_system ("tc_generic_method.pl -mode=blast2prf -infile=%s -seqfile=%s -outfile=%s -minid=%d -maxid=%d -mincov=%d -trim=%d", S->blastfile[a], seqF, prfF, prot_min_sim, prot_max_sim, prot_min_cov, psitrim);
    }
  
  return S;
}




Sequence * seq2unique_name_seq ( Sequence *S)
{
  int a;
 
  if ((a=name_list2unique_name_list (S->nseq, S->name)))
    {
      add_warning ( stderr, "Sequence %s is duplicated in file %s. The sequence will be renamed", S->name[a-1], S->file[a-1]);
    }
  return S;
}
Alignment * aln2unique_name_aln ( Alignment *S)
{
  int a;
  if ((a=name_list2unique_name_list (S->nseq, S->name)))
    {
      add_warning ( stderr, "Sequence %s is duplicated in file %s. The sequence will be renamed", S->name[a-1], S->file[a-1]);
    }
  return S;
}


int name_list2unique_name_list (int n, char **name)
{
  int duplicate=0;
  int a, b;

  for (a=0; a<n-1; a++)
    for (b=a+1; b<n; b++)
      {
	if ( strm (name[a], name[b]))
	  {duplicate=a+1;b=a=n;}
      }

  if (duplicate)
    {
      char *tmp1, *tmp2;
      Sequence *S;
      FILE *fp;

      tmp1=vtmpnam (NULL);
      tmp2=vtmpnam (NULL);
      fp=vfopen (tmp1, "w");
      for (a=0; a< n; a++)fprintf ( fp, ">%s\naggggg\n", name[a]);
      vfclose (fp);
      printf_system ("fasta_aln2fasta_aln_unique_name.pl %s > %s", tmp1, tmp2);
      S=get_fasta_sequence (tmp2, NULL);
      for (a=0; a<n; a++)
	{
	  name[a]=(char*)vrealloc (name [a], sizeof (int)*(strlen (S->name[a])+1));
	  sprintf ( name[a], "%s", S->name [a]);
	}
      free_sequence(S, -1);
    }
  return duplicate;
}
char**gene2exons    (char **seq, int nseq)
{

  int a, b, c,r;
  for (a=0; a<nseq; a++)
    {
      int in_exon=0, flag=0,l;
      l=strlen (seq[a]);
      for ( b=0; b<l; b++)
	{
	  r=seq[a][b];
	  if (isupper (r))
	    {
	      in_exon=1;
	      seq[a][b]=(flag)?r:tolower(r);
	    }
	  else if (in_exon)
	    {
	      in_exon=0;
	      flag=1-flag;
	      seq[a][b]='-';
	    }
	  else seq[a][b]='-';
	}
    }
  return seq;
}
Sequence* seq2clean_seq (Sequence *S, char *alp)
{
  int a, b, c, d, l;

  for (a=0; a< S->nseq; a++)
    {
      l=strlen (S->seq[a]);
      for (d=0,b=0; b<l; b++)
	{
	  c=S->seq[a][b];
	  if ( alp==NULL && !strchr (AA_ALPHABET, c) && !strchr (DNA_ALPHABET, c));
	  else if (alp && strchr (alp, c));
	  else S->seq[a][d++]=c;
	}
      S->seq[a][d]='\0';
      S->len[a]=strlen (S->seq[a]);
    }
  return S;
}
int ** seq2aln_pos      (Alignment *A, int *ns, int **l_s)
    {
    int **code;
    int a, b,c, d,l, p , g;


    l=MAX(strlen (A->seq_al[l_s[0][0]]), strlen (A->seq_al[l_s[1][0]]));
    code=declare_int ((A->S)->nseq,l+1);

    for (c=0; c<2; c++)
        {
	l=strlen (A->seq_al[l_s[c][0]]);
	for (d=0; d<ns[c]; d++)
	    {
	    a=A->order[l_s[c][d]][0];
	    for (p=0, b=0; b<l; b++)
	        {
		    g=is_gap (A->seq_al[l_s[c][d]][b]);
		    if (!g){p++; code[a][p]=b+1;}
		}
	    }
	}
    return code;
    }

Alignment *local_maln2global_maln (char *seq, Alignment *A)
    {
      /*inputs a BLAST alignmnent where the master sequence may be partila
	outputs the same alignment, while amkeing sure the profile is perfectly in sink with its master sequence
      */

      int a, b, c;
      int start, end, rend;
      char qname[100], *p;
      Alignment *B=NULL;

      sprintf ( qname, "%s", A->name[0]);
      p=strtok (qname, "_");
      if ( !strm (p, "QUERY"))
	   {
	     fprintf ( stderr, "\nUnappropriate format for the alignment [%s:FATAL]", PROGRAM);
	     myexit (EXIT_FAILURE);
	   }

      start=atoi(strtok (NULL, "_"));
      end=atoi(strtok (NULL, "_"));
      rend=strlen (seq);

      B=copy_aln (A,NULL);
      if ( start>1 || end<rend )A=realloc_aln (A,rend+1);

      for (a=0; a<start-1; a++)
	{
	  A->seq_al[0][a]=seq[a];
	  for ( b=1; b< A->nseq; b++)A->seq_al[b][a]='-';
	}

      for (c=0,a=start-1; a< end; a++, c++)
	{
	  A->seq_al[0][a]=seq[a];
	  for ( b=1; b< A->nseq; b++)
	    {
	      A->seq_al[b][a]=B->seq_al[b][c];
	    }
	}
      for ( a=end; a<rend; a++)
	{
	  A->seq_al[0][a]=seq[a];
	  for ( b=1; b< A->nseq; b++)A->seq_al[b][a]='-';
	}
      for ( a=0; a< A->nseq; a++) A->seq_al[a][rend]='\0';
      free_aln (B);

      A->len_aln=rend;
      return A;
    }


int ** aln2inv_pos ( Alignment *A)
{
	int **pos,a;
	pos=(int**)vcalloc (A->nseq, sizeof (char*));
	for (a=0; a< A->nseq; a++)pos[a]=seq2inv_pos (A->seq_al[a]);
	return pos;
}


int *  seq2inv_pos ( char *seq)
{
  /*returns a list where each value gives the index of the corresponding residue in seq*/
  /*Numbering: 1 to L : Analogy to the aln2pos*/

  int a,l1, l2;
  int *pos;

  l1=strlen ( seq);
  for ( l2=a=0; a< l1; a++)l2+=1-is_gap(seq[a]);
  pos=(int*)vcalloc (l2+1, sizeof (int));
  for ( l2=a=0; a< l1; a++)if (!is_gap(seq[a]))pos[++l2]=a+1;
  return pos;
}


int ** aln2pos_simple_2 (Alignment *A)
    {
    int **pos1;
    int **pos2;
    pos1=aln2pos_simple (A, A->nseq);
    pos2=duplicate_int  (pos1, A->nseq,read_size_int (pos1[0],sizeof (int)));
    pos1=aln2pos_simple (NULL, 0);
    return pos2;
    }
int ** aln2pos_simple (Alignment *A, int n_nseq, ...)
    {
    /*
    function documentation: start
    int ** aln2pos_simple (Alignment *A, int n_nseq, ...)

####with two parameter only: Alignment *A, int n_nseq

    this function turns A into pos, a matrix where each residue is replaced by its index according to the complete sequence.
    the indices in pos are computed using A->order[x][1] that contains the indice of the first residue of seq x of A

    n_nseq MUST not be null

####with more than two param:
     int ** aln2pos_simple (Alignment *A, int n_nseq, int *ns, int **ls)
     n_nseq must be set to 0 for the param 3 and four to be read

     ns[x]=number seq in group
     ls[x]=list of the sequences in group x ( size=ns[x])

    The computation of the indices is only carried out on the scpecified residues

####IMPORTANT
      in pos, the numbering of the residues goes from 1 to L:
        pos[seq=0][col=0]=1, means that the position #1  sequence #1
	in the alignmnet contains residue #1 from sequence A->order[0][0];
	gaps are negative values:
	in other words
	
	r=pos[seq][col]-1;
	if (r>0) the index is correct
    function documentation: end
    */

    int a, b,c, p, g,l;
    int **T;

    int max_nseq;
    int n_len=0;

    int *list=NULL;
    int *ns=NULL;
    int **ls=NULL;



    va_list ap;


    if ( A==NULL)
       {
	 return NULL;
       }
    else
       {
       if ( n_nseq>0)
          {
	  list=(int*)vcalloc(n_nseq, sizeof (int));
	  for ( a=0; a< n_nseq; a++)list[a]=a;
	  }
       else
          {
	  va_start (ap, n_nseq);
	  ns=va_arg(ap, int * );
	  ls=va_arg(ap, int **);
	  va_end(ap);
	  list=(int*)vcalloc ( ns[0]+ns[1], sizeof (int));
	  n_nseq=0;
	  for ( a=0; a< ns[0]; a++)list[n_nseq++]=ls[0][a];
	  for ( a=0; a< ns[1]; a++)list[n_nseq++]=ls[1][a];

	  }
       max_nseq=MAX(read_size_int(A->order,sizeof (int*)),return_max_int (A->order, read_size_int(A->order,sizeof (int*)),0))+1;
       n_len=get_longest_string ( A->seq_al,A->max_n_seq, NULL, NULL)+1;


       T=declare_int (max_nseq, n_len);
       for ( c=0; c< n_nseq; c++)
           {
	   a=list[c];
	   l=strlen ( A->seq_al[a]);

	   for ( p=A->order[a][1],b=0; b<l; b++)
	       {
	       g=1-is_gap(A->seq_al[a][b]);
	       p+=g;
	       T[a][b]=(g==1)?p:-(1+p);
	       if ( A->seq_al[a][b]==UNDEFINED_RESIDUE)T[a][b]=0;
	       if ( A->seq_cache && T[a][b]>0)T[a][b]=A->seq_cache[A->order[a][0]][T[a][b]];
	       }
	   }
       vfree (list);
       }

   return T;
   }
Alignment ** split_seq_in_aln_list ( Alignment **aln, Sequence *S, int n_seq, char **seq_list)
        {
	int a, b, c;
	char * long_seq=NULL;
	int    len,l;
	int  **translation;
	int  **table;




	if ( aln==NULL)return NULL;
	translation=declare_int ( S->nseq,2);

	for (len=0,a=0; a< S->nseq; a++)
	    {
	    if((b=name_is_in_list (S->name[a],seq_list, n_seq, 100))!=-1)
	       {
	       l=strlen(S->seq[a])+1;
	       long_seq=(char*)vrealloc(long_seq,(len+l+1)*sizeof(char));
	       long_seq=strcat(long_seq, S->seq[a]);
	       long_seq=strcat(long_seq, "*");

	       translation[a][0]=b;
	       translation[a][1]=len;
	       len+=l;
	       }
	    else translation[a][0]=-1;
	    }

	long_seq[len-1]='\0';
	len--;

	table=declare_int ( len+1, 2);

	for ( b=0,a=0; a< S->nseq; a++)
	    {
	    if ( translation[a][0]!=-1)
	       {
	       c=1;
	       while (long_seq[b]!='\0' && long_seq[b]!='*')
		   {
		   table[b+1][1]=c++;
		   table[b+1][0]=translation[a][0];
		   b++;
		   }
	       table[b][1]=c++;
	       table[b][0]=translation[a][0];
	       b++;
	       }
	    }

	for ( a=0; a< (aln[-1])->nseq; a++)
	    {
	    for ( b=0; b< (aln[a])->nseq; b++)
	        {

		(aln[a])->order[b][0]=table[(aln[a])->order[b][1]][0];
		(aln[a])->order[b][1]=table[(aln[a])->order[b][1]][1];
		sprintf ( (aln[a])->name[b],"%s_%d_%d", S->name[(aln[a])->order[b][0]],a+1,b+1);
		}
	    }
	free_int (translation, -1);
	free_int (table,       -1);
	return aln;
	}


Sequence  *  fill_sequence_struc ( int nseq, char **sequences, char **seq_name, Genomic_info *genome_co)
{
	int a;
	Sequence *S;
	int shortest, longuest;

	if (!sequences)
	{
		shortest=longuest=0;
	}
	else if ( nseq>1)
	{
		shortest=get_shortest_string( sequences, nseq, NULL, NULL);
		longuest=get_longest_string (sequences, nseq, NULL, NULL);
	}
	else if ( nseq==1)
	{
		shortest=longuest=strlen (sequences[0]);
	}
	else
	{
		return NULL;
	}


	S=declare_sequence (shortest, longuest,nseq);
	S->nseq=nseq;
	
	if (sequences)
		S->seq=copy_char ( sequences, S->seq);
	else
		S->seq=declare_char (S->nseq, 1);

	S->name=copy_char ( seq_name, S->name);

	ungap_array (S->seq,nseq);
	for ( a=0; a< S->nseq; a++)
		S->len[a]=strlen(S->seq[a]);

	if (genome_co != NULL)
	{
		unsigned int len;
		S->genome_co =(Genomic_info*)vcalloc(nseq, sizeof(Genomic_info));
		Genomic_info *tmp, *tmp_ori;
		for ( a=0; a< S->nseq; a++)
		{
			tmp = &(S->genome_co[a]);
			tmp_ori = &(genome_co[a]);
			tmp->strand = tmp_ori->strand;
			tmp->start = tmp_ori->start;
			tmp->end = tmp_ori->end;
			tmp->seg_len = tmp_ori->seg_len;
			tmp->seg_name =(char*) vcalloc(strlen(tmp_ori->seg_name)+1, sizeof(char));
			strcpy(tmp->seg_name, tmp_ori->seg_name);

// 			printf("-%s %i-\n", tmp_ori->seg_name, tmp_ori->start);
		}
	}
	else
		S->genome_co = NULL;


	return S;
}

Alignment * thread_profile_files2aln (Alignment *A, char *template_file, Fname *F)
{

  Alignment *P;
  int a;

  if (!A->S)A->S=aln2seq (A);
  if (template_file)A->S=seq2template_seq (A->S, template_file,F);
  for ( a=0; a< A->nseq; a++)
    {
      P=seq2R_template_profile (A->S, a);
      if ( P)
	{
	  P->expand=1;
	  sprintf ( P->name[0], "%s", A->name[a]);
	}
    }

  return expand_aln (A);
}




Alignment * expand_aln (Alignment *A)
  {
  /*This function expands the profiles within an alignment*/


  int a, b, d, e;
  Alignment *MAIN=NULL, *SUB=NULL;
  int n_sub_seq=0;
  int new_nseq=0;
  int *list;
  Alignment *Profile;

  if ( !A)return A;



  list=(int*)vcalloc (A->nseq, sizeof (int));
  for ( a=0; a< A->nseq; a++)
    {
      Profile=seq2R_template_profile (A->S, A->order[a][0]);
      if (Profile && Profile->expand)
	{
	  new_nseq+=Profile->nseq;
	}
      else
	{
	  new_nseq++;
	  list[n_sub_seq++]=a;
	}
    }

  if ( n_sub_seq==A->nseq){vfree(list);return A;}
  else if (n_sub_seq==0){MAIN=copy_aln (A, MAIN);MAIN->nseq=0;}
  else
    {
      MAIN=extract_sub_aln (A, n_sub_seq, list);
    }
  vfree(list);


  for ( a=0; a< A->nseq; a++)
    {
      Profile=seq2R_template_profile (A->S, A->order[a][0]);
      if ( Profile && Profile->expand)
	{

	  SUB=copy_aln (Profile,SUB);

	  SUB=realloc_aln2(SUB, SUB->nseq, A->len_aln+1);

	  for ( e=0,b=0; b< A->len_aln; b++)
	    {
	      if ( is_gap(A->seq_al[a][b]))
		{for (d=0; d< SUB->nseq; d++)SUB->seq_al[d][b]='-';}
	      else
	         {
		   for(d=0; d<SUB->nseq; d++)SUB->seq_al[d][b]=Profile->seq_al[d][e];
		   e++;
		 }

	    }
	  MAIN=stack_aln(MAIN, SUB);
	}
    }
  free_aln (A);
  free_aln (SUB);
  return MAIN;
  }
Alignment * expand_number_aln (Alignment *A,Alignment *EA)
  {
  /*This function expands the profiles within an alignment*/


  int a, b, d, e;
  Alignment *MAIN=NULL, *SUB=NULL, *C=NULL;
  int n_sub_seq=0;
  int new_nseq=0;
  int *list;
  Alignment *Profile;

  if ( !EA || !A)return EA;

  if ( EA->nseq<A->nseq)
    {
      fprintf (stderr, "\n[ERROR:expand_number_aln] Using as a master an expanded aln (%d %d) [FATAL:%s]", EA->nseq, A->nseq,PROGRAM);
      EA->A=A->A=NULL;
      print_aln (EA);
      print_aln (A);
      myexit (EXIT_FAILURE);
    }


  list=(int*)vcalloc (EA->nseq, sizeof (int));
  for ( a=0; a< EA->nseq; a++)
    {
      Profile=seq2R_template_profile (EA->S, EA->order[a][0]);
      if (Profile && Profile->expand)new_nseq+=Profile->nseq;
      else
	{
	  new_nseq++;
	  list[n_sub_seq++]=a;
	}
    }

  if ( n_sub_seq==EA->nseq){vfree(list);return EA;}
  else if (n_sub_seq==0){MAIN=copy_aln (EA, MAIN);MAIN->nseq=0;}
  else
    {
      MAIN=extract_sub_aln (EA, n_sub_seq, list);
    }


  list[0]=EA->nseq;
  C=extract_sub_aln (EA,1, list);
  vfree(list);



  for ( a=0; a< EA->nseq; a++)
    {
      Profile=seq2R_template_profile (EA->S, EA->order[a][0]);
      if ( Profile && Profile->expand)
	{
	  SUB=copy_aln (Profile,SUB);
	  SUB=realloc_aln2(SUB, SUB->nseq, EA->len_aln+1);

	  for ( e=0,b=0; b<= EA->len_aln; b++)
	    {
	      if (is_gap(A->seq_al[a][b]))
		{
		for ( d=0; d<SUB->nseq; d++)
		  SUB->seq_al[d][b]=NO_COLOR_RESIDUE;
		}
	      else
		{
		  for ( d=0; d<SUB->nseq; d++)
		    {

		      if ( is_gap (Profile->seq_al[d][e]))
			{
			  SUB->seq_al[d][b]=NO_COLOR_RESIDUE;
			}
		      else SUB->seq_al[d][b]=EA->seq_al[a][b];
		    }
		  e++;
		}
	    }
	  for (d=0; d< SUB->nseq; d++)SUB->score_seq[d]=EA->score_seq[a];

	  MAIN=stack_aln(MAIN, SUB);
	}
    }

  MAIN=stack_aln(MAIN, C);
  MAIN->nseq--;
  MAIN->score=MAIN->score_aln=EA->score_aln;

  free_aln (SUB);
  free_aln (EA);

  free_aln (C);

  return MAIN;
  }

Alignment * probabilistic_rm_aa ( Alignment *A, int pos, int len)
{
  int random_len=0;
  int a, b;
  int left, right;

  if ( len<0)
    {
      random_len=1;
      len=-len;
    }

  vsrand(0);

  if (pos==0)pos= (rand()%(A->len_aln-(2*len+len))) +len;


  for ( a=0; a< A->nseq; a++)
	{
	  if (random_len)left =rand()%len;
	  else left=len;
	  if (random_len)right=rand()%len;
	  else right=len;
	  if ( (pos-right)<0 || (pos+left)>A->len_aln)
	    {
	      add_warning ( stderr, "probabilistic_rm_aa, pos out of range [%s]\n", PROGRAM);
	    }
	  else
	    for ( b=pos-right; b<pos+left; b++)A->seq_al[a][b]=(b==pos)?'~':'*';
	}

  ungap_aln (A);
  free_sequence ( A->S, A->nseq);
  A->S=aln2seq (A);
  return A;

}

Alignment * remove_gap_column ( Alignment *A, char *mode)
  {
    int   a, b;
    char *p;
    int  *seq_list;
    int   nseq=0;
    int keep_col, cl;


    seq_list =(int*)vcalloc ( A->nseq, sizeof (int));
    while (  (p=strtok(mode, ":")))
      {
	mode=NULL;
	if (p[0]=='#')
	  {
	    seq_list[nseq++]=atoi(p+1)-1;
	  }
	else if ( (a=name_is_in_list (p, A->name, A->nseq, 100))!=-1)
	  {
	    seq_list[nseq++]=a;
	  }
      }

    if ( nseq==0)
      {
	for ( a=0; a< A->nseq; a++)seq_list[a]=a;
	nseq=A->nseq;
      }

    for ( cl=0,a=0; a<=A->len_aln; a++)
      {
	for (keep_col=1, b=0; b< nseq && keep_col; b++)
	  {
	    keep_col=(is_gap(A->seq_al[seq_list[b]][a]))?0:keep_col;
	  }

	if ( keep_col)
	  {
	    for ( b=0; b< A->nseq; b++)
	      {
		A->seq_al[b][cl]=A->seq_al[b][a];
	      }
	    cl++;
	  }
	else
	  {
	    for ( b=0; b< A->nseq; b++)
	      {
		A->seq_al[b][cl]='-';
	      }
	    cl++;
	  }
      }
    A->len_aln=cl;
    vfree (seq_list);

    return A;
  }


Alignment * ungap_sub_aln (Alignment *A, int ns, int *ls)
        {

	int a, b, c,t;
	int len;

	len=strlen ( A->seq_al[ls[0]]);

	for ( c=0,a=0; a<len; a++)
		{
		for ( t=0,b=0; b<ns; b++)
			t+=is_gap(A->seq_al[ls[b]][a]);
		if (t==ns);
		else
		    {
		    for ( b=0; b<ns; b++)
		    	A->seq_al[ls[b]][c]=A->seq_al[ls[b]][a];
		    c++;
		    }
	 	}
	 for ( b=0; b<ns; b++)A->seq_al[ls[b]][c]='\0';
	 return A;
	}

Sequence * ungap_seq ( Sequence *S)
        {
	  int a;

	  if ( !S)return NULL;
	  ungap(S->seq[0]);
	  S->max_len=S->min_len=strlen (S->seq[0]);
	  for ( a=0; a< S->nseq; a++)
	    {
	      ungap(S->seq[a]);
	      S->len[a]=strlen (S->seq[a]);
	      S->max_len=MAX(S->max_len,S->len[a]);
	      S->min_len=MAX(S->min_len,S->len[a]);
	    }
	  return S;

	}
Alignment* shift_column (Alignment *A, int from, int to);
int max_shift (Alignment *A, int p);
int column_is_lower (Alignment *A, int p);

Alignment * unalign_aln_2 (Alignment *A, Alignment *C, int t)
{
  int a, b, pos, len;
  Sequence *S;
  int n, insert;
  if (C)
    {
      for (a=0; a<A->nseq; a++)
	for (b=0; b<A->len_aln; b++)
	  {
	    int res=C->seq_al[a][b];
	    A->seq_al[a][b]=toupper(A->seq_al[a][b]);
	    if ((isdigit (res) && (res-'0')<=t))
	      A->seq_al[a][b]=tolower(A->seq_al[a][b]);
	  }
    }

  n=0;
  while ( A->seq_al[0][n])
    {
      insert=0;
      for (b=0; b<A->nseq; b++)if (islower (A->seq_al[b][n]))insert=1;
      if (insert)
	{
	  insert_gap_col (A,n,1);
	  for (b=0; b<A->nseq; b++)
	    {
	      if ( islower (A->seq_al[b][n+1]))
		{
		  A->seq_al[b][n]=A->seq_al[b][n+1];
		  A->seq_al[b][n+1]='-';
		}
	    }
	}
      n++;
    }
  for (a=A->len_aln-1; a>=0; a--)
    {
      if (column_is_lower (A,a))
	{
	  int s;
	  s=max_shift (A,a);
	  shift_column (A,a, a+s);
	}
    }
  return A;
}
Alignment* shift_column (Alignment *A, int from, int to)
{
  char *buf;
  int a;

  buf=(char*)vcalloc (A->nseq, sizeof (char));
  for (a=0; a<A->nseq; a++)
    {
      buf[a]=A->seq_al[a][from];
      A->seq_al[a][from]='-';
    }
  to++;
  insert_gap_col (A, to, 1);
  for ( a=0; a<A->nseq; a++)A->seq_al[a][to]=buf[a];
  vfree (buf);
  ungap_aln (A);
  return A;
}
int max_shift (Alignment *A, int p)
{
  int shift, max_shift, a;
  for (max_shift=A->len_aln,a=0; a< A->nseq; a++)
    {
      shift=0;

      if (!islower (A->seq_al[a][p]) || A->seq_al[a][p]=='-')continue;
      while (A->seq_al[a][p+shift+1]=='-')shift++;
      max_shift=MIN(shift,max_shift);
    }
  return max_shift;
}
int column_is_lower (Alignment *A, int p)
{
  int a;

  for ( a=0; a<A->nseq; a++)
    if ( !is_gap (A->seq_al[a][p]) && !islower(A->seq_al[a][p]))return 0;
  return 1;
}

Alignment * unalign_aln (Alignment *A, Alignment *C, int t)
{
  int a, b, pos, len;
  Sequence *S;

  for (a=0; a<A->nseq; a++)
    for (b=0; b<A->len_aln; b++)
      {
	int res=C->seq_al[a][b];
	A->seq_al[a][b]=toupper(A->seq_al[a][b]);
	if ((isdigit (res) && (res-'0')<=t))
	  A->seq_al[a][b]=tolower(A->seq_al[a][b]);
      }


  for (pos=-1, a=0; a<C->nseq; a++)
    {
      b=0;
      while ( C->seq_al[a][b])
	{
	  int res=C->seq_al[a][b];
	  if ((isdigit (res) && (res-'0')<=t))
	    {
	      if (pos==-1){pos=b;len=1;}
	      else len++;
	    }
	  else if (pos!=-1)
	    {

	      C=unalign_aln_pos(C,a,pos, len);
	      pos=-1;
	    }
	  b++;
	}
      if ( pos!=-1){C=unalign_aln_pos(C,a,pos, len);pos=-1;}
    }
  S=aln2seq (A);
  thread_seq_struc2aln (C, S);
  A=realloc_aln2 (A, A->nseq, C->len_aln+1);
  A->len_aln=C->len_aln;
  for (a=0; a<A->nseq; a++)sprintf ( A->seq_al[a], "%s", C->seq_al[a]);
  ungap_aln (A);

  free_sequence (S, -1);
  return A;
}
Alignment * unalign_aln_pos (Alignment *A, int s, int p, int l)
{
  int a;
  char *buf;
  int unalign=0;


  buf=(char*)vcalloc (l+1, sizeof (char));
  for (a=0; a<l; a++)
    {
      buf[a]=A->seq_al[s][p+a];
      A->seq_al[s][p+a]='-';
    }


  A=insert_gap_col (A,p, l);
  for (a=0; a<l; a++)
    {
      A->seq_al[s][p+a]=buf[a];
    }
  vfree (buf);
  return A;
}
Alignment * insert_gap_col (Alignment *A, int p, int l)
{
  int a, c;
  char *buf;
  char *gap;

  gap=generate_null(l);
  if ( !A || p>=A->len_aln || p<0)return A;

  buf=(char*)vcalloc (A->len_aln+l+1, sizeof (char));
  A=realloc_aln2(A,A->nseq, A->len_aln+l+1);
  for (a=0; a<A->nseq; a++)
    {
      c=A->seq_al[a][p];
      A->seq_al[a][p]='\0';
      sprintf ( buf, "%s%s%c%s", A->seq_al[a],gap,c,A->seq_al[a]+p+1);
      sprintf (A->seq_al[a], "%s", buf);
    }
  vfree (buf);
  A->len_aln+=l;
  return A;
}
Alignment * unalign_residues (Alignment *A, int si1, int si2)
{
  char *s1, *s2, *ns1, *ns2;
  int l, a, b,r1, r2;

  s1=A->seq_al[si1];s2=A->seq_al[si2];
  l=strlen (s1);

  ns1=(char*)vcalloc (2*l+1, sizeof (char));
  ns2=(char*)vcalloc (2*l+1, sizeof (char));

  for (b=a=0; a< l; a++)
    {
      r1=s1[a]; r2=s2[a];
      if (is_gap(r1) || is_gap(r2) || isupper (r1) || isupper(r2))
	{
	  ns1[b]=(r1=='.')?'-':r1;
	  ns2[b]=(r2=='.')?'-':r2;
	  b++;
	}
      else
	{
	  ns1[b]=r1;
	  ns2[b]='-';
	  b++;
	  ns2[b]=r2;
	  ns1[b]='-';
	  b++;
	}
    }
  ns1[b]='\0';
  ns2[b]='\0';
  A->seq_al[si1]=ns1;
  A->seq_al[si2]=ns2;


  A->len_aln=strlen (ns1);
  return A;
}
Alignment *degap_aln (Alignment *A)
{
  //Reomove all the gaps
  int a;
  for ( a=0; a< A->nseq; a++)ungap (A->seq_al[a]);
  return A;
}
Alignment *RmLowerInAln(Alignment *A, char *gap)
{
  char g;
  int a, b;

  g=(gap)?gap[0]:'-';
  for (a=0; a<A->nseq; a++)
    for (b=0; b<A->len_aln; b++)
      {
	char c=A->seq_al[a][b];
	if (is_gap(c));
	else if (islower(c))A->seq_al[a][b]=g;
      }
  return A;
}

char *ungap_fastaF_big (char *in, char *out, int p)
{
  FILE *fp;
  char *tmp=vtmpnam (NULL);
  int *cache=NULL;
  long *map;
  int l, ml=0;
  int a,b,nseq;
  char *rec, *name, *header, *seq;
  if (!out) out=vtmpnam(NULL);

  
  map=fasta2map(in);
  nseq=read_array_size_new (map)-1;

  if (p==-1)
    {
       fp=vfopen (tmp, "w");  
       for (a=0; a<nseq; a++)
	 {
	   rec=file2record_it (in,a,map);
	   seq=seq2clean(FastaRecord2seq(rec));
	   header=FastaRecord2header(rec);
	   fprintf (fp, "%s\n",header);
	   l=strlen (seq);
	   for (b=0; b<l; b++)
	     {
	       if (seq[b]!='-')fprintf ( fp, "%c", seq[b]);
	     }
	 }
      fprintf ( fp, "\n");
      vfclose (fp);
    }
  else
    {
      for (a=0; a<nseq; a++)
	{
	  rec=file2record_it (in,a, map);
	  seq=seq2clean(FastaRecord2seq(rec));
	  l=strlen (seq);
	  if (!cache){cache=(int*)vcalloc (l, sizeof (int));ml=l;}
	  else if (ml<l){cache=(int*)vrealloc (cache, ml*sizeof (int));ml=l;}
	  
	  for ( b=0; b<l; b++)if(seq[b]=='-')cache[b]++;
	}
      for ( a=0; a<ml; a++)cache[a]=(cache[a]*100)/nseq;
      
      
      fp=vfopen (tmp, "w");  
      for (a=0; a<nseq; a++)
	{
	  rec=file2record_it (in,a, map);
	  seq=seq2clean(FastaRecord2seq(rec));
	  header=FastaRecord2header(rec);
	  fprintf (fp, "%s\n",header);
	  l=strlen (seq);
	  for (b=0; b<l; b++)
	    {
	      if (cache[b]<p)fprintf ( fp, "%c", seq[b]);
	    }
	  fprintf ( fp, "\n");
	}
      vfclose (fp);
      vfree (cache);
    }
  
  printf_system ("mv %s %s", tmp, out);
  vfree (map);
  return out;
}



Alignment *ungap_aln_n ( Alignment *A, int p)
	{
/*remove all the columns of gap-only within an alignment*/
	int a, b, c;
	int t;
	int gp;

	if ( A->nseq==0)return A;

	for ( c=0,a=0; a< A->len_aln; a++)
		{
		for ( t=0,b=0; b<A->nseq; b++)
			t+=is_gap(A->seq_al[b][a]);
		gp=(t*100)/A->nseq;
		if (p>0 && (gp>=p || (t==A->nseq && p==100) || (t && p==1)));//Remove columns containing more than p% gaps
		else if (p<0 && (gp<=p || (t==0 && p==-100) ||(t && p==-1)));//remove columns containing less than p% gaps
		else
		  {
		    for ( b=0; b<A->nseq; b++)
		      A->seq_al[b][c]=A->seq_al[b][a];
		    c++;
		  }
	 	}
	for ( b=0; b<A->nseq; b++)A->seq_al[b][c]='\0';
	A->len_aln=c;
	return A;
	}

Alignment *ungap_aln ( Alignment *A)
{
  return ungap_aln_n (A, 100);
}
/*
Alignment *ungap_aln ( Alignment *A)
	{
	int a, b, c,t;

	for ( c=0,a=0; a< A->len_aln; a++)
		{
		for ( t=0,b=0; b<A->nseq; b++)
			t+=is_gap(A->seq_al[b][a]);
		if (t==A->nseq);
		else
		    {
		    for ( b=0; b<A->nseq; b++)
		    	A->seq_al[b][c]=A->seq_al[b][a];
		    c++;
		    }
	 	}
	 for ( b=0; b<A->nseq; b++)A->seq_al[b][c]='\0';
	 A->len_aln=c;
	 return A;

	 }
*/


Alignment *remove_end (Alignment *A)
        {
	int a, b, d;
	int left, right;

	for (a=0; a< A->len_aln; a++)
	    {
	    for ( b=0, d=0; b< A->nseq; b++)
		if ( !is_gap(A->seq_al[b][a]))d++;
	    if ( d>1)break;
	    }
	left=a;
	for (a=A->len_aln-1; a>0; a--)
	    {
	    for ( b=0, d=0; b< A->nseq; b++)
		if ( !is_gap(A->seq_al[b][a]))d++;
	    if ( d>1)break;
	    }
	right=a;

	return extract_aln(A, left, right+1);
	}

Alignment* condense_aln (Alignment *A)
{
  /* condense complementarz columns:
     X-       X
     -X  ....>X
     X-       X

  */
  int a, b, plen, n,m, r1, r2;

  plen=0;
  while ( A->len_aln !=plen)
    {
      plen=A->len_aln;
      for ( a=0; a< A->len_aln-1; a++)
	{
	  for ( n=m=b=0; b< A->nseq; b++)
	    {
	      r1=is_gap(A->seq_al[b][a]);
	      r2=is_gap(A->seq_al[b][a+1]);
	      n+=(r1 || r2);
	      m+=r1;
	    }

	  if ( n==A->nseq && m!=A->nseq)
	    {
	      for (b=0; b< A->nseq; b++)
		{
		  if (!is_gap(A->seq_al[b][a+1]))
		      {
			A->seq_al[b][a]=A->seq_al[b][a+1];
			A->seq_al[b][a+1]='-';
			}
		}
	      a++;
	    }
	}
    }
  A=ungap_aln(A);
  return A;
}




void compress_aln ( Alignment *A)
        {

	  /*remove all the columns of gap-only within an alignment*/
	int a, b, c, d;



	for (c=0, a=0; a< A->len_aln; a++)
	  {
	    for ( b=0, d=0; b< A->nseq; b++)
		if ( A->seq_al[b][a]!='-'){d=1; break;}
	    if ( d==0);
	    else
	      {
		for (b=0; b< A->nseq; b++)
		  A->seq_al[b][c]=A->seq_al[b][a];
	      c++;
	      }
	  }
	 A->len_aln=c;

	for ( a=0; a< A->nseq; a++)
	  A->seq_al[a][c]='\0';
	}

Alignment *seq_coor2aln ( Sequence *S, Alignment *A, int **coor, int nseq)
        {
	int a;
	char *buf;

	A=realloc_alignment2(A, nseq, return_maxlen ( S->seq, S->nseq)+1);
	for ( a=0; a< S->nseq; a++)sprintf ( A->file[a], "%s", S->file[a]);
	for ( a=0; a< nseq; a++)
	    {
	    sprintf (A->name[a], "Repeat_%d_%d", a, coor[a][0]);
	    buf=extract_char ( S->seq[coor[a][0]], coor[a][1]-1, coor[a][2]);
	    sprintf ( A->seq_al[a],"%s", buf);
	    vfree(buf);
	    A->order[a][0]=0;
	    A->order[a][1]=coor[a][1]-1;
	    }
	A->nseq=nseq;
	return A;
	}

Alignment *strings2aln (int nseq,...)
        {
	  /*strings2aln(nseq, <name1>, <seq1>, <name2>, <seq2>....)*/
	  va_list ap;
	  char **list, **list2;
	  char **name, **name2;
	  Sequence *S;
	  Alignment *A;
	  int a, max;

	  va_start(ap, nseq);
	  list=(char**)vcalloc (nseq, sizeof (char*));
	  name=(char**)vcalloc (nseq, sizeof (char*));
	  for ( a=0; a< nseq; a++)
	    {
	      name[a]=va_arg(ap,char*);
	      list[a]=va_arg(ap,char*);

	    }
	  va_end(ap);

	  for ( max=0,a=0; a< nseq; a++)
	    {
	      max=(strlen (list[a])>max)?strlen(list[a]):max;
	    }
	  list2=declare_char (nseq, max+1);
	  name2=declare_char (nseq, MAXNAMES+1);

	  for ( a=0; a< nseq; a++)
	    {
	      sprintf ( list2[a], "%s", list[a]);
	      sprintf ( name2[a], "%s", name[a]);
	    }


	  S=fill_sequence_struc(nseq,list2,name2, NULL);

	  free_char (list2, -1);
	  free_char (name2, -1);
	  vfree (list);
	  vfree(name);
	  A=seq2aln(S,NULL, 1);
	  return A;
	}

/**
 * Transform a ::Sequence into an ::Alignment.
 */
Alignment *seq2aln (Sequence *S, Alignment *A,int rm_gap)
	{
	int a;

	
	A=realloc_alignment2(A, S->nseq, S->max_len+1);
	for ( a=0; a< S->nseq; a++)sprintf ( A->file[a], "%s", S->file[a]);
	A->nseq=S->nseq;
	A->max_len=S->max_len;
	A->min_len=S->min_len;

	
	for ( a=0; a< S->nseq; a++)
		{
		A->order[a][0]=a;
		A->order[a][1]=0;

		A->seq_comment[a]=csprintf ( A->seq_comment[a], "%s", S->seq_comment[a]);
		A->aln_comment[a]=csprintf ( A->aln_comment[a], "%s", S->aln_comment[a]);

		A->name[a]=csprintf ( A->name[a], "%s", S->name[a]);
		A->seq_al[a]=csprintf ( A->seq_al[a], "%s", S->seq[a]);


		ungap ( A->seq_al[a]);
		A->len[a]=strlen ( A->seq_al[a]);

		if ( rm_gap==0 || rm_gap==NO_PAD)A->seq_al[a]=csprintf ( A->seq_al[a], "%s", S->seq[a]);

		}
	if (rm_gap!=NO_PAD)padd_aln (A);
	A->S=S;
	return A;
	}




Alignment *padd_aln ( Alignment *A)
{
  A->seq_al=padd_string (A->seq_al, A->nseq, '-');
  A->len_aln=strlen (A->seq_al[0]);
  return A;
}

char **padd_string ( char **string, int n,char pad)
{
  /*Pads a the strings so that they all have the same length*/

  int max_len, a;
  char *buf;

  max_len=get_longest_string  (string,n, NULL, NULL);
  for (a=0; a<n; a++)
	    {
	    buf=generate_null (max_len-strlen (string[a]));
	    string[a]=vcat ( string[a], buf);
	    vfree (buf);
	    }
  return string;
}


char *msaF2fastaF(char *file)
{
  Alignment *A;
  Sequence *S;
  char *tmp=vtmpnam(NULL);
  int a;
  FILE *fp;
  
  A=main_read_aln (file, NULL);
  fp=vfopen (tmp, "w");
  for (a=0; a<A->nseq; a++)
    {
      fprintf ( fp, ">%s %s\n%s\n", A->name [a], A->aln_comment[a], A->seq_al[a]);
    }
  vfclose (fp);
  S=free_aln (A);
  if (S)free_sequence (S,S->nseq);  
  return tmp;
}
  
int trim_fastaF_big(char *in1, char*in2, char *nout1, char *nout2, long **nmap1, long **nmap2)
{
  //Same as trim_aln_file, but works on big files: no need to read the entire msa
  //returns a file map
  long *map1, *map2;
  int nseq1, nseq2;
  int tot=0, pos=0;
  FILE *fp;
  char **name;
  int i1, i2, a, l;
  char *s;
  int *order;
  char *out1=vtmpnam (NULL);
  char *out2=vtmpnam (NULL);
  
  if (!format_is_fasta(in1))in1=msaF2fastaF(in1);
  if (!format_is_fasta(in2))in2=msaF2fastaF(in2);
  

  map1=fasta2map(in1);
  nseq1=read_array_size_new(map1)-1;
    
  map2=fasta2map(in2);
  nseq2=read_array_size_new(map2)-1;
  
  if (nmap1)nmap1[0]=(long*)vcalloc (nseq1+1, sizeof (long));
  if (nmap2)nmap2[0]=(long*)vcalloc (nseq2+1, sizeof (long));
  order=(int*)vcalloc ( nseq1, sizeof (int));  
  name=(char **)vcalloc (nseq1, sizeof (char*));
  
  for (a=0; a<nseq1; a++)
    {
      name[a]=csprintf(name[a],"%s",FastaRecord2name(file2record_it(in1,a,map1))); 
    }
 
  fp=(out2)?vfopen (out2, "w"):NULL;
  for (tot=0,pos=0, i2=0; i2<nseq2; i2++)
    {
      s=file2record_it(in2,i2, map2);
      if ((i1=name_is_in_hlist (FastaRecord2name(s), name, nseq1))!=-1)
	{
	  if (fp)fprintf ( fp, "%s",s);
	  if (nmap2)nmap2[0][tot]=pos; 
	  pos+=strlen (s);
	  order[tot]=i1;
	  tot++;
	}
    }
  if (nmap2)nmap2[0][tot]=pos;
  if (fp)vfclose (fp);
  
  fp=(out1)?vfopen (out1, "w"):NULL;
  for (pos=0,i1=0; i1<tot; i1++)
    {
      s=file2record_it(in1,order[i1], map1);
      if (fp)fprintf (fp, "%s", s);
      if (nmap1)nmap1[0][i1]=pos; pos+=strlen (s);
    }
  if (nmap1)nmap1[0][tot]=pos;
  if (fp)vfclose (fp);

  //must be set right because it is used to determine nseq
  if (nmap1)nmap1[0]=(long*)vrealloc (nmap1[0],(tot+1)*sizeof (long));
  if (nmap2)nmap2[0]=(long*)vrealloc  (nmap2[0],(tot+1)*sizeof (long));
  
  free_char (name, -1);
  vfree (order);
  vfree (map1);
  vfree (map2);
  
  if (nout1)printf_system ("mv %s %s", out1, nout1);
  if (nout2)printf_system ("mv %s %s", out2, nout2);

  return tot;
}
  
Alignment * trim_aln_with_seq ( Alignment *S, Alignment *P)
{
  Alignment *A, *R;
  int a, b, c;
  static int seqindex;
  P=aln2profile (P);
  S=aln2profile (S);

  A=align_two_aln (S,P, "blosum62mt",-8,-1, "myers_miller_pair_wise");
  for (a=0; a<A->nseq; a++) sprintf (A->name[a], "tmpname_%d", seqindex++);

  R=copy_aln (A, NULL);
  for (c=0, a=0; a< A->len_aln; a++)
    {
      if ( is_gap (A->seq_al[0][a]));
      else
	{
	  for ( b=0; b<A->nseq; b++)
	    R->seq_al[b][c]=A->seq_al[b][a];
	  c++;
	}
    }
  for ( a=0; a< A->nseq; a++)R->seq_al[a][c]='\0';
  R->len_aln=c;
  R->S=aln2seq (R);

  free_aln (S);
  free_aln (P);
  free_aln (A);

  return R;
}

Alignment * add_align_seq2aln ( Alignment *A, char *seq, char *seq_name)
        {
	  if ( !A)
	    {
	      A=declare_aln (NULL);
	      A=realloc_aln2 ( A, 1, strlen (seq)+1);
	      A->nseq=0;
	      sprintf ( A->name[A->nseq], "%s", seq_name);
	      sprintf ( A->seq_al[A->nseq], "%s", seq);
	      A->nseq++;

	    }
	  else if ( strlen (seq)!=A->len_aln)
	    {
	      fprintf ( stderr, "\nError: Attempt to stack incompatible aln and aligned sequence[FATAL]\n");
	      myexit (EXIT_FAILURE);
	      A=NULL;
	    }
	  else
	    {

	      A=realloc_aln2 ( A, A->nseq+1, A->len_aln+1);
	      sprintf ( A->name[A->nseq], "%s", seq_name);
	      sprintf ( A->seq_al[A->nseq], "%s", seq);
	      A->nseq++;
	    }
	  return A;
	}


Alignment *aln2number (Alignment *A)
        {
	A->seq_al=char_array2number(A->seq_al, A->nseq);
	return A;
	}
Sequence *seq2number (Sequence *A)
        {
	A->seq=char_array2number(A->seq, A->nseq);
	return A;
	}
Alignment * aln2X (Alignment *A, int x)
{
  char *buf;
  int a,b,c,d;
  
  if (!A) return NULL;
  buf=(char*)vcalloc ((A->len_aln*x)+1, sizeof (char));
  A=realloc_aln (A, (A->len_aln)*3+1);
  for (a=0; a<=A->nseq; a++)
    {
      for (d=0,b=0; b<A->len_aln; b++)
	{
	  for (c=0; c<x; c++)
	    buf[d++]=A->seq_al[a][b];
	}
      buf[d]='\0';
      sprintf (A->seq_al[a], "%s", buf);
    }
  A->len_aln*=x;
  vfree(buf);
  return A;
}


Sequence * aln2seq (Alignment *A)
{
  return aln2seq_main(A, RM_GAP);
}
Sequence * aln2seq_main (Alignment *A, int mode)
	{
	Sequence *S;
	int a;
	
	
	if ( !A) return NULL;
	S=declare_sequence(0, 0, A->nseq);
	for (a=0; a<A->nseq; a++)
	  {
	    S->seq[a]=csprintf (S->seq[a], "%s", A->seq_al[a]);
	    S->name[a]=csprintf (S->name[a], "%s", A->name[a]);
	    S->seq_comment[a]=csprintf (S->seq_comment[a], "%s", A->aln_comment[a]);
	    S->aln_comment[a]=csprintf (S->aln_comment[a], "%s", A->aln_comment[a]);
	    
	    if (mode==RM_GAP)ungap(S->seq[a]);
	    S->len[a]=strlen (S->seq[a]);	    
	    S->file[a]=csprintf(S->file[a],"%s", A->file[a]);
	  }
	S->nseq=A->nseq;
	return S;
	}

Sequence  *keep_residues_in_seq ( Sequence *S, char *list, char replacement)
{
  Alignment *A=NULL;
  int a;

  A=seq2aln (S, A,1);
  A=keep_residues_in_aln ( A, list, replacement);
  for ( a=0; a< A->nseq; a++)
    {
      ungap (A->seq_al[a]);
      sprintf ( S->seq[a], "%s", A->seq_al[a]);
    }
  free_aln (A);
  return S;
}


Alignment *aln2short_aln ( Alignment *A, char *list, char *nnew, int spacer)
{
  int a, b, r, cl, l;
  char *buf;

  for ( a=0; a< A->nseq; a++)
    {
      buf=(char*)vcalloc ( strlen (A->seq_al[a])+1, sizeof (char));

      for (l=0,cl=0, b=0; b< A->len_aln; b++)
	{
	  r=A->seq_al[a][b];
	  if ( is_gap(r));
	  else if ( is_in_set (r, list))
	    {
	      if (cl){cl=0; buf[l++]=nnew[0];}
	      buf[l++]=r;
	    }
	  else
	    {
	      if ( cl==spacer){buf[l++]=nnew[0];cl=0;}
	      cl++;
	    }

	}

      buf[l]='\0';
      sprintf (A->seq_al[a], "%s", buf);
      vfree (buf);
    }
  return A;
}

Alignment *keep_residues_in_aln ( Alignment *A, char *list, char replacement)
{
  return filter_keep_residues_in_aln (A,NULL, 0, -1, list, replacement);
}
Alignment *filter_keep_residues_in_aln ( Alignment *A,Alignment *ST, int use_cons, int value, char *list, char replacement)
{
  char **sl;
  int n, a;

  n=strlen (list);
  sl=declare_char (n+1, 256);
  for (a=0; a< n; a++)
    sprintf ( sl[a], "%c%c", list[a], list[a]);
  sprintf ( sl[a],"#%c", replacement);
  A=filter_aln_convert (A, ST,use_cons,value, n+1, sl);
  free_char (sl, -1);
  return A;
}


Alignment *filter_convert_aln ( Alignment *A,Alignment *ST, int use_cons, int value, int n, ...)
{
  va_list ap;
  char **sl;
  int a;
  va_start (ap, n);
  sl=(char**)vcalloc ( n,sizeof(char*));
  for ( a=0; a< n; a++)
    {
      sl[a]=va_arg(ap, char * );
    }
  va_end(ap);
  A=filter_aln_convert (A,ST,use_cons,value, n,sl);
  vfree(sl);
  return A;
}

Alignment * filter_aln ( Alignment *A, Alignment *ST, int value)
        {
	  return filter_aln_convert (A, ST,0,value,DELETE, NULL);
	}
Alignment * filter_aln_switchcase ( Alignment *A, Alignment *ST,int use_cons, int value)
        {
	  return filter_aln_convert (A, ST,0,value,SWITCHCASE, NULL);
	}
Alignment * filter_aln_upper_lower ( Alignment *A, Alignment *ST,int use_cons, int value)
        {
	  return filter_aln_convert (A, ST,use_cons,value, LOWER, NULL);
	}
Alignment * filter_aln_lower_upper ( Alignment *A, Alignment *ST,int use_cons, int value)
        {

	  return filter_aln_convert (A, ST,use_cons,value, UPPER, NULL);
	}
Alignment * STseq2STaln ( Alignment *A, Alignment *ST)
        {
	  int a, i=0;

	  if  (ST && ST->len_aln !=A->len_aln)
		{
		  Sequence *S_T, *S_A;

		  S_T=aln2seq (ST);
		  S_A=aln2seq (A);

		  for (a=0; a< A->nseq; a++)
		    {
		      i=name_is_in_list (A->name[a], S_T->name,S_T->nseq, 100);
		      if (i!=-1)
			{
			  char *s1, *s2;
			  s1=(S_T)->seq[i];ungap(s1);
			  s2=(S_A)->seq[a];ungap(s2);

			  if ( strlen (s1)!=strlen(s2))
			    {
			      fprintf ( stderr, "%s\n%s\n", s1, s2);
			      printf_exit (EXIT_FAILURE, stderr, "ERROR: Sequence %s has different length in the alignment and in the  structure Alignment [FATAL:%s]\n", A->name[a], PROGRAM);
			    }
			}
		    }
		   ST=copy_aln (A, ST);
		   thread_seq_struc2aln (ST,S_T);
		}

	  return ST;
	}
Alignment * merge_annotation   ( Alignment *A, Alignment *ST, char *seq)
{
  int s, a, b;

  ST=STseq2STaln (A, ST);
  if ( seq==NULL)s=0;
  else
    s=name_is_in_list ( seq, A->name, A->nseq, 100);

  if (s==-1)
    {
      printf_exit(EXIT_FAILURE,stderr,"%s is not in your MSA [FATAL: %s]", PROGRAM);
    }

  for (a=0; a<A->len_aln; a++)
    {
      int t, r;

      t=A->seq_al[s][a];
      if (is_gap (t))continue;
      for (b=0; b<A->nseq; b++)
	{
	  t=A->seq_al[s][a];
	  r=ST->seq_al[b][a];
	  if ( isdigit (r))
	    {
	      if (!isdigit(t) || (isdigit (t) && t<r))
		A->seq_al[s][a]=r;
	    }
	}
    }
  return A;
}



Alignment * filter_aln_convert ( Alignment *A, Alignment *ST,int use_cons, int value, int n_symbol,char **symbol_list)
        {
	  int a, b, c;
	  int st;
	  int cons=0;


	  ST=STseq2STaln (A, ST);
	  if ( ST && use_cons)
	    {
	      cons=name_is_in_list ("con", ST->name,ST->nseq+1, 100);
	      if ( cons==-1)cons=name_is_in_list ("cons", ST->name,ST->nseq+1, 100);
	      if ( cons==-1)cons=name_is_in_list ("Cons", ST->name,ST->nseq+1, 100);
	      if ( cons==-1)
		{
		  use_cons=0;
		  fprintf (stderr, "WARNING: Could not Use the Consensus Sequence [WARNING:%s]\n", PROGRAM);
		}
	    }

	  A->residue_case=KEEP_CASE;
	  for ( a=0; a< A->nseq; a++)
	    {
	      if(value!=10 && ST && !use_cons)
		{
		  c=name_is_in_list (A->name[a], ST->name, ST->nseq,100);
		  if (c==-1)st=11;
		}

	      for ( b=0; b< A->len_aln; b++)
		{
		  if ( value==10 || !ST)st=11;
		  else if ( ST && use_cons)
		    {
		      st=(isdigit(ST->seq_al[cons][b]))?ST->seq_al[cons][b]-'0':ST->seq_al[cons][b];
		    }
		  else st=(isdigit(ST->seq_al[c][b]))?ST->seq_al[c][b]-'0':ST->seq_al[c][b];


		  if ( st==value || value==-1 || st==NO_COLOR_RESIDUE)
		    {
		      if      ( n_symbol==UPPER  && !symbol_list)A->seq_al[a][b]=toupper (A->seq_al[a][b]);
		      else if ( n_symbol==LOWER  && !symbol_list)A->seq_al[a][b]=tolower (A->seq_al[a][b]);
		      else if ( n_symbol==SWITCHCASE && !symbol_list)
			{
			  if ( !isalpha(A->seq_al[a][b]));
			  else if (isupper (A->seq_al[a][b]))A->seq_al[a][b]=tolower (A->seq_al[a][b]);
			  else if (islower (A->seq_al[a][b]))A->seq_al[a][b]=toupper (A->seq_al[a][b]);
			}
		      else if ( n_symbol==DELETE && !symbol_list)A->seq_al[a][b]='-';
		      else
			{
			  A->seq_al[a][b]=convert(A->seq_al[a][b],n_symbol,symbol_list);
			}
		    }

		}
	    }
	  return A;
	}


char ** sar_aln2motif (Alignment *A, Alignment *B, int *pos, int c);
char ** sar_aln2motif (Alignment *A, Alignment *B, int *pos, int c)
{
  static Alignment *I;
  static Alignment *O;
  int a, b, o, i;

  float tp,tn,fp,fn,best, sp, sn, sen2;
  float best_pred=-1;
  int best_motif=0;


  int n1;
  static char ***alp;
  static int *alp_size;

  char ***motif_list;
  int n;


  if (!I)
    {
      I=copy_aln(A, NULL);
      O=copy_aln(A, NULL);
    }



  I->nseq=O->nseq=I->len_aln=O->len_aln=0;
  for (a=0; a<A->len_aln; a++)
    {
      if (pos[a])
	{
	  for (i=o=0,b=0; b<A->nseq; b++)
	    {

	      if ( is_gap(A->seq_al[b][a]))return 0;
	      if (B->seq_al[b][c]=='I')I->seq_al[i++][I->len_aln]=A->seq_al[b][a];
	      else O->seq_al[o++][O->len_aln]=A->seq_al[b][a];
	    }
	  I->len_aln++;
	  O->len_aln++;
	}
    }

  if (O->len_aln==0 || I->len_aln==0) return 0;
  O->nseq=o;
  I->nseq=i;
  for (a=0; a<o; a++)O->seq_al[a][O->len_aln]='\0';
  for (a=0; a<i; a++)I->seq_al[a][I->len_aln]='\0';

  alp=(char***)vcalloc ( sizeof (char**), I->len_aln);
  alp_size= (int*)vcalloc ( I->len_aln, sizeof (int));
  for (a=0; a<I->len_aln; a++)
    {
      char *col;
      alp[a]=string2alphabet ( (col=aln_column2string (I,a)),2, &alp_size[a]);
      vfree (col);
    }



  motif_list=generate_array_string_list (I->len_aln, alp, alp_size, &n, NULL, OVERLAP);
  best_pred=best_motif=0;
  for (a=0; a<n; a++)
    {

      tp=tn=fp=fn=0;

      for (b=0; b<I->nseq; b++)
	{
	  if (match_motif (I->seq_al[b], motif_list[a]))tp++;
	  else fn++;
	}
      for (b=0; b<O->nseq; b++)
	{
	  if (match_motif (O->seq_al[b], motif_list[a]))fp++;
	  else tn++;
	}
      rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);

      if (best> best_pred)
	{
	  best_pred=best;
	  best_motif=a;
	}
    }

  output_Alignment_without_header ( I, stdout);
  fprintf ( stdout, "\n");
  output_Alignment_without_header ( O, stdout);


  fprintf ( stdout, "\nMotifCompound %d pred: %.2f motif: ", c, best_pred);
  for (n1=0, a=0; a<I->len_aln; a++)
    {
      char *m;
      int l;
      m=motif_list[best_motif][a];
      fprintf ( stdout, "[%s]-", m);
      l=strlen (m);
      n1+=(l==1 && !strm ("*",m) )?1:0;
    }
  fprintf (stdout, "SCORE: %d", n1);

  for (a=0; a<n; a++)vfree (motif_list[a]);
  vfree (motif_list);
  free_arrayN((void ***) alp, 3);
  vfree (alp_size);

  return NULL;
}




void explore_weight_matrix (Alignment *A, Alignment *B, int range, int n, int *array);
void explore_weight_matrix (Alignment *A, Alignment *B, int range, int n, int *array)
{
  int a;
  if ( n==A->len_aln)
    {
      fprintf ( stdout, "\n W:");
      for (a=0; a<A->len_aln; a++)fprintf ( stdout, "%d", array[a]);
      fprintf ( stdout, " %.4f",(float)sar_aln2r(A,B,array,0));
      return;
    }
  else
    {
      for ( a=0; a<range; a++)
	{
	  array[n]=a;
	  explore_weight_matrix (A, B, range, n+1, array);
	}
    }
}
float search_best_combo(Alignment *A, Alignment *B);
void search_best_combo_sar_aln(Alignment *A, Alignment *B);
void search_best_combo_sar_aln(Alignment *A, Alignment *B)
{
  int a,b,c;
  Alignment *S;
  float s;
  int w=5;

  S=copy_aln (B, NULL);
  S->len_aln=w;
  for ( a=0; a<B->len_aln-w;a++)
    {
      for (b=0; b<B->nseq; b++)
	{
	  for (c=0; c<w; c++)
	    {
	      S->seq_al[b][c]=B->seq_al[b][a+c];
	    }
	  S->seq_al[b][c]='\0';
	    }

      s=search_best_combo (A, S);
      fprintf ( stdout,"\nP: XXXX \nP: XXXXX A=%d / %d", a, B->len_aln);

    }

}

float search_best_combo(Alignment *A, Alignment *B)
{
  int a, b, c, d, best_pos,nl, max;
  float best_score, score;
  int *list, *pos;

  int  w;
  int combo_mode=1; //1: greedy 2: consider all thw w combinations;
  FILE *fp2;
  static int **M;
  max=2;
  int delta=0;
  w=1;

  pos=(int*)vcalloc ( A->len_aln, sizeof (int));
  list=(int*)vcalloc (A->len_aln, sizeof (int));
  nl=0;

  if ( combo_mode==1)
    {
      for (a=0; a< max; a++)
	{
	  for (best_score=-9999,best_pos=0,b=0; b< A->len_aln-w; b++)
	    {
	      for (c=0; c<nl; c++)pos[list[c]]=1;
	      for (c=0; c<w; c++)pos[b+c]=1;
	      score=sar_aln2r(A,B,pos,0);
	      if ( score>best_score)
		{
		  best_score=score;
		  best_pos=b;
		}
	      for (c=0; c<w; c++)pos[b+c]=0;
	    }
	  if (best_pos==list[nl-1])break;
	  list[nl++]=best_pos;
	  for (b=0; b<nl; b++) pos[list[b]]=1;
	  fprintf ( stdout, "\n%2d P: %d S:%.3f Delta= %d", nl,best_pos, best_score, (int)sar_aln2delta(A,B, pos,0));
	  for (b=0; b<nl; b++) pos[list[b]]=0;


	}
      for (a=0; a<nl; a++) pos[list[a]]=1;
      fprintf ( stdout, "\nR: %3f " ,(float)sar_aln2r(A,B,pos,1));

    }
  else if ( combo_mode==2)
    {
      int  *array;
      char *tmpf;
      FILE *fp;
      char *buf=NULL;
      int *preset,  n_preset;

      tmpf=vtmpnam (NULL);
      max=1;
      generate_array_int_list (max, 0,A->len_aln-1, 1,NULL, tmpf);
      printf_system ( "cp %s testfile", tmpf);
      buf=(char*)vcalloc ( 1000, sizeof (char));
      fp=vfopen (tmpf, "r");
      best_score=-99999;

      n_preset=0;
      preset=(int*)vcalloc (A->len_aln, sizeof (int));
      preset[n_preset++]=353;
      preset[n_preset++]=361;
      //preset[n_preset++]=365;
      //preset[n_preset++]=187;
      //preset[n_preset++]=397;
      //preset[n_preset++]=492;


      while ( (buf=vfgets ( buf, fp))!=NULL)
	{

	  array=string2num_list (buf);

	  for (a=1; a<=max; a++)
	    {
	      pos[array[a]]=1;
	    }
	  for ( a=0; a<n_preset; a++)pos[preset[a]]=1;

	  score=sar_aln2r(A,B,pos,0);

	  if ( score>best_score)
	    {
	      best_score=score;
	      fprintf ( stdout, "\n");
	      for (a=0; a<n_preset; a++)fprintf (stdout, "%2d ", preset[a]);
	      for (a=1; a<=max; a++)fprintf (stdout, "%2d ", array[a]);
	      fprintf ( stdout, " R: %.3f", best_score);
	      for (nl=0,a=0; a<n_preset; a++)list[nl++]=preset[a];
	      for (a=1; a<=max; a++)list[nl++]=array[a];
	    }
	  //if ( score!=0)HERE ("R=%.2f", score);
	  for (b=1; b<=max; b++)
	    pos[array[b]]=0;
	  vfree (array);
	}
      fprintf ( stdout, "\n");
      vfclose (fp);
      //for (a=0; a<max; a++)fprintf (stdout, "%2d ", array[best_pos][a]);
      //fprintf ( stdout, " R: %.3f", best_score);
    }
  for (c=0; c<B->len_aln; c++)
    {
      sar_aln2motif (A,B,pos, c);

    }
  myexit (EXIT_FAILURE);
  HERE ("***************");
  fp2=vfopen ("aln.aln", "w");
  for (a=0; a<A->nseq; a++)
    {
      fprintf (fp2, ">%s\n", A->name[a]);
      for ( b=0; b<nl; b++)fprintf (fp2, "%c", A->seq_al[a][list[b]]);
      fprintf ( fp2, "\n");
    }
  vfclose (fp2);
  HERE ("Output aln.aln");
  if (1)
    {
      float tp=0, tn=0, fp=0, fn=0, pp2=0,pp=0, sn,sn2, sp;
      int **result,**result2,**compound_score, *ref_score,n2,n, s, p, c;
      Alignment *AI, *AO;
      int simI, simO;

      compound_score=declare_int (B->len_aln, 2);
      ref_score=(int*)vcalloc (nl, sizeof (int));

      result=declare_int (B->len_aln*A->nseq*A->nseq, 2);
      result2=declare_int (B->len_aln*A->nseq*A->nseq, 2);

      for (n2=c=0; c< B->len_aln; c++)
	{

	  int sar1, sar2;
	  pp=tp=tn=fp=fn=0;
	  if (!M)M=read_matrice ("blosum62mt");
	  for (n=0,a=0; a<A->nseq-1; a++)
	    {
	      for (b=a+1; b<A->nseq;b++)
		{
		  for (s=0,p=0; p<nl; p++)
		    {
		      char r1, r2;

		      r1=A->seq_al[a][list[p]];
		      r2=A->seq_al[b][list[p]];
		      if ( !is_gap (r1) && !is_gap(r2))s+=M[r1-'A'][r2-'A'];
		    }
		  result2[n2][0]=result[n][0]=s;

		  sar1=B->seq_al[a][c];sar2=B->seq_al[b][c];

		  if (sar1=='I' && sar1==sar2)
		    {
		      result2[n2][1]=result[n][1]=1;
		      pp++;pp2++;
		      n++;n2++;
		    }
		  else if ( sar1==sar2 && sar1=='O')
		    {
		      ;
		    }
		  else
		    {
		      result2[n2][1]=result[n][1]=0;
		      n++;n2++;
		    }
		  //else if ( s1==s2=='O')result[n][1]=-1;
		}
	    }

	  if (pp==0)continue;
	  sort_int_inv (result, 2, 0, 0, n-1);


	  for (tp=0,a=0; a<n; a++)
	    {
	      tp+=result[a][1];
	      if ((pp-tp) == (a-tp))break;
	    }
	  fp=a-tp;
	  fn=pp-tp;
	  tn=n-pp;

	  sn=(tp/(tp+fn));
	  sn2=(tp/(tp+fp));
	  sp=(tn/(tn+fp));
	  fprintf ( stdout, "\nCompound %3d sn: %.3f sn2: %.3f sp: %.3f MIN: %.3f",c,sn, sn2,sp, MIN((MIN(sn,sn2)),sp));
	  compound_score[c][0]=c;
	  compound_score[c][1]=1000*MIN((MIN(sn,sn2)),sp);
	}

      sort_int_inv (compound_score,2, 1, 0, B->len_aln-1);

      fp2=vfopen ("compound.fasta", "w");
      for (d=0; d<nl; d++)
	{
	  int r1, r2;
	  for (n=0,a=0;a<A->nseq; a++)
	    for (b=0; b<A->nseq; b++)
	      {
		r1= A->seq_al[b][list[d]];
		r2= A->seq_al[b][list[d]];
		if (is_gap(r1) || is_gap(r2))continue;
		else
		  {
		    ref_score[d]+=M[r1-'A'][r2-'A'];
		    n++;
		  }
	      }
	  ref_score[d]/=n;
	}
      AO=copy_aln (A, NULL);
      AI=copy_aln (A,NULL);
      AO->len_aln=AI->len_aln=nl;
      for (a=0; a<A->nseq; a++)AO->seq_al[a][nl]=AI->seq_al[a][nl]='\0';

      for (a=0; a<B->len_aln; a++)
	{
	  fprintf (stdout, "\n>%4d %4d ", compound_score[a][0], compound_score[a][1]);
	  for (b=0; b<B->nseq; b++) fprintf (stdout, "%c", B->seq_al[b][compound_score[a][0]]);
	  fprintf ( stdout, "\n");

	  for (AI->nseq=0,b=0; b<B->nseq; b++)
	    {
	      if (B->seq_al[b][compound_score[a][0]]=='O')continue;
	      fprintf ( stdout, "\n\t");
	      for (c=0; c<nl; c++)
		{
		  fprintf ( stdout, "%c", A->seq_al[b][list[c]]);
		  AI->seq_al[AI->nseq][c]=A->seq_al[b][list[c]];
		}
	      AI->nseq++;
	    }
	  fprintf ( stdout, "\n\t");
	  for (d=0; d<nl; d++)
	    {
	      for (score=0,n=0,b=0; b<B->nseq; b++)
		{
		  if (B->seq_al[b][compound_score[a][0]]=='O')continue;
		  for (c=0; c<B->nseq; c++)
		    {
		      if (B->seq_al[c][compound_score[a][0]]=='O')continue;
		      {
			int r1, r2;

			r1= A->seq_al[b][list[d]];
			r2= A->seq_al[b][list[d]];
			if (is_gap(r1) || is_gap(r2))continue;
			else score+=M[r1-'A'][r2-'A'];
			n++;
		      }
		    }
		}
	      score/=n;
	      if ((float)score/(float)ref_score[d]>1.2)fprintf ( stdout, "*");
	      else fprintf ( stdout, " ");
	    }
	  for (AO->nseq=0,b=0; b<B->nseq; b++)
	    {
	      if (B->seq_al[b][compound_score[a][0]]=='I')continue;
	      fprintf ( stdout, "\n\t");
	      for (c=0; c<nl; c++)
		{
		  AO->seq_al[AO->nseq][c]=A->seq_al[b][list[c]];
		  fprintf ( stdout, "%c", A->seq_al[b][list[c]]);
		}
	      AO->nseq++;
	    }
	  simI=aln2sim (AI, "blosum62mt"); simO=aln2sim (AO, "blosum62mt");
	  fprintf ( stdout, "\nDELTA: I: %d O: %d %d",simI,simO, simI-simO);
	  delta+=simI-simO;
	}

      for ( a=0; a<B->nseq; a++)
	{

	  fprintf ( fp2, ">%s\n", B->name[a]);
	  for (b=0; b<B->len_aln/2; b++)
	    fprintf ( fp2, "%c", B->seq_al[a][compound_score[b][0]]);
	  fprintf (fp2, "\n");
	}
      vfclose (fp2);
      HERE ("OUTPUT compound.fasta");
      result=result2;
      n=n2;
      pp=pp2;

      sort_int_inv (result, 2, 0, 0, n-1);


      for (tp=0,a=0; a<n; a++)
	{
	  tp+=result[a][1];
	  if ((pp-tp) == (a-tp))break;
	}
      fp=a-tp;
      fn=pp-tp;
      tn=n-pp;

      sn=(tp/(tp+fn));
      sn2=(tp/(tp+fp));
      sp=(tn/(tn+fp));
      fprintf ( stdout, "\nTOT:  sn: %.3f sn2: %.3f sp: %.3f MIN: %.3f",sn, sn2,sp, MIN((MIN(sn,sn2)),sp));

    }
  HERE ("Delta= %d", delta);


  /*
  C=copy_aln(A, NULL);
  for (a=0; a< nl; a++)
    for (b=0; b<A->nseq; b++)
      C->seq_al[b][a]=A->seq_al[b][list[a]];
  C->len_aln=nl;
  array=vcalloc (C->len_aln, sizeof (int));
  explore_weight_matrix (C, B, 6,0, array);
  */

  return best_score;
}


void count_misc (Alignment *A, Alignment *B)
{
  int **done, a, b, c, d, e,f, g, *list, n, score;
  double **slist, *r;
  int *pos;
  int w=1;

  search_best_combo (A,B);
  myexit (EXIT_FAILURE);
  pos=(int*)vcalloc (A->len_aln+1, sizeof (int));
  /*
  pos[354]=1;
  pos[362]=1;
  pos[366]=1;
  pos[398]=1;
  pos[476]=1;


  fprintf ( stdout, "\nR: %3f " ,(float)sar_aln2r(A,B,pos,1));myexit (EXIT_FAILURE);
  */
  for (a=0; a< A->len_aln-w; a++)
    {
      for (c=0; c<w; c++)
	{
	  pos[a+c]=1;
	}
      pos[398]=1;
      pos[362]=1;
      pos[354]=1;
      pos[366]=1;
      pos[419]=1;
      pos[494]=1;
      pos[476]=1;
      pos[337]=1;
      fprintf ( stdout, "\nP: %3d  W:2 R: %3f ",a+1, (float)sar_aln2r(A,B,pos,0));
      for (c=0; c<w; c++)
	  {
	    pos[a+c]=0;
	  }
    }

  myexit (EXIT_FAILURE);
  for (a=0; a<w; a++) pos[a]=1;
  for (a=w; a< A->len_aln-1; a++)
    {
      pos[a-w]=0;
      pos[a]=1;
      fprintf ( stdout, "\nP: %3d W:2 R: %3f ",a, (float)sar_aln2r(A,B,pos,0));
    }

  myexit (EXIT_FAILURE);
  pos[2]=1;
  pos[3]=1;



  explore_weight_matrix (A, B,3, 0,pos);
  myexit (EXIT_FAILURE);

  for (a=0; a<A->len_aln; a++)
    for ( b=0; b<A->len_aln; b++)
      for (c=0; c<A->len_aln; c++)
	for (d=0; d<A->len_aln; d++)
	  for (f=0; f<A->len_aln; f++)
	    for (g=0; g<A->len_aln; g++)
	    {
	      e=0;
	      pos[e++]=a;
	      pos[e++]=b;
	      pos[e++]=c;
	      pos[e++]=d;
	      pos[e++]=f;
	      pos[e++]=g;
	      pos[e++]=-1;
	      fprintf ( stdout, "\n%d %d %d %d %d %d %.3f", a, b,c,d,f, g, sar_aln2r(A,B, pos,0));

	    }

  myexit (EXIT_FAILURE);


  slist=declare_double (A->nseq*A->nseq*10, 2);
  done=declare_int (256, 256);
  list=(int*)vcalloc ( A->nseq, sizeof (int));

  for (a=0; a<A->len_aln-1; a++)
    {
      for (b =0; b<256; b++)for (c=0; c<256; c++)done[b][c]=0;

      for (b=0; b<A->nseq-1; b++)
	{
	  int r1, r2;
	  r1=A->seq_al[b][a];
	  r2=A->seq_al[b][a+1];
	  if (done[r1][r2])continue;
	  n=0;
	  done[r1][r2]=1;
	  list[n++]=b;
	  fprintf ( stdout, "\n%3d %c%c: %s ",a+1, r1, r2, A->name[b]);
	  for ( c=b+1; c<A->nseq; c++)
	    {
	      if (r1==A->seq_al[c][a] && r2==A->seq_al[c][a+1])
		{
		  fprintf ( stdout, "%s ", A->name[c]);
		  list[n++]=c;
		}

	    }
	  if (B && n>1)
	    {
	      for (e=0,score=0,c=0; c<n-1; c++)
		for (d=c+1; d<n; d++,e++)
		  score+=get_sar_sim2(B->seq_al[list[c]], B->seq_al[list[d]]);
	      fprintf ( stdout, " Score=%d", score/e);
	    }
	}
    }
  for (score=0,e=0,a=0; a<A->nseq-1; a++)
    for (b=a+1; b<A->nseq; b++,e++)
      {
	score+=get_sar_sim2(B->seq_al[a], B->seq_al[b]);
      }
  fprintf  (stdout,"AVG=%d", score/e);
  for (n=0,a=0; a< A->nseq-1; a++)
    {
      static int **M;
      int sim;
      if (!M)M=read_matrice ("blosum62mt");


      for (b=a+1; b<A->nseq; b++)
	{
	  int n11, n01, n10, n00, n1;

	  for (sim=d=0;d<A->len_aln; d++)
	    {
	      int r1, r2;
	      r1=A->seq_al[a][d];
	      r2=A->seq_al[b][d];
	      sim+=(r1==r2)?1:0;
	      //sim +=(M[r1-'A'][r2-'A']>0)?1:0;
	    }

	  sim=(100*sim)/(A->len_aln);//+rand()%10;
	  for (n1=n00=n11=n10=n01=score=0, d=0; d<B->len_aln; d++)
	    {
	      int r1, r2;
	      r1=B->seq_al[a][d];
	      r2=B->seq_al[b][d];
	      n11+=(r1=='I' && r2=='I');
	      n00+=(r1=='O' && r2=='O');
	      n10+=(r1=='I' && r2=='0');
	      n01+=(r1=='O' && r2=='I');
	      n1+=(r1=='I' || r2=='I');
	    }
	  score =((n11+n00)*100)/B->len_aln;

	  //score=get_sar_sim2(B->seq_al[a], B->seq_al[b]);

	  fprintf ( stdout, "\nSIM: %d SC: %d", sim, score);
	  slist[n][0]=(double)sim;
	  slist[n][1]=(double)score;
	  n++;
	}
    }
  r=return_r(slist, n);
  fprintf ( stdout, "\nR= %.4f", (float)r[0]);
  myexit (EXIT_FAILURE);
}

int aln2ngap ( Alignment *A)
{
  int ngap=0, a, b;
  for (a=0; a< A->len_aln; a++)
    for (b=0; b<A->nseq; b++) ngap+=is_gap (A->seq_al[b][a]);
  return ngap;
}
int  * count_in_aln ( Alignment *A, Alignment *ST, int value, int n_symbol,char **symbol_list, int *table)
        {
	  int a, b, c=0, d;
	  int st;

	  if (!table)table=(int*)vcalloc (n_symbol, sizeof (int));

	  A->residue_case=KEEP_CASE;
	  for ( a=0; a< A->nseq; a++)
	        {
		if(value!=10 && ST)for ( c=0; c< ST->nseq; c++)if ( strm(ST->name[c], A->name[a]))break;
		for ( b=0; b< A->len_aln; b++)
		    {
		      if ( value==10 || !ST)st=11;
		      else st=(isdigit(ST->seq_al[c][b]))?ST->seq_al[c][b]-'0':ST->seq_al[c][b];
		      if ( st==value || value==-1)
			{
			  for ( d=0; d<n_symbol; d++)table[d]+=is_in_set ( A->seq_al[a][b], symbol_list[d]);
			}
		    }
		}
	  return table;
	}

char *dna_aln2cons_seq ( Alignment *A)
        {
	int a, b, best;
	static int **column_count;
	static int **old_tot_count;
	static int **new_tot_count;
	static char *string1, *string2;
	int **count_buf;
	char r1, r2,*seq;
	int NA=0, NG=1, NC=2, NT=3, IGAP=4;
	static int   MAX_EST_SIZE=10000;
	static int   size_increment=1000;
	static int first;
	int overlap=0, best_overlap=0;


	seq=(char*)vcalloc ( A->len_aln+1, sizeof (char));

	if (!column_count )
	  {
	    column_count=(int**)vcalloc(MAX_EST_SIZE, sizeof (int*));
	    for ( a=0; a< MAX_EST_SIZE; a++)
	      column_count[a]=(int*)vcalloc (5, sizeof (int));

	    old_tot_count=(int**)vcalloc(MAX_EST_SIZE, sizeof (int*));
	    new_tot_count=(int**)vcalloc(MAX_EST_SIZE, sizeof (int*));
	    A->P=declare_profile( "agct-",MAX_EST_SIZE);
	    string1=(char*)vcalloc (MAX_EST_SIZE, sizeof (char));
	    string2=(char*)vcalloc (MAX_EST_SIZE, sizeof (char));
	  }
	else if (A->len_aln>MAX_EST_SIZE)
	  {
	    if ( column_count)
	      {
		for ( a=0; a< MAX_EST_SIZE; a++)
		  vfree(column_count[a]);
		vfree(column_count);
		vfree(old_tot_count);
		vfree(new_tot_count);
		vfree(string1);
		vfree(string2);
	      }

	  column_count=(int**)vcalloc(MAX_EST_SIZE+ size_increment, sizeof (int*));
	  for ( a=0; a< MAX_EST_SIZE+ size_increment; a++)
	      column_count[a]=(int*)vcalloc (5, sizeof (int));

	  old_tot_count=(int**)vcalloc(MAX_EST_SIZE+ size_increment, sizeof (int*));
	  new_tot_count=(int**)vcalloc(MAX_EST_SIZE+ size_increment, sizeof (int*));

	  for (a=0; a< MAX_EST_SIZE; a++)
	    {
	      old_tot_count[a]=*(column_count++);
	      for ( b=0; b<5; b++)old_tot_count[a][b]=(A->P)->count[b][a];
	    }
	  free_int ( (A->P)->count, -1);

	  (A->P)->count=declare_int (5, MAX_EST_SIZE+ size_increment);
	  (A->P)->max_len=MAX_EST_SIZE+ size_increment;
	  MAX_EST_SIZE+= size_increment;
	  string1=(char*)vcalloc (MAX_EST_SIZE, sizeof (char));
	  string2=(char*)vcalloc (MAX_EST_SIZE, sizeof (char));
	  }


	sprintf ( string1, "%s",A->seq_al[0]);
	sprintf ( string2, "%s",A->seq_al[1]);


	string1=mark_internal_gaps(string1,'.');
	string2=mark_internal_gaps(string2,'.');



	for (b=0,a=0; a< A->len_aln; a++)
	  {
	    r1=string1[a];
	    r2=string2[a];

	    if ( r1==r2)
	      {
		overlap++;
	      }
	    else
	      {
		best_overlap=MAX(overlap, best_overlap);
		overlap=0;
	      }


	    if (!is_gap(r1) && first==1)new_tot_count[a]=old_tot_count[b++];
	    else if (is_gap(r1) || first==0){new_tot_count[a]=*column_count;column_count++;};

	    if ( first==0)
	      {
		if(r1=='a')       new_tot_count[a][NA]++;
		else if ( r1=='g')new_tot_count[a][NG]++;
		else if ( r1=='c')new_tot_count[a][NC]++;
		else if ( r1=='t')new_tot_count[a][NT]++;
		else if (is_gap(r1));
		else
		  {
		   new_tot_count[a][NA]++;
		   new_tot_count[a][NG]++;
		   new_tot_count[a][NC]++;
		   new_tot_count[a][NT]++;
		  }
	      }
	    if ( a> 0 && a<A->len_aln-1 && r1=='.')
	      {
		new_tot_count[a][IGAP]+=((new_tot_count[a-1][NA]+new_tot_count[a-1][NG]+new_tot_count[a-1][NC]+new_tot_count[a-1][NT]));
	      }


	    if(r2=='a')       new_tot_count[a][NA]++;
	    else if ( r2=='g')new_tot_count[a][NG]++;
	    else if ( r2=='c')new_tot_count[a][NC]++;
	    else if ( r2=='t')new_tot_count[a][NT]++;
	    else if ( r2=='.')new_tot_count[a][IGAP]++;
	    else if ( r2=='-');
	    else
	      {
		new_tot_count[a][NA]++;
		new_tot_count[a][NG]++;
		new_tot_count[a][NC]++;
		new_tot_count[a][NT]++;
	      }
	    (A->P)->count[0][a]=new_tot_count[a][NA];
	    (A->P)->count[1][a]=new_tot_count[a][NG];
	    (A->P)->count[2][a]=new_tot_count[a][NC];
	    (A->P)->count[3][a]=new_tot_count[a][NT];
	    (A->P)->count[4][a]=new_tot_count[a][IGAP];

	    best_int(4,1, &best,new_tot_count[a][NA], new_tot_count[a][NG],new_tot_count[a][NC],new_tot_count[a][NT]);
	    if( best==0)      seq[a]='a';
	    else if ( best==1)seq[a]='g';
	    else if ( best==2)seq[a]='c';
	    else if ( best==3)seq[a]='t';
	  }

	first=1;

	seq[a]='\0';
	fprintf ( stderr, "[Best Overlap: %d Residues]", best_overlap);
	count_buf=old_tot_count;
	old_tot_count=new_tot_count;
	new_tot_count=count_buf;

	return seq;

	}

char *aln2cons_maj ( Alignment *A, int ns, int *ls, int n_groups, char **group_list)
        {
	char *seq;
	int a, b;
	int len;
	int clean_ls=0;
	static int *aa;

	if ( !aa) aa=(int*)vcalloc (1000, sizeof (int));

	len=strlen  (A->seq_al[ls[0]]);
	seq=(char*)vcalloc (len+1, sizeof (char));

	if ( ns==0)
	  {
	    ns=A->nseq;
	    ls=(int*)vcalloc ( A->nseq, sizeof (int));
	    for ( a=0; a< A->nseq; a++)ls[a]=a;
	    clean_ls=1;
	  }

	for ( a=0; a<len; a++)
	    {
	      int best_s=0, best_aa=0, r;
	      for (b=0; b< ns; b++)
		    {
		      r=tolower(A->seq_al[ls[b]][a]);
		      aa[r]++;
		      if (!is_gap(r) && aa[r]>best_s)
			{
			  best_s=aa[r];
			  best_aa=r;
			}
		      seq[a]=best_aa;
		    }
		for (best_s=0, best_aa=0,b=0; b< ns; b++)
		  {
		    aa[tolower(A->seq_al[ls[b]][a])]=0;
		  }
	    }
	if ( clean_ls)vfree(ls);
	seq[a]='\0';

	return seq;
	}





char *aln2cons_seq ( Alignment *A, int ns, int *ls, int n_groups, char **group_list)
        {
	char *seq;
	int a, b, c;
	int best_group=0;
	int aa_group=0;
	int *group;
	int len;
	int clean_ls=0;

	len=strlen  (A->seq_al[ls[0]]);
	seq=(char*)vcalloc (len+1, sizeof (char));

	if ( ns==0)
	  {
	    ns=A->nseq;
	    ls=(int*)vcalloc ( A->nseq, sizeof (int));
	    for ( a=0; a< A->nseq; a++)ls[a]=a;
	    clean_ls=1;
	  }


	if ( !group_list)
	   {
	       group_list=declare_char ( 26, 2);
	       for ( a=0; a<26; a++)group_list[a][0]=a+'a';
	       n_groups=26;
	       aa_group=1;
	   }


	for ( a=0; a<len; a++)
	    {
		group=(int*)vcalloc (n_groups+1, sizeof (int));
		for (best_group=0,b=0; b< ns; b++)
		    {
		    if ( !is_gap(A->seq_al[ls[b]][a]))
			 {
			 for (c=0; c< n_groups; c++)
			     if ( is_in_set (tolower(A->seq_al[ls[b]][a]), group_list[c]))
			                 {group[c]++;
					  best_group=(group[c]>group[best_group])?c:best_group;
					 }
			 }
		    seq[a]=group_list[best_group][0];
		    }
		vfree (group);
	    }
	seq[a]='\0';
	if ( aa_group) free_char (group_list, -1);

	if ( clean_ls)vfree(ls);

	return seq;
	}

Alignment *aln2conservation ( Alignment *A, int threshold,char *seq)
{
  int a, b, c, d, i, c1, c2;
  int   *pos;
  float *eval;
  float tot=0;
  float tn=0;
  int **sim;
  int w=0;

  pos =(int*)vcalloc (A->len_aln, sizeof (int));
  eval=(float*)vcalloc (A->len_aln, sizeof (int));
  sim=aln2sim_mat (A, "idmat");
  if (seq)i=name_is_in_list (seq, A->name, A->nseq, 100);
  else i=0;

  if ( i==-1) {HERE ("%s is an unknown:sequence [FATAL]"); myexit (EXIT_FAILURE);}

  for (a=0; a<A->len_aln; a++)
    {
      double s;
      int e;
      for (c=0,e=a-w; e<=a+w; e++)
	{
	  if (e<0 || e==A->len_aln)continue;
	  c1=toupper (A->seq_al[i][e]);
	  for (b=0; b<A->nseq; b++)
	    {
	      c2=toupper (A->seq_al[b][a]);
	      if (c1==c2)
		{
		  c++;
		  s=(double)((double)sim[i][b]/(double)(100));

		}
	      else
		{
		  s=(double)(((double)100-(double)sim[i][b])/(double)(100));
		}
	      eval[a]+=(s==0)?0:log(s);
	    }
	}
      pos[a]=(c*100)/A->nseq;
      if (!is_gap(c1)){tot+=pos[a]; tn++;}

      if (pos[a]>=threshold)A->seq_al[i][a]=toupper (A->seq_al[i][a]);
      else A->seq_al[i][a]=tolower (A->seq_al[i][a]);
    }
  fprintf (stdout, ">%s %s [i=%d]\n%s\n", A->name[i],A->aln_comment[i],i, A->seq_al[i]);
  tot=(tn>0)?(float)tot/(float)tn:0;

  for (d=0,a=0; a<A->len_aln; a++)
    {
      fprintf (stdout, "# %c %4d", A->seq_al[i][a],pos[a]);


      if ( !is_gap (A->seq_al[i][a]))
	{
	  fprintf (stdout, " LogOdd: %6.2f ", (tot==0 || pos[a]==0)?0:(float)log((float)pos[a]/tot));
	  fprintf ( stdout, " Pos: %5d E-Val: %9.2f", ++d, eval[a]/(A->nseq));
	}
      fprintf ( stdout, "\n");
    }
  fprintf ( stdout, "#average conservation: %.2f", tot);
  myexit (EXIT_SUCCESS);
}


char *aln2cons_seq_cov (Alignment *A)
{
  int *score;
  char *cons;
  int pleft,n,left, a, b;
  

  cons=(char*)vcalloc (A->len_aln+1, sizeof (int));
  score=(int*)vcalloc (A->nseq, sizeof (int));
  left=A->len_aln;
  

  //1 deal with empty columns --- there should not be any
  for (a=0; a<A->len_aln; a++)
    {
      n=0;
      for (b=0; b<A->nseq; b++)
	if (!is_gap(A->seq_al[b][a]))n++;
      if (n==0)
	{
	  cons[a]='-';
	  left--;
	}
    }
  //add residues from the sequences covering max.
  pleft=left;
  while (left>0)
    {
      int best_score=0;
      int best_seq=0;
      
      for (a=0; a<A->nseq; a++)score[a]=0;
      for (a=0; a<A->len_aln; a++)
	{
	  if (!cons[a])
	    {
	      n=0;
	      for (b=0; b<A->nseq; b++)
		{
		  if (!is_gap(A->seq_al[b][a]))n++;
		}
	      for (b=0; b<A->nseq; b++)
		{
		  if (!is_gap(A->seq_al[b][a]))
		    {
		      score[b]+=n;
		      if (score[b]>best_score){best_score=score[b];best_seq=b;}
		    }
		}
	    }
	}
      for (a=0; a<A->len_aln; a++)
	{
	  if (!cons[a] && !is_gap(A->seq_al[best_seq][a]))
	    {
	      cons[a]=A->seq_al[best_seq][a];
	      left--;
	    }
	}
      
      if (left==pleft)
	{
	  fprintf ( stderr, "\nCONS: %s\n", cons);
	  print_aln (A);
	  printf_exit (EXIT_FAILURE, stderr, "\nCould not Build profile consensus::aln2cons_seq_cov [FATAL:%s]", PROGRAM);
	}
	  
      pleft=left;
    }
  vfree(score);
  return cons;
}
		

char *aln2cons_seq_mat ( Alignment *A, char *mat_name)
{
  return sub_aln2cons_seq_mat (A, A->nseq, NULL, mat_name);
}
char *sub_aln2cons_seq_mat2 ( Alignment *A,int ns, char **ls, char *mat_name)
{
  char *cons;
  int *list;
  list=name_array2index_array(ls, ns, A->name, A->nseq);
  cons=sub_aln2cons_seq_mat  ( A,ns, list, mat_name);
  vfree (list);
  return cons;
}

char *sub_aln2cons_seq_mat  ( Alignment *A,int ns, int *ls, char *mat_name)
{
 int a, b, c, s;
 char *seq, r1, r2;
 int **mat;
 int score=0, best_score=0, best_r=0;
 int len;
 int naa;

 mat=read_matrice (mat_name);
 len=strlen ( A->seq_al[(ls==NULL)?0:ls[0]]);
 seq=(char*)vcalloc (len+1, sizeof (char));
 for ( a=0; a<len; a++)
   {
     for (b=0; b<20; b++)
       {
	 r1=tolower(AA_ALPHABET[b]);
	 for ( naa=0,score=0,c=0; c<ns; c++)
	   {
	     s=(ls==NULL)?c:ls[c];
	     if ( ls && ls[c]==-1) continue;
	     else if (is_gap(A->seq_al[s][a]))continue;
	     else
	       {
		 naa++;
		 r2=tolower(A->seq_al[s][a]);
		 if (!is_gap(r2))score+=mat[r1-'A'][r2-'A'];
	       }
	   }
	 if (naa==0)best_r='-';
	 if ( b==0 || score>best_score){best_score=score; best_r=r1;}
       }
     seq[a]=best_r;
   }
 free_int (mat, -1);
 return seq;
}

int  seq_list2in_file ( TC_method *M, Sequence *S, char *list, char *file)
{
  X_template *T=NULL;

  if ( !S)return 0;
  else
    {
      int t;
      t=tolower(M->seq_type[0]);

      if ( t=='s')
	{
	  return seq_list2fasta_file ( S, list, file, M->out_mode);

	}
      else
	{
	  FILE *fp, *fp2;
	  int a, n, s, c;
	  int *slist;



	  fp=vfopen ( file, "w");
	  slist=string2num_list (list);
	  n=slist[0];

	  if (strlen (M->seq_type) >1)
	    {
	      printf_exit( EXIT_FAILURE,stderr, "\Mixed seq_type not supported for external methods\n[FATAL:%s]", PROGRAM);
	    }

	  for ( a=2; a<n; a++)
	    {
	      s=slist[a];
	      if (t=='p')T=(S->T[s])->P;
	      else if (t=='r')T=(S->T[s])->R;
	      else if (t=='g')T=(S->T[s])->G;

	      if (!T && t=='r')
		{
		  fprintf ( fp, ">%s\n%s%s", S->name[s], S->seq[s], LINE_SEPARATOR);
		}
	      else if ( T && T->template_file && T->template_file[0])
		{
		  fp2=vfopen (T->template_file, "r");
		  while ( (c=fgetc (fp2))!=EOF)
		    {
		      fprintf ( fp, "%c", c);
		    }
		  fprintf (fp, "%s", LINE_SEPARATOR);
		  vfclose (fp2);
		}
	    }

	  fprintf (fp, "TARGET_SEQ_NAME: ");
	  for (a=2; a<n; a++)fprintf ( fp, "%s ", (S->name[slist[a]]));
	  fprintf ( fp, "%s", LINE_SEPARATOR);

	  vfclose (fp); vfree (slist);

	}

      return 1;
    }
}

int  seq_list2fasta_file( Sequence *S,  char *list, char *file, char *outmode)
        {
	FILE *fp;
	int n, a, s;
	static char *buf;
	static int blen;
	int l;
	//out_mode: names can only be re-converted when out mode is aln

	/*Buf is used because cmalloced functions cannot go through strtok*/
	if ( !S)return 0;
	else
	  {
	    fp=vfopen ( file, "w");
	    if ( !list)
	      {
		for ( a=0; a<S->nseq; a++)
		  {
		    if (outmode && strm (outmode, "aln"))fprintf ( fp, ">%s %s\n%s\n", decode_name (S->name[a], CODE),S->name[a], S->seq[a]);
		    else fprintf ( fp, ">%s %s\n%s\n", S->name[a],S->name[a], S->seq[a]);
		  }
	      }
	    else
	      {
		int **list2;
		int max;

		l=strlen (list);
		if ( l>blen)
		  {
		    if (buf)vfree(buf);
		    buf=(char*)vcalloc ( strlen (list)+1, sizeof (char));
		    sprintf ( buf, "%s", list);
		    blen=l;
		  }
		n=atoi(strtok (list,SEPARATORS));

		list2=declare_int (n, 2);
		max=n*1000;
		for ( a=0; a<n; a++)
		  {
		    list2[a][0]=atoi(strtok (NULL, SEPARATORS));
		    list2[a][1]=rand()%max;
		  }
		if ( atoigetenv ("HoT_4_TCOFFEE"))sort_int ( list2,2, 1, 0, n-1);
		for ( a=0; a< n; a++)
		  {
		    int i=list2[a][0];
		    if (outmode && strm (outmode, "aln"))fprintf ( fp, ">%s %s\n%s\n", decode_name (S->name[i], CODE), S->name[a],S->seq[i]);
		    else fprintf ( fp, ">%s %s\n%s\n", S->name[a], S->name[a],S->seq[i]);
		  }
	      }
	    vfclose (fp);
	  }
	return 1;
	}
Structure * seq2struc ( Sequence *S, Structure *ST)
        {
	int a, b;

	for ( a=0; a< S->nseq; a++)
	    for ( b=0; b< S->len[a]; b++)
		ST->struc[a][b+1][ST->n_fields-1]=S->seq[a][b];
	return ST;
	}

void aln2struc (Alignment *A, Structure *ST)
        {
	int a, b, c;

	for ( a=0; a< A->nseq; a++)
	    for (c=0, b=0; b< A->len_aln; b++)
	        {
		if ( !is_gap (A->seq_al[a][b]))
		     {
		     ST->struc[a][c][ST->n_fields-1]=A->seq_al[a][b];
		     c++;
		     }
		}
	}
Alignment *stack_aln (Alignment *A, Alignment *B)
        {
	int a,b;
	int max_len=0, max_nseq=0;
	if ( B==NULL)return A;
	if ( A==NULL)return B;

	max_nseq=A->nseq+B->nseq;
	for (a=0; a< A->nseq; a++)max_len=MAX(strlen(A->seq_al[a]),max_len);
	for (a=0; a< B->nseq; a++)max_len=MAX(strlen(B->seq_al[a]),max_len);

	A=realloc_aln2 ( A,max_nseq,max_len+1);

	for (a=A->nseq,b=0; b< B->nseq; b++, a++)
	    {
	    sprintf ( A->seq_comment[a] , "%s", B->seq_comment[b]);
	    sprintf ( A->aln_comment[a] , "%s", B->aln_comment[b]);

	    sprintf ( A->seq_al [a] , "%s", B->seq_al [b]);
	    sprintf ( A->name   [a] , "%s", B->name[b]);
	    sprintf ( A->file   [a], "%s" , B->file[b]);
	    A->order[a][0]=B->order[b][0];
	    A->order[a][1]=B->order[b][1];
	    A->score_seq[a]=B->score_seq[b];
	    A->len[a]=B->len[b];
	    }

	A->len_aln=MAX(A->len_aln, B->len_aln);
	A->nseq=A->nseq+B->nseq;
	A->score_aln=A->score_aln+B->score_aln;

	A->finished=A->finished+B->finished;
	return A;
	}

Alignment *chseqIaln(char *name, int seq_n, int start,int len,Sequence *S, int seqIaln, Alignment *A)
        {
	char *seq;

	seq=extract_char ( S->seq[seq_n], start, len);
	A=realloc_aln2 (A, (A==NULL)?(seqIaln+1):MAX(A->nseq,seqIaln+1), ((A==NULL)?(strlen (seq)):MAX(strlen (seq),A->len_aln))+1);


	sprintf ( A->seq_al[seqIaln], "%s",seq);


	A->order[seqIaln][0]=seq_n;
	A->order[seqIaln][1]=start;
	sprintf ( A->name[seqIaln], "%s", name);
	A->nseq=MAX(A->nseq, seqIaln+1);
	A->len_aln=return_maxlen(A->seq_al, A->nseq);
	A->S=S;
	vfree (seq);
	return A;
	}

Alignment * aln_gap2random_aa(Alignment *A)
        {
	 int a, b,l;
	 char alp[200];

	 if (strm ( (A->S)->type, "PROTEIN"))
	   sprintf ( alp, "acefghiklmnpqrstuvwy");
	 else if ( strm ( (A->S)->type, "DNA") ||strm ( (A->S)->type, "RNA") )
	   sprintf ( alp, "agct");
	 l=strlen (alp);


	 for (a=0; a<A->nseq; a++)
	    for ( b=0; b<A->len_aln; b++)
	      if ( is_gap (A->seq_al[a][b]))A->seq_al[a][b]=alp[(int)rand()%(l)];
	  return A;
	}

Alignment * make_random_aln(Alignment *A,int nseq, int len, char *alphabet)
        {
	int a;


	A=realloc_aln2(A, nseq, len+1);

	A->nseq=0;
	A->len_aln=len;
	for ( a=0; a< A->nseq; a++)sprintf ( A->file[a], "random alignment");
	for ( a=0; a< nseq; a++)
	    A=add_random_sequence2aln(A,alphabet);
	return A;
	}
Alignment * add_random_sequence2aln( Alignment *A, char *alphabet)
        {
	int a, n;

	vsrand(0);

	n=strlen(alphabet);
	A=realloc_alignment2 (A, A->nseq+1, A->len_aln+1);

	for ( a=0; a< A->len_aln; a++)A->seq_al[A->nseq][a]=alphabet[rand()%n];
	if (! A->name[A->nseq][0])
	  {
	    for ( a=0; a<10; a++)A->name[A->nseq][a]=alphabet[rand()%n];
	    A->name[A->nseq][a]='\0';
	  }

	A->nseq++;
	return A;
	}

Sequence *get_defined_residues( Alignment *A)
        {
	    char *buf;
	    Sequence *S;
	    int a, b, s, l, r;
	    if ( !A || !A->S) return NULL;

	    S=duplicate_sequence (A->S);
	    for ( a=0; a< S->nseq; a++)
		for ( b=0; b< S->len[a]; b++)S->seq[a][b]=UNDEFINED_RESIDUE;
	    buf=(char*)vcalloc(A->len_aln+1,sizeof (char));
	    for ( a=0; a< A->nseq; a++)
	        {
		    sprintf ( buf, "%s",A->seq_al[a]);
		    ungap(buf);
		    l=strlen (buf);
		    s=A->order[a][0];

		    for ( b=1; b<= l; b++)
			{
			r=A->seq_cache[s][b];

			if ( r>=0)S->seq[s][r-1]=(A->S)->seq[s][r-1];
			}
		}
	    vfree(buf);
	    return S;
	}
Alignment *thread_defined_residues_on_aln ( Alignment *A, Sequence *S1)
	{
	int a, b;
	int gap, r,s, r2;
	for ( a=0; a< A->nseq; a++)
	    {
		s=A->order[a][0];
		r=A->order[a][1];
		for (b=0;b< A->len_aln; b++)
		    {
		    gap=is_gap(A->seq_al[a][b]);

		    if (!gap)
			{
			r+=!gap;
			r2=A->seq_cache[s][r]-1;

			if (r2>=0 && S1->seq[s][r2]==UNDEFINED_RESIDUE)
			    A->seq_al[a][b]=UNDEFINED_RESIDUE;
			}
		    }
	    }
	return A;
	}

int ** trim_aln_borders (char **seq1, char **seq2, int nseq)
{
	int a, b, c,l1,l2;
	char *buf1;
	char *buf2;
	int max;




	max=MAX(get_longest_string (seq1,-1, NULL, NULL),get_longest_string (seq2,-1, NULL, NULL))+1;
	buf1=(char*)vcalloc ( max, sizeof(char));
	buf2=(char*)vcalloc ( max, sizeof(char));

	for ( a=0; a< nseq; a++)
	{
		sprintf ( buf1, "%s", seq1[a]);
		sprintf ( buf2, "%s", seq2[a]);



		ungap (buf1);
		ungap (buf2);

		if (str_overlap ( buf1, buf2,'*')!=0)
		{
			l1=strlen ( seq1[a]);
			l2=strlen ( seq2[a]);
			for ( b=0,c=0; c< l1; c++)
				if ( !is_gap(seq1[a][c]))seq1[a][c]=buf1[b++];
				seq1[a][c]='\0';
			for ( b=0,c=0; c< l2; c++)
				if ( !is_gap(seq2[a][c]))seq2[a][c]=buf2[b++];
				seq2[a][c]='\0';
		}
	}
	vfree (buf1);
	vfree (buf2);
	return NULL;

}



Sequence * merge_seq( Sequence *IN, Sequence *OUT)
{
	int a;


	if ( OUT==NULL)
	  {
	    return duplicate_sequence (IN);
	  }
	else
	  {
	    char *dup;
	    if ( IN && (dup=check_hlist_for_dup( IN->name, IN->nseq)))
	      {
		fprintf ( stderr, "\nERROR: %s is duplicated in file %s[FATAL]\n", dup, IN->file[0]);
		myexit (EXIT_FAILURE);
	      }
	    for ( a=0; a< IN->nseq; a++)
	      if ((OUT=add_sequence ( IN, OUT, a))==NULL)return NULL;
	    return OUT;
	  }
}

Alignment *seq_name2removed_seq_name(Sequence *S, Alignment *NA, float **diff)
{
  int a, b, rb, s;
  float min_diff;
  for (a=0; a< S->nseq; a++)
    {
      if (name_is_in_list( S->name[a], NA->name, NA->nseq, 100)!=-1) continue;
      for ( min_diff=100, s=0, b=0; b< NA->nseq; b++)
	{
	  rb=name_is_in_list ( NA->name[b], S->name, S->nseq, 100);
	  if ( diff[a][rb]<min_diff)
	    {
	      s=b;
	      min_diff=diff[a][rb];

	    }
	}
      strcat ( NA->seq_comment[s], " ");
      strcat ( NA->seq_comment[s], S->name[a]);
    }
  return NA;
}




int seq_name2index (char *name, Sequence *S)
{
  if ( !S) return -1;
  else return name_is_in_list ( name, S->name, S->nseq, MAXNAMES+1);
}
char * seq_name2coor ( char *s, int *start, int *end, char sep)
{
  /*name|start|end */
  char n1[100], n2[100];
  int a=0, b=0, c=0;

  n1[0]=n2[0]='\0';
  start[0]=end[0]=0;

  while ( s[a]!=sep && s[a]!='\0')a++;
  if ( s[a]=='\0')return s;
  else
    s[a++]='\0';



  while ( s[a]!=sep && s[a]!='\0')n1[b++]=s[a++];

  if ( s[a]=='\0'){n1[b]='\0';if ( n1[0])start[0]=atoi(n1);return s;}
  else s[a++]=n1[b]='\0';


  while ( s[a]!=sep && s[a]!='\0')n2[c++]=s[a++];
  n2[c]='\0';


  if ( n1[0])start[0]=atoi(n1);
  if ( n2[0])end[0]=atoi(n2);


  return s;
}

Sequence *extract_one_seq(char *n,int start, int end, Alignment *S, int keep_name)
       {

	 int seq, a;
	 FILE*fp;
	 char *name;
	 Sequence *OUT_S;


	 if ( n[0]=='#')seq=S->nseq;
	 else if ( (seq=name_is_in_list (n, S->name, S->nseq, 100)+1)!=0);
	 else if (is_number (n) && (seq=atoi(n))!=0) seq=atoi(n);
	 else
	   {
	     fprintf ( stderr, "\nCould not find Sequence %s [FATAL]", n);
	     myexit (EXIT_FAILURE);
	   }
	 seq--;

	 name=vtmpnam ( NULL);
	 fp=vfopen ( name, "w");
	 if ( start && end &&!keep_name)fprintf (fp, ">%s_%d_%d\n",S->name[seq],start, end);
	 else if ( start && end==0 && !keep_name)fprintf (fp, ">%s_%d_%d\n",S->name[seq],start,(int)strlen ( S->seq_al[seq]));
	 else fprintf (fp, ">%s\n", S->name[seq]);

	 if ( start==0 && end==0){fprintf (fp, "%s\n", S->seq_al[seq]);}
	 else if (end==0){fprintf (fp, "%s\n", S->seq_al[seq]+start-1);}
	 else
	   {
	     for ( a=start-1; a<end; a++){fprintf ( fp, "%c", S->seq_al[seq][a]);}
	     fprintf ( fp, "\n");
	   }


	 vfclose (fp);
	 OUT_S=get_fasta_sequence_num (name, NULL);

	 return OUT_S;
       }



Sequence * extract_sub_seq( Sequence  *COOR, Sequence *S)
{
	int a, b, c,s;
	int start, end;

	for ( a=0; a< S->nseq; a++)
	{
		if ( (s=name_is_in_list ( S->name[a], COOR->name, COOR->nseq, 100))!=-1)
		{

			sscanf ( COOR->seq_comment[s], "%d %d", &start, &end);
			for (c=0,b=start-1; b< end; b++, c++)S->seq[a][c]=S->seq[a][b];
			S->seq[a][c]='\0';
			sprintf ( S->seq_comment[a], "%s",COOR->seq_comment[s]);

		}
	}
	S=reorder_seq ( S, COOR->name, COOR->nseq);
	return S;
}



char *    aln_column2string (Alignment *A, int p)
  {
    char *s;
    int a;
    if (p>=A->len_aln)
      {
	HERE ("ERROR: index (p=%d) loger than aln (l=%d) [FATAL]", p, A->len_aln);
	myexit (EXIT_FAILURE);
      }
    else
      {
	s=(char*)vcalloc (A->nseq+1, sizeof (char));
	for (a=0; a< A->nseq; a++)s[a]=A->seq_al[a][p];
      }
    return s;
  }


int **fix_seq_aln (Sequence *S, Alignment*A, int **cache)
{
  int s, b,i,nr;

  if (!cache)cache=(int**)vcalloc (S->nseq, sizeof (int*));

  for (s=0; s<A->nseq; s++)
    {
      if ((i=name_is_in_list (A->name[s], S->name, S->nseq, 100)==-1))continue;
      for (nr=0,b=0; b<A->len_aln; b++)
	{
	  if (!is_gap(A->seq_al[s][b]))
	    cache[i][++nr]=b+1;
	}
    }
  return cache;
}

int **fix_seq_seq (Sequence *S0, Sequence *Sx)
{
  //Expresses seq1 in terms of s2
  //sequences 0-N
  //residues  1-N+1
  int s0, r0,i;
  int **index;

  index=(int**)vcalloc ( S0->nseq, sizeof (int*));
  for (s0=0; s0<S0->nseq; s0++)
    {
      int l=S0->len[s0];
      index[s0]=(int*)vcalloc (l+1, sizeof (int));
      i=index[s0][0]=name_is_in_list (S0->name[s0], Sx->name, Sx->nseq, 100);

      if (i==-1);
      else if (strim(S0->seq[s0], Sx->seq[i]))
	{
	  for (r0=1; r0<=l; r0++)
	    {
	      index [s0][r0]=r0;
	    }
	}
      else
	{
	  int c;
	  int nr0=0;
	  int nr1=0;
	  Alignment *B;
	  
	  B=align_two_sequences (S0->seq[s0], Sx->seq[i], const_cast<char*>( (strm(S0->type, "PROTEIN"))?"blosum62mt":"idmat"), -4,-1, "myers_miller_pair_wise");
	  for (c=0; c<B->len_aln; c++)
	    {

	      int g0=is_gap(B->seq_al[0][c]);
	      int g1=is_gap(B->seq_al[1][c]);
	      nr0+=1-g0;
	      nr1+=1-g1;
	      if (!g0 && !g1)index[s0][nr0]=nr1;
	    }
	  if (aln2sim(B, "idmat")<20) add_warning (stderr,"Unreliable reconciliation for sequence %s. If it is a PDB, check source file", S0->name[s0]);
	  
	}
    }
  return index;
}
int **fix_aln_seq_new (Alignment *A, Sequence *Sx)
{
  Sequence *S;
  int **f;

  S=aln2seq (A);
  f=fix_seq_seq(S, Sx);
  free_sequence (S, S->nseq);
  return f;
}
Alignment * fix_aln_seq  ( Alignment *A, Sequence *S)
        {
	int a, b, c;
	char *buf1, *buf2;
	int g0, g1, nr0, nr1;
	int id, tot;
	Alignment *B;


	/*This function establishes the correspondance between every (1..N+1) residue of each aligned sequence
	  and its correspondance in S:
	  A->seq_cache[a][b]=x means that residue b of aligned sequence a corresponds to residue x of the sequence with tye same index in S
	  A->seq_cache[a][b]=0 means there is no correspondance.
	  a is the index of the sequence
	  Applying this function is needed for turning an alignment into a constraint list
	*/


	if ( S==NULL)return A;
	A=reorder_aln (A, S->name,S->nseq);
	if (A->seq_cache)free_int (A->seq_cache, -1);
	A->seq_cache=declare_int ( S->nseq, MAX((A->len_aln+1), S->max_len+1));

	for (a=0; a< S->nseq; a++)
	  for ( b=0; b< A->len_aln; b++)A->seq_cache[a][b]=-1;

	buf1=buf2=NULL;
	for ( a=0; a< S->nseq; a++)
	    {
	    for (b=0; b< A->nseq; b++)
	        {
		if (strm ( S->name[a], A->name[b]))
		   {
		    
		   A->order[b][0]=a;

		   vfree (buf1);
		   
		   buf1=(char*)vcalloc ( A->len_aln+1, sizeof (char));
		   sprintf (buf1, "%s", A->seq_al[b]);
		   ungap (buf1);
		   upper_string (buf1);
		   
		   vfree(buf2);
		   buf2=(char*)vcalloc (strlen(S->seq[a])+1, sizeof (char));
		   sprintf (buf2, "%s",S->seq[a]);
		   ungap (buf2);
		   upper_string (buf2);



		   if ( strm (buf1,buf2))
		       {

			   for ( c=0; c<S->len[a]; c++)A->seq_cache[a][c+1]=c+1;
		       }
		   else
		       {

			   B=align_two_sequences (buf2,buf1,"blosum62mt",-4,-1, "myers_miller_pair_wise");
			   if ( getenv ("DEBUG_RECONCILIATION"))
			     {
			       fprintf (stderr, "\n[DEBUG_RECONCILIATION:fix_aln_seq]\nReconciliation of %s\nA=Ref_sequence\nB=New_seq", S->name[a]);
			       print_aln (B);
			     }

			   for (id=0, tot=0,nr0=0,nr1=0,c=0; c<B->len_aln; c++)
			     {
			       g0=is_gap(B->seq_al[0][c]);
			       g1=is_gap(B->seq_al[1][c]);
			       nr0+=1-g0;
			       nr1+=1-g1;
			       if ( !g0 && !g1)
				 {
				   tot++;
				   id+=(B->seq_al[0][c]==B->seq_al[1][c])?1:0;
				   A->seq_cache[a][nr1]=nr0;
				 }
			       else if (g0 && !g1)
				 {
				   A->seq_cache[a][nr1]=0;
				 }
			     }
			   if ( ((id*100)/tot)<20)
			     {
			       print_aln (B);
			       fprintf ( stderr, "\nTwo different sequences have the same name: %s", S->name[a]);
			       fprintf ( stderr, "\nIf %s is a PDBID, Make sure it identifies the right chain (A, B, 1, 2...)", S->name[a]);
			       fprintf ( stderr, "\nChain number or index must be added to the PDB id (i.e. 1gowA)");
			       fprintf ( stderr, "\nIf You want to use %s anyway, rename it with a non-PDB identifier such as seq_%s\n",S->name[a],S->name[a]);
			       myexit (EXIT_FAILURE);
			     }

			   free_sequence ( B->S, -1);
			   free_aln (B);
		       }

		   }
		}
	    }
	    vfree(buf1);vfree(buf2);
	    return A;
	}

Sequence * add_prf2seq  ( char *file, Sequence *S)
    {
      static int n;
      static char** prf_name;
      char **new_seq;
      Sequence *NS;

      if (!prf_name)
	{
	  prf_name=declare_char (1,100);
	}
      if (file && !atoigetenv("KM_COFFEE_PRF"))
	sprintf (prf_name[0], "%s", file);
      else
	sprintf (prf_name[0], "prf_%d", ++n);
      
      if ( !is_aln (file)&& !is_seq (file))return S;
      else
	{
	  X_template *R;
	  Alignment *A;


	  R=fill_R_template(file,file, S);

	  A=(R->VR)->A;
	  ((R->VR)->A)->expand=1;
	  new_seq=declare_char (1,A->len_aln+1);
	  
	  


	  if (atoigetenv("KM_COFFEE_CONS_COV"))
	    sprintf ( new_seq[0], "%s",aln2cons_seq_cov(A));
	  else
	    sprintf ( new_seq[0], "%s",aln2cons_seq_mat(A, "blosum62mt"));
	  
	  
	  NS=fill_sequence_struc(1, new_seq,prf_name, NULL);

	  S=add_sequence (NS, S, 0);
	  (S->T[S->nseq-1])->R=R;

	  free_sequence (NS, NS->nseq);
	  free_char( new_seq, -1);

	  return S;
	}
    }
int prf_in_seq ( Sequence *S)
{
  int a;

  if ( !S) return 0;
  else
    {
      for ( a=0; a< S->nseq; a++)
	if (seq2R_template_profile(S, a)) return 1;
    }
  return 0;
}



Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i)
{
	int s, a;


	char *buf;
	if (OUT==NULL)
	  {
	    OUT=duplicate_sequence (IN);
	    return OUT;
	  }
	for (a=0; a<OUT->nseq; a++)
	  {
	    Alignment *P;
	    P=seq2R_template_profile (OUT, a);
	    if (!P)
	      continue;
	    else if (name_is_in_list (IN->name[i], P->name, P->nseq, 100)!=-1)
	      return OUT;
	  }

	/*Adds sequence i of IN at the end of OUT*/

	if ((s=name_is_in_list ( IN->name[i], OUT->name, OUT->nseq,STRING))==-1 )
	{
		OUT=realloc_sequence (OUT, OUT->nseq+1, IN->len[i]);
		sprintf ( OUT->name[OUT->nseq],"%s",IN->name[i]);
		sprintf ( OUT->file[OUT->nseq],"%s",IN->file[i]);
		sprintf ( OUT->seq_comment[OUT->nseq],"%s",IN->seq_comment[i]);
		sprintf ( OUT->aln_comment[OUT->nseq],"%s",IN->aln_comment[i]);

		sprintf ( OUT->seq[OUT->nseq],"%s",IN->seq[i]);
		if (IN -> genome_co != NULL)
		{
			Genomic_info *tmp_in = &(IN->genome_co[i]);
			Genomic_info *tmp_out = &(OUT->genome_co[OUT->nseq]);
			tmp_out->strand = tmp_in->strand;
			tmp_out->start = tmp_in->start;
			tmp_out->end = tmp_in->end;
			tmp_out->seg_len = tmp_in->seg_len;
			tmp_out->seg_name =(char*) vcalloc(strlen(tmp_in->seg_name)+1, sizeof(char));
			strcpy(tmp_out->seg_name, tmp_in->seg_name);

		}

		OUT->len[OUT->nseq]=IN->len[i];
		OUT->T[OUT->nseq][0]=IN->T[i][0];
		OUT->nseq++;
		return OUT;
	}
	else if ( s!=-1 && !case_insensitive_strcmp ( IN->seq[i], OUT->seq[s]))
	{
	  fprintf ( stderr,"[DEBUG_RECONCILIATION:add_sequence]\n%s\n%s\n", IN->seq[i], OUT->seq[s]);
	  if ( getenv4debug("DEBUG_RECONCILIATION"))fprintf ( stderr,"[DEBUG_RECONCILIATION:add_sequence]\n%s\n%s\n", IN->seq[i], OUT->seq[s]);

	  add_warning (stderr, "DISCREPANCY:%s in [%s] and  [%s]\n", IN->name[i], IN->file[i], OUT->file[s]);


	  if (((buf=build_consensus(IN->seq[i], OUT->seq[s],"cfasta_pair_wise" ))!=NULL) || ((buf=build_consensus(IN->seq[i], OUT->seq[s],"myers_miller_pair_wise" ))!=NULL))
	    {

	      OUT->max_len=MAX(OUT->max_len, strlen(buf));
	      OUT->min_len=MIN(OUT->min_len, strlen(buf));
	      OUT->seq    =realloc_char ( OUT->seq, -1, -1,OUT->nseq,OUT->max_len+1);

	      sprintf ( OUT->seq[s],"%s",buf);
	      OUT->len[s]=strlen (buf);
	      vfree (buf);
	      return OUT;
	    }
	  else
	    {
	      fprintf ( stderr, "IMPOSSIBLE TO RECONCILIATE SOME SEQUENCES[FATAL:%s]\n", PROGRAM);
	      print_aln ( align_two_sequences (IN->seq[i], OUT->seq[s], "idmat", 0, 0, "fasta_pair_wise"));
	      myexit (EXIT_FAILURE);
	      return NULL;
	    }

	}
	else
	  {
	    return OUT;
	  }
}


Sequence  * trim_seq       ( Sequence  *A, Sequence  *B)
        {
	int a;
	Sequence *R;

	if (A->nseq>B->nseq)
	  {
	    Sequence *I;
	    I=A;A=B;B=I;
	  }

	R=declare_sequence (MIN(A->min_len,B->min_len), MAX(A->max_len, B->max_len), MIN(A->nseq, B->nseq));
	R->nseq=0;

	for (a=0; a< A->nseq; a++)
	    {
	    if ( name_is_in_list ( A->name[a], B->name, B->nseq,STRING+1)!=-1)
	        {
		  sprintf ( R->name[R->nseq], "%s", A->name[a]);
		  sprintf ( R->seq[R->nseq], "%s", A->seq[a]);
		  sprintf ( R->file[R->nseq], "%s", A->file[a]);
		  sprintf ( R->aln_comment[R->nseq], "%s", A->aln_comment[a]);
		  sprintf ( R->seq_comment[R->nseq], "%s", A->seq_comment[a]);

		  R->len[R->nseq]=A->len[a];
		  R->nseq++;
		}
	    }
	return R;
	}

Sequence * trim_aln_seq ( Alignment *A, Alignment *B)
        {
	int a;
	static char **name_list;
	int n=0;
	Sequence *SA, *SB;
	int **cache_A=NULL;
	int **cache_B=NULL;
	int * p;

	/*This function inputs two alignments A and B
	  It removes sequences that are not common to both of them
	  It rearange the sequences so that they are in the same order
	  A decides on the order
	  The Sequences (A->S) and (B->S) are treated the same way
	  Sequences are also merged in order to detects discrepencies.
	  A pointer to S is returned
	*/
	if (name_list)free_char (name_list, -1);
	name_list=declare_char (MAX(A->nseq, B->nseq), STRING+1);

	for ( a=0; a< A->nseq; a++)
	    {
	    if ( name_is_in_list ( A->name[a], B->name, B->nseq,STRING)!=-1)
	        {
		sprintf ( name_list[n++], "%s", A->name[a]);
		}
	    }



	reorder_aln ( A, name_list, n);
	if (A->seq_cache)cache_A=duplicate_int (A->seq_cache, -1, -1);
	if (B->seq_cache)cache_B=duplicate_int (B->seq_cache, -1, -1);
	reorder_aln ( B, name_list, n);
	for ( a=0; a< n; a++)
	  {
	    if ( cache_A)
	      {
		p=A->seq_cache[A->order[a][0]];
		A->seq_cache[A->order[a][0]]=cache_A[a];
		cache_A[a]=p;
	      }
	   if ( cache_B)
	      {
		p=B->seq_cache[B->order[a][0]];
		B->seq_cache[B->order[a][0]]=cache_B[a];
		cache_B[a]=p;
	      }
	   A->order[a][0]=B->order[a][0]=a;
	  }
        free_int(A->seq_cache, -1);
	free_int(B->seq_cache, -1);

	A->seq_cache=cache_A;
	B->seq_cache=cache_B;



	SA=aln2seq(A);
	SB=aln2seq(B);

	A->S=B->S=merge_seq (SA, SB);
	return A->S;
	}



Sequence * trim_aln_seq_name ( Alignment *A, Alignment *B)
        {
	int a;
	Sequence *S;

	/*This function inputs two alignments A and B
	  It removes sequences that are not common to both of them
	  It rearange the sequences so that they are in the same order
	  A decides on the order
	*/
	S=declare_sequence ( 1, 1, A->nseq+B->nseq);
	S->nseq=0;
	for ( a=0; a< A->nseq; a++)
	    {
	    if ( name_is_in_list ( A->name[a], B->name, B->nseq,STRING)!=-1)
	        {
		sprintf ( S->name[S->nseq++], "%s", A->name[a]);
		}
	    }
	return S;
	}



char ** rm_name_tag (char **name, int nseq, char *tag)
{
  int a , b, ntag;
  char **tag_list;
  char *s;
  char **template_list;
  if ( !name )return NULL;

  tag_list=declare_char (10, 4);

  if ( tag)
    {
      ntag=1; sprintf ( tag_list[0], "%s", tag);
    }
  else
    {
      ntag=0;
      sprintf ( tag_list[ntag++], "_S_");
      sprintf ( tag_list[ntag++], "_G_");
    }
  template_list=declare_char (nseq, 100);
  for ( a=0; a<nseq ; a++)
    {
      for ( b=0; b<ntag; b++)
	{
	s=strstr(name[a], tag_list[b]);
	if ( s)
	  {
	    s[0]='\0';
            s[2]='\0';
	    sprintf ( template_list[a], ">%s _%s_ %s", name[a], s+1, s+3);
	    break;
	  }
	}
    }

  free_char (tag_list, -1);
  return template_list;
}
Sequence * swap_header ( Sequence *S, Sequence *H)
{
  int a, b, n;

  for ( a=0; a< S->nseq; a++)
    {
      if ( (n=name_is_in_list (S->name[a],H->name, H->nseq, 1000))!=-1)
	   {
	     char **list;


	     list=string2list (H->seq_comment[n]);
	     if ( list==NULL || atoi(list[0])==1)continue;
	     sprintf (S->name[a], "%s%s%s",H->name[n], list[1], list[2]);
	     vfree ( S->seq_comment[a]);S->seq_comment[a]=(char*)vcalloc ( strlen (H->seq_comment[n])+1, sizeof (char));
	     for (b=3; b< atoi(list[0]); b++)S->seq_comment[a]=strcat (S->seq_comment[a], list[b]);
	     free_char (list, -1);
	   }
    }
  return S;
}

char *template_file2abs_template_file(char *name)
{
  FILE *out;
  FILE *in;
  char *outF;
  int  a;
  char ***l;
  char *pdb=NULL;

  
  
  if (!name) return NULL;
  if (!file_exists (NULL, name)) return name;
 
  
  l=file2list(name, " ");
  a=0;
  out=vfopen (outF=vtmpnam (NULL), "w");
  while (l[a])
    {
      //if pdb does not exist, check if its *.pdb version exists and seek the abs path
      //if no file to be found, assume it is an identifier and leave it untouched
      
      if (check_file_exists(l[a][3]))pdb=csprintf (pdb, "%s", l[a][3]);
      else pdb=csprintf (pdb, "%s.pdb", l[a][3]);
      
      //if (check_file_exists (pdb))pdb=csprintf (pdb, "%s",fname2abs(pdb));
      if (check_file_exists (pdb)){;}
      else pdb=csprintf (pdb, "%s", l[a][3]);
      
      fprintf (out,"%s %s %s\n", l[a][1], l[a][2],pdb);
      a++;
    }

  vfclose (out);
  free_arrayN ((void ***)l, 3);
  vfree(pdb);
  return outF;
}
	
	       

Sequence * profile_seq2template_seq ( Sequence *S, char *template_file, Fname *F)
{
  /*This function fetches potential templates associated with sequences within a profile*/
  int i,b;
  Alignment *A;
  
  
  for ( i=0; i< S->nseq; i++)
    {
      if ( (A=seq2R_template_profile (S, i)))
	{
	  A->S=aln2seq (A);
	  A->S=seq2template_seq (A->S, template_file, F);
	}
    }
  return S;
}

/**
 * Specify types of used templates for each sequence.
 *
 * Writes the types of templates existing for each Seqeuence into Template::seq_type of
 * the Sequence::T template object. The seq_type is a string like "P..S.......", for example,
 * where dots are actually white spaces.
 *
 */
Sequence * seq2template_type(Sequence *Seq)
{
  //add template
  int a, e;
  int s;
  struct X_template *S=NULL;
  struct X_template *P=NULL;
  struct X_template *R=NULL;
  struct X_template *G=NULL;
  struct X_template *F=NULL;
  struct X_template *T=NULL;
  struct X_template *E=NULL;
  struct X_template *U=NULL;
  Alignment *A;


  e=' ';
  for (a=0; a< Seq->nseq; a++)
    {
      if (!Seq->T[a])continue;
      //HERE ADD a Template
      P=seq_has_template (Seq, a, "_P_");
      S=seq_has_template (Seq, a, "_S_");
      R=seq_has_template (Seq, a, "_R_");
      G=seq_has_template (Seq, a, "_G_");
      F=seq_has_template (Seq, a, "_F_");
      T=seq_has_template (Seq, a, "_T_");
      E=seq_has_template (Seq, a, "_E_");
      U=seq_has_template (Seq, a, "_U_");

      s=(!P)?1:0;
      sprintf ( (Seq->T[a])->seq_type, "%c%c%c%c%c%c%c%c", (P)?'P':e, (S)?'S':e, (S &&!P)?'s':e,(R)?'R':e, (G)?'G':e,(T)?'T':e,(E)?'E':e,(U)?'U':e);

      if (R && (A=seq2R_template_profile (Seq,a)) && A->S)
	{

	  A->S=seq2template_type ( A->S);
	}
    }
  return Seq;
}

char * string_contains_template_tag (char *string_in)
{
  char string[100];

  if ( strstr (string, "_P_"))return "_P_";
  if ( strstr (string, "_S_"))return "_S_";
  if ( strstr (string, "_R_"))return "_R_";
  if ( strstr (string, "_G_"))return "_G_";
  if ( strstr (string, "_F_"))return "_F_";
  if ( strstr (string, "_T_"))return "_T_";
  if ( strstr (string, "_E_"))return "_E_";
  if ( strstr (string, "_U_"))return "_U_";

  return NULL;
}
static int check_blast_is_installed (char *server);



static int check_blast_is_installed (char *server)
{
  if (strm (server, "EBI"));
  else if ( strm (server, "NCBI"))
	// The BLAST client to access the NCBI services it is embedded in the BLAST+ tools by adding the '-remote' command line option
    return check_program_is_installed ("blastp",NULL, NULL,NCBIBLAST_ADDRESS, INSTALL_OR_DIE);
  else if ( strm (server, "LOCAL"))
    return check_program_is_installed (NCBIBLAST_4_TCOFFEE,NULL, NULL,NCBIBLAST_ADDRESS, INSTALL_OR_DIE);
  return 1;
}


Sequence * vremove_seq_template_files(Sequence *S)
{
    return handle_seq_template_file (S, "remove");
}
Sequence * display_seq_template_files(Sequence *S)
{
  return handle_seq_template_file (S, "display");
}
Sequence * handle_seq_template_file (Sequence *S, char *mode)
{
  int a;
  Template *T;

  for (a=0; a< S->nseq; a++)
    {
      T=S->T[a];
      if (T)
	{
	  handle_X_template_files (T->P, mode);
	  handle_X_template_files (T->F, mode);
	  handle_X_template_files (T->R, mode);
	  handle_X_template_files (T->T, mode);
	  handle_X_template_files (T->E, mode);
	}
    }

  return S;
}
int handle_X_template_files ( X_template *T, char *mode)
  {
    if (!T)return 0;

    if ( strm (mode, "remove"))
      {
	vremove (T->template_file);
	vremove (T->template_name);
      }
    else if (strm (mode, "display"))
      {
	static char *buf;
	//do not diplay the nameof template files that are in the cache
	if ( !strstr (T->template_name, get_cache_dir()))
	  {
	    buf=csprintf (buf,"Template %s",  template_type2type_name (T->template_type));
	    if (check_file_exists (T->template_name))display_output_filename ( stdout,buf,T->template_format,T->template_name, STORE);
	  }
      }
    else
      {
	printf_exit (EXIT_FAILURE, stderr, "\nUnkonwn mode %s for template handling [FATAL:%s]", mode, PROGRAM);
      }
    return 1;
  }

char *trim_template_file (char *file, Sequence *S);//Remove from template file all sequences that cannot be used


/**
 * Adds templates to a Sequence object.
 *
 * This function is recursive: Depending on the type of \c template_file,
 * it performs some preparations and calles itselt again with a different template_file.
 *     - \c no_template and the sequence will be returned, unchanged.
 *     - \c MODE_ and the prefix \c MODE_ will be cut off before calling the same function again.
 *     - \c PSIBLAST,\c BLAST, \c EXPRESSO, \c RCOFFEE and secveral others:
 *          The script \b tc_generic_method.pl is initiated (that means, the name of the script is given as new
 *          parameter for the recursive call with an additional prefix \c SCRIPT_ ). @ work here!!!
 *          In this case
 *
 * \todo @ work here!!! Understand this function and then Document it.
 *
 * \param[in,out] S Sequence object, will be modified.
 * \param[in] template_list String containing the template file or commands how to get one
 */
static int ntemp;
Sequence * seq2mmseqs_template_seq(Sequence *S,  Fname *F)
{
  /*Expected format for the template file:
    >seq_name _X_ Target_template
    X: S for Structures
    G for genomes (Exoset)
    When alternative templates are given for a sequence, the first one superseeds all the others
  */
  
  
  /*Fill the sequences*/
  /*1: No template*/
  char buf[1000];
 
  int PmC,PmI,PMI;
  int BmC,BmI,BMI, Trim;
  char *server;
  char *pdb_db,*prot_db;
  char pdb_type[100];
  char *p;
  int remove_template_file=0;
  static char *seqdb;
 
  static char *seq=vtmpnam (NULL);
  static char *outfile=vtmpnam(NULL);
  static char *tf=NULL;
  static char *command;
  char *cache=get_cache_dir();
  ntemp++;

  remove_template_file=get_int_variable ("remove_template_file");
  server=get_string_variable ("blast_server");
  pdb_db=get_string_variable ("pdb_db");
  prot_db=get_string_variable ("prot_db");         

  PmI=get_int_variable ("pdb_min_sim");
  PMI=get_int_variable ("pdb_max_sim");
  PmC=get_int_variable ("pdb_min_cov");

  BmI=get_int_variable ("prot_min_sim");
  BMI=get_int_variable ("prot_max_sim");
  BmC=get_int_variable ("prot_min_cov");
  Trim=get_int_variable("psitrim");
  
  output_fasta_seqS(seq,S);
  if (!F)F=parse_fname (S->file[0]);


  
  tf=csprintf (tf, "%s%s_R_%d.template_list", F->path,F->name,ntemp);
  fprintf ( stderr, "\n! Running MMSEQS against %s -- This may take a while...\n", prot_db);
command=csprintf ( command, "t_coffee -other_pg mmseqs2prf.pl -q %s -db %s -o %s  -template_file %s  -cachedb %s -prot_min_sin %d -prot_max_sim %d -prot_min_cov %d -psitrim %d -quiet", seq, prot_db,cache, tf,cache, BmI, BMI, BmC, Trim);
  printf_system (command);
  if ( check_file_exists (tf) && format_is_fasta(tf))
	{
	  S=seq2template_seq (S,tf, F);
	  trim_template_file (tf,S);
	}
  else
    {
      
      add_warning (stderr, "Could not Run %s to find templates[%s](unforked mode)\n",command, PROGRAM);
      return NULL;
    }
  
  vfree (command);
  return S;
}
Sequence * seq2template_seq ( Sequence *S, char *template_list, Fname *F)
{
  /*Expected format for the template file:
    >seq_name _X_ Target_template
    X: S for Structures
    G for genomes (Exoset)
    When alternative templates are given for a sequence, the first one superseeds all the others
  */
  
  
  /*Fill the sequences*/
  /*1: No template*/
  char buf[1000];
 
  int PmC,PmI,PMI;
  int BmC,BmI,BMI, Trim;
  char *server;
  char *pdb_db,*prot_db;
  char pdb_type[100];
  char *p;
  int remove_template_file=0;
  static char *seqdb;

  remove_template_file=get_int_variable ("remove_template_file");
  server=get_string_variable ("blast_server");
  pdb_db=get_string_variable ("pdb_db");
  prot_db=get_string_variable ("prot_db");         

  PmI=get_int_variable ("pdb_min_sim");
  PMI=get_int_variable ("pdb_max_sim");
  PmC=get_int_variable ("pdb_min_cov");

  BmI=get_int_variable ("prot_min_sim");
  BMI=get_int_variable ("prot_max_sim");
  BmC=get_int_variable ("prot_min_cov");
  Trim=get_int_variable("psitrim");
  
  if (template_list && strm(template_list, "MMSEQS"))return seq2mmseqs_template_seq(S,F);
      

  if (strm (prot_db, "dataset") || strm (prot_db, "self"))
    {
      if (!seqdb)
	{
	  seqdb=vtmpnam(NULL);
	  seq2blastdb (seqdb,S->blastdbS);

	}
      prot_db=seqdb;
      strcpy(server,"LOCAL");
    }

  
  //Set the type of the PDB structure
  if ((p=get_string_variable ("pdb_type")))
    {
      sprintf ( pdb_type, "%s",p);
    }
  else
    {
      sprintf (pdb_type, "dmn");
    }

  if ( (template_list && template_list[0]=='\0') || strm ( template_list, "no_template"))
    {
      return S;
    }
  else if ( strstr (template_list, "MODE_"))//pre_set mode
    {
      return seq2template_seq ( S,template_list+strlen ("MODE_"),F);
    }
  else if ( strm ( template_list, "SSP")|| strm ( template_list, "GOR"))
    {

      /*use GOR to Predict the secondary structure*/
      check_program_is_installed (GOR4_4_TCOFFEE,NULL, NULL,GOR4_ADDRESS, INSTALL_OR_DIE);
      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#ssp_template@seq#%s/%s@obs#%s/%s@cache#%s@type#_E_",get_mcoffee_4_tcoffee(), "New_KS.267.seq", get_mcoffee_4_tcoffee(), "New_KS.267.obs", get_cache_dir());
      S=seq2template_seq (S,buf, F);
      return S;
    }
  else if ( strm ( template_list, "PSISSP") || strm (template_list, "PSIGOR"))
    {

      /*Computes a GOR consensus on a psi-blast output*/
      check_program_is_installed (GOR4_4_TCOFFEE,NULL, NULL,GOR4_ADDRESS, INSTALL_OR_DIE);
      check_blast_is_installed(server);

      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#psissp_template@seq#%s/%s@obs#%s/%s@cache#%s@minid#%d@maxid#%d@mincov#%d@server#%s@type#_E_",get_mcoffee_4_tcoffee(), "New_KS.267.seq", get_mcoffee_4_tcoffee(), "New_KS.267.obs", get_cache_dir(), BmI,BMI,BmC,server);
      S=seq2template_seq (S,buf, F);
      return S;
    }
  else if ( strm ( template_list, "TM"))
    {
      
      /*predict transmembrane structure*/
      check_program_is_installed (HMMTOP_4_TCOFFEE,NULL, NULL,HMMTOP_ADDRESS, INSTALL_OR_DIE);
      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#tm_template@arch#%s/%s@psv#%s/%s@type#_T_",get_mcoffee_4_tcoffee(), "hmmtop.arch", get_mcoffee_4_tcoffee(), "hmmtop.psv");
      
      S=seq2template_seq (S,buf, F);
      return S;
    }
  else if ( strm ( template_list, "PSITM"))
    {
      cputenv ("psiJ_4_TCOFFEE=5");
      /*predict transmembrane structure*/
      check_program_is_installed (HMMTOP_4_TCOFFEE,NULL, NULL,HMMTOP_ADDRESS, INSTALL_OR_DIE);
      check_blast_is_installed(server);

      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#psitm_template@database#%s@arch#%s/%s@psv#%s/%s@cache#%s@minid#%d@maxid#%d@mincov#%d@server#%s@type#_T_", prot_db, get_mcoffee_4_tcoffee(), "hmmtop.arch", get_mcoffee_4_tcoffee(), "hmmtop.psv",get_cache_dir(), BmI,BMI,BmC,server);
      S=seq2template_seq (S,buf, F);
      return S;
    }
  else if (strm ( template_list, "HH"))
    {
      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#psiprofile_template@database#%s@method#hh@cache#%s@minid#%d@maxid#%d@mincov#%d@trim#%d@server#%s@type#_R_", prot_db,get_cache_dir(),BmI,BMI,BmC,Trim,server);
      S=seq2template_seq (S,buf, F);

      return S;
    }
  else if (strm ( template_list, "PSIBLAST"))
    {
      //Note PSIBLAST ans BLAST are NOW THE SAME, but BLAST=psiblast -num_iterations 1
      check_blast_is_installed(server);
    
      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#psiprofile_template@database#%s@method#blastp@cache#%s@minid#%d@maxid#%d@mincov#%d@trim#%d@server#%s@type#_R_", prot_db,get_cache_dir(),BmI,BMI,BmC,Trim,server);
      S=seq2template_seq (S,buf, F);

      return S;
    }
  else if (strm ( template_list, "BLAST") )
    {

      check_blast_is_installed(server);
      
      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#profile_template@database#%s@method#blastp@cache#%s@minid#%d@maxid#%d@mincov#%d@trim#%d@server#%s@type#_R_", prot_db,get_cache_dir(),BmI,BMI,BmC,Trim,server);
      S=seq2template_seq (S,buf, F);

      return S;
    }
  else if ( strm ( template_list, "EXPRESSO") || strm (template_list, "PDB"))
    {
      check_blast_is_installed(server);
      
      int isRNA = 0;
      int i;
      for (i= 0; i < S->len[0]; ++i)
	{
	   isRNA =  (isRNA || is_rna(S->seq[0][i]));
	}

      if (isRNA)
	{
	  sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#pdb_template@database#%s@method#blastn@cache#%s@minid#%d@maxid#%d@mincov#%d@server#%s@type#_P_@pdb_type#%s",pdb_db, get_cache_dir(),PmI,PMI,PmC, server,pdb_type);
	}
      else
	{

	  sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#pdb_template@database#%s@method#blastp@cache#%s@minid#%d@maxid#%d@mincov#%d@server#%s@type#_P_@pdb_type#%s",pdb_db, get_cache_dir(),PmI,PMI,PmC, server,pdb_type);
	}
      
       return seq2template_seq (S,buf, F);
    }

  else if ( strm (template_list, "RCOFFEE") || strm (template_list, "RNA"))
    {
      
      //extract structure from sequences if possible otherwise use RNAPlfold
      char *file_struc_calc = vtmpnam (NULL);
      FILE* struc_calc_f =vfopen(file_struc_calc,"w");
      int i;
      
      for (i = 0; i< S->nseq; ++i)
	{
	  if (S->T[i]->P)
	    {
	      fprintf(struc_calc_f,"%s %s\n",S->name[i],S->T[i]->P->template_file);
	    }
	  else
	    {
	      fprintf(struc_calc_f,"%s\n",S->name[i]);
	    }
	}
      vfclose(struc_calc_f);
      sprintf ( buf, "SCRIPT_tc_generic_method.pl@mode#RNA_template@pdbfile#%s@cache#%s@type#_F_", file_struc_calc,get_cache_dir());
      
      return seq2template_seq (S,buf,F);
    }


  /*2: Templates from seqnames (SELF) or named like the sequences (SEQFILE)*/
  else if ( strstr (template_list, "SELF_") ||strstr (template_list, "SEQFILE_") )
    {
      int a;
      char *p;
      
      //add template
      for (a=0; a< S->nseq; a++)
	{

	  if ( (p=strstr (template_list,"SELF_")))p=S->name[a];
	  else if ( strstr (template_list, "SEQFILE_"))p=template_list;
	  else
	    {
	      fprintf ( stderr, "\nUnkown mode for Template [FATAL:%s]\n", PROGRAM);
	      myexit (EXIT_FAILURE);
	    }

	  if (      strstr (template_list, "_P_") && !(S->T[a])->P)(S->T[a])->P  =fill_P_template  ( S->name[a], p,S);//PDB
	  else if ( strstr (template_list, "_S_") && !(S->T[a])->S)(S->T[a])->S  =fill_S_template  ( S->name[a], p,S);//Sequence
	  else if ( strstr (template_list, "_R_" )&& !(S->T[a])->R)(S->T[a])->R  =fill_R_template  ( S->name[a], p,S);//pRofile
	  else if ( strstr (template_list, "_G_" )&& !(S->T[a])->G)(S->T[a])->G  =fill_G_template  ( S->name[a], p,S);//Genomic
	  else if ( strstr (template_list, "_F_" )&& !(S->T[a])->F)(S->T[a])->F  =fill_F_template  ( S->name[a], p,S);//Fold
	  else if ( strstr (template_list, "_T_" )&& !(S->T[a])->T)(S->T[a])->T  =fill_T_template  ( S->name[a], p,S);//Trans Membrane
	  else if ( strstr (template_list, "_E_" )&& !(S->T[a])->E)(S->T[a])->E  =fill_E_template  ( S->name[a], p,S);//Secondary Structure
	  else if ( strstr (template_list, "_U_" )&& !(S->T[a])->U)(S->T[a])->U  =fill_U_template  ( S->name[a], p,S);//unicode, list template

	}
      return S;
    }

  /*2: Templates comes in a template_file*/
  else if ( template_list==NULL || format_is_fasta (template_list))
    {
      Sequence *T;
      int a, i;
      int ntemp=0;
      
      T=(template_list!=NULL)?get_fasta_sequence (template_list, NULL):S;
      for (a=0; a< T->nseq; a++)
	{
	  
	  
	  char *p;
	  if ((i=name_is_in_list(T->name[a], S->name, S->nseq, MAXNAMES))!=-1)
	    {
	      //This thing because the new parsing does not keep the space between name and comments in FASTA
	      if (T->seq_comment[a] && T->seq_comment[a][0]=='_')
		{
		  char *buf=csprintf (NULL, " %s", T->seq_comment[a]);
		  T->seq_comment[a]=csprintf ( T->seq_comment[a], "%s",buf);
		  vfree(buf);
		}
	      if (       (p=strstr (T->seq_comment[a], " _P_ ")) && !(S->T[i])->P &&( (S->T[i])->P=fill_P_template (S->name[i],p,S)))ntemp++;
	      else if (  (p=strstr (T->seq_comment[a], " _F_ ")) && !(S->T[i])->F &&( (S->T[i])->F=fill_F_template (S->name[i],p,S)))ntemp++;
	      else if (  (p=strstr (T->seq_comment[a], " _S_ ")) && !(S->T[i])->S &&( (S->T[i])->S=fill_S_template (S->name[i],p,S)))ntemp++;

	      else if (  (p=strstr (T->seq_comment[a], " _R_ ")) && !(S->T[i])->R &&( (S->T[i])->R=fill_R_template (S->name[i],p,S)))ntemp++;
	      else if (  (p=strstr (T->seq_comment[a], " _G_ ")) && !(S->T[i])->G &&( (S->T[i])->G=fill_G_template (S->name[i],p,S)))ntemp++;
	      else if (  (p=strstr (T->seq_comment[a], " _T_ ")) && !(S->T[i])->T &&( (S->T[i])->T=fill_T_template (S->name[i],p,S)))ntemp++;
	      else if (  (p=strstr (T->seq_comment[a], " _E_ ")) && !(S->T[i])->E &&( (S->T[i])->E=fill_E_template (S->name[i],p,S)))ntemp++;
	      else if (  (p=strstr (T->seq_comment[a], " _U_ ")) && !(S->T[i])->U &&( (S->T[i])->E=fill_U_template (S->name[i],p,S)))ntemp++;

	      if (T!=S)strcat (S->seq_comment[i], T->seq_comment[a]);

	    }
	}

      if (T!=S)free_sequence (T, -1);

      if ( remove_template_file==2)
	{
	  vremove (template_list);
	}
      else
	if (template_list)display_output_filename ( stdout, "Template_List","fasta_seq", template_list, STORE);
      return S;
    }

  /*3 Templates are generated with a script*/
  else if (strstr (template_list, "SCRIPT_") && get_string_variable ("multi_core") && strstr (get_string_variable ("multi_core"), "templates") && get_nproc()>1)
    {
      char *tmp1,*command;
      Alignment *A;
      char **temp_file,**seq_file;
      int  * pid_list, pid, npid, submited;
      int nproc, max_nproc;
      int num=0;

      char outfile[1000];
      static char *script;
      static int ntemp;
      char *p;
      int z, i;
      int freeF=0;
     
      if (!script)script=(char*)vcalloc ( 1000, sizeof(char));

      ntemp++;

      command=(char*)vcalloc ( 1000, sizeof (char));
      tmp1=vtmpnam (NULL);

      A=seq2aln (S,NULL, 0);
      string_array_upper(A->seq_al, A->nseq);
      output_fasta_seq (tmp1, A);
      sprintf ( script, "%s", after_strstr (template_list, "SCRIPT_"));

      if ((p=strstr (template_list, "@type#")))
	p+=strlen ("@type#");

      if (!F){F=parse_fname (S->file[0]);freeF=1;}
      sprintf (outfile, "%s%s_%s%d.template_list", F->path,F->name,template_type2short_type_name(p),ntemp);
      while ( check_file_exists (outfile))
	{
	  sprintf (outfile, "%s%s_%s%d.%d.template_list",F->path, F->name,template_type2short_type_name(p),ntemp, ++num);
	}
      if (freeF)free_fname(F);

      nproc=get_nproc();
      //max_nproc=2*nproc;
      max_nproc=10; //EBI recommended maximum
      script=substitute(script, "@", " -");
      script=substitute(script, "#", "=");

      temp_file=(char**)vcalloc ( A->nseq, sizeof (char*));
      seq_file =(char**)vcalloc (A->nseq, sizeof (char*));
      pid_list =(int *)vcalloc (MAX_N_PID, sizeof (int *));

      fprintf ( stderr, "\n\t------ Fetch Templates [Multi Core Mode %d CPUs]\n",get_nproc());
      for (npid=0, submited=0,i=0; i<S->nseq; i++)
	{
	  FILE *fp2;
	  seq_file[i]=vtmpnam (NULL);
	  temp_file[i]=vtmpnam (NULL);
	  fp2=vfopen (seq_file[i], "w");
	  fprintf ( fp2, ">%s\n%s\n", S->name[i], S->seq[i]);
	  vfclose (fp2);

	  pid=vvfork(NULL);
	  if (pid==0)
	    {
	      initiate_vtmpnam (NULL);
	      if  ( strstr (script, "tc_generic_method"))
		{
		  //sprintf ( command, "%s -other_pg %s -infile=%s -outfile=%s -tmpdir=%s",get_string_variable ("t_coffee"),script,seq_file[i],temp_file[i],get_tmp_4_tcoffee());
		  sprintf ( command, "%s -infile=%s -outfile=%s -tmpdir=%s",script,seq_file[i],temp_file[i],get_tmp_4_tcoffee());
		}
	      else
		//sprintf ( command, "%s -other_pg %s -infile=%s -outfile=%s",get_string_variable("t_coffee"),script,seq_file[i],temp_file[i]);
		sprintf ( command, "%s -infile=%s -outfile=%s",script,seq_file[i],temp_file[i]);
	      command=substitute(command, "@", " ");
	      //my_system ( command);
	      myexit (my_system(command));
	    }
	  else
	    {
	      pid_list[pid]=npid;
	      //set_pid(pid);
	      npid++;
	      submited++;
	      submited=vwait_npid(submited,max_nproc,nproc);
	    }
	}

      submited=vwait_npid(submited,0,0);
      //Concatenate all the files
      vremove (outfile);
      for (i=0; i<npid; i++)  file_cat (temp_file[i],outfile);

      //Free the process table
      vfree (temp_file);
      vfree (pid_list);
      vfree (seq_file);
      
      free_aln (A);
      
      if ( check_file_exists (outfile) && format_is_fasta(outfile))
	{
	  S=seq2template_seq (S, outfile, F);
	  trim_template_file (outfile,S);
	}
      else if (strstr (command, "webblast.pl"))return S;
      else
	{

	  add_warning (stderr, "Could not Run %s to find templates[%s](Forked mode)\n",command, PROGRAM);
	  return NULL;
	}
      
      vfree (command);
      return S;
    }

  else if (strstr (template_list, "SCRIPT_"))
    {
      char x[299];
      char *tmp1,*command;
      Alignment *A;
      char outfile[1000];
      static char *script;
      static int ntemp;
      char *p;
      int z;
      if (!script)script=(char*)vcalloc ( 1000, sizeof(char));
     
      ntemp++;

      command=(char*)vcalloc ( 1000, sizeof (char));
      tmp1=vtmpnam (NULL);

      A=seq2aln (S,NULL, 0);
      string_array_upper(A->seq_al, A->nseq);
      output_fasta_seq (tmp1, A);
      sprintf ( script, "%s", after_strstr (template_list, "SCRIPT_"));
      fprintf ( stderr, "\n");
      if ((p=strstr (template_list, "@type#")))
	p+=strlen ("@type#");
      if (F)
	{
	  sprintf (outfile, "%s%s_%s%d.template_list", F->path,F->name,template_type2short_type_name(p),ntemp);
	}
      else
	{
	  F=parse_fname (S->file[0]);
	  sprintf (outfile, "%s%s_%s%d.template_list",F->path, F->name,template_type2short_type_name(p),ntemp);
	  free_fname (F);
	}

      script=substitute(script, "@", " -");
      script=substitute(script, "#", "=");

      if  ( strstr (script, "tc_generic_method"))
	{
	  sprintf ( command, "%s -other_pg %s -infile=%s -outfile=%s -tmpdir=%s",get_string_variable ("t_coffee"),script, tmp1,outfile,get_tmp_4_tcoffee());
	}
      else sprintf ( command, "%s -other_pg %s -infile=%s -outfile=%s",get_string_variable("t_coffee"),script, tmp1, outfile);

      vremove (outfile);
      command=substitute(command, "@", " ");

      my_system ( command);

      free_aln (A);
      
      if ( check_file_exists (outfile) && format_is_fasta(outfile))
	{
	  S=seq2template_seq (S, outfile, F);
	  trim_template_file (outfile,S);
	}
      else if (strstr (command, "webblast.pl"))return S;
      else
	{

	  add_warning (stderr, "Could not Run %s to find templates[%s](unforked mode)\n",command, PROGRAM);
	  return NULL;
	}

      vfree (command);
      return S;
    }

  return S;
}
char *trim_template_file (char *file, Sequence *S)//Remove from template file all sequences that cannot be used
{
  Sequence *T;
  int a;
  FILE *fp;

  T=get_fasta_sequence (file,NULL);
  fp=vfopen (file, "w");
  for (a=0; a<T->nseq; a++)
    {
      if (!strstr (T->seq_comment[a], "_P_"))fprintf (fp, ">%s %s\n", T->name[a], T->seq_comment[a]);
      else if (strstr (T->seq_comment[a], "_P_"))
	{
	  int n=name_is_in_list (T->name[a], S->name, S->nseq, 100);
	  if (n!=-1 && seq2P_template_file (S, n))
	    fprintf (fp, ">%s %s\n", T->name[a], T->seq_comment[a]);
	  else
	    HERE ("SKIPPED:%s %s", T->name[a], T->seq_comment [a]);
	}
    }
  vfclose (fp);
  return file;
}
char* seq2template_file (Sequence *S, char *file)
{
  Alignment *A;
  int i;
  if (!S)return file;
  if (file==NULL)file=vtmpnam (NULL);

  seq2template_file2 (S, file, "w");

  for (i=0; i<S->nseq; i++)
    {
      if ( (A=seq2R_template_profile (S, i)))
	{
	  Sequence *S;
	  S=A->S;
	  if (S)seq2template_file2 (A->S, file, "a");
	}
    }
  return file;
}

int seq2template_file2 (Sequence *S, char *file, char *mode)
{
  FILE *fp;
  int i;
  char buf1[10000];
  char buf2[10000];
  struct X_template *X;

  fp=vfopen ( file, mode);
  for ( i=0; i< S-> nseq; i++)
    {
      buf1[0]=0;
      if ( S->T)
	{
	  if (S->T[i])
	    {
	      if ( (X=(S->T[i])->P)){sprintf (buf2, " %s %s ", X->template_type, X->template_file);strcat (buf1, buf2);}
	      /*if ( (X=(S->T[i])->S)){sprintf (buf2, " %s %s ", X->template_type, X->template_file);strcat (buf1, buf2);}*/
	      if ( (X=(S->T[i])->R)){sprintf (buf2, " %s %s ", X->template_type, X->template_file);strcat (buf1, buf2);}
	      if ( (X=(S->T[i])->G)){sprintf (buf2, " %s %s ", X->template_type, X->template_file);strcat (buf1, buf2);}
	      if (buf1[0])fprintf ( fp, ">%s %s\n", S->name[i], buf1);
	    }
	}
    }
  vfclose (fp);
  return EXIT_SUCCESS;
}




int seq2n_X_template ( Sequence *S, char *type)
{
  int a, n;

  for (n=0,a=0; a< S->nseq; a++)
    {
      if ( strm2 (type, "_P_","_*_") && (S->T[a])->P)n++;
      if ( strm2 (type, "_F_","_*_") && (S->T[a])->F)n++;
      if ( strm2 (type, "_S_","_*_") && (S->T[a])->S)n++;
      if ( strm2 (type, "_R_","_*_") && (S->T[a])->R)n++;
      if ( strm2 (type, "_G_","_*_") && (S->T[a])->G)n++;
    }
  return n;
}
struct X_template *fill_X_template ( char *name, char *p, char *token)
{
  struct X_template *X;




  char *k;

  X=(X_template*)vcalloc (1, sizeof (X_template));
  sprintf ( X->seq_name, "%s", name);
  if ( (k=strstr (p, token)))sscanf (k+strlen(token), "%s",X->template_name);
  else sprintf (X->template_name, "%s", p);


  /*Add a Structure HERE*/
  sprintf ( X->template_type, "%s", token);
  if ( strm (token, "_P_"))X->VP=(P_template*)vcalloc (1, sizeof (P_template));
  if ( strm (token, "_F_"))X->VF=(F_template*)vcalloc (1, sizeof (F_template));

  if ( strm (token, "_S_"))X->VS=(S_template*)vcalloc (1, sizeof (S_template));
  if ( strm (token, "_R_"))X->VR=(R_template*)vcalloc (1, sizeof (R_template));
  if ( strm (token, "_G_"))X->VG=(G_template*)vcalloc (1, sizeof (G_template));
  if ( strm (token, "_T_"))X->VT=(T_template*)vcalloc (1, sizeof (T_template));
  if ( strm (token, "_E_"))X->VE=(E_template*)vcalloc (1, sizeof (E_template));
  if ( strm (token, "_U_"))X->VU=(U_template*)vcalloc (1, sizeof (U_template));

  return X;
}

struct X_template* free_X_template ( struct X_template *X)
{
  if (X->VP)
    {
      vfree (X->VP);
    }
  if (X->VF)
    {
      vfree (X->VF);
    }
  if ( X->VS)
    {
      free_sequence ((X->VS)->S, -1);
      vfree (X->VS);
    }
  if ( X->VR)
    {
      free_aln ((X->VR)->A);
      vfree (X->VR);
    }
  if ( X->VG)
    {
      free_sequence ((X->VG)->S, -1);
      vfree (X->VG);
    }

  vfree (X);
  return NULL;
}

FILE * display_sequence_templates (Sequence *S,int i, FILE *io)
{


  io=display_X_template ( (S->T[i])->P, io);

  io=display_X_template ( (S->T[i])->F, io);

  io=display_X_template ( (S->T[i])->S, io);

  io=display_X_template ( (S->T[i])->R, io);
  io=display_X_template ( (S->T[i])->G, io);
  io=display_X_template ( (S->T[i])->T, io);
  io=display_X_template ( (S->T[i])->E, io);

  return io;
}

FILE * display_X_template (struct X_template *X, FILE *io)
{

  if ( !X) return io;
  if ( !strm (X->template_type, "_S_"))
    fprintf (io, "\n\t%s: Template=%s, File=%s",template_type2type_name (X->template_type), X->template_name,X->template_file);
  return io;
}
char *template_type2short_type_name (char *type)
{
  //add_template
  if (!type)return "";
  else if ( strstr (type, "_P_"))       return "pdb";
  else if ( strstr (type, "_F_")) return "rfold";
  else if ( strstr (type, "_S_")) return "seq";
  else if ( strstr (type, "_R_")) return "prf";
  else if ( strstr (type, "_G_")) return "genome";
  else if ( strstr (type, "_E_")) return "ssp";
  else if ( strstr (type, "_T_")) return "tmp";
  else if ( strstr (type, "_U_")) return "unicode";
  else return type;
}
char *template_type2type_name (char *type)
{
  //add_template
  if ( strstr (type, "_P_"))      return "PDB struc";
  else if ( strstr (type, "_F_")) return "RNA Fold";
  else if ( strstr (type, "_S_")) return "Sequeence";
  else if ( strstr (type, "_R_")) return "Profile";
  else if ( strstr (type, "_G_")) return "Genomic";
  else if ( strstr (type, "_E_")) return "Protein Secondary Structure";
  else if ( strstr (type, "_T_")) return "Protein Transmembrane Topology";
  else if ( strstr (type, "_U_")) return "Unicode and strings";

  else return type;
}
struct X_template *fill_F_template ( char *name,char *p, Sequence *S)
{
  /*Profile template*/
  struct X_template *F;

  F=fill_X_template ( name, p, "_F_");
  sprintf (F->template_format , "TCOFFEE_LIBRARY");
  if (!F || !check_file_exists (F->template_name))
    {
      fprintf ( stderr, "Could not fill _F_ (Fold) template for sequence |%s|", name);
      free_X_template (F);
      return NULL;
    }
  else if ( check_file_exists (F->template_name))
    {
      sprintf ( F->template_file, "%s", F->template_name);
    }

  return F;

}
struct X_template *fill_P_template ( char *name,char *p, Sequence *S)
{
  struct X_template *P;
  Sequence *PS;
  Alignment *A;
  int sim, cov, i;
  char *buf;
  char *template_name=NULL;
  char *template_file=NULL;
  char *k;
  
  if (!name || !S) return NULL;
  else if (!p)
    {
      add_warning(stderr, "_P_ Template | %s | Could not be found",name);
      return NULL;
    }
 
  if ( (k=strstr (p, "_P_")))template_name=csprintf (template_name, "%s",k+4);
  else template_name=csprintf (template_name, "%s", p); 
  
  if (!(template_name=is_pdb_struc(template_name)))//make sure a valid PDB exists
    {
      add_warning(stderr, "_P_ Template | %s | Could not be found",name);
      return NULL;
    }
  if (!(template_file=fix_pdb_file(template_name)))//make sure the PDB can be used
    {
      add_warning(stderr, "_P_ Template | %s | Could not be used",name);
      vfree(template_name);
      return NULL;
    }
 
  //We now have a valid PDB file 
  P=fill_X_template (name,p, "_P_");
  sprintf (P->template_format , "pdb");
  sprintf (P->template_file, "%s",template_file);
  sprintf (P->template_name, "%s",name);
  vfree(template_file); vfree(template_name);
  
  /*Check the target sequence is similar enough*/

  PS=get_pdb_sequence (P->template_file);



  if ( PS==NULL)
    {
      add_warning( stderr, "_P_  Template |%s| Could not be used for Sequence |%s|: Structure Not Found", P->template_name, name);
      free_X_template (P);P=NULL;
    }
  else
    {
      int minsim=get_int_variable ("pdb_min_sim");
      int mincov=get_int_variable ("pdb_min_cov");


      i=name_is_in_list (name, S->name, S->nseq, 100);

      A=align_two_sequences (S->seq[i], PS->seq[0],"idmat",-3,0, "fasta_pair_wise");

      sprintf ( A->name[0], "seq");
      sprintf ( A->name[1], "pdb");
      cov=aln2coverage (A, 0);
      sim=aln2sim (A, "idmat");

      if (sim<=minsim)
	{
	  add_information( stderr, "_P_  Template %s Could not be used for Sequence %s: Similarity too low [%d, Min=%d]",P->template_name,name,sim,minsim);
	  add_information( stderr, "If you want to include %s in anycase,add -pdb_min_sim=%d to the command line",name,sim);
	  print_aln (A);
	  free_X_template (P);
	  P=NULL;
	}
      else if ( cov<=mincov)
	{
	  add_information(stderr, "_P_  Template |%s| Could not be used for Sequence |%s|: Coverage too low [%d, Min=%d]",P->template_name,name, cov, mincov);
	  add_information( stderr, "If you want to include this sequence in anycase add -pdb_min_cov=%d to the command line", cov);
	  print_aln (A);
	  free_X_template (P);P=NULL;
	}
      free_aln(A);
      free_sequence (PS, -1);
    }

  return P;
}

struct X_template *fill_P_template_new ( char *name,char *p, Sequence *S)
{
  struct X_template *P;
  Sequence *PS;
  Alignment *A;
  int sim, cov, i;
  char *buf;
  
 
  P=fill_X_template ( name, p, "_P_");
 
  sprintf (P->template_format , "pdb");

  if (!P ||(check_file_exists (P->template_name) && !is_pdb_file (P->template_name) ))
    {
      
      fprintf ( stderr, "Could not fill _P_ template for sequence |%s|", name);
      free_X_template (P);
      return NULL;
    }
  else if ( check_file_exists (P->template_name))
    {
      
      sprintf ( P->template_file, "%s", P->template_name);
      buf=path2filename (P->template_name);
      if (P->template_name!=buf)
	{
	  sprintf ( P->template_name, "%s",buf );
	  vfree (buf);
	}
    }
   else
     {
       char *st;
      
       st=is_pdb_struc (P->template_name);

       if (st)
	 {
	   if (st!=P->template_file)sprintf ( P->template_file, "%s", st);
	 }
     }
  
  /*Make a first run to fix relaxed PDB files*/
  buf=fix_pdb_file (P->template_file);
  
  if ( buf!=P->template_file)
  {

    sprintf ( P->template_file, "%s",buf);

  }

  /*Check the PDB FILE EXISTS*/

  if (!is_pdb_file (P->template_file))
    {

      if (p)add_warning(stderr, "_P_ Template | %s | Could not be found\n",p);
      else if (name)add_warning(stderr, "_P_ Template | %s | Could not be found\n",name);
      free_X_template (P);
      return NULL;
    }
  else
    {
      buf= get_pdb_id (P->template_file);
      if (buf!=(P->VP)->pdb_id)
	{
	  sprintf ((P->VP)->pdb_id, "%s", buf);
	  vfree (buf);
	}
    }

  /*Check the target sequence is similar enough*/

  PS=get_pdb_sequence (P->template_file);



  if ( PS==NULL)
    {
      add_warning( stderr, "_P_  Template |%s| Could not be used for Sequence |%s|: Structure Not Found", P->template_name, name);
      free_X_template (P);P=NULL;
    }
  else
    {
      int minsim=get_int_variable ("pdb_min_sim");
      int mincov=get_int_variable ("pdb_min_cov");


      i=name_is_in_list (name, S->name, S->nseq, 100);

      A=align_two_sequences (S->seq[i], PS->seq[0],"idmat",-3,0, "fasta_pair_wise");

      sprintf ( A->name[0], "seq");
      sprintf ( A->name[1], "pdb");
      cov=aln2coverage (A, 0);
      sim=aln2sim (A, "idmat");

      if (sim<=minsim)
	{
	  add_information( stderr, "_P_  Template %s Could not be used for Sequence %s: Similarity too low [%d, Min=%d]",P->template_name,name,sim,minsim);
	  add_information( stderr, "If you want to include %s in anycase,add -pdb_min_sim=%d to the command line",name,sim);
	  print_aln (A);
	  free_X_template (P);
	  P=NULL;
	}
      else if ( cov<=mincov)
	{
	  add_information(stderr, "_P_  Template |%s| Could not be used for Sequence |%s|: Coverage too low [%d, Min=%d]",P->template_name,name, cov, mincov);
	  add_information( stderr, "If you want to include this sequence in anycase add -pdb_min_cov=%d to the command line", cov);
	  print_aln (A);
	  free_X_template (P);P=NULL;
	}
      free_aln(A);
      free_sequence (PS, -1);
    }

  return P;
}

int **seq2pdb_index (Sequence *S)
{
  //index contains for each residue of S, the Pos value of the ATOM lines in the original PDB file
  int a, b;
  int **index=(int**)vcalloc (S->nseq, sizeof (int*));
  for (a=0; a<S->nseq; a++)
    {
      char *file=seq2P_template_file (S, a);
      if (file)
	{
	  int *plist;
	  int l=strlen (S->seq[a]);
	  Sequence  *PS=get_pdb_sequence_from_field(file, "ATOM");
	  plist=pdb2atom_pos_list(file);
	  Alignment *A=align_two_sequences (S->seq[a], PS->seq[0],"idmat",-3,0, "fasta_pair_wise");
	  int seq=0;
	  int pdb=0;
	  //print_aln (A);
	  index[a]=(int*)vcalloc(l, sizeof (int));
	  for (b=0; b<l; b++)index[a][b]=-1;
	  
	  
	  for (b=0; b<A->len_aln; b++)
	    {
	      
	      int seq_r=1-is_gap(A->seq_al[0][b]);
	      int pdb_r=1-is_gap(A->seq_al[1][b]);
	      seq+=seq_r;
	      pdb+=pdb_r;
	      if (seq_r)
		{
		  if (pdb_r)index[a][seq-1]=plist[pdb-1];
		  else index[a][seq-1]=-1;
		}
	    }
	  free_aln(A);
	  free_sequence (PS, -1);
	  vfree (plist);
	}
    }
  return index;
}			     

struct X_template *fill_S_template ( char *name,char *p, Sequence *Seq)
{
  struct X_template *S;
  S=fill_X_template ( name, p, "_S_");
  if ( strm (name, p))sprintf ( S->template_file, "%s",output_fasta_seqX (NULL,"w",Seq,NULL, seq_name2index (name, Seq)));
  (S->VS)->S=get_fasta_sequence (S->template_file, NULL);
  return S;
}
struct X_template *fill_R_template ( char *name,char *p, Sequence *S)
{
  /*Profile template*/
  struct X_template *R;


  R=fill_X_template ( name, p, "_R_");
  sprintf (R->template_format , "fasta_aln");
  
  if (!isfile(R->template_name))
    {
      static char *buf;
      buf=csprintf (buf, "%s/%s", get_cache_dir(),R->template_name);
      sprintf (R->template_name, "%s", buf);
    }

  if (!is_aln(R->template_name) && !is_seq (R->template_name))
    {

      add_information ( stderr, "_R_ Template %s Could not be found\n",R->template_name);
      free_X_template (R);
      return NULL;
    }
  else
    {
      int s;
      Sequence *S1;
      Alignment *A1;

      (R->VR)->A=main_read_aln (R->template_name, NULL);

      if ( !S)
	sprintf ( R->template_file, "%s", R->template_name);
      else
	{
	  s=name_is_in_list(name, S->name, S->nseq, 100);
	  if ( s!=-1)
	    {
	      char **lseq =(char**)vcalloc (1, sizeof (char*));
	      char **lname=(char**)vcalloc (1, sizeof (char*));
	      lseq [0]=S->seq [s];
	      lname[0]=S->name[s];
	      
	      S1=fill_sequence_struc (1,lseq,lname, NULL);
	      vfree (lseq); vfree(lname);
	      A1=seq2aln (S1,NULL, RM_GAP);
	      
	      (R->VR)->A=trim_aln_with_seq (A1, (R->VR)->A);
	      
	      sprintf ( R->template_file, "%s", vtmpnam (NULL));
	      output_clustal_aln (R->template_file, (R->VR)->A);
	    }
	  else
	    sprintf ( R->template_file, "%s", R->template_name);
	}
      (R->VR)->A=aln2profile ((R->VR)->A);

      //free_data_in_aln ((R->VR)->A);

    }
  return R;
}

struct X_template *fill_T_template ( char *name,char *p, Sequence *S)
{
  /*Profile template*/
  struct X_template *T;

  T=fill_X_template ( name, p, "_T_");
  sprintf (T->template_format , "fasta_seq");

  if (!is_aln(T->template_name) && !is_seq (T->template_name))
    {

      add_information ( stderr, "_T_ Template %s Could not be found\n",T->template_name);
      free_X_template (T);
      return NULL;
    }
  else
    {

      (T->VT)->S=main_read_seq(T->template_name);
      sprintf ( T->template_file, "%s", T->template_name);
    }
  return T;
}
//add template
struct X_template *fill_U_template ( char *name,char *p, Sequence *S)
{
  /*Profile template*/
  struct X_template *U;

  U=fill_X_template ( name, p, "_U_");
  sprintf (U->template_format , "string list");

  if (!check_file_exists(U->template_name))
    {
      add_information ( stderr, "_U_ Template %s Could not be found\n",U->template_name);
      free_X_template (U);
      return NULL;
    }
  else
    {
      //(U->VU)->list=file2string(U->template_name);
      sprintf ( U->template_file, "%s", U->template_name);
    }
  return U;
}
struct X_template *fill_E_template ( char *name,char *p, Sequence *S)
{
  /*Profile template*/
  struct X_template *E;


  E=fill_X_template ( name, p, "_E_");
  sprintf (E->template_format , "fasta_seq");

  if (!is_aln(E->template_name) && !is_seq (E->template_name))
    {

      add_information ( stderr, "_E_ Template %s Could not be found\n",E->template_name);
      free_X_template (E);
      return NULL;
    }
  else
    {
      (E->VE)->S=main_read_seq (E->template_name);
      sprintf ( E->template_file, "%s", E->template_name);
    }
  return E;
}
struct X_template *fill_G_template ( char *name,char *p, Sequence *S)
{
  struct X_template *G;
  G=fill_X_template ( name, p, "_G_");
  sprintf (G->template_format , "fasta_seq");

  /*1: Get the sequence from another file if needed*/
  if ( strm (name, p))sprintf ( G->template_file, "%s",output_fasta_seqX (NULL,"w",S,NULL, seq_name2index (name, S)));
  else if ( strstr (p, "SEQFILE_"))
    {
      Sequence *ST;
      int i2;


      ST=main_read_seq (after_strstr ( p,"SEQFILE_G_"));

      i2=seq_name2index (name, ST);
      if ( i2!=-1)
	{
	  sprintf ( G->template_file, "%s",output_fasta_seqX (NULL,"w",ST,NULL, i2));
	  sprintf ( G->template_name, "%s", name);
	}
      free_sequence (ST, -1);
    }
  else sprintf (G->template_file, "%s", G->template_name);


  /*2: Put the template in VG->S*/
  if (!is_seq (G->template_file))
    {
      add_information ( stderr, "_G_ Template %s Could not be found \n",p);

      free_X_template (G);
      return NULL;
    }
  else
    {
      (G->VG)->S=get_fasta_sequence (G->template_file, NULL);
    }
  return G;
}


char *seq2T_value ( Sequence *S, int n, char *value, char *type)
{
  static char *rv_buf;
  X_template *X;

  if ( !rv_buf)rv_buf=(char*)vcalloc (100, sizeof(char));
  if (!(X=seq_has_template (S, n, type)))return NULL;
  else
    {
      if (strm (value, "template_file"))return X->template_file;
      else if ( strm (value, "template_name"))return X->template_name;
      else if ( strm (value, "seq_name"))return X->seq_name;
      else if (strm (type, "_P_"))
	{
	  if ( strm (value, "pdb_id"))return (X->VP)->pdb_id;
	}
      else if ( strm (type, "_R_"))
	{
	  if ( strm (value, "A"))
	    {
	      if ((X->VR)->A)
		{sprintf ( rv_buf, "%ld", (long)(X->VR)->A);return rv_buf;}
	      else return NULL;
	    }
	}

    }
  return NULL;
}
char *seq2P_pdb_id (Sequence *S, int n)
{
  if (!S->T || !S->T[n] || !(S->T[n])->P ) return NULL;
  else return ((S->T[n])->P)->template_name;
}


char *seq2P_template_file(Sequence *S, int n)
{

  return seq2T_value (S, n, "template_file", "_P_");
}

char *profile2P_template_file (Sequence *S, int n)
{
  Alignment *A;
  int a;
  char *p;

  if ( !(A=seq2R_template_profile (S, n)))return NULL;
  for (a=0; a<A->nseq; a++)
    {
      if ((p=seq2P_template_file (A->S, a))!=NULL)return p;
    }
  return NULL;
}



Alignment * seq2R_template_profile (Sequence *S, int n)
{
  X_template *X;

  return (Alignment *)atop(seq2T_value (S, n, "A", "_R_"));

  if (!(X=seq_has_template (S, n, "_R_")))return NULL;
  else
    {
      if (!(X->VR))return NULL;
      else return (X->VR)->A;
    }
  return NULL;



}
char * seq2E_template_string (Sequence *S, int n)
{
  struct X_template *T;

  if ( (T=seq_has_template (S, n, "_E_"))!=NULL)
    return  ((T->VE)->S)->seq[0];
  else
    return NULL;
}
//add template
int* seq2U_template (Sequence *S, int n)
{
   struct X_template *T;

   if ( (T=seq_has_template (S, n, "_U_"))!=NULL)
     return  (T->VU)->list;
   else
     return NULL;
}
char * seq2T_template_string (Sequence *S, int n)
{
  struct X_template *T;

  if ( (T=seq_has_template (S, n, "_T_"))!=NULL)
    return  ((T->VT)->S)->seq[0];
  else
    return NULL;
}
int seq2n_template (Sequence *S, char *mode)
{
  int a, n;
  if (!S || !S->nseq || !mode || !mode)return 0;
  for (a=0, n=0; a<S->nseq; a++)
    {
      if (seq_has_template(S,a,mode))n++;
    }
  return n;
}
struct X_template* seq_has_template ( Sequence *S, int n, char *mode)
{
  Template *T;

  if ( !S || !mode) return NULL;
  else if ( n<0 || n>=S->nseq)return NULL;
  else if ( !(S->T)) return NULL;
  else if ( !(S->T[n]))return NULL;

  T=S->T[n];
  //ADD STRUCTURE
  //add template
  if      ( strm (mode, "_P_"))return T->P;
  else if ( strm (mode, "_F_"))return T->F;
  else if ( strm (mode, "_S_"))return T->S;
  else if ( strm (mode, "_R_"))return T->R;
  else if ( strm (mode, "_T_"))return T->T;
  else if ( strm (mode, "_E_"))return T->E;
  else if ( strm (mode, "_U_"))return T->U;
  else if ( strm (mode, "_G_"))return T->G;
  else return NULL;
}

char ** name2random_subset (char **in_name, int n_in, int n_out)
{
  char **out_name;

  int **list;
  int a,max;


  vsrand (0);
  max=n_in*10000;
  out_name=declare_char (n_out,MAXNAMES+1 );
  list=declare_int (n_in, 2);

  for (a=0; a<n_in; a++)
      {
	list[a][0]=a;
	list[a][1]=rand ()%max;
      }
  sort_int ( list,2, 1, 0, n_in-1);

  for ( a=0; a<n_in; a++)
  {
    sprintf ( out_name[a], "%s", in_name[list[a][0]]);
  }
  free_int (list, -1);
  return out_name;
}

Alignment * aln2random_order (Alignment *A)
{

  char **name_list;

  name_list=name2random_subset (A->name, A->nseq, A->nseq);

  A=reorder_aln (A, name_list, A->nseq);

  free_char (name_list, -1);
  return A;
}
Alignment *aln2jacknife (Alignment *A, int nseq, int len)
{
  int a, b;

  if (nseq!=0 && nseq<A->nseq)
    {
      char **name;

      name=name2random_subset (A->name, A->nseq, nseq);
      A=reorder_aln (A, name, nseq);
      free_char (name, -1);
    }

  if (len!=0 && len<A->len_aln)
    {
      int **l;
      Alignment *B;

      l=declare_int (A->len_aln, 2);
      for (a=0; a< A->len_aln; a++)
	{
	  l[a][0]=a;
	  l[a][1]=rand()%(A->len_aln*1000);
	}
      sort_int ( l,2, 1, 0, A->len_aln-1);
      B=copy_aln (A, NULL);
      for ( a=0; a< len; a++)
	{
	  for ( b=0; b<A->nseq; b++)
	    {
	      A->seq_al[b][a]=B->seq_al[b][l[a][0]];
	    }
	}
      for (b=0; b<A->nseq; b++)A->seq_al[b][len]='\0';
      free_aln (B);
      free_int (l, -1);
    }
  return A;
}
Alignment * aln2scramble_seq (Alignment *A)
{
  int **list;
  char **name_list;
  int a,max;

  max=100*A->nseq;
  vsrand (0);

  list=declare_int (A->nseq, 2);
  name_list=(char**)vcalloc (A->nseq, sizeof (char*));


  for (a=0; a<A->nseq; a++)
      {
	list[a][0]=a;
	list[a][1]=rand ()%max;
      }
  sort_int ( list,2, 1, 0, A->nseq-1);

  for ( a=0; a< A->nseq; a++)
    name_list[a]=A->seq_al[a];
  for (a=0; a<A->nseq; a++)
    {
      A->seq_al[a]=name_list[list[a][0]];
    }
  vfree (name_list);
  free_int (list, -1);
  return aln2random_order (A);
}


Alignment * shuffle_aln ( Alignment *A,int N, char *name_i, char *mode)
{
  int **index;
  int a, b;
  FILE *fp;
  FILE *fp2;
  char *fname;
  char *name;
  char *tmp;

  
  tmp=vtmpnam (NULL);
  fp2=vfopen (tmp, "w");
  index=declare_int (A->nseq, 2);
  
  if (!name_i){name= (char*)vcalloc (100, sizeof (char)); sprintf (name, "replicate");}
  else 
    {
      name=(char*)vcalloc (strlen (name_i)+1, sizeof(char));
      sprintf (name, "%s", name_i);
    }
  fname=(char*)vcalloc ( strlen (name)+100, sizeof (char));
  

  if (!mode)
    for (a=0; a<A->nseq; a++)
      {
	ungap (A->seq_al[a]);
      }
  
  
  for (a=0; a<N; a++)
    {
      for (b=0; b<A->nseq; b++)
	{
	  index[b][0]=b;
	  index[b][1]=rand()%(A->nseq*100);
	}
      sort_int (index,2,1, 0, A->nseq-1);
      sprintf (fname, "%s.%d.shuffled.fa",name,a+1);
      fprintf (fp2,">%s\n", fname);
      
      fp=vfopen (fname, "w");
      for (b=0;b<A->nseq; b++)
	{
	  int i=index[b][0];
	  fprintf (fp,">%s %s\n%s\n", A->name[i], A->seq_comment[i], A->seq_al[i]);
	}
      vfclose (fp);
    }
  
  vfclose (fp2);
  vfree (name); vfree (fname);
  free_aln (A);
  A=main_read_aln (tmp, NULL);
  return A;
}
	  
Alignment * reorder_aln ( Alignment *A, char **iname, int nseq)
{
 
  FILE*fp;
  int a, i, n;
  char **seq_al=(char**)vcalloc (nseq, sizeof (char*));
  char **name=(char**)vcalloc (nseq, sizeof (char*));
  char **aln_comment=(char**)vcalloc(nseq, sizeof (char*));
  int  **order=(int**)vcalloc(nseq, sizeof (int*));
  
  
  for ( n=0,a=0; a<nseq; a++)
    {
      if ((i =name_is_in_hlist ( iname[a],A->name, A->nseq))!=-1)
	{
	  seq_al[n]=A->seq_al[i];
	  name[n]=A->name[i];
	  aln_comment[n]=A->aln_comment[i];
	  order[n]=A->order[i];
	  n++;
	}
    }
  A->nseq=n;
  for (a=0; a<A->nseq; a++)
    {
      A->seq_al[a]=seq_al[a];
      A->name[a]=name[a];
      aln_comment[a]=A->aln_comment[a];
      A->order[a]=order[a];
    }
  vfree (seq_al);vfree(name); vfree(aln_comment); vfree(order);
  
  if ( A->A)A->A=reorder_aln(A->A, name, nseq);
  
  return A;
}

Sequence * reorder_seq_2 ( Sequence *A, int **order,int field, int nseq)
	{
	  char **name;
	  int a;

	  if (!A || !order) return A;
	  name=declare_char (A->nseq, 100);
	  for (a=0; a<nseq; a++)
	    sprintf ( name[a], "%s", A->name[order[a][field]]);
	  A=reorder_seq (A, name,nseq);
	  free_char (name, -1);
	  return A;
	}
Sequence * reorder_seq ( Sequence *A, char **name, int nseq)
	{
	int a,sn;
	Sequence *nA;


	nA=duplicate_sequence (A);


	for ( a=0; a< nseq; a++)
	  {
	    sn=name_is_in_list (name[a] ,nA->name, nA->nseq, 100);
	    if (sn==-1)continue;

	    if ( nA->file)       sprintf ( A->file[a], "%s", nA->file[sn]);

	    //if ( nA->seq_comment)sprintf ( A->seq_comment[a], "%s", nA->seq_comment[sn]);
	    //if ( nA->aln_comment)sprintf ( A->aln_comment[a], "%s", nA->aln_comment[sn]);
	    //sprintf ( A->seq[a], "%s", nA->seq[sn]);
	    //sprintf ( A->name[a], "%s", nA->name[sn]);
	    A->seq[a]=csprintf (A->seq[a], "%s", nA->seq[sn]);
	    A->name[a]=csprintf (A->name[a], "%s", nA->name[sn]);
	    if ( nA->seq_comment)A->seq_comment[a]=csprintf (A->seq_comment[a], "%s", nA->seq_comment[sn]);
	    if ( nA->aln_comment)A->aln_comment[a]=csprintf (A->aln_comment[a], "%s", nA->aln_comment[sn]);
	    
	    
	    A->len[a]=nA->len[sn];
	    
	    A->T[a][0]=nA->T[sn][0];
	  }
	A->nseq=nseq;
	free_sequence (nA, nA->nseq);

	return A;
}

char * concatenate_seq ( Sequence *S, char *conc, int *order)
        {
	    int a;

	    vfree (conc);
	    conc=(char*)vcalloc ( S->nseq*S->max_len, sizeof (char));

	    for ( a=0; a< S->nseq; a++)
	        {
		    conc=strcat ( conc, S->seq[order[a]]);
		}
	    return conc;

	}




Alignment * rotate_aln ( Alignment *A, char *name)
{
  Alignment *B;
  int a, b;

  B=declare_aln2 (A->len_aln, A->nseq+1);
  for ( a=0; a< A->nseq; a++)
    for ( b=0; b< A->len_aln; b++)
      {
	B->seq_al[b][a]=A->seq_al[a][b];
      }
  for (a=0; a< A->len_aln; a++)
    if (name && name[0])sprintf ( B->name[a], "%s_%s%d", name, (a<9)?"0":"",a+1);
    else
      sprintf ( B->name[a], "%d", a+1);


  for (a=0; a< A->len_aln; a++)B->seq_al[a][A->nseq]='\0';
  B->len_aln=A->nseq;
  B->nseq=A->len_aln;
  /*free_aln (A);*/
  return B;
}
int invert_seq_file    (char *file)
{
  Sequence *A;
  A=main_read_seq (file);
  A=invert_seq2(A);
  output_fasta_seqS (file,A);
  free_sequence(A, A->nseq);
  return 1;
}

int invert_aln_file    (char *file)
{
  Alignment *A;
  A=main_read_aln (file,NULL);
  A=invert_aln(A);
  output_fasta_aln (file,A);
  free_aln(A);
  return 1;
}

Alignment * invert_aln ( Alignment *A)
{
  char *buf;
  int l, a, b, c;

  for ( a=0; a< A->nseq; a++)
    {
        l=strlen ( A->seq_al[a]);
	buf=(char*)vcalloc ( l+1,sizeof (char) );

	for ( c=l-1,b=0; b< l; b++, c--)
	  {
	    buf[c]=A->seq_al[a][b];
	  }
	buf[l]='\0';
	sprintf ( A->seq_al[a], "%s", buf);
    }
  vfree(buf);
  return A;
}
Sequence * invert_seq2 (Sequence *A)
{
   char *buf;
  int l, a, b, c;

  for ( a=0; a< A->nseq; a++)
    {
        l=strlen ( A->seq[a]);
	buf=(char*)vcalloc ( l+1,sizeof (char) );

	for ( c=l-1,b=0; b< l; b++, c--)
	  {
	    buf[c]=A->seq[a][b];
	  }
	buf[l]='\0';
	sprintf ( A->seq[a], "%s", buf);
    }
  vfree(buf);
  return A;
}
char * complement_string (char *s)
{
  char *buf;
  int l, a, b, c;

  l=strlen (s);
  for ( b=0; b< l; b++)
    {
      char r;
      r=s[b];
      if ( r=='a')r='t';
      else if (r=='A')r='T';
      else if (r=='t')r='a';
      else if (r=='T')r='A';
      else if (r=='g')r='c';
      else if (r=='G')r='C';
      else if (r=='c')r='g';
      else if (r=='C')r='G';
      s[b]=r;
    }

  return invert_string (s);
}
Alignment * complement_aln ( Alignment *A)
{
  char *buf;
  int l, a, b, c;

  for ( a=0; a< A->nseq; a++)
    {
      A->seq_al[a]=complement_string (A->seq_al[a]);
    }

  return A;
}

Alignment * extract_nol_local_aln(Alignment *A, int start, int max_end)
     {
     A=extract_aln ( A, start, max_end);
     A=trunkate_local_aln (A);
     return A;
     }

Alignment * alnpos_list2block (Alignment *A, int n, char **in_list)
{
  int *pos;
  int a;
  char **list;
  int list_declared=0;
  Alignment *B;

  if (check_file_exists (in_list[0]))
    {
      int mn;
      char ***tmp_list;

      mn=count_n_line_in_file (in_list[0]);
      list=declare_char (mn, 100);
      list_declared=1;
      tmp_list=file2list (in_list[0], " ");
      a=0;
      n=0;
      while (tmp_list[a])
	{
	  if (tmp_list[a][1][0]!='!')
	    {
	      sprintf (list[n++], "%s", tmp_list[a][1]);
	    }
	  a++;
	}
      free_arrayN ((void **)tmp_list, 3);
    }
  else
    {
      list=in_list;
    }


  pos=(int*)vcalloc (A->len_aln, sizeof (int));
  for (a=0; a<n; a++)
    {

      if (strstr (list[a], "-"))
	{
	  int start, end, x;
	  x=sscanf (list[a], "%d-%d", &start, &end);
	  if (x!=2 || !A || start<=0 || start>=end || end>A->len_aln+1)
	    {
	      add_warning ( stderr, "Illegal coordinates in extract_pos_list [%s]", list[a]);
	      return A;
	    }
	  start--; end--;
	  for (a=start; a<end; a++)pos[a]=1;
	}
      else
	{
	  int p;
	  p=atoi (list[a]);
	  if (p<1 || p>A->len_aln)
	    {
	      add_warning ( stderr, "Illegal coordinates in extract_pos_list [%s]", list[a]);
	    }
	  p--;
	  pos[p]=1;
	}
    }
  B=alnpos2block(A, pos, NULL);
  vfree (pos);
  if ( list_declared)free_char (list, -1);

  return B;
}
Alignment * aln2block   (Alignment  *A, int start, int end, Alignment *B)
{
  if ( !A || start<=0 || start>=end || end>A->len_aln+1)
    {
      add_warning ( stderr, "Illegal coordinates in extract_block start=%d end=%d len=%d [Note : [start-end[, with [1...n] ** Block Ingored", start, end, A->len_aln);
      return A;
    }
  else
    {
      int *pos, p;
      start--;
      end--;
      pos=(int*)vcalloc (A->len_aln, sizeof (int));
      for (p=start;p<end;p++)
	    {
	      pos[p]=1;
	    }
      B=alnpos2block (A, pos, B);
      vfree (pos);
      return B;
    }
}
Alignment * alnpos2block   (Alignment  *A, int *pos, Alignment *B)
{

  //extract a subset of B without over-writing A
  int a, b;

  B=copy_aln (A, B);
  B->len_aln=0;
  for (a=0; a<=A->len_aln; a++)
    {
      if ( pos[a]!=0 || a==A->len_aln)
	{
	  for ( b=0; b<A->nseq; b++)
	    B->seq_al[b][B->len_aln]=A->seq_al[b][a];
	  if ( a!=A->len_aln)B->len_aln++;
	}
    }

  return B;
}
Alignment * extract_aln ( Alignment *A, int start, int end)
{
  return extract_aln2 ( A, start, end, "cons");
}

Alignment * extract_aln2 ( Alignment *A, int in_start, int in_end, char *seq)
     {
       char *tmp;
       FILE *fp;


       tmp=vtmpnam (NULL);
       fp=vfopen (tmp, "w");
       fprintf ( fp, "%s %d %d\n", seq, in_start, in_end);
       vfclose (fp);
       return extract_aln3 (A,tmp);
     }
Alignment * extract_aln3 ( Alignment *B, char *file)
{
  char *tmp=vtmpnam(NULL);
  char ***list=file2list(file, " ");
  FILE*fp;
  int start, end, i, n, c;
  
  if (!B || !file ||!(list=file2list (file, " ")))return B;
  fp=vfopen ( tmp, "w");
  n=0;    
  while (list[n])
    {
      if ( list[n][1][0]=='#'){;}
      else
	{
	  int start=atoi (list[n][2])-1;
	  if (atoi(list[n][0])==4)end=atoi (list[n][3]);
	  else end=B->len_aln;
	  
	  if ( end>B->len_aln)
	    add_warning ( stderr, "Illegal coordinates [%s offset: %d] (line %d) ** Line ignored",list[n][1],start, n+1);
	  else if ( strm (list[n][1], "cons"))
	    {
	      for ( i=0; i<B->nseq; i++)
		{
		  c=B->seq_al[i][end];
		  B->seq_al[i][end]='\0';
		  
		  fprintf ( fp, ">%s\n%s\n", B->name[i], B->seq_al[i]+start);
		  B->seq_al[i][end]=c;
		}
	    }
	  else if ((i=name_is_in_hlist (list[n][1],B->name, B->nseq))!=-1)
	    {
	      c=B->seq_al[i][end];
	      fprintf ( fp, ">%s\n%s\n", B->name[i], B->seq_al[i]+start);
	      B->seq_al[i][end]=c;
	    }
	  else 
	    {
	      add_warning ( stderr, "Seq %s does not belong to the alignment (line %d) ** Line ignored",list[n][1],n+1);
	    }
	  n++;
	}
    }
  vfclose (fp);
  free_arrayN(list, 3);
  
  
  return quick_read_fasta_aln (B,tmp);
}
	      

Alignment * trunkate_local_aln ( Alignment *A)
     {
     int a, b;
     int **pos;
     int **cache;
     int seq;


     cache=declare_int (return_max_int (A->order,read_size_int ( A->order,sizeof (int*)),0)+1,return_max_int (A->order,read_size_int ( A->order,sizeof (int*)),1)+A->len_aln+1);
     pos=aln2pos_simple(A,A->nseq);

     for ( b=0; b<A->len_aln; b++)
	 for ( a=0; a< A->nseq; a++)
	     {
	     seq=A->order[a][0];
	     if ( pos[a][b]<=0);
	     else if ( pos[a][b]>0)
		 {

		 if (cache[seq][pos[a][b]]==0)cache[seq][pos[a][b]]++;
		 else if ( cache[seq][pos[a][b]]>=1)
		      {
		      cache[seq][pos[a][b]]++;
		      A->seq_al[a][b]='\0';
		      }
		 }
	     }

     A->len_aln=get_shortest_string ( A->seq_al, A->nseq, NULL, NULL);
     pad_string_array ( A->seq_al, A->nseq, A->len_aln, '-');

     free_int (pos, -1);
     free_int ( cache,-1);


     return A;
     }

int get_nol_aln_border ( Alignment *A, int start, int direction)
     {
     int a, b;
     int **pos;
     int **cache;
     int seq,end;

     /*This Function Returns the limit position for a non overlaping alignment*/

     cache=declare_int (return_max_int (A->order,read_size_int ( A->order,sizeof (int*)),0)+1,return_max_int (A->order,read_size_int ( A->order,sizeof (int)),1)+A->len_aln+1);
     pos=aln2pos_simple(A,A->nseq);
     end=(direction==GO_RIGHT)?A->len_aln:-1;


     for ( b=start; b!=end;b+=direction)
	 for ( a=0; a< A->nseq; a++)
	     {
	     seq=A->order[a][0];
	     if ( pos[a][b]<=0);
	     else if ( pos[a][b]>0)
		 {

		 if (cache[seq][pos[a][b]]==0)cache[seq][pos[a][b]]++;
		 else if ( cache[seq][pos[a][b]]>=1)
		      {
		      cache[seq][pos[a][b]]++;
		      free_int(cache, -1);
		      return b-direction;
		      }
		 }
	     }

     free_int ( cache,-1);
     free_int (pos, -1);
     return end-direction;
     }





char * extract_defined_seq ( char *in, int in_of, int in_start, int *aa_def, int dir, int *out_start, char *out)
     {
     int start=0, end,l;
     int b, c, d;



     if ( dir==GO_LEFT){start=in_start-1;}
     else if ( dir==GO_RIGHT){start=in_start+1;}

     end=start;
     while (aa_def[end]!=UNDEFINED)
	 {
	 end+=dir;
	 }
     end-=dir;

     if (end<start)SWAP(end,start);

     l=strlen ( in);
     out_start[0]=-1;
     for (b=0,d=0,c=in_of;b<l; b++)
         {
	 c+=1-is_gap(in[b]);
	 if ( c>=start && c<=end)
	     {
	     if ( out_start[0]==-1)out_start[0]=c-!is_gap(in[b]);
	     out[d++]=in[b];
	     }
	 }
     out[d]='\0';


     return out;
     }

Alignment * aln2N_replicate (Alignment *A,char *nn, char *name)
{
  int a, n;
  char *fname;

  fname=(char*)vcalloc (100, sizeof (char));
  if (nn)n=atoi(nn);
  else n=100;
  if (!name){name=(char*)vcalloc (100, sizeof (char)); sprintf (name, "replicate");}


  for (a=0; a< n;a++)
    {
      FILE *fp;
      sprintf (fname, "%s.%d.rep",name, a+1);
      fp=vfopen (fname, "w");

      vfclose(aln2replicate (A, fp));
      fprintf ( stdout, ">%s Alignment Replicate #%d\n",fname, a+1);
    }
  myexit (EXIT_SUCCESS);
}
FILE *aln2replicate (Alignment *A, FILE *fp)
{
  int a, b;
  int *p;
  float tot=0;
  float corr;
  if (A->col_weight)for (a=0; a<A->len_aln; a++)tot+=A->col_weight[a];
  else tot=A->len_aln;

  p=(int*)vcalloc (A->len_aln, sizeof (int));
  corr=(float)A->len_aln/tot;

  for (a=0; a<A->len_aln; a++)
    {
      int x;
      x=rand()%(int)tot;

      p[a]=(int)(x*corr);
    }

  for (a=0; a<A->nseq; a++)
    {
      fprintf ( fp, ">%s\n", A->name[a]);
      //for (b=0;b<A->len_aln; b++)fprintf ( stdout, "%d ", (int)p[b]);
      for (b=0;b<A->len_aln; b++)fprintf ( fp, "%c", A->seq_al[a][p[b]]);
      fprintf ( fp, "\n");
    }

  vfree (p);
  return fp;
}

Alignment * orthologous_concatenate_aln (Alignment *A, Sequence *S, char *mode)
{
  Alignment *C;
  char **name, *cname;
  int nname=0;
  int a, b,c, i;

  if (mode && strm (mode, "voronoi"))seq_weight2species_weight (A, S);


  cname=(char*)vcalloc ( 100, sizeof (char));
  name=declare_char (A->nseq, 100);
  for (a=0; a<A->nseq; a++)
    {
      char *p=strstr (A->name[a], "_");
      if (!p)
	{
	  fprintf ( stderr, "\nWARNING: Seq %s could not be included.", A->name[a]);
	}
      p+=1;
      if ( name_is_in_list (p, name,nname, 100)==-1)
	{
	  sprintf ( name[nname++], "%s", p);
	}
    }

  C=declare_aln2 (nname, (A->len_aln*S->nseq)+1);
  free_char (C->name,-1); C->name=name;
  C->nseq=nname;
  C->col_weight=(float*)vcalloc ( A->len_aln*S->nseq, sizeof(float));

  C->len_aln=0;
    for (a=0; a<S->nseq; a++)
    {
      for (b=0; b<C->nseq; b++)
	{
	  sprintf (cname, "%s_%s", S->name[a],C->name[b]);
	  if ((i=name_is_in_list (cname, A->name, A->nseq, 100))==-1)
	    {
	      char *s=generate_null (A->len_aln);
	      strcat (C->seq_al[b], s);
	      vfree (s);
	    }
	  else
	    strcat (C->seq_al[b], A->seq_al[i]);
	}
      for (c=C->len_aln, b=0;b<A->len_aln;b++, c++)
	{
	  C->col_weight[c]=(S->W)->SEQ_W[a];
	}
      C->len_aln+=A->len_aln;
    }
  return C;
}


Alignment * concatenate_aln ( Alignment *A1, Alignment *A2, char *spacer)
{
  Alignment *A;
  int a, i;

  A=declare_aln2( A1->nseq+A2->nseq , A1->len_aln+A2->len_aln+1);
  for ( a=0; a< A1->nseq; a++)
    {
      if ((i=name_is_in_list ( A1->name[a], A2->name, A2->nseq, 100))!=-1)
	{
	  sprintf ( A->name[A->nseq], "%s", A1->name[a]);
	  sprintf (A->seq_al[A->nseq], "%s%s%s", A1->seq_al[a],(spacer)?spacer:"", A2->seq_al[i]);
	  A->nseq++;
	}
      else
	{
	  char *buf;
	  buf=generate_string (A2->len_aln, '-');
	  sprintf ( A->name[A->nseq], "%s", A1->name[a]);
	  sprintf (A->seq_al[A->nseq], "%s%s",  A1->seq_al[a], buf);
	  A->nseq++;
	  vfree (buf);
	}
    }
  for ( a=0; a< A2->nseq; a++)
    {
      if ((i=name_is_in_list ( A2->name[a], A1->name, A1->nseq, 100))==-1)
	{
	  char *buf;
	  buf=generate_string (A1->len_aln, '-');
	  sprintf ( A->name[A->nseq], "%s", A2->name[a]);
	  sprintf (A->seq_al[A->nseq], "%s%s",  buf, A2->seq_al[a]);
	  A->nseq++;
	  vfree (buf);
	}
    }
  A->len_aln=A1->len_aln+A2->len_aln;
  return A;
}
Alignment * aln_cat ( Alignment *A, Alignment *B)
     {
     int a;

     if ( A->nseq!=B->nseq)
	 {
	 fprintf ( stderr, "\nERROR IN ALN CAT: DIFFERENT NSEQ\n");
	 myexit(EXIT_FAILURE);
	 }

     A=realloc_alignment2(A, A->nseq,A->len_aln+B->len_aln+1);

     for ( a=0;a< A->nseq; a++)
         {
	 strcat ( A->seq_al[a], B->seq_al[a]);
	 }
     A->len_aln+=B->len_aln;
     return A;
     }
int verify_aln ( Alignment *A, Sequence *S, char *message)
     {
     int a, b, c,s,r;


     for ( a=0;a< A->nseq; a++)
         {
	 s=A->order[a][0];
	 r=A->order[a][1];
	 for ( b=0, c=0; b< A->len_aln; b++)
	     {
	     if ( !is_gap(A->seq_al[a][b]))
		  {
		  if (tolower(A->seq_al[a][b])!=tolower(S->seq[s][c+r]))
		      {
		      fprintf ( stderr, "\n%s\nResidue [%c %d, %c %d] line %d seq %d",message,A->seq_al[a][b], b,S->seq[s][c+r], c+r,a,s);
		      output_Alignment_with_res_number(A, stderr);
		      myexit(EXIT_FAILURE);
		      return 0;
		      }
		  c++;
		  }
	     }
	 }
     return 1;
     }

Alignment *adjust_est_aln ( Alignment *PW, Alignment *M, int s)
{
  /*This function reajusts M, threading M onto PW
   two seqences in PW
   s+1 seq in M

   seq 0 PW ----> 0->s-1 in M
   seq 1 PW ----> 1->s   in M;

   */
  int a, b;
  static char **array;


  int top_M=0;
  int bottom_M=0;


  if ( array==NULL)
    {
      array=declare_char (500, 100000);
    }

  for ( a=0; a< PW->len_aln; a++)
    {
      if ( is_gap(PW->seq_al[0][a]))
	{
	  for ( b=0; b< s; b++)
	    array[b][a]='-';
	}
      else
	{
	  for ( b=0; b< s; b++)
	    array[b][a]=M->seq_al[b][top_M];
	top_M++;
	}

      if ( is_gap(PW->seq_al[1][a]))
	{
	  array[s][a]='-';
	}
      else
	{

	  array[s][a]=M->seq_al[s][bottom_M];
	bottom_M++;
	}
    }

  M->len_aln=PW->len_aln;
  for (a=0; a<s; a++)
    {
      for (b=0; b<PW->len_aln; b++)
	M->seq_al[a][b]=array[a][b];
      M->seq_al[a][b]='\0';
    }


  M->nseq=s+1;

  return M;
}


Alignment * rename_seq_in_aln (Alignment *A, char ***list)
{
  int n, i;
  if ( !A)return A;



  n=0;
  while ( list[n][0][0])
    {
      if ( (i=name_is_in_list (list[n][0], A->name, A->nseq, 100))!=-1)
	{
	  sprintf ( A->name[i], "%s", list[n][1]);
	}
      n++;
    }

  A->S=rename_seq_in_seq (A->S, list);
  return A;
}
Sequence * rename_seq_in_seq (Sequence *A, char ***list)
{
  int n, i;
  if ( !A || !list)return A;

  n=0;
  while ( list[n][0][0])
    {
      if ( (i=name_is_in_list (list[n][0], A->name, A->nseq, 100))!=-1)
	{
	  sprintf ( A->name[i], "%s", list[n][1]);
	}
      n++;
    }
  return A;
}
/********************************************************************/
/*                                                                  */
/*                   FLOAT SIMILARITIES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
float get_seq_fsim ( char *string1, char *string2, char *ignore, char *similarity_set,int **matrix, int MODE )
	{
	int len, a, r1, r2, nr1=0, nr2=0;
	float pos=0, sim=0;


	len=MIN((strlen (string1)),(strlen (string2)));
	if ( len==0)return 0;

	for ( a=0; a< len; a++)
		{

		  r1=string1[a];
		  r2=string2[a];
		  nr1+=!is_gap(r1);
		  nr2+=!is_gap(r2);

		  if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
			if ( matrix)sim+=matrix[r1-'A'][r2-'A'];
			else if (is_in_same_group_aa(r1,r2,0, NULL,similarity_set))
				{
				sim++;
				}
			}
		}
	if ( MODE==UNGAPED_POSITIONS)return ( sim*100)/pos;
	else if ( MODE==ALIGNED_POSITIONS)return (sim*100)/len;
	else if ( MODE==AVERAGE_POSITIONS)return (sim*200)/(nr1+nr2);
	else
	  {
	    return 0;
	  }

	}
float get_seq_fsim2 ( char *string1, char *string2, char *ignore, char *in_mode)
	{
	int len1;
	int a;
	int p1, p2;
	int r1=0,r2=0;
	char *p;
	char mode[1000];
	float r=0, pos1, pos2, pos0, gap, sim;


	sprintf ( mode, "%s", in_mode);

	/*mode: <mat>__<sim_mode>
	  mat: idscore to get the alignment done
	       any legal cw matrix
	  sim_mode: sim1->identities/matches
                    sim2->identities/min len
	*/


	if ( (p=strstr (mode, "_"))!=NULL)
	  {
	    p[0]='\0';
	    p++;
	  }

	if (strstr (mode, "idscore"))
	  {
	    static int **mat;
	    if (!mat) mat=read_matrice ("blosum62mt");
	    return idscore_pairseq (string1, string2, -12, -1, mat,mode);

	  }

	len1=strlen (string1);
	for ( sim=pos1=pos2=pos0=gap=0,a=0; a< len1; a++)
		{
		  r1=string1[a];
		  r2=string2[a];
		  p1=1-is_in_set (r1, ignore);
		  p2=1-is_in_set (r2, ignore);
		  pos1+=p1; pos2+=p2;
		  if (p1 && p2)
			{
			  pos0++;
			  if (is_in_same_group_aa(r1,r2,0, NULL,p))
			    {
			      sim++;
			    }
			}
		  else if (p1+p2==1)
		    {
		      gap++;
		    }
		}

	if ( p==NULL || strm (p, "sim1") || strm (p, "sim"))
	  {
	    r=(pos0==0)?0:(sim*MAXID)/pos0;
	  }
	else if ( strm (p, "sim2"))
	  {
	    r=(pos1==0 || pos2==0)?0:(sim*MAXID)/MIN(pos1,pos2);
	  }
	else if ( strm (p, "sim3"))
	  {
	    r=(pos1==0 || pos2==0)?0:(sim*MAXID)/MAX(pos1,pos2);
	  }
	else if ( strm (p, "gap1"))
	  {
	    r=(len1==0)?MAXID:(gap*MAXID)/len1;
	    r=MAXID-r;
	  }
	else if ( strm (p, "logid"))
	  {
	    r=logid_score (pos0, sim);
	  }
	else
	  {
	  r=(pos0==0)?0:(sim*MAXID)/pos0;
	  }

	return r;

	}

/********************************************************************/
/*                                                                  */
/*                   ALIGNMENT ANALYSES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
int **dist_array2sim_array ( int **p, int max)
{
  int s1, s2, a, b;
  s1=read_array_size ((void *)p, sizeof (void *));
  s2=read_array_size ((void*)p[0],sizeof (int));
  /*  s2=read_array_size ((void*)p[0],sizeof (void *)); OLD before 64 BITS*/
  for ( a=0; a< s1; a++)
    for ( b=0; b< s2; b++)
      {
	p[a][b]=max-p[a][b];
      }
  return p;
}

int **sim_array2dist_array ( int **p, int max)
{
  int s1, s2, a, b;
  s1=read_array_size ((void *)p, sizeof (void *));
  s2=read_array_size ((void*)p[0],sizeof (int));

  /*s2=read_array_size ((void*)p[0],sizeof (void *)); OLD before 64 Bits stuff*/
 for ( a=0; a< s1; a++)
    for ( b=0; b< s2; b++)
      {
	p[a][b]=max-(int)p[a][b];
      }
  return p;
}

int **normalize_array (int **p, int max, int norm)
{
int s1, s2, a, b;
 s1=read_array_size ((void *)p, sizeof (void *));
 s2=read_array_size ((void*)p[0],sizeof (int));

 /*s2=read_array_size ((void*)p[0],sizeof (void *)); OLD before 64 Bits stuff*/
 for ( a=0; a< s1; a++)
   for ( b=0; b< s2; b++)
     {
       p[a][b]=(p[a][b]*norm)/max;
      }
 return p;
}

int aln2most_similar_sequence ( Alignment *A, char *mode)
{
  int **w;
  int a, b;
  int avg, best_avg=0, best_seq=0;
  char *buf;
  int coverage;


  if ( !A) return -1;
  else if ( A->nseq==1)return 0;
  else
    {
      buf=(char*)vcalloc ( A->len_aln+1, sizeof (char));
      w=get_sim_aln_array ( A, mode);

      for ( a=0; a< A->nseq; a++)
	{
	  sprintf ( buf, "%s", A->seq_al[a]);
	  ungap(buf);
	  coverage=(strlen(buf)*MAXID)/A->len_aln;

	  for ( avg=0,b=0; b< A->nseq; b++)avg+=w[a][b]*coverage;
	  if ( avg>best_avg){best_avg=avg; best_seq=a;}
	}
      free_int (w, -1);
      vfree (buf);
      return best_seq;
    }

}

int aln2coverage ( Alignment *A, int ref_seq)
{
  int a,b;
  int cov_pos=0, npos=0;


  for ( a=0; a< A->len_aln; a++)
    {
      if ( !is_gap ( A->seq_al[ref_seq][a]))
	{
	  npos++;
	  for ( b=0; b< A->nseq; b++)
	    {
	      if ( b!=ref_seq && !is_gap ( A->seq_al[b][a])){cov_pos++;break;}
	    }
	}
    }

  return  (int) (npos==0)?0:(( MAXID*cov_pos)/npos);
}


int sub_aln2sim ( Alignment *A, int *ns, int **ls, char *mode)
{
  int a, b, n;
  float avg;

  n=0; avg=0;
  if (!A || (ns==NULL && A->nseq<2))return -1;
  else if (ns==NULL)
    {
      for (a=0; a< A->nseq-1; a++)
	for ( b=a+1; b< A->nseq;b++, n++)
	  avg+=generic_get_seq_sim (A->seq_al[a], A->seq_al[b], NULL, mode);
    }
  else
    {
      for (a=0; a<ns[0]; a++)
	for (b=0; b< ns[1]; b++, n++)
	  {
	    avg+=generic_get_seq_sim (A->seq_al[ls[0][a]], A->seq_al[ls[1][b]], NULL, mode);
	  }
    }
  return (int)(n==0)?0:((float)avg/(float)n);
}
int sub_aln2max_sim ( Alignment *A, int *ns, int **ls, char *mode)
{
  int a, b, n;
  float avg;

  n=0; avg=0;
  if (!A || (ns==NULL && A->nseq<2))return -1;
  else if (ns==NULL)
    {
      for (a=0; a< A->nseq-1; a++)
	for ( b=a+1; b< A->nseq;b++, n++)
	  avg=MAX(avg,generic_get_seq_sim (A->seq_al[a], A->seq_al[b], NULL, mode));
    }
  else
    {
      for (a=0; a<ns[0]; a++)
	for (b=0; b< ns[1]; b++, n++)
	  {
	    avg=MAX(avg,generic_get_seq_sim (A->seq_al[ls[0][a]], A->seq_al[ls[1][b]], NULL, mode));
	  }
    }
  return avg;
}
double* aln2column_normalized_entropy (Alignment *A)
{
  int a;
  double max=0;
  double *e=aln2column_entropy (A);
  for (a=0; a<A->len_aln; a++)
    if (e[a]>max) max=e[a];
  for (a=0; a<A->len_aln; a++)
    {
      e[a]=(double)100-((double)100*(e[a]/max));
    }
  return e;
}
double* aln2column_entropy (Alignment *A)
{
  double*e, *p;
  double n;
  int a, b;
  
  if (!A || A->len_aln==0) return NULL;

  e=(double*)vcalloc (A->len_aln, sizeof (double));
  p=(double*)vcalloc (256, sizeof (double));


  for (a=0; a<A->len_aln; a++)
    {
      n=0;
      for (b=0; b<A->nseq; b++)
	{
	  char r=A->seq_al[b][a];
	  if ( !is_gap(r))
	    {
	      r=tolower(r);
	      p[r]++;
	      n++;
	    }
	}
      for (b=0; b<256; b++)
	{
	  if (n<0.00000001 || p[b]<0.0000001)continue;
	  p[b]/=n;
	  e[a]+=p[b]*(log(p[b])/log(2));
	  p[b]=0;
	}
      e[a]*=(double)-1;
      if (e[a]<0.0000001)e[a]=0;
    }
  vfree (p);
  return e;
}

double aln2entropy (Alignment *A, int *in_ls, int in_ns, float gap_threshold)
{
  int ns, a, s, col, r,ncol;
  int *ls;
  double *count;
  double entropy=0;
  float ng;

  ls=(int*)vcalloc ( A->nseq, sizeof (int));
  count=(double*)vcalloc ( 26, sizeof (double));


  if ( in_ls)
    {
      ns=in_ns;
      for ( a=0; a< ns; a++)ls[a]=in_ls[a];
    }
  else
    {
      ns=A->nseq;
      for ( a=0; a< ns; a++)ls[a]=a;
    }

  if ( ns==0)
    {
      vfree(ls);vfree(count);return 0;
    }
  for (ncol=0,col=0; col<A->len_aln; col++)
    {
      for (ng=0,a=0; a< ns; a++)
	{
	  s=ls[a];
	  ng+=is_gap(A->seq_al[s][col]);
	}
      ng/=ns;
      if ( ng>gap_threshold)continue;

      ncol++;

      for ( a=0; a<ns; a++)
	{
	  s=ls[a];
	  r=tolower(A->seq_al[s][col]);
	  if (!is_gap(r))count[r-'a']++;
	}
      for (a=0; a<26; a++)
	{
	  if ( count[a]==0);
	  else
	    {
	      count[a]/=(double)ns;

	      entropy+=count[a]*log(count[a]);
	      count[a]=0;
	    }
	}
    }
  entropy/=-ncol;
  vfree (ls); vfree(count);

  return entropy;
}

int aln2sim2 (Alignment *A)
{
	int a, b, c;
	double *score;
	double tscore=0;
	score=(double*)vcalloc ( 256, sizeof (double));

	for (a =0; a<A->len_aln; a++)
	{
		for (b=0; b<A->nseq; b++)
		{
			c=tolower(A->seq_al[b][a]);
			if (c!='-')score[c]++;
		}
		for (b=0; b<256; b++)
		{
			tscore+=(score[b]*score[b]);
			score[b]=0;
		}
	}
	tscore/=A->nseq;
	vfree (score);
	return (int)tscore;
}


int aln2sim ( Alignment *A, char *mode)
{
  return sub_aln2sim ( A, NULL, NULL, mode);
  /*
  if ( !A || A->nseq<2) return -1;
  w=get_sim_aln_array ( A, mode);

  for (c=0, a=0; a< A->nseq-1; a++)
    for ( b=a+1; b< A->nseq; b++, c++)
      {
	avg+=(float)w[a][b];
      }
  free_int (w, -1);
  return (int)((float)avg/(float)c);
  */
}

int aln_is_aligned ( Alignment *A)
{
  int a, b;

  if ( !A)return 0;
  for (a=0; a< A->nseq; a++)
    for ( b=A->len_aln-1; b>0; b--)
      {
	if (!is_gap(A->seq_al[a][b]) && is_gap(A->seq_al[a][b-1]))return 1;
      }
  return 0;
}


int seq2aln2sim_old ( char *seq1, char *seq2, char *mode_aln, char *mode_id)
{
  Alignment *A;
  int sim;

  A=align_two_sequences (seq1, seq2, "pam250mt", -10, -1, mode_aln);
  sim=aln2sim (A, mode_id);
  free_aln (A);
  return sim;
}
int seq2aln2sim ( char *seq1, char *seq2, char *mode_aln, char *mode_id)
{
  Alignment *A;
  int sim;
  static int gop;

  if (!gop)
    {
      int **m;
      m=read_matrice ("blosum62mt");
      gop=get_avg_matrix_mm(m, AA_ALPHABET)*10;
      free_int (m, -1);
    }

  A=align_two_sequences (seq1, seq2, "blosum62mt",gop,-1, mode_aln);
  sim=aln2sim (A, mode_id);
  free_aln (A);
  return sim;
}
int* get_cdna_seq_winsim ( int *cache, char *string1, char *string2, char *ignore, char *mode,int *w )
	{
	int len1, len2;
	int a, x;


	len1=strlen (string1);
	len2=strlen (string2);

	if ( len1!=len2)
	  {
	    fatal_exit( stderr,EXIT_FAILURE, "\nTHE TWO cDNAs DO NOT HAVE THE SAME LENGTH [FATAL:get_cdna_seq_sim:%s", PROGRAM);
	  }

	x=get_cdna_seq_sim(cache, string1, string2, ignore, "");
	for ( a=0; a< len1; a++)
	  w[a]=x;

	add_warning (stderr, "winsim not implemented for cDNA");
	return w;
	}

int get_cdna_seq_sim ( int *cache, char *string1, char *string2, char *ignore, char *mode)
	{
	int len1;
	int len2;
	int a;
	int pos=0;
	int sim=0;
	char r1=0, r2=0;

	len1=strlen (string1);
	len2=strlen (string2);



	if ( len1!=len2)
	  {
	    fprintf ( stderr, "\nTHE TWO cDNAs DO NOT HAVE THE SAME LENGTH [FATAL:get_cdna_seq_sim:%s", PROGRAM);
	    crash("");
	  }

	for ( a=0; a< len1;)
		{

		  if ( cache[a]==0){a++;continue;}
		  else if ( cache[a]==1)
		    {

		      r1=translate_dna_codon (string1+a, 'x');
		      r2=translate_dna_codon (string2+a, 'x');
		      a+=3;
		    }

		  if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
			if (is_in_same_group_aa(r1,r2,0, NULL,mode+4))
				{
				sim++;
				}
			}
		}



	if (pos==0)
		 return 0;
	else
		return (int) (sim*MAXID)/pos;

	}

int* get_seq_winsim ( char *string1, char *string2, char *ignore, char *mode, int*w)
	{
	int len1, len2, len;
	int left, right;
	int a,b;
	int sim=0;
	int window;
	int r1, r2;

	len1=strlen (string1);
	len2=strlen (string2);
	window=atoi(mode);
	len=2*window+1;

	if ( len1!=len2)return 0;
	if (window==0 || (window*2+1)>=len1)
	  {
	    sim=get_seq_sim (string1, string2, ignore, "");
	    for (a=0; a<len1; a++)w[a]=sim;
	    return w;
	  }


	for ( a=0; a< len1; a++)
		{

		  left =MAX(0, a-window);
		  right=MIN(len1, left+len);
		  for (sim=0,b=left; b<right; b++)
		    {
		      r1=string1[b];
		      r2=string2[b];
		      if (  !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			  if (r1==r2)sim++;
			}
		    }
		  w[a]=(sim*MAXID)/len;
		}
	return w;
	}


int get_seq_sim ( char *string1, char *string2, char *ignore, char *in_mode)
	{
	int len1;
	int a;
	int pos1, pos2, pos0,gap=0, sim;
	int p1, p2;
	int r=0,r1=0,r2=0;
	char *p;
	static char *mode;

	if (!mode)mode=(char*)vcalloc (100, sizeof (char));
	else mode[0]='\0';
	if (in_mode)
	  {
	    while (in_mode[0]=='_')in_mode++;
	    sprintf ( mode, "%s", in_mode);
	  }

	/*mode: <mat>__<sim_mode>
	  mat: idscore to get the alignment done
	       any legal cw matrix
	  sim_mode: sim1->identities/matches
                    sim2->identities/min len
	*/


	if ( (p=strstr (mode, "_"))!=NULL)
	  {
	    p[0]='\0';
	    p++;
	  }


	if (strstr (mode, "idscore"))
	  {
	    static int **mat;
	    if (!mat) mat=read_matrice ("blosum62mt");
	    return idscore_pairseq (string1, string2, -12, -1, mat,mode);

	  }
	len1=strlen (string1);
	for ( sim=pos1=pos2=pos0=0,a=0; a< len1; a++)
		{
		  r1=string1[a];
		  r2=string2[a];
		  p1=1-is_in_set (r1, ignore);
		  p2=1-is_in_set (r2, ignore);

		  pos1+=p1; pos2+=p2;
		  if (p1 && p2)
			{
			  pos0++;
			  if (is_in_same_group_aa(r1,r2,0, NULL, mode))
			    {
			      sim++;
			    }
			}
		  else if (p1+p2==1)
		    {
		      gap++;
		    }
		}

	if ( strstr (mode, "cov"))
	  {
	    r=(pos0+gap==0)?0:(pos0*MAXID)/(pos0+gap);
	  }
	else if ( p==NULL || strm (p, "sim1") || strm (p, "sim"))
	  {
	    r=(pos0==0)?0:(sim*MAXID)/pos0;
	  }
	else if ( strm (p, "sim2"))
	  {
	    r=(pos1==0 || pos2==0)?0:(sim*MAXID)/MIN(pos1,pos2);
	  }
	else if ( strm (p, "sim3"))
	  {
	    r=(pos1==0 || pos2==0)?0:(sim*MAXID)/MAX(pos1,pos2);
	  }
	else if ( strm (p, "gap1"))
	  {
	    r=(len1==0)?MAXID:(gap*MAXID)/len1;
	    r=MAXID-r;
	  }
	else if ( strm (p, "logid"))
	  {
	    r=logid_score (pos0, sim);
	  }
	else if ( strstr (mode, "sim"))
	  {
	    r=(pos0==0)?0:(sim*MAXID)/pos0;
	  }


	return r;

	}
int get_seq_sim_2 ( char *string1, char *string2, char *ignore, char **gr, int ng)
	{
	int len1;
	int len2;
	int a;
	int pos=0;
	int sim=0;
	char r1, r2;


	len1=strlen (string1);
	len2=strlen (string2);

	if ( len1!=len2)return 0;

	for ( a=0; a< len1; a++)
		{
		r1=string1[a];
		r2=string2[a];
		if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
		    	if (is_in_same_group_aa(r1,r2,ng, gr, NULL)==1)
				{
				  sim++;
				}
			}
		}

	if (pos==0)
		 return 0;
	else
		return (int) (sim*MAXID)/pos;

	}

int get_seq_sim_3 ( char *string1, char *string2, char *ignore, int **mat)
	{
	int len1;
	int len2;
	int a;

	int sim=0;
	char r1, r2;


	len1=strlen (string1);
	len2=strlen (string2);

	if ( len1!=len2)return 0;

	for ( a=0; a< len1; a++)
		{
		r1=string1[a];
		r2=string2[a];
		if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			sim+=mat[r1-'A'][r2-'A'];
			}
		}
	return sim;

	}
int * get_aln_col_weight ( Alignment *A, char *mode)
	{
	int a, b;
	char *col;
	int *weight;

	col=(char*)vcalloc ( A->nseq, sizeof (int));
	weight=(int*)vcalloc (A->len_aln, sizeof (int));

	for (a=0; a< A->len_aln; a++)
		{
		for ( b=0; b< A->nseq; b++)
			col[b]=A->seq_al[b][a];
		weight[a]=(find_group_aa_distribution (col, A->nseq,0,NULL,NULL, mode )*MAXID)/A->nseq;
		}
	vfree (col);
	return weight;

	}

int analyse_aln_column ( Alignment *B, int col)
    {

    char r=' ';
    int a, b, c=0;
    static char *mat;
    static int ng_cw_star;
    static char **cw_star;
    int *cw_star_count;

    static int ng_cw_col;
    static char **cw_col;
    int *cw_col_count;

    static int ng_cw_dot;
    static char **cw_dot;
    int *cw_dot_count;






    if ( !B->S || !(B->S)->type)B= get_aln_type (B);

    if ( !mat)mat=(char*)vcalloc ( STRING, sizeof (char));

    if ( !ng_cw_star)
       {
	   cw_star=make_group_aa ( &ng_cw_star, strcpy ( mat,"idmat"));
	   cw_col=make_group_aa ( &ng_cw_col, strcpy (mat,"clustalw_col"));
	   cw_dot=make_group_aa ( &ng_cw_dot, strcpy (mat, "clustalw_dot"));
       }

    cw_star_count=(int*)vcalloc (ng_cw_star, sizeof (int));
    cw_col_count=(int*)vcalloc ( ng_cw_col, sizeof (int));
    cw_dot_count=(int*)vcalloc (ng_cw_dot, sizeof (int));

    for ( a=0; a< B->nseq; a++)
        {
	c=tolower (B->seq_al[a][col]);
	if (is_gap(c)){r=' ';break;}

	for ( b=0; b< ng_cw_star; b++)
	    cw_star_count[b]+=is_in_set (c, cw_star[b]);
	for ( b=0; b< ng_cw_col; b++)
	    cw_col_count[b]+=is_in_set  (c, cw_col[b]);
	for ( b=0; b< ng_cw_dot; b++)
	    cw_dot_count[b]+=is_in_set  (c, cw_dot[b]);
	}





    if ( !is_gap(c) && r==' ')
	for ( b=0; b< ng_cw_star; b++)if ( cw_star_count[b]==B->nseq){r='*'; break;}
    if ( !is_gap(c) && r==' ' && !(strm((B->S)->type, "DNA")||strm ((B->S)->type,"RNA")))
	for ( b=0; b< ng_cw_col ; b++)if ( cw_col_count [b]==B->nseq){r=':'; break;}
    if ( !is_gap(c) && r==' ' && !(strm((B->S)->type, "DNA")||strm ((B->S)->type,"RNA")))
	for ( b=0; b< ng_cw_dot ; b++)if ( cw_dot_count [b]==B->nseq){r='.'; break;}



    vfree(cw_star_count);
    vfree(cw_col_count);
    vfree(cw_dot_count);

    return r;
    }


int ** get_cov_aln_array ( Alignment *A, char *mode)
{
  int **w;
  int a, b, c, t;

  w=declare_int ( A->nseq, A->nseq);


  for ( a=0; a< A->nseq-1; a++)
    {
      w[a][a]=100;
      for ( t=0,b=a+1; b< A->nseq; b++)
	{
	  for ( c=0; c< A->len_aln; c++)
	    {
	      t+=(!is_gap(A->seq_al[a][c]) &&!is_gap(A->seq_al[b][c]));
	    }
	  w[a][b]=w[b][a]=(t*100)/A->len_aln;
	}
    }
  return w;
}

int ** get_cov_master_aln_array ( Alignment *A,int n, char *mode)
{
  int **w;
  int    b, c, t;

  w=declare_int ( A->nseq, A->nseq);


  for (b=0; b< A->nseq; b++)
	{

	  for (t=0, c=0; c< A->len_aln; c++)
	    {
	      t+=(!is_gap(A->seq_al[n][c]) &&!is_gap(A->seq_al[n][c]));
	    }
	  w[n][b]=w[b][n]=(t*100)/A->len_aln;
	}

  return w;
}
int ** get_sim_master_aln_array ( Alignment *A,int n, char *mode)
	{
	int **w;
	int a;

	w=declare_int ( A->nseq, A->nseq);


	for ( a=0; a< A->nseq; a++)
	  {
	    if ( strm (mode, "cdna"))
	    w[n][a]=w[a][n]=get_cdna_seq_sim ( A->cdna_cache[0], A->seq_al[a], A->seq_al[n],GAP_LIST, mode);
	  else
	    w[n][a]=w[a][n]=get_seq_sim ( A->seq_al[n], A->seq_al[a],GAP_LIST, mode);
	  }
	return w;
	}
int ** get_dist_aln_array ( Alignment *A, char *mode)
{

  int **w;

  w=get_sim_aln_array ( A, mode);
  return sim_array2dist_array(w,MAXID);
}
Sequence * seq2filter (Sequence *Sin, int min, int max)
{
  int *keep;
  char *tmpfile;
  Sequence *S, *Sout;
  int a, b, sim;
  int **M;
  FILE *fp;
  int n;

  S=duplicate_sequence (Sin);
  for (a=0; a<S->nseq; a++)ungap(S->seq[a]);
  keep=(int*)vcalloc (S->nseq, sizeof (int));
  M=read_matrice ("blossum62mt");
  for (a=0; a<S->nseq; a++)
    {
      output_completion ( stderr, a, S->nseq, 100, "Distance Matrix Computation: ");
      for ( b=a+1; b<S->nseq; b++)
	{

	  sim=idscore_pairseq(S->seq[a], S->seq[b],-10, -2,M, "sim");
	  if ( sim>min && sim<max)keep[a]=keep[b]=1;
	  fprintf ( stderr, "\nSim %d Min %d Max %d", sim, min, max);
	}
    }

  tmpfile=vtmpnam (NULL);
  fp=vfopen (tmpfile, "w");
  for (n=0,a=0; a< S->nseq; a++)
    if ( keep[a])
      {
      fprintf ( fp, ">%s %s\n%s", S->name[a], S->seq_comment[a], S->seq[a]);
      n++;
      }
  vfclose (fp);
  if (n==0) return NULL;
  Sout=main_read_seq(tmpfile);
  free_int (M, -1); vfree (keep); free_sequence (S, -1);
  return Sout;
}

Alignment * grep_seq (Alignment *S,char *field, char *mode, char *string)
{
  int a;
  FILE *fp;
  char *tmp;
  int n=0;

  tmp=vtmpnam (NULL);
  fp=vfopen (tmp, "w");

  if ( !strm(mode, "KEEP") && ! strm (mode, "REMOVE"))
    {
      printf_exit(EXIT_FAILURE,stderr, "+grep <field> <KEEP|REMOVE> <string> [FATAL: %s]", PROGRAM);
    }
  else if ( !strm(field, "SEQ") && ! strm (field, "COMMENT") && ! strm(field, "NAME"))
    {
      printf_exit(EXIT_FAILURE,stderr,"ERROR: +grep <NAME|COMMENT|SEQ> <mode> <string> [FATAL: %s]", PROGRAM);
    }


  for (n=0, a=0; a< S->nseq; a++)
    {
      int found=0;

      if (strm(field, "NAME") && perl_strstr (S->name[a], string))found=1;
      else if (strm(field, "COMMENT") && S->seq_comment[a][0] && perl_strstr (S->seq_comment[a], string) )found=1;
      else if (strm(field, "SEQ") && perl_strstr (S->seq_al[a], string))found=1;

      if ( (strm (mode, "KEEP") && found) || (strm (mode, "REMOVE") && !found))
	{
	  n++;
	  fprintf (fp, ">%s", S->name[a]);
	  if (S->seq_comment[a][0])fprintf (fp, " %s", S->seq_comment[a]);
	  fprintf (fp, "\n%s\n", S->seq_al[a]);
	}
    }

  vfclose (fp);

  free_aln (S);
  if ( n==0) return NULL;
  else
    return main_read_aln (tmp, NULL);
}

Alignment * modify_seq (Alignment *S, char *field, char *string1, char *string2)
{
  int a;
  FILE *fp;
  char *tmp;

  tmp=vtmpnam (NULL);
  fp=vfopen (tmp, "w");
  for ( a=0; a< S->nseq; a++)
    {
      if (strm(field, "NAME"))S->name[a]=substitute ( S->name[a], string1, string2);
      else if (strm(field, "COMMENT"))S->seq_comment[a]=substitute ( S->seq_comment[a], string1, string2);
      else if (strm(field, "SEQ"))S->seq_al[a]=substitute ( S->seq_al[a], string1, string2);
      fprintf (fp, ">%s", S->name[a]);
      if (S->aln_comment[a][0])fprintf (fp, " %s", S->aln_comment[a]);
      fprintf (fp, "\n%s\n", S->seq_al[a]);
    }
  vfclose (fp);
  free_aln (S);
  S=main_read_aln (tmp, NULL);
  return S;
}

int ** seq2sim_mat (Sequence *S, char *mode)
{
  return seq2comp_mat ( S,mode, "sim");
}
int ** seq2cov_mat (Sequence *S, char *mode)
{
  return seq2comp_mat ( S,mode, "cov");
}

int ** seq2comp_mat (Sequence *S, char *mode, char *comp_mode)
{
  int a, b;
  int **sim;
  char file[1000];
  Alignment *A;
  char *name;


  /*Use pre_computed value if available in the current dir*/

  name=path2filename(S->file[0]);
  sprintf ( file, "%s%s.%s.%s_file", get_cache_dir(),name, mode, comp_mode);
  A=seq2aln(S,NULL, RM_GAP);
  if ( check_file_exists (file) && is_distance_matrix_file (file) && (sim=input_similarities(file, A, NULL))!=NULL)
    {
      display_input_filename (stderr, "SIMILARITY_MATRIX", "SIMILARITY_MATRIX_FORMAT_01", file, CHECK);
      fprintf ( stderr, "\n");
    }
  else
    {
      char mode2[1000];
      int **M;

      M=read_matrice (mode);
      sim=declare_int ( S->nseq, S->nseq);
      for ( a=0; a< S->nseq; a++)
	{
	  ungap (S->seq[a]);
	  sim[a][a]=100;
	}

      for ( a=0; a<S->nseq-1; a++)
	{

	  output_completion4halfmat ( stderr, a, S->nseq, 100, "Similarity Matrix Computation: ");
	  for ( b=a+1; b< S->nseq; b++)
	    {
	      sim[a][b]=sim[b][a]=idscore_pairseq(S->seq[a], S->seq[b],-12, -1,M, comp_mode);
	    }
	}
      free_int (M,-1);
      sprintf ( mode2, "_memory_%ld", (long int)sim);
      output_similarities( file, A, mode2);
      display_output_filename (stderr, "SIMILARITY_MATRIX", "SIMILARITY_MATRIX_FORMAT_01", file, CHECK);
      fprintf ( stderr, "\n");
    }
  free_aln (A);
  return sim;
}

int ** fast_aln2sim_list (Alignment *A,  char *mode, int *ns, int **ls)
{
  int **simm;
  int p1, p2, p3, r1, r2;
  int gap,pos0,pos1,pos2,len,sim;
  int a, b, c, m, s=0,s1, s2, n;
  int free_ns=0;

  if (ns==NULL)
    {
      free_ns=1;
      ns=(int*)vcalloc (2, sizeof (int));
      ns[0]=ns[1]=A->nseq;
      ls=declare_int (2, A->nseq);
      for ( a=0; a< 2; a++)
	for (b=0; b<A->nseq; b++)
	  ls[a][b]=b;
    }


  simm=declare_int (ns[0]*ns[1]+1, 3);

  if (strstr (mode, "sim1"))m=0;
  else if (strstr (mode, "sim2"))m=1;
  else if (strstr (mode, "sim3"))m=2;
  else if (strstr (mode, "gap1"))m=3;
  else if (strstr (mode, "cov1"))m=4;
  else if (strstr (mode, "logid"))m=5;
  else m=0;



  for (n=0,a=0; a<ns[0]; a++)
    {
      s1=ls[0][a];
      for ( b=0; b<ns[1]; b++, n++)
	{
	  s2=ls[1][b];
	  gap=pos0=pos1=pos2=len=sim=0;

	  for ( c=0; c< A->len_aln; c++)
	    {
	      r1=tolower (A->seq_al[s1][c]);
	      r2=tolower (A->seq_al[s2][c]);
	      p1=(r1!='-')?1:0;
	      p2=(r2!='-')?1:0;
	      p3=p1+p2;
	      if ( p3==0)continue;
	      if ( p3==1)gap++;
	      if ( r1==r2)sim++;
	      pos1+=p1;
	      pos2+=p2;
	      pos0+=(p3==2)?1:0;
	      len++;
	    }

	  if (m==0)s=(pos0==0)?0:(sim*MAXID)/pos0; //sim1
	  else if (m==1)  s=(MIN(pos1,pos2)==0)?0:(sim*MAXID)/MIN(pos1,pos2);//sim2
	  else if (m==2)  s=(MAX(pos1,pos2)==0)?0:(sim*MAXID)/MAX(pos1,pos2);//sim3
	  else if (m==3)  s=(len==0) ?0:((len-gap)*MAXID)/len;//gap1
	  else if (m==4)  s=(len==0) ?0:((pos0)*MAXID)/len; //cov
	  else if (m==5)
	    {
	      s=logid_score ( sim, len);
	    }
	  simm[n][0]=s1;
	  simm[n][1]=s2;
	  simm[n][2]=s;
	}
    }

  if ( free_ns) {vfree(ns); free_int (ls, -1);}
  simm[n][0]=-1;
  return simm;
}

int ** fast_aln2sim_mat (Alignment *A,  char *mode)
{
  int **simm;
  int p1, p2, p3, r1, r2;
  int gap,pos0,pos1,pos2,len,sim;
  int a, b, c, m;

  simm=declare_int (A->nseq, A->nseq);



  if (strstr (mode, "sim1"))m=0;
  else if (strstr (mode, "sim2"))m=1;
  else if (strstr (mode, "sim3"))m=2;
  else if (strstr (mode, "gap1"))m=3;
  else if (strstr (mode, "cov1"))m=4;
  else if (strstr (mode, "logid"))m=5;
  else m=0;



  for ( a=0; a< A->nseq-1; a++)
    {
      simm[a][a]=MAXID;
      for ( b=a+1; b< A->nseq; b++)
	{
	  gap=pos0=pos1=pos2=len=sim=0;

	  for ( c=0; c< A->len_aln; c++)
	    {
	      r1=tolower (A->seq_al[a][c]);
	      r2=tolower (A->seq_al[b][c]);
	      p1=(r1!='-')?1:0;
	      p2=(r2!='-')?1:0;
	      p3=p1+p2;
	      if ( p3==0)continue;
	      if ( p3==1)gap++;
	      if ( r1==r2)sim++;
	      pos1+=p1;
	      pos2+=p2;
	      pos0+=(p3==2)?1:0;
	      len++;
	    }

	  if (m==0)simm[a][b]=simm[b][a]=(pos0==0)?0:(sim*MAXID)/pos0; //sim1
	  else if (m==1)  simm[a][b]=simm[b][a]=(MIN(pos1,pos2)==0)?0:(sim*MAXID)/MIN(pos1,pos2);//sim2
	  else if (m==2)   simm[a][b]=simm[b][a]=(MAX(pos1,pos2)==0)?0:(sim*MAXID)/MAX(pos1,pos2);//sim3
	  else if (m==3)   simm[a][b]=simm[b][a]=(len==0) ?0:((len-gap)*MAXID)/len;//gap1
	  else if (m==4)   simm[a][b]=simm[b][a]=(len==0) ?0:((pos0)*MAXID)/len; //cov
	  else if (m==5)
	    {

	      //Inspired from Muscle +mafft 5
	      simm[a][b]=simm[b][a]=logid_score ( sim, len);
	    }
	}
    }
  return simm;
}
int logid_score ( int sim, int len)
{
  float score;

  if ( len==0)return (int)(0.33*(float)MAXID);

  score=(float)sim/(float)len;
  if (score>0.9) score=1.0;
  else score=-log10 (1.0-score);

  score=(score*MAXID);
  return score;
}
int ** aln2dist_mat(Alignment *A)
{
  int a, b, c;

  int **d=declare_int (A->nseq, A->nseq);
  for (a=0; a<A->nseq-1; a++)
    for (b=a+1; b<A->nseq; b++)
      {
	int rr,rg,l1,l2;
	rg=rr=l1=l2=0;

	for (c=0; c<A->len_aln; c++)
	  {
	    int r1=A->seq_al[a][c];
	    int r2=A->seq_al[b][c];

	    if (r1=='-' || r2=='-');
	    else if ( r1==r2)rr++;
	    else rg++;
	  }

	d[a][b]=100-((rr*100)/(rr+rg));
	d[b][a]=d[a][b];
      }
  return d;
}
int ** aln2dist_mat_gap(Alignment *A)
{
  int a, b, c;

  int **d=declare_int (A->nseq, A->nseq);
  for (a=0; a<A->nseq-1; a++)
    for (b=a+1; b<A->nseq; b++)
      {
	int rr,rg,l1,l2;
	rg=rr=l1=l2=0;

	for (c=0; c<A->len_aln; c++)
	  {
	    int r1=A->seq_al[a][c];
	    int r2=A->seq_al[b][c];
	    if (r1!='-')l1++;
	    if (r2!='-')l2++;
	    if (r1=='-' && r2=='-');
	    else if (r1!='-' && r2!='-')rr++;
	    else if (r1=='-' || r2=='-')rg++;
	  }
	d[a][b]=((rg-fabs((l1-l2)))*100)/MIN(l1,l2);
	d[b][a]=d[a][b];
      }
  return d;
}

int ** aln2sim_mat_km (Alignment *A, char *mode)
{
  int a, b, c, d, dim;
  int n=A->nseq;
  double max_d=0;
  double **v, **di;
  int **si;
  int max=1000;
  dim=60;
  v=aln2km_vector  (A, mode, &dim);
  di=declare_double (n,n);
  si=declare_int (n,n);


  for (a=0; a<n; a++)
    for (b=0; b<n; b++)
      {
	for (c=0; c<dim; c++)
	  di[a][b]+=(v[a][c]-v[b][c])*(v[a][c]-v[b][c]);
	di[a][b]=sqrt(di[a][b]);
	if (di[a][b]>max_d)max_d=di[a][b];
      }

  for (a=0; a<n; a++)
    for (b=0; b<n; b++)
      si[a][b]=(int)(((max_d-di[a][b])/max_d)*(double)max);

  free_double (di, -1);
  free_double (v, -1);
  return si;
}
int ** aln2sim_mat (Alignment *A, char*mode)
{


  if ( strstr (mode, "idmat"))return fast_aln2sim_mat(A, mode);
  return get_sim_aln_array(A, mode);
}
int ** aln2cov (Alignment *A)
{
  int a, b, c;
  int r1, r2, gr1, gr2, pos0, gap;
  int **cov;
  cov=declare_int (A->nseq, A->nseq);

  for (a=0; a< A->nseq-1; a++)
    {
      cov[a][a]=100;
      for ( b=a+1; b<A->nseq; b++)
	{
	  for (gap=0,pos0=0,c=0;c<A->len_aln; c++)
	    {
	      r1=A->seq_al[a][c];
	      r2=A->seq_al[b][c];
	      gr1=is_gap(r1); gr2=is_gap(r2);
	      if ( gr1+gr2==0)pos0++;
	      else if ( gr1+gr2<2)gap++;
	    }
	  cov[a][b]=cov[b][a]=((gap+pos0)==0)?0:((pos0*100)/(gap+pos0));
	}
    }
  return cov;
}
int ** get_raw_sim_aln_array (Alignment *A, char *mode)
{
  int **w;
  int **M;
  int a, b, c, r1, r2, set, max, min;

  w=declare_int (A->nseq, A->nseq);
  if (strstr(mode, "sar"))M=NULL;
  else M=read_matrice (mode);

  HERE ("RAW STUFF");

  for ( set=0,a=0; a< A->nseq; a++)
    for (b=a; b<A->nseq; b++)
      {
	if (M)
	  {
	    for (c=0; c<A->len_aln; c++)
	      {
		r1=A->seq_al[a][c];
		r2=A->seq_al[b][c];

		if ( !is_gap(r1) && !is_gap(r2))
		  w[a][b]+=M[r1-'A'][r2-'A'];
	      }
	  }
	else if ( strm (mode, "sarmat2"))
	  {
	    w[a][b]=get_sar_sim2 (A->seq_al[a], A->seq_al[b]);
	  }
	else
	  {
	    HERE ("ERROR: %s is an unknown mode of raw_sim\n", mode); myexit (EXIT_FAILURE);
	  }

	w[b][a]=w[a][b];
	if (!set){min=max=w[a][b];set=1;}
	min=MIN(min,w[a][b]);
	max=MAX(max,w[a][b]);
      }
  for (a=0; a<A->nseq; a++)
    for (b=a; b<A->nseq; b++)
      {
	w[b][a]=((max-min)==0)?0:((w[b][a]-min)*100)/(max-min);
	w[a][b]=w[b][a];
      }
  free_int (M, -1);
  return w;
}
int ** array2sim (char **seq_al,int nseq, char *mode)
	{
	int **w;
	int a, b;

	
	w=declare_int ( nseq, nseq);
	
	for ( a=0; a<nseq-1; a++)
	  {
	    for ( b=a+1; b<nseq; b++)
	      {

		w[a][b]=w[b][a]=generic_get_seq_sim ( seq_al[a], seq_al[b],NULL, mode);
	      }
	  }
	return w;
	}
int ** get_sim_aln_array ( Alignment *A, char *mode)
	{
	int **w;
	int a, b;


	w=declare_int ( A->nseq, A->nseq);

	for ( a=0; a< A->nseq-1; a++)
	  {
	    for ( b=a+1; b< A->nseq; b++)
	      {

		w[a][b]=w[b][a]=generic_get_seq_sim ( A->seq_al[a], A->seq_al[b], (A->cdna_cache)?A->cdna_cache[0]:NULL, mode);
	      }
	  }
	return w;
	}
int generic_get_seq_sim ( char *seq1, char *seq2, int*cache, char *mode)
{


   if ( strm (mode, "cdna"))
     return get_cdna_seq_sim ( cache, seq1, seq2,GAP_LIST, mode);
   else if ( strnm (mode, "ktup",4))
     return ktup_comparison (seq1, seq2,atoi(mode+4));
   else if ( strstr (mode, "sarmat2"))
     {

       return get_sar_sim2 (seq1, seq2);
     }
   else if ( strstr (mode, "sarmat"))
     return (int) get_sar_sim (seq1,seq2);
   else
     {
       return get_seq_sim ( seq1,seq2,GAP_LIST, mode);
     }
}
int *** get_winsim_aln_array ( Alignment *A,char *mode, int ***w)
	{
	int a, b;
	for ( a=0; a< A->nseq; a++)
		for ( b=0; b< A->nseq; b++)
			{
			  if ( strm (mode, "cdna"))
			    w[a][b]=get_cdna_seq_winsim ( A->cdna_cache[0], A->seq_al[a], A->seq_al[b],GAP_LIST, mode, w[a][b]);
			  else
			    w[a][b]=get_seq_winsim ( A->seq_al[a], A->seq_al[b],GAP_LIST, mode, w[a][b]);
			}
	return w;
	}

Alignment * seq2profile (Sequence *S, int i)
{
  Alignment *A;

  if ((A=seq2R_template_profile (S, i)))
    {
      return A;
    }
  else
    {
      char *tmp;
      FILE *fp;
      tmp=vtmpnam (NULL);
      fp=vfopen ( tmp, "w");
      fprintf (fp, ">%s\n%s\n", S->name[i], S->seq[i]);
      vfclose (fp);

      (S->T[i])->R=fill_R_template (S->name[i], tmp, S);

      return  seq2R_template_profile (S, i);
    }
}
Alignment* remove_seq_from_aln (Alignment *A, char *seq)
{
  int a, n;
  for (n=0,a=0; a<A->nseq; a++)
    {
      if ( strm (seq, A->name[a]))continue;
      else if ( n==a);
      else
	{
	  sprintf (A->name[n], "%s",A->name[a]);
	  sprintf (A->seq_al[n], "%s",A->seq_al[a]);
	  if (A->seq_comment[a])sprintf (A->seq_comment[n], "%s", A->seq_comment[a]);
	  if (A->aln_comment[a])sprintf (A->aln_comment[n], "%s", A->aln_comment[a]);
	  A->order[n][0]=A->order[a][0];
	  A->order[n][1]=A->order[a][1];
	}
      n++;
    }
  A->nseq=n;
  return A;
}


Alignment* aln2sub_aln_file (Alignment *A, int n, char **string)
{
  char ***list;
  int a;

  list=(char***)vcalloc (A->nseq, sizeof (char***));
  if ( n==0)return A;
  else if (n>1)
    {
      int l;
      char *buf;

      for (l=0,a=0; a< n; a++)l+=strlen (string[a]);
      buf=(char*)vcalloc ( 2*n+l+1, sizeof (char));
      for (a=0; a< n; a++){buf=strcat (buf,string[a]), buf=strcat ( buf, " ");}
      list[0]=string2list (buf);
      vfree (buf);
    }
  else if ( file_exists (NULL,string[0]))
    {
      list=read_group (string[0]);

    }
  else
    {
      fprintf (stderr, "\nERROR: file <%s> does not exist [FATAL:%s]\n",string[0], PROGRAM);
      myexit (EXIT_FAILURE);
    }


  a=0;
  while (list[a])
    {
      int i, b;
      FILE *fp;
      n=atoi (list[a][0]);
      fp=vfopen (list[a][1], "w");
      for (b=2; b<n; b++)
	{
	  i=name_is_in_list (list[a][b], A->name, A->nseq, MAXNAMES);
	  if (n==3)ungap (A->seq_al[i]);
	  fprintf (fp, ">%s\n%s\n", A->name[i], A->seq_al[i]);
	}
      vfclose (fp);
      free_char (list[a], -1);
      a++;
    }
  vfree(list);
  return A;
}
Sequence *remove_empty_sequence (Sequence *S)
{
  int a, b;
  static char *c;
  Sequence *NS;
  if (!S) return S;

  
  for (a=0, b=0; a< S->nseq; a++)
    {
      c=csprintf (c, "%s", S->seq[a]);
      ungap (c);
      if ( strlen (c)==0)
	{
	  S->seq[a]=NULL;
	  add_warning ( stderr, "Sequence %s does not contain any residue: automatically removed from the set [WARNING:%s]",S->name[a], PROGRAM);
	}
    }
  NS=duplicate_sequence (S);
  free_sequence (S, S->nseq);
  
  return NS;
}
Alignment* aln2sub_seq (Alignment *A, int n, char **string)
{
  char ***list;
  int a;
  Sequence *S=NULL;

  list=(char***)vcalloc (A->nseq, sizeof (char***));
  if ( n==0)return A;
  else if (n>1)
    {
      int l;
      char *buf;

      for (l=0,a=0; a< n; a++)l+=strlen (string[a]);
      buf=(char*)vcalloc ( 2*n+l+1, sizeof (char));
      for (a=0; a< n; a++){buf=strcat (buf,string[a]), buf=strcat ( buf, " ");}
      list[0]=string2list (buf);
      vfree (buf);
    }
  else if ( file_exists (NULL,string[0]))
    {
      list=read_group (string[0]);

    }
  else
    {
      fprintf (stderr, "\nERROR: file <%s> does not exist [FATAL:%s]\n",string[0], PROGRAM);
      myexit (EXIT_FAILURE);
    }



  a=0;
  while (list[a])
    {
      int t;
      Alignment *B;
      Sequence *subS;


      B=main_read_aln (list[a][1], NULL);
      t=aln2most_similar_sequence(B, "idmat");
      subS=extract_one_seq(B->name[t],0,0,B,KEEP_NAME);
      S=add_sequence (subS,S,0);
      free_aln (B);free_sequence (subS, -1);
      vremove (list[a][1]);
      a++;
    }
  vfree(list);
  return seq2aln (S, NULL, RM_GAP);
}

Alignment * aln2collapsed_aln (Alignment * A, int n, char **string)
{
  Alignment *B;
  char ***list;
  char **list2;
  char *buf=NULL;
  FILE *fp;
  int a, b,c, ns, m, l;
  int *collapsed;

  list=(char***)vcalloc (A->nseq, sizeof (char***));
  ns=0;
  if ( n==0)return A;
  else if (n>1)
    {
      for (l=0,a=0; a< n; a++)l+=strlen (string[a]);
      buf=(char*)vcalloc ( 2*n+l+1, sizeof (char));
      for (a=0; a< n; a++){buf=strcat (buf,string[a]), buf=strcat ( buf, " ");}

      list[0]=string2list (buf);ns=1;

    }
  else if ( file_exists (NULL,string[0]))
    {
      /*Format: Fasta like, the name fo the group followed with the name of the sequences
	><Group name> <First Seq> <second seq> ....
	Groups must NOT be overlaping
      */
      l=measure_longest_line_in_file (string[0])+1;
      buf=(char*)vcalloc (l, sizeof (char));
      ns=0;
      fp=vfopen (string[0], "r");
      while ((c=fgetc(fp))!=EOF)
	{
	  buf=fgets (buf,l-1, fp);
	  if ( c=='>')list[ns++]=string2list (buf);
	}
      vfclose (fp);
    }
  else
    {
      fprintf (stderr, "\nERROR: file <%s> does not exist [FATAL:%s]\n",string[0], PROGRAM);
      myexit (EXIT_FAILURE);
    }

  vfree (buf); buf=NULL;

  /*Identify lost sequences*/
  collapsed=(int*)vcalloc (A->nseq, sizeof (int));
  for ( a=0; a< ns; a++)
      {
	m=atoi (list[a][0]);
	for (b=2; b<m ; b++)
	  {
	    c=name_is_in_list (list[a][b], A->name, A->nseq, MAXNAMES);
	    if ( c>=0)collapsed[c]=1;
	  }
      }
  for ( a=0; a< A->nseq; a++)
    {
      if ( collapsed[a]==0)
	{
	  list[ns]=declare_char (3, MAXNAMES);
	  sprintf ( list[ns][0], "3");
	  sprintf ( list[ns][1], "%s", A->name[a]);
	  sprintf ( list[ns][2], "%s", A->name[a]);
	  ns++;
	}
    }
  vfree (collapsed);





  list2=declare_char (A->nseq, 100);
  /*1 Collapse the alignment*/
  for ( a=0; a< ns; a++)
    {
      sprintf ( list2[a], "%s", list[a][2]);
    }
   B=extract_sub_aln2 ( A, ns, list2);
  /*2 Rename the sequences*/
  for ( a=0; a< ns; a++)
    {
      sprintf ( B->name[a], "%s", list[a][1]);
    }
  /*replace sequence with consensus*/

  for ( a=0; a< ns; a++)
    {
      m=atoi (list[a][0]);
      for (c=0, b=2; b<m;c++, b++)
	{
	  sprintf ( list2[c], "%s", list[a][b]);
	}
      buf=sub_aln2cons_seq_mat2 ( A,m-2,list2, "blosum62mt");
      sprintf (B->seq_al[a], "%s", buf);
    }
  vfree (buf);

  free_aln (A);
  B->S=aln2seq(B);
  return B;
}
Alignment * aln2profile (Alignment * A)
    {
      Alignment *B=NULL;
      char *cons;

      if (!A->P)
	{
	  A->P=declare_profile (AA_ALPHABET,A->len_aln+1);
	}
      B=copy_aln (A, B);
      free_int ((A->P)->count, -1);
      free_int ((A->P)->count2, -1);
      free_int ((A->P)->count3, -1);
      (A->P)->count=aln2count_mat (A);
      (A->P)->count2=aln2count_mat2 (A);

      cons=aln2cons_seq_mat (A, "blosum62mt");

      sprintf (B->seq_al[0], "%s", cons);
      B->nseq=1;
      (A->P)->count3=aln2count_mat2 (B);
      vfree (cons);
      free_aln (B);



      return A;

    }

int** aln2count_mat2 ( Alignment *A)
{
  return sub_aln2count_mat2 (A, 0, NULL);
}

int sub_aln2nseq_prf ( Alignment *A, int ns, int *ls)
{


  int a, c, s;
  Alignment *R;
  int n;
  int free_ls=0;


  if ( ns==0)
    {
      n=ns=A->nseq;
      ls=(int*)vcalloc (n, sizeof (int));
      for ( a=0; a<A->nseq; a++)ls[a]=a;
      free_ls=1;
    }
  else
    {
      n=ns;
    }

  for (c=0,a=0; a<ns; a++)
    {
      s=ls[a];
      if ( A->S && (R=seq2R_template_profile (A->S, A->order[s][0]))!=NULL)
	{
	  n+=R->nseq;
	}
      else
	{
	  ;
	}
    }

  if ( free_ls) vfree (ls);
  return n;
}

int** sub_aln2count_mat2 ( Alignment *A, int ns, int *ls)
{
  char **p;
  int **count;
  int a, b, c, s;
  Alignment *R;
  int n;
  int free_ls=0;

  if ( ns==0)
    {
      n=ns=A->nseq;
      p=(char**)vcalloc ( n, sizeof (char*));
      ls=(int*)vcalloc (n, sizeof (int));
      for ( a=0; a<A->nseq; a++)ls[a]=a;
      free_ls=1;
    }
  else
    {
      n=ns;
      p=(char**)vcalloc (n, sizeof (char*));
    }

  for (c=0,a=0; a<ns; a++)
    {
      s=ls[a];
      if ( A->S && (R=seq2R_template_profile (A->S, A->order[s][0]))!=NULL)
	{
	  n+=R->nseq;
	  p=(char**)vrealloc (p, n*sizeof (char*));
	  for (b=0; b<R->nseq; b++)
	    {
	      p[c++]=R->seq_al[b];
	    }
	}
      else
	{
	  int w;
	  w=A->order[s][4]+1;

	  for (b=0; b<w; b++)
	    p[c++]=A->seq_al[s];
	}
    }
  count=sub_aln2count_mat3 (p,c);
  vfree (p);
  if ( free_ls) vfree (ls);
  return count;
}
int** sub_aln2count_mat3 (char **al, int ns)
{
  int **count;
  int used[1000];
  int a, b;
  int r;

  int len;
  int us;


  /*count[x][0]=n symbols in column
    count[x][1]=total_size of line
    count[x][2]=Gap frequency

    count[x][n]=symbol n
    count[x][n+1]=N occurence symbol n;
    count[x][n+2]=N frequence symbol n*100;

    special multi-channeling
    count[x][count[x][1]]=Nseq
    count[x][count[x][1]+s]=residue col x, sequence s
  */


  for (a=0; a< 1000; a++)used[a]=0;
  len=strlen (al[0]);

  count=declare_int (len+2,100+ns+3);
  count[len][0]=END_ARRAY;
  count[len][1]=ns;
  count[len][2]=len;



  for (a=0; a<len; a++)
    {
      for (us=ns, b=0; b<ns; b++)
	{
	  r=tolower (al[b][a]);

	  if (is_gap(r))
	    {
	      int op=0;
	      int cl=0;
	      us--;

	      if (a>0     && al[b][a-1]=='-')op=1;
	      if (a<len-1 && al[b][a+1]=='-')cl=1;
	      if      (!cl && !op)count[a][100]++;
	      else if (op)count[a][101]++;
	      else count[a][102]++;
	    }
	  else if (used[r])
	    {
	      count[a][used[r]*3+1]++;
	    }
	  else
	    {
	      used[r]=++count[a][0];
	      count[a][used[r]*3]=r;
	      count[a][used[r]*3+1]++;
	    }
	}
      count[a][1]=count[a][0]*3+2;
      /*count[a][2]=(A->nseq-us)*100/A->nseq;*/
      count[a][2]=ns-us;

      for (b=3; b<count[a][1]; b+=3)
	{
	  count[a][b+2]=(count[a][b+1]*100)/us;
	  used[count[a][b]]=0;
	}


      /*Option for multi channeling*/

      /*
      count[a][count[a][1]]=A->nseq;
      for (b=1; b<=A->nseq; b++)
	count [a][count[a][1]+b]=(is_gap(A->seq_al[b-1][a]))?0:A->seq_al[b-1][a];
      */
    }
#ifdef XXXXXX
  HERE ("Display ");
  for (a=0; a< 5; a++)
    {
      fprintf ( stderr, "\n");
      for ( b=3; b< count[a][1]; b+=3)
	{
	  fprintf ( stderr, "[%c %d]", count[a][b], count[a][b+1]);
	}
      fprintf ( stderr, "\n");
      for ( b=0; b<ns; b++)
	{
	  fprintf ( stderr, "%c", al[b][a]);
	}
    }
  HERE ("End of Display");
#endif
  return count;
}

int** aln2count_mat ( Alignment *A)
    { /*
	function documentation: start

	int output_freq_mat ( char *outfile, Aligmnent *A)

	This function counts the number of residues in each column of an alignment (Prot/NA)
	It outputs these values in the following format

	This format can be piped into:
	The routine used for computing the p-value  gmat-inf-gc-v2c

	function documentation: end
      */

    int a, b,x;
    int **freq_mat;
    int alp_size;

    alp_size=sizeof (AA_ALPHABET);
    freq_mat=declare_int (alp_size+2, A->len_aln);


    for ( a=0; a<A->len_aln; a++)
      {
	for ( b=0; b< A->nseq; b++)
	  {
	    if ( is_gap ( A->seq_al[b][a]))freq_mat[alp_size][a]++;
	    else
	      {
		x=tolower(A->seq_al[b][a]);
		freq_mat[x-'a'][a]++;
		freq_mat[alp_size+1][a]++;

	      }
	  }
      }

    return freq_mat;
    }
char *aln2random_seq (Alignment *A, int pn1, int pn2, int pn3, int gn)
    {

      /*


	 Given the frequencies in A ( read as total counts of each Residue in
	 freq[A->nseq][A->len_aln], and pn1, pn2 and pn3:

	                  1-Generate a new amino-acid at each position
			  2-Insert Gaps, using a HMM.


	 pn3=Weight of the noise induced with sub mat.

	 pn1=% noise type 1 ( Varies with entropi)
	  n1=Ratio noise type 1

	 T =Nseq
	 t1=Noise 1 expressed in Nseq
	 al=alphabet size;
	 ncat=number of non 0 cat for a given position
	 ICi initial count for residue i

	 Ci=freq[seq][AA]
	 t1=T*n1*(1-1/ncat);
	 t2=T*n2;

	 Ci= ICi*(T-(t1+t2))/T +(t1)/al+(t2)/al

      */

      int **freq;
      int **count;
      float T, tot_t1, tot_t2,tot_t3, n1, n2, n3;
      float ncat;

      double gf;
      double *init_freq;
      double *blur_freq;
      double *t1, *t2,*t3;
      int a, b, c, x;
      char *seq;
      int tot;
      /*Viterbi  Parameters */

      int p;
      int AL=0;        /*Allowed Transition*/
      int F=-100000; /*Forbiden Transition*/

      int GAP_TRANSITION;
      int IGAP=0, IAA=1;

      int state,best_state=0, score, best_score=0;
      int p_state;
      int e=0;
      int **score_tab;
      int **state_tab;
      int nstate=2;
      int **transitions;

      int max;

      seq=(char*)vcalloc ( A->len_aln+1, sizeof (char));
      count=aln2count_mat(A);
      freq=aln2count_mat(A);

      T=100;

      n1=(float)pn1/100;
      n2=(float)pn2/100;
      n3=(float)pn3/100;

      for ( a=0; a< A->len_aln; a++)
	{
	  for ( b=0; b<26; b++)
	    freq[b][a]=freq[b][a]*((T)/(A->nseq-freq[26][a]));
	  freq[26][a]= (freq[26][a]*T)/A->nseq;
	}


      init_freq=(double*)vcalloc ( 26, sizeof (double));
      blur_freq=(double*)vcalloc ( 26, sizeof (double));

      tot_t1=tot_t2=tot_t3=0;

      t1=(double*)vcalloc ( 27, sizeof (double));
      t2=(double*)vcalloc ( 27, sizeof (double));
      t3=(double*)vcalloc ( 27, sizeof (double));
      for (a=0; a< A->len_aln; a++)
	{

        /*Compute Frequencies*/
	  for (tot=0, b=0; b<26; b++)
		{
		  if ( is_aa(b+'A'))
		    {
		      init_freq[b]=freq[b][a];
		      tot+=freq[b][a];
		    }
		}
        /*Count the number of  different amino acids*/
	for ( ncat=0, b=0; b<=26; b++)
	  {
	    ncat+=(freq[b][a]!=0)?1:0;
	  }
	/*Blurr the distribution using */
       	blur_freq=compute_matrix_p (init_freq);


	/*compute noise 1: biased with blurred content * enthropy--> keeps prosite motifs*/
	tot_t1=T*n1*(1-1/ncat);
	for (  b=0; b< 26; b++)if ( is_aa(b+'A')){t1[b]=blur_freq[b]*(1-1/ncat)*n1;}

	/*Compute noise 2: completely random*/
	tot_t2=T*n2;
	for (  b=0; b< 26; b++)if ( is_aa(b+'A')){t2[b]=tot_t2/21;}

	/*compute noise 3: biased with the sole content(pam250mt)*/
	tot_t3=T*n3;
	for (  b=0; b<26; b++)if ( is_aa(b+'A')){t3[b]=blur_freq[b]*n3;}

	for ( b=0; b<26; b++)
	  {
	    if ( is_aa('A'+b))
	      freq[b][a]=freq[b][a]*(T-(tot_t1+tot_t2+(tot_t3)))/T+t1[b]+t2[b]+t3[b];
	  }

	/*end of the loop that mutates position a*/
	}

        vfree (blur_freq);
	vfree (init_freq);
	vfree ( t3);

      /*1-Generate the amino acids of the new sequence new*/


      vsrand (0);

      for ( a=0; a< A->len_aln; a++)
	{

	  for (T=0,b=0; b<26; b++)T+=freq[b][a];
	  x=rand ()%((int)T);
	  for (c=0,b=0; b<26; b++)
	    {
	     c+=freq[b][a];
	     if ( c>=x)
	       {
		 seq[a]='A'+b;
		 c=-1;
		 break;
	       }
	    }
	  if ( c!=-1)seq[a]='-';
	}
      seq[a]='\0';


      /*2 Generate the gaps in the new sequence*/



      if ( gn<0);
      else
	{

	  transitions=declare_int ( nstate, nstate);
	  score_tab=declare_int ( A->len_aln+2, nstate       );
	  state_tab=declare_int ( A->len_aln+2, nstate       );



	  for (a=0; a<nstate;a++)
	    for (b=0; b<nstate;b++)
	      {transitions[a][b]=F;}

	  GAP_TRANSITION=AL-gn;

	  transitions[IGAP ][IGAP ]=AL;
	  transitions[IAA][IAA]=AL;
	  transitions[IAA ][IGAP]=GAP_TRANSITION;
	  transitions[IGAP][IAA ]=GAP_TRANSITION;


	  for ( p=1; p<=A->len_aln; p++){for (state=0; state< nstate; state++){score_tab[p][state]=F;state_tab[p][state]=-1;} }

	  for (p=1; p<= A->len_aln; p++)
	    {
	      for (max=0,a=0; a<26; a++)max=MAX(max, freq[a][p-1]);
	      max=(max*(A->nseq-count[26][p-1]))/A->nseq;

	      for (state=0; state< nstate; state++)
		{


		  gf=freq[26][p-1];
		  if      ( state==IGAP)  e=gf-50;
		  else if ( state==IAA )  e=max-50;
		  for (p_state=0; p_state<nstate; p_state++)
		    {
		      score=(score_tab[p-1][p_state]==F)?F:(e+transitions[p_state][state]+score_tab[p-1][p_state]);
		      if(p_state==0 || score>best_score){ best_score=score;best_state=p_state;}
		    }
		  score_tab[p][state]=best_score;
		  state_tab[p][state]=best_state;
		}
	    }

	  for (state=0; state<nstate; state++)
	    {
	      if (state==0 || score_tab[p-1][state]>best_score){best_score=score_tab[p-1][state]; best_state=state;}
	    }

	  for (p=A->len_aln; p>0;)
	    {
	      if ( best_state==IGAP)
		{
		  seq[p-1]='-';
		}
	      else if ( best_state==IAA)
		{
		  seq[p-1]=seq[p-1];
		}
	      best_state=state_tab[p][best_state];
	      p--;
	    }
        }

      free_int (freq, -1);
      return seq;
    }

/********************************************************************/
/*                                                                  */
/*			Weighting functions                         */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
Alignment * master_trimseq( Alignment *A, Sequence *S,char *mode)
     {
       Alignment *NA;
       char *p;
       int a, b;
       int use_aln=0, upper_sim=0, min_nseq=0, lower_sim=0;
       float f_upper_sim, f_lower_sim;
       char weight_mode[1000];
       char method[1000];
       int statistics=0;
       int trim_direction=TOP;
       float **sim_weight;
       int *seq_list;
       int table=0;




     /*
       mode:
           (trim)_<seq or aln>_%<percentage of tot weight to keep>_n<number of seq to keep>_w<weight mode>
     */



     seq_list=(int*)vcalloc ( S->nseq, sizeof (int));
     for ( a=0; a< A->nseq; a++)
       {
	 seq_list[a]=1;
       }


     use_aln=aln_is_aligned(A);

     if ( mode[0]=='\0')
       {

	 upper_sim=50;
	 lower_sim=0;
	 min_nseq=0;
	 sprintf (weight_mode, "pwsim");
	 sprintf ( method, "clustering2");
       }
     else
       {

	 upper_sim=lower_sim=min_nseq;
	 sprintf (weight_mode, "pwsim");
	 sprintf ( method, "clustering2");
       }

     /*
      U or % (deprecated) Upper bound for pairwise similarity
      L or m (depercated) Lower  bound for pairwise similarity
      n max number of sequences
      N max number of sequences as a fraction of thet total
      S print Statistics
      T print Table of distances
     */



     while ( (p=strtok(mode, "_")))
	   {
	     mode=NULL;
	     if (strm (p, "seq"))use_aln=0;
	     else if ( strm(p,"aln"))use_aln=1;
	     else if  (p[0]=='s')statistics=1;
	     else if  (p[0]=='t')table=1;
	     else if  (p[0]=='U')upper_sim=atoi(p+1);
	     else if  (p[0]=='L')lower_sim=atoi(p+1);
	     else if  (p[0]=='n')min_nseq=atoi(p+1);
	     else if  (p[0]=='N')min_nseq=atoi(p+1)*-1;
	     else if  (p[0]=='B')trim_direction=BOTTOM;
	     else if  (p[0]=='T')trim_direction=TOP;
	     else if  (p[0]=='W')sprintf (weight_mode, "%s", p+1);
	     else if  (p[0]=='M')sprintf (method, "%s", p+1);
	     else if  (p[0]=='K')
	       {

		 while ((p=strtok(NULL, ":")))
		   {

		     if ( p[0]=='#')
		       {
			 seq_list[atoi(p+1)-1]=2;
		       }
		     else if ( (a=name_is_in_list (p, A->name, A->nseq, 100))!=-1)

		       {
			 seq_list[a]=2;
		       }
		   }
	       }
	   }

     if ( !upper_sim && !min_nseq && !lower_sim)upper_sim=50;



     if  (!S)
       {
	 fprintf ( stderr, "\ntrimseq requires a set of sequences[FATAL:%s]\n", PROGRAM);
	 crash("");
       }

     else if ( min_nseq> S->nseq)
       {
	 min_nseq=S->nseq;
       }
     else if ( min_nseq<0)
       {
	 if ( min_nseq<-100)
	   {
	     add_warning ( stderr, "trimseq: Nseq(N)  max_val=100%% [Automatic reset]\n");
	     min_nseq=-100;
	   }

	 min_nseq=(int)((float)S->nseq*((float)min_nseq/100)*-1);
       }


     NA=seq2subseq3 (A, S,use_aln,lower_sim,upper_sim,min_nseq,trim_direction, weight_mode,&sim_weight, seq_list );

     if ( table)
       {
	 fprintf ( stderr, "\nSIMILARITY MATRIX\n");
	 for ( a=0; a< A->nseq-1; a++)
	   for ( b=a+1; b< A->nseq; b++)
	     {
	       fprintf ( stderr, "%15s Vs %15s : %3.2f %% id\n", A->name[a], A->name[b], 100-sim_weight[a][b]);
	     }
       }
     if ( statistics)
       {
	 f_upper_sim=(upper_sim>100)?((float)upper_sim/(float)100):upper_sim;
	 f_lower_sim=(upper_sim>100)?((float)lower_sim/(float)100):lower_sim;

	 fprintf ( stderr, "\nTRIM Informations:\n");
	 fprintf ( stderr, "\tUse...........: %s\n",(use_aln)?"multiple_aln":"pairwise_aln");
	 fprintf ( stderr, "\tcluster_mode..: %s\n"  ,method);
	 fprintf ( stderr, "\tsim_mode......: %s\n"  ,weight_mode);
	 fprintf ( stderr, "\tlower_id_bound: %.2f%%\n"  ,(f_lower_sim==0)?-1:f_lower_sim);
	 fprintf ( stderr, "\tupper_id_bound: %.2f%%\n",(f_upper_sim==0)?-1:f_upper_sim);
	 fprintf ( stderr, "\tnseq_kept.....: %d (out of %d)\n"  ,NA->nseq, S->nseq);
	 fprintf ( stderr, "\treduction.....: %d%% of original set\n"  ,(NA->nseq*100)/S->nseq);
	 fprintf ( stderr, "\tTrim_direction: From %s \n"  ,(trim_direction==BOTTOM)?"Bottom":"Top");
       }

     return NA;
   }

Alignment *sim_filter (Alignment *A, char *in_mode, char *seq)
{
  int **sim, **cov;
  int *list;
  int *keep;
  int maxnseq, nseq_ratio, nc;
  int new_nseq;
  int a, s, n, k;
  Alignment *R;
  char *mode;
  int outlayers;
  int direction=1;//remove the higher than
  int coverage=0; //remove based on coverage
  static char *field;
  int maxsim, minsim, maxcov, mincov;

  if ( !field) field=(char*)vcalloc (1000, sizeof (char));

  mode=(char*)vcalloc ( strlen (in_mode)+10, sizeof (char));
  sprintf ( mode, "_%s_", in_mode);

  strget_param ( mode, "_I", "100", "%d", &maxsim);
  strget_param ( mode, "_i", "0", "%d",  &minsim);
  strget_param ( mode, "_C", "100", "%d",  &maxcov);
  strget_param ( mode, "_c", "0", "%d",  &mincov);





  keep=(int*)vcalloc ( A->nseq, sizeof (int));
  list=(int*)vcalloc ( A->nseq, sizeof (int));






  if (!seq)s=0;
  else s=name_is_in_list (seq, A->name, A->nseq, 100);
  if (s==-1)
    {

      if ( s==-1)printf_exit (EXIT_FAILURE, stderr, "ERROR: %s is not a valid sequence", seq);
    }
  else
    keep[s]=1;

  //get the distances
  if ( strstr (mode, "_seq_"))
    {
      char **seq;
      int **M;

      M=read_matrice ("blosum62mt");
      seq=declare_char (A->nseq, A->len_aln+1);
      for (a=0; a<A->nseq; a++)
	{
	  sprintf ( seq[a], "%s", A->seq_al[a]);
	  ungap (seq[a]);
	}

      sim=declare_int (A->nseq, A->nseq);
      cov=declare_int (A->nseq, A->nseq);

      for (a=0; a<A->nseq; a++)
	{
	  if ( s!=a)
	    {
	      sim[s][a]=sim[a][s]=idscore_pairseq(seq[s], seq[a],-12, -1,M,"sim");
	      cov[s][a]=cov[a][s]=idscore_pairseq(seq[s], seq[a],-12, -1,M,"cov");

	    }
	}
      free_char (seq, -1);
      free_int (M,-1);
    }
  else
    {
      sim=aln2sim_mat (A, "idmat");
      cov=aln2cov (A);
    }

  for (a=0; a< A->nseq; a++)
    {
      if (a==s)continue;
      else
	{
	  if ( sim[s][a]>maxsim || sim[s][a]<minsim|| cov[s][a]<mincov||cov[s][a]>maxcov)keep[a]=-1;
	  else keep[a]=1;
	}
    }

  for ( n=0, a=0; a< A->nseq; a++)
    {
      if ( keep[a]!=-1)
	{
	  list[n++]=a;
	}
    }

  R=extract_sub_aln (A, n, list);
  free_int (sim, -1); free_int (cov, -1);vfree (list);

  return R;
}
int km_group2centroid  (int g,double **v, int n, int dim, int *size);
int * km2centroids (Alignment *A, int k, char *mode, int *keep)
{

  double **v;
  int a,b;
  int ndim=60;
  int **group;
  int g;
  Sequence *S;

  S=aln2seq(A);


  if (k>A->nseq)k=A->nseq;
  group=declare_int (k, 3);

  
  v=aln2km_vector(A,mode,&ndim);
 
  km_kmeans (v,A->nseq,ndim,k,0.0001,NULL);

  for (a=0; a<A->nseq; a++)if (!keep[a])keep[a]=-1;
  for (a=0; a<k; a++)
    {
      int s=km_group2centroid (a,v, A->nseq, ndim, &g);
      if ( s!=-1)
	{
	  keep[s]=0;
	  sprintf (A->seq_comment[s], "Kmeans Cluster Size: %4d",g);
	}
    }

  free_double (v, -1);
  return keep;
}

int km_group2centroid  (int g,double **v, int n, int dim, int *size)
{
  double *avg;
  double dist, bdist, val;
  int a, b, gs,ba;

  size[0]=0;
  avg=(double*)vcalloc ( dim, sizeof (double));

  for (a=0; a<n; a++)
    {
      if ((int)v[a][dim+1]==g)
	{
	  size[0]++;
	  for (b=0; b<dim; b++)avg[b]+=v[a][b];
	}
    }

  if (size[0]==0)return -1;
  for (b=0; b<dim; b++)avg[b]/=(double)size[0];

  for (bdist=-1,a=0; a<n; a++)
    {
      if ((int)v[a][dim+1]==g)
	{
	  for (dist=0,b=0; b<dim; b++)
	    {
	      val=v[a][b]-avg[b];
	      dist+=(val*val);
	    }

	  if ( dist<bdist || bdist==-1)
	    {
	      bdist=dist;
	      ba=a;
	    }
	}
    }
  vfree (avg);
  return ba;
}


int *seq2kmeans_class (Alignment *A, int k, char *mode)
{
  double **v;
  int dim=60;
  int *classs;
  int a;

  v=aln2km_vector(A,mode,&dim);

  km_kmeans (v,A->nseq,dim,k,0.0001,NULL);

  classs=(int*)vcalloc (A->nseq, sizeof (int));
  for (a=0; a<A->nseq; a++)classs[a]=(int)v[a][dim+1];
  free_double (v,-1);
  return classs;
}

Alignment** seq2kmeans_subset (Alignment*A, int k, int *n, char *mode)
{
  Alignment **AL;
  n[0]=0;

  if (A->nseq<=k)
    {
      AL=(Alignment**)vcalloc (1, sizeof (Alignment*));
      AL[n[0]++]=A;
    }
  else
    {
      int a,b, dim;
      double **v;
      static char *tfile;
      int *gn;
      int **gl;

      //Declare memory
      if (!tfile)tfile=vtmpnam(NULL);
      gn=(int*)vcalloc    (k, sizeof (int));
      gl=declare_int(k, A->nseq);
      AL=(Alignment**)vcalloc  (k, sizeof (Alignment*));

      //Run KM
      dim=0;
      v=aln2km_vector(A,mode,&dim);
      km_kmeans (v,A->nseq,dim,k,0.0001,NULL);

      //distribut sequences. Warning::some groups may be empty
      for (a=0; a<A->nseq; a++)
	{
	  int g=(int)v[a][dim+1];
	  int s=(int)v[a][dim+2];
	  gl[g][gn[g]++]=s;
	}
	for (a=0; a<A->nseq; ++a)
		vfree(v[a]);
	vfree(v);

      for (a=0; a<k; a++)
	{
	  FILE *fp;
	  if (gn[a])
	    {
	      fp=vfopen (tfile, "w");
	      for (b=0; b<gn[a]; b++)
		{
		  int s=gl[a][b];
		  fprintf (fp, ">%s\n%s\n", A->name[s], A->seq_al[s]);
		}
	      vfclose (fp);
	      AL[n[0]++]=main_read_aln(tfile, NULL);

	    }
	}
      vfree (gn);
      free_int (gl, -1);
    }

  return AL;

}
Alignment** seq2id_subset (Alignment*A,int k,int *ng, char *mode)
{
  int n, nn, a,b;
  int *l, *nl;
  Alignment **AL;
  static char *file1;
  static char *file2;
  int minid;
  int count=0;
  Alignment *B;

  if (!A || !A->nseq)return NULL;

  if (!file1)file1=vtmpnam (NULL);
  minid=atoi (mode);

  if (minid<100)
    myexit(fprintf_error (stderr, "minid<100 Not supported in seq2id_subset"));

  ng[0]=0;
  n=A->nseq;


  nn=0;
  l=(int*)vcalloc  ( n*2, sizeof (int));
  nl=(int*)vcalloc ( n*2, sizeof (int));
  for (a=0; a<n; a++)l[a]=a;
  AL=(Alignment**)vcalloc (n, sizeof (Alignment *));

  while (n)
    {
      int n2=0;
      int l1;
      FILE *fp2;
      char tmpfile[100];
      l1=strlen(A->seq_al[l[0]]);
      fp2=vfopen (file1, "w");
      fprintf (fp2, ">%s\n%s\n", A->name[l[0]], A->seq_al[l[0]]);
      for (a=1; a<n; a++)
	{
	  int l2=strlen(A->seq_al[l[a]]);
	  if ( l1!=l2)
	    nl[nn++]=l[a];
	  else
	    {
	      int mm=0;
	      for (b=0; b<l1; b++)if (A->seq_al[l[a]][b]!=A->seq_al[l[0]][b])mm++;
	      if (mm>5)nl[nn++]=l[a];
	      else fprintf (fp2, ">%s\n%s\n", A->name[l[a]], A->seq_al[l[a]]);
	    }


	}

      vfclose (fp2);
      fp2=vfopen (file1, "r");vfclose (fp2);
      for (a=0; a<nn; a++)l[a]=nl[a];
      B=AL[ng[0]++]=main_read_aln(file1, NULL);



      if (B->nseq==0)
	{
	  HERE ("%s", B->seq_al[0]);
	  printf_system ( "cp %s error_file", file1);exit (0);
	}
      n=nn;
      nn=0;
    }

  vfree (l);
  vfree (nl);
  return AL;
}
Alignment* km_seq (Alignment *A, int k, char *mode, char *name)
{
  FILE **f;
  int ndim=60;
  double **data2;
  int a,b, *gs;
  int **sd;
  Sequence *S;

  if (k>A->nseq)k=A->nseq;

  S=aln2seq(A);
  gs=(int*)vcalloc ( k, sizeof (int));
  f=(FILE**)vmalloc (k*sizeof (FILE*));
  data2=aln2km_vector(A,mode,&ndim);
  km_kmeans (data2,A->nseq,ndim,k,0.0001,NULL);

  for (a=0;a<k; a++)
    {
      char buf[1000];
      sprintf ( buf, "%s.cluster_%d.fasta", (name)?name:"kmeans", a+1);
      f[a]=vfopen (buf, "w");
    }

  for (a=0; a<A->nseq; a++)
    {
      int g=(int)data2[a][ndim+1];
      gs[g]++;
      fprintf (f[g], ">%s\n%s\n", S->name[a], S->seq[a]);
    }
  for (a=0; a<k; a++)
    {
      fprintf ( stderr, "\nGroup %d: %d seq ---> %s.cluster_%d.fasta", a+1, gs[a],(name)?name:"kmeans", a+1);
      vfclose (f[a]);
    }
  fprintf ( stderr, "\n");

  exit (EXIT_SUCCESS);
}

int aln2gap_trimmed (Alignment *A, int n, char *alnf, char *seqf)
{
  int **v;
  int a,b,c,t;
  int ng, nr;
  FILE *aln;
  FILE *seq;

  if (A->nseq<=n)return 0;
  if (n<0){n*=-1; n=(A->nseq*n)/100;}

  v=declare_int (A->nseq, 2);
  for (a=0; a< A->nseq; a++)v[a][0]=a;

  for (a=0; a<A->len_aln; a++)
    {
      for (ng=nr=0,b=0; b<A->nseq; b++)
	{

	  ng+=(A->seq_al[b][a]=='-')?1:0;
	  nr+=(A->seq_al[b][a]!='-')?1:0;
	}

      for (b=0; b<A->nseq; b++)
	{
	  if (A->seq_al[b][a]!='-') v[b][1]+=(nr==0)?0:(ng*100)/nr;
	}
    }
  sort_int (v, 2, 1, 0, A->nseq-1);

  seq=vfopen (seqf, "w");
  for (a=A->nseq-n; a< A->nseq; a++)
    {
      int s=v[a][0];
      ungap (A->seq_al[s]);
      fprintf ( seq, ">%s\n%s\n", A->name[s], A->seq_al[s]);
    }
  vfclose (seq);

  A->nseq-=n;
  for (c=0,a=0; a< A->len_aln; a++)
    {
      for (t=0,b=0; b<A->nseq; b++)t+=(A->seq_al[v[b][0]][a]=='-');
      if (t!=A->nseq)
	{
	  for (b=0; b<A->nseq; b++)A->seq_al[v[b][0]][c]=A->seq_al[v[b][0]][a];
	  c++;
	}
    }
  A->len_aln=c;

  aln=vfopen (alnf, "w");
  for (a=0; a<A->nseq; a++)
    {
      int s=v[a][0];
      A->seq_al[s][c]='\0';
      fprintf ( aln, ">%s\n%s\n", A->name[s], A->seq_al[s]);
    }
  vfclose (aln);


  free_aln (A);
  free_int (v, -1);
  return 1;
}



Alignment *gap_trim (Alignment *A, int f)
{
	int **v, *list;
	Alignment *R;
	size_t a,b, n;
	size_t cmax, ng, nr,sc;
	double max=0;
	if (!f) f=50;

	list=(int*)vcalloc (A->nseq,sizeof (int));
	v=declare_int (A->nseq, 2);
	for (a=0; a< A->nseq; a++)v[a][0]=a;

	for (a=0; a<A->len_aln; a++)
	{
		for (ng=nr=0,b=0; b<A->nseq; b++)
		{

			ng+=(A->seq_al[b][a]=='-')?1:0;
			nr+=(A->seq_al[b][a]!='-')?1:0;
		}

		for (b=0; b<A->nseq; b++)
		{
			if (A->seq_al[b][a]!='-')
			{
				nr=(nr==0)?1:nr;
				sc=(ng/nr)*100;
				v[b][1]+=sc;
				max+=sc;
			}
		}
	}

	max=(max*(100-f))/100;
// 	fprintf(stderr, "max %f %i", max,f);
	sort_int (v, 2, 1, 0, A->nseq-1);
	for (n=0,cmax=0,a=0;a<A->nseq; a++)
	{

		cmax+=v[a][1];
 	//	fprintf(stderr,"%i\n", cmax);
		if (cmax<max)
		{
// 			fprintf (stderr, ">%s GapScore: %d \n", A->name[v[a][0]], v[a][1]);
			list[n++]=v[a][0];
		}
		else
		{
// 			fprintf (stderr, ">%s GapScore: %d \n", A->name[v[a][0]], v[a][1]);
			fprintf (stderr, ">%s\n", A->name[v[a][0]]);
		}
	}

	HERE ("Removed %d Sequences\n", A->nseq-n);

	R=extract_sub_aln (A, n, list);
	vfree (list); free_int (v, -1);
	return R;
}
int*aln2subset_cov (Alignment *A, char *mode, int *n);
int *aln2subset    (Alignment *A, char *mode, int *n)
{
  int a,*list;


  if (strm (mode, "cov"))list=aln2subset_cov (A, "cov", n);
  else if (strm (mode, "all"))
    {
      
      n[0]=A->nseq;
      list=(int*)vcalloc ( A->nseq, sizeof (int));
      for (a=0; a<n[0]; a++) list[a]=a;
    }
  else
    {
      printf_exit (EXIT_FAILURE, stderr, "ERROR: [%s] is an unknown mode of aln2subset [FATAL:%s]\n",mode,PROGRAM);
    }
  return list;
}
int*aln2subset_cov (Alignment *A, char *mode, int *ns)
{
  int *score, *list, *rlist;
  char *cons;
  int pleft,n,left, a, b;
  int *seql;
  ns[0]=0;
  list=(int*)vcalloc (A->nseq, sizeof (int));
  cons=(char*)vcalloc (A->len_aln+1, sizeof (int));
  score=(int*)vcalloc (A->nseq, sizeof (int));
  left=A->len_aln;
  seql=(int*)vcalloc (A->nseq, sizeof(int));
  
  
  //1 deal with empty columns --- there should not be any
  for (a=0; a<A->len_aln; a++)
    {
      n=0;
      for (b=0; b<A->nseq; b++)
	if (!is_gap(A->seq_al[b][a]))n++;
      if (n==0)
	{
	  cons[a]='-';
	  left--;
	}
    }
  //add residues from the sequences covering max.
  
  pleft=left;
  while (left>0)
    {
      int best_score=0;
      int best_seq=0;
      
      for (a=0; a<A->nseq; a++)score[a]=0;
      for (a=0; a<A->len_aln; a++)
	{
	  if (!cons[a])
	    {
	      n=0;
	      for (b=0; b<A->nseq; b++)
		{
		  if (!seql[b] && !is_gap(A->seq_al[b][a]))n++;
		}
	      for (b=0; b<A->nseq; b++)
		{
		  if (!seql[b] && !is_gap(A->seq_al[b][a]))
		    {
		      score[b]+=n;
		      if (score[b]>best_score){best_score=score[b];best_seq=b;}
		    }
		}
	    }
	}
      seql[best_seq]=1;
      list[ns[0]++]=best_seq;
      for (a=0; a<A->len_aln; a++)
	{
	  if (!cons[a] && !is_gap(A->seq_al[best_seq][a]))
	    {
	      cons[a]=A->seq_al[best_seq][a];
	      left--;
	    }
	}
      if (left && left==pleft)
	{
	  fprintf ( stderr, "\nCONS: %s\n", cons);
	  print_aln (A);
	  printf_exit (EXIT_FAILURE, stderr, "\nCould not Build profile consensus::aln2cons_seq_cov [FATAL:%s]", PROGRAM);
	}
      
      pleft=left;
    }
  
  rlist=(int*)vcalloc (ns[0], sizeof (int));
  for (a=0; a<ns[0]; a++)rlist[a]=list[a];
  vfree (list);
  vfree(score);
  vfree(cons);
  vfree (seql);
  return rlist;
}



static int find_worst_seq ( int **sim, int n, int *keep, int max, int direction);
Alignment *simple_trimseq (Alignment *A, Alignment *K, char *in_mode, char *seq_list, int **sim)
{
  int *list;
  int *keep;
  int maxnseq, maxsim, nseq_ratio, nc;
  int new_nseq;
  int a,b, s, n, k;
  Alignment *R;
  char *mode;
  int outlayers;
  int direction=1;//remove the higher than
  int coverage=0; //remove based on coverage
  static char *field;
  int *tot_avg;
  int KeepN=0;
  int Print=0;

  if ( !field) field=(char*)vcalloc (1000, sizeof (char));

  mode=(char*)vcalloc ( strlen (in_mode)+10, sizeof (char));
  sprintf ( mode, "_%s_", in_mode);

  strget_param ( mode, "_%%", "0", "%d", &maxsim);
  strget_param ( mode, "_n", "0", "%d",  &maxnseq);
  strget_param ( mode, "_N", "0", "%d",  &nseq_ratio);
  strget_param ( mode, "_F", "0", "%d",  &nc);
  strget_param ( mode, "_O", "0", "%d",  &outlayers);
  strget_param ( mode, "_K", "0", "%d",  &KeepN);

  strget_param ( mode, "_f", "NAME", "%s", field);

  if ( strstr (mode, "_P_"))Print=1;

  if ( strstr (mode, "_min"))direction=-1;
  else direction=1;

  if ( strstr (mode, "_cov"))coverage=1;
  else coverage=0;


  if ( nseq_ratio)
    {
      maxnseq=(A->nseq*nseq_ratio)/100;
      maxsim=0;
    }
  else if ( maxnseq)
    {
      maxsim=0;
    }
  else if ( !maxsim)
    {
      maxsim=100;
    }


  keep=(int*)vcalloc ( A->nseq, sizeof (int));
  list=(int*)vcalloc ( A->nseq, sizeof (int));




  /*Remove Sequences that do not have at least one residue in the first and last nc columns*/
  if ( nc)
    {
      int left, right, full_n,x, y;
      int *full_list;

      Alignment *F;

      full_list=(int*)vcalloc ( A->nseq, sizeof (int));
      full_n=0;
      for (x=0; x< A->nseq; x++)
	{
	  for ( left=0,y=0; y<MIN(A->len_aln,nc); y++)
	    if (!is_gap(A->seq_al[x][y]))left=1;

	  for ( right=0,y=MAX(0,(A->len_aln-nc)); y<A->len_aln; y++)
	    if (!is_gap(A->seq_al[x][y]))right=1;

	  if ( left && right)full_list[full_n++]=x;
	}
      F=extract_sub_aln (A, full_n, full_list);
      free_aln (A);
      vfree (full_list);
      A=F;
    }

  /*Reorder the sequences according to the tree order: hopefully better phylogenetic coverage after trim*/
  if (strstr (mode, "_T") && !strstr (mode, "_kmeans"))
    {
      NT_node **T;
      Sequence *O;

      if (!sim)sim=sim_array2dist_array ( NULL, MAXID);
      T=int_dist2nj_tree (sim, A->name, A->nseq, NULL);
      O=tree2seq (T[3][0], NULL);
      A=reorder_aln (A, O->name, O->nseq);

      free_int (sim, -1);
      free_sequence (O, -1);
    }

  if (strstr (mode, "_kmeans_"))
    {
      if ( coverage)myexit (fprintf_error (stderr, "_kmeans_ does not support coverage"));
    }
  else if ( coverage==0)
    {
      if ( strstr (mode, "seq_") && !sim)sim=seq2comp_mat (aln2seq(A), "blosum62mt", "sim");
      else sim=aln2sim_mat (A, "idmat");
    }
  else
    {
      int b;
      if ( strstr (mode, "seq_") && !sim)sim=seq2comp_mat (aln2seq(A), "blosum62mt", "cov");
      else sim=aln2cov (A);

    }


  if ( K && K->nseq>0)
    {
      for ( a=0; a< K->nseq; a++)
	if ( (k=name_is_in_list (K->name[a], A->name, A->nseq, MAXNAMES+1))!=-1)
	  {
	    keep[k]=1;
	  }
    }
  if ( seq_list)
    {
      for ( a=0; a< A->nseq; a++)
	{
	  if (strstr (field, "NAME") && perl_strstr (A->name[a], seq_list)){keep[a]=1;}
	  else if (strstr (field, "COMMENT") && A->seq_comment && perl_strstr(A->seq_comment[a], seq_list)){keep[a]=1;}
	  else if (strstr (field, "SEQ") && perl_strstr((A->S)->seq[a], seq_list)){keep[a]=1;}
	}

    }
  for (a=0; a<KeepN; a++)keep[a]=1;

  if (Print)
    {
      for ( a=0; a< A->nseq; a++)
	if ( keep[a]) fprintf ( stderr, "\nFORCED KEEP %s", A->name[a]);
    }

  new_nseq=A->nseq;

  if (strstr (mode, "_kmeans_"))
    {
      //get master sequences via kmeans
      if (outlayers)myexit (fprintf_error (stderr, "_kmeans_ does not support outlayers"));
      keep=km2centroids (A, maxnseq, mode,keep);
    }
  else
    {
      while ( (s=find_worst_seq (sim, A->nseq, keep, maxsim, direction))!=-1 && new_nseq>maxnseq)
	{
	  for ( a=0; a< A->nseq; a++)sim[a][s]=sim[s][a]=-1;
	  keep[s]=-1;
	  new_nseq--;
	}
      /*Trim Outlayers*/
      if (outlayers!=0)
	{
	  int nn, b;
	  tot_avg=(int*)vcalloc ( A->nseq, sizeof (int));

	  for (a=0; a<A->nseq; a++)
	    {
	      if ( keep[a]==-1)tot_avg[a]=-1;
	      else
		{
		  for (nn=0, b=0; b< A->nseq; b++)
		    {
		      if (a==b || keep[b]==-1)continue;
		      else
			{
			  tot_avg[a]+=sim[a][b];
			  nn++;
			}
		    }
		  tot_avg[a]=(nn==0)?-1:(tot_avg[a])/nn;
		}
	    }
	  for ( a=0; a<A->nseq; a++)
	    {
	      if (tot_avg[a]!=-1 && tot_avg[a]<outlayers)
		{
		  fprintf ( stderr, "\nREMOVED OUTLAYER: %3d %% avg similarity with remaining sequences [Seq %s]", tot_avg[a],A->name[a]);
		  keep[a]=-1;
		}
	    }
	  vfree ( tot_avg);
	}
    }
  for ( n=0, a=0; a< A->nseq; a++)
    {
      if ( keep[a]!=-1)
	{
	  list[n++]=a;
	}
    }

  R=extract_sub_aln (A, n, list);
  free_int (sim, -1); vfree (list);

  return R;
}

int find_worst_seq ( int **sim, int n, int *keep,int max,int direction)
{
  int **sc;
  int a, b, r=0;
  int si;

  sc=declare_int (n, 2);
  if (direction==-1)max=100-max;

  for ( a=0; a< n; a++) sc[a][0]=a;
  for ( a=0; a< n-1; a++)
    {
      for ( b=a+1; b<n; b++)
	{

	  if (sim[a][b]>=0)si=(direction==-1)?100-sim[a][b]:sim[a][b];
	  else si=sim[a][b];
	  if ( si>max)
		{
		  if ( keep[a]!=1)sc[a][1]+=si;
		  if ( keep[b]!=1)sc[b][1]+=si;
		}
	}
    }

  sort_int_inv ( sc, 2, 1, 0, n-1);
  if ( sc[0][1]>0)r=sc[0][0];
  else r=-1;

  free_int (sc, -1);
  if (r!=-1 && keep && keep[r])return -1;
  else return r;
}

int find_worst_seq_old ( int **sim, int n, int *keep,int max,int direction)
{
  int **sc;
  int a, b, r=0;

  sc=declare_int (n, 2);

  for ( a=0; a< n; a++) sc[a][0]=a;
  for ( a=0; a< n-1; a++)
    {
      for ( b=a+1; b<n; b++)
	{
	  if ( direction==1)
	    {
	      if ( sim[a][b]>max)
		{
		  if ( keep[a]!=1)sc[a][1]+=sim[a][b];
		  if ( keep[b]!=1)sc[b][1]+=sim[a][b];
		}
	    }
	  else if ( direction == -1)
	    {
	      if ( sim[a][b]<max && sim[a][b]>=0)
		{
		  if ( keep[a]!=1)sc[a][1]+=sim[a][b];
		  if ( keep[b]!=1)sc[b][1]+=sim[a][b];
		}
	    }
	}
    }

  if ( direction ==1) //remove max
    {
      sort_int_inv ( sc, 2, 1, 0, n-1);
      if ( sc[0][1]>0)r=sc[0][0];
      else r=-1;

    }
  else if ( direction ==-1)//remove min
    {
      sort_int_inv ( sc, 2, 1, 0, n-1);
      if ( sc[0][1]>=0)r=sc[0][0];
      else r=-1;
      HERE ("** %d %d\n", r,sc[0][1]);
    }
  free_int (sc, -1);
  if (r!=-1 && keep && keep[r])return -1;
  else return r;
}


Alignment * trimseq( Alignment *A, Sequence *S,char *mode)
   {
     Alignment *NA;
     char *p;
     int a, b;
     int use_aln=0, upper_sim=0, min_nseq=0, lower_sim=0;
     char weight_mode[1000];
     char method[1000];
     int statistics=0;
     int trim_direction=TOP;
     float **sim_weight;
     int *seq_list;
     int table=0;
     int print_name=0;
     float f_lower_sim, f_upper_sim;



     /*
       mode:
           (trim)_<seq or aln>_%<percentage of tot weight to keep>_n<number of seq to keep>_w<weight mode>
     */



     seq_list=(int*)vcalloc ( S->nseq, sizeof (int));
     for ( a=0; a< A->nseq; a++)
       {
	 seq_list[a]=1;
       }


     use_aln=aln_is_aligned(A);


     if ( mode[0]=='\0')
       {

	 upper_sim=50;
	 lower_sim=0;
	 min_nseq=0;
	 sprintf (weight_mode, "pwsim_fragment");
	 sprintf ( method, "clustering2");
       }
     else
       {

	 upper_sim=lower_sim=min_nseq;
	 sprintf (weight_mode, "pwsim_fragment");
	 sprintf ( method, "clustering2");
       }

     /*
      U or % (deprecated) Upper bound for pairwise similarity
      L or m (depercated) Lower  bound for pairwise similarity
      n max number of sequences
      N max number of sequences as a fraction of thet total
      S print Statistics
      T print Table of distances
     */



     while ( (p=strtok(mode, "_")))
	   {
	     mode=NULL;
	     if (strm (p, "seq"))use_aln=0;
	     else if ( strm(p,"aln"))use_aln=1;
	     else if  (p[0]=='s')statistics=1;
	     else if  (p[0]=='t')table=1;
	     else if  (p[0]=='p')print_name=1;
	     else if  (p[0]=='U')upper_sim=atoi(p+1);
	     else if  (p[0]=='L')lower_sim=atoi(p+1);
	     else if  (p[0]=='n')min_nseq=atoi(p+1);
	     else if  (p[0]=='N')min_nseq=atoi(p+1)*-1;
	     else if  (p[0]=='B')trim_direction=BOTTOM;
	     else if  (p[0]=='T')trim_direction=TOP;
	     else if  (p[0]=='W')sprintf (weight_mode, "%s", p+1);
	     else if  (p[0]=='M')sprintf (method, "%s", p+1);
	     else if  (p[0]=='K')
	       {

		 while ((p=strtok(NULL, ":")))
		   {

		     if ( (a=name_is_in_list (p, A->name, A->nseq, 100))!=-1)
		       {
			 seq_list[a]=2;
		       }
		   }
	       }
	   }

     if ( !upper_sim && !min_nseq && !lower_sim)upper_sim=50;



     if  (!S)
       {
	 fprintf ( stderr, "\ntrimseq requires a set of sequences[FATAL:%s]\n", PROGRAM);
	 crash("");
       }

     else if ( min_nseq> S->nseq)
       {
	 min_nseq=S->nseq;
       }
     else if ( min_nseq<0)
       {
	 if ( min_nseq<-100)
	   {
	     add_warning ( stderr, "trimseq: Nseq(N)  max_val=100%% [Automatic reset]\n");
	     min_nseq=-100;
	   }

	 min_nseq=(int)((float)S->nseq*((float)min_nseq/100)*-1);
       }


     NA=seq2subseq2 (A, S,use_aln,lower_sim,upper_sim,min_nseq,trim_direction, weight_mode,&sim_weight, seq_list );

     if ( table)
       {
	 fprintf ( stderr, "\nSIMILARITY MATRIX\n");
	 for ( a=0; a< A->nseq-1; a++)
	   for ( b=a+1; b< A->nseq; b++)
	     {
	       fprintf ( stderr, "%15s Vs %15s : %3.2f %% id\n", A->name[a], A->name[b], 100-sim_weight[a][b]);
	     }
       }

     NA=seq_name2removed_seq_name(S, NA,sim_weight);

     if ( print_name)
       {
	 fprintf ( stderr, "\nList of sequences with their closest removed neighbors\n");
	 for ( a=0; a< NA->nseq; a++)fprintf ( stderr, "\n%s: %s\n", NA->name[a], NA->seq_comment[a]);
       }

     if ( statistics)
       {
	 f_lower_sim=(lower_sim>100)?(float)lower_sim/100:lower_sim;
	 f_upper_sim=(upper_sim>100)?(float)upper_sim/100:upper_sim;

	 fprintf ( stderr, "\nTRIM seq Informations:\n");
	 fprintf ( stderr, "\tUse...........: %s\n",(use_aln)?"multiple_aln":"pairwise_aln");
	 fprintf ( stderr, "\tcluster_mode..: %s\n"  ,method);
	 fprintf ( stderr, "\tsim_mode......: %s\n"  ,weight_mode);
	 fprintf ( stderr, "\tlower_id_bound: %.2f%%\n"  ,(f_lower_sim==0)?-1:f_lower_sim);
	 fprintf ( stderr, "\tupper_id_bound: %.2f%%\n",(f_upper_sim==0)?-1:f_upper_sim);
	 fprintf ( stderr, "\tnseq_kept.....: %d (out of %d)\n"  ,NA->nseq, S->nseq);
	 fprintf ( stderr, "\treduction.....: %d%% of original set\n"  ,(NA->nseq*100)/S->nseq);
	 fprintf ( stderr, "\tTrim_direction: From %s \n"  ,(trim_direction==BOTTOM)?"Bottom":"Top");
       }

     return NA;
   }

Alignment * tc_trimseq( Alignment *A, Sequence *S,char *mode)
   {
     Alignment *NA;
     Sequence  *TS;
     char *trimfile, *alnfile;
     int *seq_list;
     int a, nseq=0, sim=0;
     char *p;
     char command[100000];
     char keep_list[10000];

     int top, bottom, middle, pmiddle;

     keep_list[0]='\0';

     seq_list=(int*)vcalloc ( S->nseq, sizeof (int));
     for ( a=0; a< A->nseq; a++)
       {
	 seq_list[a]=1;
       }

     trimfile=vtmpnam (NULL);
     alnfile=vtmpnam (NULL);
     if ( !aln_is_aligned (A))
       {
	 fprintf ( stderr, "\ntrimTC: computation of an Approximate MSA  [");
	 A=compute_tcoffee_aln_quick ( A, NULL);
	 fprintf ( stderr, "DONE]\n");
       }
     output_clustal_aln (alnfile, A);


     while ( (p=strtok(mode, "#")))
	   {
	     mode=NULL;


	     if (p[0]=='%' || p[0]=='S')sim=(p[1]=='%')?atoi(p+2):atoi(p+1);
	     else if  (p[0]=='n' || p[0]=='N')nseq=atoi(p+1);
	     else if  (p[0]=='K')
	       {
		 if ( (a=name_is_in_list (p+1, A->name, A->nseq, 100))!=-1)
		   {
		     seq_list[a]=2;
		   }

	       }
	   }
     if ( nseq ==0 && sim ==0)
       {
	 fprintf ( stderr, "\nERROR: trimTC\nIndicate the maximum number of sequences Nnseq\nOR the maximum average similarity of the chosen sequencesSx\nEX: +trimTC S20 OR +trimTC N5");
	 fprintf ( stderr, "\n[FATAL:%s]", PROGRAM);
	 myexit (EXIT_FAILURE);
       }

     for ( a=0; a<A->nseq; a++)if (seq_list[a]==2){strcat ( keep_list, A->name[a]);strcat ( keep_list," ");}

     if ( sim)
       {
	 sprintf ( command , "%s -infile %s -trim  -trimfile=%s  -split_score_thres %d -convert -iterate 0 ",get_string_variable("t_coffee"), alnfile, trimfile,sim);
	 if ( keep_list[0]){strcat ( command, " -seq_to_keep ");strcat ( command, keep_list);}
	 my_system ( command);
	 TS=read_sequences (trimfile);
       }
     else if ( nseq && A->nseq>nseq)
       {

	 top=100;bottom=0;
	 pmiddle=0;middle=50;

	 sprintf ( command , "%s -infile %s -trim  -trimfile=%s  -split_score_thres %d -convert -iterate 0",get_string_variable("t_coffee"), alnfile, trimfile,middle);
	 if ( keep_list[0]){strcat ( command, " -seq_to_keep ");strcat ( command, keep_list);}
	 my_system ( command);

	 TS=read_sequences (trimfile);
	 fprintf ( stderr, "\n\tTrimTC: Sim %d Nseq %d\t",middle, TS->nseq);

	 if ( TS->nseq>nseq)top=middle;
	 else if ( TS->nseq<nseq)bottom=middle;
	 pmiddle=middle;
	 middle=(top-bottom)/2+bottom;

	 while (TS->nseq!=nseq && pmiddle!=middle)
	   {

	     sprintf ( command , "%s -infile %s -trim  -trimfile=%s  -split_score_thres %d -convert -iterate 0 ",get_string_variable("t_coffee"), alnfile, trimfile,middle);
	     if ( keep_list[0]){strcat ( command, " -seq_to_keep ");strcat ( command, keep_list);}
	     my_system ( command);
	     free_sequence (TS, -1);
	     TS=read_sequences (trimfile);
	     fprintf ( stderr, "\n\tTrimTC: Sim %d Nseq %d\t", middle, TS->nseq);

	     if ( TS->nseq>nseq)top=middle;
	     else if ( TS->nseq<nseq)bottom=middle;
	     pmiddle=middle;
	     middle=(top-bottom)/2+bottom;
	   }
       }
     else
       {
	 TS=aln2seq (A);
       }
     NA=seq2aln (TS, NULL, 1);
     vremove ( alnfile);
     fprintf ( stderr, "\n");

     return NA;
   }

Alignment* seq2subseq3( Alignment *A, Sequence *S,int use_aln, int int_lower_sim,int int_upper_sim, int min_nseq, int trim_direction, char *weight_mode, float ***sim_weight, int *seq_list)
{
  int a, b;
  int new_nseq;

  /*OUTPUT*/
  char **seq, **name;
  Sequence *NS;
  Alignment *NA;
  float sim, lower_sim, upper_sim;

  lower_sim=(int_lower_sim>100)?(float)int_lower_sim/100:int_lower_sim;
  upper_sim=(int_upper_sim>100)?(float)int_upper_sim/100:int_upper_sim;

  sim_weight[0]=get_weight   ((use_aln)?A:NULL, S, weight_mode);

  name=declare_char (S->nseq, (MAXNAMES+1));
  seq= declare_char (S->nseq, S->max_len+1);

  /*
    Remove every sequence that is more than upper_sim and less than lower_sim similar to the master sequences
    the master sequence(s) are those for which seq_list[x]==2
  */




  new_nseq=A->nseq;


  for (a=0; a< A->nseq; a++)
    {
      if ( seq_list[a]==2)
	{

	  for ( b=0; b< A->nseq;b++)
	    {
	      sim=100-sim_weight[0][a][b];
	      if (seq_list[b]==1 && (sim>upper_sim || sim<lower_sim))
		{
		 seq_list[b]=0;
		 new_nseq--;
		}
	    }

	}
    }

  /*Prepare the new sequence List*/

  for (b=0, a=0; a<S->nseq; a++)
	  {
	    if ( seq_list[a])
	      {
		sprintf ( name[b], "%s", S->name[a]);
		sprintf ( seq[b] , "%s",(use_aln)?A->seq_al[a]: S->seq[a] );
		b++;
	      }
	  }


  NS=fill_sequence_struc (new_nseq,seq,name, NULL);
  NA=seq2aln(NS,NULL,1);

  if ( use_aln && A)
    {
      NA=realloc_aln2  ( NA,A->max_n_seq,A->len_aln+1);

      for (b=0, a=0; a<S->nseq; a++)
	{
	  if ( seq_list[a])
	    {
	      sprintf ( NA->seq_al[b] , "%s",A->seq_al[a]);
	      b++;
	      }
	}

      NA->len_aln=A->len_aln;
      ungap_aln(NA);
    }


  return NA;
}
Alignment* seq2subseq2( Alignment *A, Sequence *S,int use_aln, int int_lower_sim,int int_upper_sim, int min_nseq, int trim_direction, char *weight_mode, float ***sim_weight, int *seq_list)
{
  int a, b;
  int new_nseq;
  int seq_index=0;
  /*OUTPUT*/
  char **seq, **name;
  Sequence *NS;
  Alignment *NA;
  float lower_sim, upper_sim;

  lower_sim=(int_lower_sim>100)?(float)int_lower_sim/100:int_lower_sim;
  upper_sim=(int_upper_sim>100)?(float)int_upper_sim/100:int_upper_sim;


  sim_weight[0]=get_weight   ((use_aln)?A:NULL, S, weight_mode);

  name=declare_char (S->nseq, (MAXNAMES+1));
  seq= declare_char (S->nseq, S->max_len+1);

  /*
    1 REMOVE OUTLAYERS
    2 REMOVE CLOSELY RELATED SEQUENCES
    3 IF STILL TOO MANY SEQUENCES:
      REMOVE THE MOST CLOSELY RELATED ONES
  */


  /*1 Remove outlayers*/

  new_nseq=A->nseq;


  /*1 Remove outlayers*/
  while ( lower_sim && (extreme_seq(BOTTOM,A,sim_weight[0],seq_list, &seq_index) <lower_sim) && ((new_nseq)>min_nseq) && seq_index!=-1)
    {

      if ( seq_list[seq_index]==1)
	{
	  seq_list[seq_index]=0;
	  new_nseq--;
	}
    }
   /*2 Remove close relative*/


  while ( upper_sim && (extreme_seq(TOP, A,sim_weight[0],seq_list, &seq_index)>upper_sim) && ((new_nseq)>min_nseq)&& seq_index!=-1)
    {

      if ( seq_list[seq_index]==1)
	{
	  seq_list[seq_index]=0;
	  new_nseq--;
	}
    }


  /*Remove extra sequences*/

  while ( min_nseq>0 && new_nseq>min_nseq && seq_index!=-1)
    {

      extreme_seq(trim_direction, A,sim_weight[0],seq_list, &seq_index);

      if ( seq_index==-1)break;
      if ( seq_list[seq_index]==1)
	{
	  seq_list[seq_index]=0;
	  new_nseq--;
	}
    }


  /*Prepare the new sequence List*/

  for (b=0, a=0; a<S->nseq; a++)
	  {
	    if ( seq_list[a])
	      {
		sprintf ( name[b], "%s", S->name[a]);
		sprintf ( seq[b] , "%s",(use_aln)?A->seq_al[a]: S->seq[a] );
		b++;
	      }
	  }


  NS=fill_sequence_struc (new_nseq,seq,name, NULL);
  NA=seq2aln(NS,NULL,1);

  if ( use_aln && A)
    {
      NA=realloc_aln2  ( NA,A->max_n_seq,A->len_aln+1);

      for (b=0, a=0; a<S->nseq; a++)
	{
	  if ( seq_list[a])
	    {
	      sprintf ( NA->seq_al[b],"%s",A->seq_al[a]);
	      b++;
	      }
	}

      NA->len_aln=A->len_aln;
      ungap_aln(NA);
    }


  return NA;
}

float extreme_seq (int direction, Alignment *A,float **sim_weight,int *seq_list, int *seq_index)
{

  /*find the closest relative of each sequence
    Return:
          Direction= BOTTOM: the sequence whose closest relative is the most distant
	  Direction= TOP:    the sequence whose closest relative is the closest
  weight: different sequences=100
          similar sequences  =0
  */
  int a, b;

  float top_sim,bottom_sim, best_sim, sim;
  int   top_seq, bottom_seq;

  bottom_seq=top_seq=seq_index[0]=-1;
  top_sim=-1;
  bottom_sim=101;

  for (a=0; a< A->nseq; a++)
    {
      if (seq_list[a]!=1)continue;

      for ( best_sim=0, b=0; b< A->nseq; b++)
	{
	  if ( a==b || !seq_list[b])continue;

	  sim=100-sim_weight[a][b];
	  if (sim>best_sim)
	    {
	      best_sim=sim;
	    }
	}

      if ( best_sim>top_sim)
	{
	  top_seq=a;
	  top_sim=best_sim;
	}

      if ( best_sim<bottom_sim)
	{
	  bottom_seq=a;
	  bottom_sim=best_sim;
	}

    }
  if ( direction==BOTTOM  ){seq_index[0]= bottom_seq; return bottom_sim;}
  else if ( direction==TOP){seq_index[0]= top_seq; return top_sim;}
  else
    {
      seq_index[0]=-1;
      return -1;
    }
}




Alignment* seq2subseq1( Alignment *A, Sequence *S,int use_aln, int percent,int max_nseq, int ms,char *weight_mode)
   {
     float **pw_weight,**sim_weight, **seq_weight;
     int a,b,c,d;
     float sum, chosen,last_chosen, last_nchosen,nchosen;
     int condition1, condition2;
     Sequence *NS;
     Alignment *NA;
     char **name, **seq;
     float score, best_score;
     int best_seq=0;
     int *seq_list, *used_seq_list;

     /*
       mode:
           (trim)_<seq or aln>_%<percentage of tot weight to keep>_n<number of seq to keep>_w<weight mode>
     */

     sim_weight=get_weight   ((use_aln)?A:NULL, S, weight_mode);
     pw_weight=declare_float (S->nseq, S->nseq);
     seq_weight=declare_float ( S->nseq, 2);


     for (best_score=0,a=0; a<S->nseq; a++)
       {
	 for ( b=0; b<S->nseq; b++)
	   {
	     if ( a==b)continue;
	     seq_weight[a][0]+=sim_weight[a][b];
	   }
	 seq_weight[a][0]=seq_weight[a][0]/(S->nseq-1);
	 score=seq_weight[a][0]=100-seq_weight[a][0];

	 if ( score>best_score)
	   {
	     best_seq=a;
	     best_score=score;
	   }

       }
      for (a=0; a<S->nseq; a++)
       {
	 for ( b=0; b<S->nseq; b++)
	   {
	     if ( a==b)continue;
	     pw_weight[a][b]=sim_weight[a][b]*seq_weight[a][0]*seq_weight[b][0]/(100*100);

	   }
       }


     seq_list=(int*)vcalloc ( S->nseq, sizeof (int));
     used_seq_list=(int*)vcalloc ( S->nseq, sizeof (int));



     name=declare_char (S->nseq, (MAXNAMES+1));
     seq= declare_char (S->nseq, S->max_len+1);

     /*compute the normalization factor*/
     for (sum=0,d=0; d< S->nseq; d++)
	   {
	     for (score=0,c=0; c<S->nseq; c++)
	       {
		 if ( c!=d)
		   score=MAX(score, 100-sim_weight[c][d]);
	       }
	     sum+=score;
	   }
     sum=sum/S->nseq;
     /*chose the first sequence */
     for ( best_score=0,a=0; a< S->nseq; a++)
       {
	 for (score=0, b=0; b< S->nseq; b++)
	   {
	     score+=100-sim_weight[a][b];
	   }
	 if ( score>best_score)
	   {
	     best_seq=a;
	     best_score=score;
	   }

       }


     last_chosen=chosen=((best_score/S->nseq)*100)/sum;
     nchosen=last_nchosen=1;
     seq_list[0]=best_seq;
     used_seq_list[best_seq]=1;

     sprintf ( name[0],"%s", S->name[seq_list[0]]);
     sprintf ( seq[0],"%s", S->seq[seq_list[0]]);
     nchosen=last_nchosen=1;


     fprintf ( stderr, "\nTRIM:\n");
     fprintf ( stderr, "\n1-Chosen Sequences\n");
     /*Assemble the list of sequences*/
     for (a=1; a< S->nseq; a++)
       {
	 for (best_score=0,b=0; b< S->nseq; b++)
	   {
	     if (used_seq_list[b]);
	     else
	       {
		 score=pw_weight[seq_list[0]][b]+1;
		 for (c=0; c<a; c++)
		   score=MIN(score,pw_weight[seq_list[c]][b]);

		 if ( score>=best_score)
		   {
		   best_seq=b;
		   best_score=score;
		   }

	       }
	   }
	 seq_list[a]=best_seq;
	 used_seq_list[best_seq]=1;



	 for ( chosen=0,d=0; d< S->nseq; d++)
	   {
	     for (score=0, c=0; c<=a; c++)
	       {
		 if ( seq_list[c]!=d)
		   score=MAX(score, 100-sim_weight[seq_list[c]][d]);
	       }
	     chosen+=score;

	   }

	 chosen=((chosen/S->nseq)*100)/sum;
	 nchosen=a+1;

	 condition1= (int)chosen<=(int)percent || !percent;
	 condition2=(nchosen)<=max_nseq     || !max_nseq;

	 if (condition1 && condition2)
	   {
	     fprintf ( stderr, "\tADD %s (set score: %.2f %%)\n", S->name[seq_list[a]], chosen);
	     sprintf ( name[a],"%s", S->name[seq_list[a]]);
	     sprintf ( seq[a],"%s", S->seq[seq_list[a]]);

	   }
	 else
	   {
	     break;
	   }
	 last_chosen=chosen;
	 last_nchosen=nchosen;
       }

     NS=fill_sequence_struc (last_nchosen,seq,name, NULL);
     NA=seq2aln(NS,NULL,1);
     fprintf ( stderr, "\n2-Informations:\n");
     fprintf ( stderr, "\tUse...........: %s\n",(use_aln)?"multiple_aln":"pairwise_aln");
     fprintf ( stderr, "\tweight_mode...: %s\n"  ,weight_mode);
     fprintf ( stderr, "\tpercent_weight: %.2f%% (max=%d%%)\n",last_chosen,percent);
     fprintf ( stderr, "\tn_seq.........: %d\n"  ,NS->nseq);
     fprintf ( stderr, "\treduction.....: %d%% of original set\n"  ,(NS->nseq*100)/S->nseq);

     return NA;
   }
Sequence  * seq_weight2species_weight (Alignment *A, Sequence *S)
{
  float *wsp;
  float *wseq;
  int a,b;

  S->W=declare_weights(S->nseq);
  if (!A->S || !(A->S)->W)aln2voronoi_weights (A);

  wseq=((A->S)->W)->SEQ_W;
  wsp=(S->W)->SEQ_W;
  for ( a=0; a< S->nseq; a++)
    {
      for (b=0; b<A->nseq; b++)
	if ( strstr (A->name[b], S->name[a]))wsp[a]+=wseq[b];
    }
  for (a=0; a<S->nseq; a++)
    fprintf ( stderr, "\nVoronoi Weights: Species %s ---> %.2f\n", S->name[a], wsp[a]);
  return S;
}
Alignment * aln2voronoi_weights (Alignment *A)
{
  int a, b, c;
  float t=0;
  int **tab;
  float *w;

  tab=declare_int (256, A->nseq+1);
  if (A->S)free_sequence (A->S, (A->S)->nseq);
  A->S=aln2seq(A);
  (A->S)->W=declare_weights (A->nseq);
  w=((A->S)->W)->SEQ_W;

  for (a=0; a<A->len_aln; a++)
    {
      for ( b=0; b<A->nseq; b++)
	{
	  c= A->seq_al[b][a];
	  if (!is_gap(c))
	    {
	      c=tolower(c);
	      tab[c][++tab[c][0]]=b;
	    }
	}
      for (c=0; c<256; c++)
	{
	  if (tab[c][0])
	    {
	      for (b=1; b<=tab[c][0]; b++)
		{
		  w[tab[c][b]]+=(float)1/(float)tab[c][0];
		  t+=(float)1/(float)tab[c][0];
		}
	    }
	  tab[c][0]=0;
	}
    }
  for (a=0; a<A->nseq; a++)
    {
      w[a]=(w[a]/t)*A->nseq;
    }

  return A;
}

float ** get_weight ( Alignment *A, Sequence *S, char *mode)
{
 char *aln_name;
 char *weight_name;
 char *seq_name;
 char command[LONG_STRING];
 char program[LONG_STRING];
 float **weight;
 FILE *fp;
 int c;

 if ( !mode || !mode[0] || strm (mode, "msa"))
      {
	if ( getenv ( "SEQ2MSA_WEIGHT")==NULL)sprintf (program, "%s",SEQ2MSA_WEIGHT);
	else sprintf ( program, "%s", (getenv ( "SEQ2MSA_WEIGHT")));
      }
 else if ( strm(mode, "pwsim") ||strm(mode, "pwsim_fragment") )
      {
	return seq2pwsim (A, S, mode);
      }
 else
     {
       if (getenv (mode))sprintf ( program, "%s", (getenv (mode)));
       else fprintf ( stderr, "\nERROR: %s is not a valid mode for weight computation [FATAL:%s]", mode, PROGRAM);
     }

 /*MSA weights*/
 seq_name=vtmpnam(NULL);
 aln_name=vtmpnam(NULL);
 weight_name=vtmpnam(NULL);
 weight=declare_float (S->nseq+1, 2);



  if (A)
    {
      output_clustal_aln (seq_name,A);
      output_fasta_seq   (aln_name,A);
      sprintf ( command, "%s %s -i %s -w %s", program, seq_name, aln_name, weight_name);
    }
  else
    {
      A=seq2aln(S,A,1);
      output_fasta_seq   (seq_name,A);
      sprintf ( command, "%s %s -w %s", program, seq_name, weight_name);
    }


  my_system ( command);

  fp=vfopen( weight_name, "r");
  while ( (c=fgetc(fp))!='$');
  c=fgetc(fp);
  c=0;
  while ( (fscanf (fp, "%*s %f\n",&(weight[c][1])))==1)
    {weight[c][0]=c;c++;}
  vfclose (fp);


  return weight;
}

float **seq2pwsim (	   Alignment *A, Sequence *S, char *mode)
{
  int a, b, c;
  float d,t;
  float  **W;
  Alignment *B;
  W=declare_float (S->nseq, S->nseq);



  for (a=0; a< S->nseq; a++)
	for ( b=a; b<S->nseq; b++)
	  {
	    if ( a==b){d=1;}
	    else if (!A)
	      {

		B=align_two_sequences ((S)->seq[a], (S)->seq[b],"pam250mt", -10, -1, "fasta_pair_wise");
		for (t=0,d=0,c=0; c<B->len_aln; c++)
		  {
		    d+=(B->seq_al[0][c]==B->seq_al[1][c] && !is_gap(B->seq_al[0][c]));
		    t+=(!is_gap(B->seq_al[0][c]) && !is_gap(B->seq_al[1][c]));
		  }
		t=(strm ( mode, "pwsim_fragment"))?B->len_aln:t;

		d=d/((t==0)?1:t);
		free_aln(B);
	      }
	    else
	      {
		for (t=0,d=0,c=0; c<A->len_aln; c++)
		  {
		    d+=(A->seq_al[a][c]==A->seq_al[b][c] && !is_gap(A->seq_al[a][c]));
		    t+=(!is_gap(A->seq_al[a][c]) && !is_gap(A->seq_al[b][c]));
		  }
		d=d/((t==0)?1:t);
	      }


	    W[a][b]=W[b][a]=(1-d)*100;
	  }


  return W;

}

float **seq2pwsim_fragment (	   Alignment *A, Sequence *S, char *mode)
{


  int a, b, c;
  float d,t;
  float  **W;
  Alignment *B;
  W=declare_float (S->nseq, S->nseq);




  for (a=0; a< S->nseq; a++)
	for ( b=a; b<S->nseq; b++)
	  {
	    if ( a==b){d=1;}
	    else if (!A)
	      {

		B=align_two_sequences ((S)->seq[a], (S)->seq[b],"pam250mt", -10, -1, "fasta_pair_wise");
		for (t=0,d=0,c=0; c<B->len_aln; c++)
		  {
		    d+=(B->seq_al[0][c]==B->seq_al[1][c] && !is_gap(B->seq_al[0][c]));
		    t+=(!is_gap(B->seq_al[0][c]) && !is_gap(B->seq_al[1][c]));
		  }

		d=d/((t==0)?1:t);
		free_aln(B);
	      }
	    else
	      {
		for (t=0,d=0,c=0; c<A->len_aln; c++)
		  {
		    d+=(A->seq_al[a][c]==A->seq_al[b][c] && !is_gap(A->seq_al[a][c]));
		    t+=(!is_gap(A->seq_al[a][c]) && !is_gap(A->seq_al[b][c]));
		  }
		d=d/((t==0)?1:t);
	      }


	    W[a][b]=W[b][a]=(1-d)*100;
	  }


  return W;

}

/********************************************************************/
/*                                                                  */
/*			AMINO ACID FUNCTIONS                        */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
//Builds an extended alphabet from a string
char** string2alphabet (char *string, int depth, int *falp_size)
{
  int max_s;
  int a, b,c, l, n;
  char buf[1000];
  char **alp;
  int alp_size;

  char ***alp2;
  int *alp2_size;

  int *array;
  char **falp;


  l=strlen (string);
  array=(int*)vcalloc ( 256, sizeof (int));


  max_s=l+1;
  falp_size[0]=0;
  falp=declare_char (l+1, 2);

  alp=declare_char(l,2);
  alp_size=0;

  array=(int*)vcalloc ( 256, sizeof (int));
  for (a=0;a<l; a++)
    {
      if (!array[(int)string[a]])
	{
	  array[(int)string[a]]=1;
	  sprintf (alp[alp_size++], "%c", string[a]);
	  sprintf (falp[falp_size[0]++], "%c", string[a]);
	}
    }
  sprintf ( falp[falp_size[0]++], "*");
  vfree (array);

  if ( depth==1)
    {
      free_char (alp, -1);
      return falp;
    }
  alp2=(char***)vcalloc ( depth, sizeof (char**));
  alp2_size=(int*)vcalloc (depth, sizeof (int));

  for (a=0; a<depth; a++)
    {
      alp2[a]=alp;
      alp2_size[a]=alp_size;
    }


  for (a=2; a<=depth; a++)
    {
      char ***result_array;

      result_array=generate_array_string_list (a, alp2, alp2_size, &n, NULL, NO_OVERLAP);
      max_s+=n+1;
      falp=(char**)vrealloc (falp, sizeof (char**)*max_s);
      for (b=0; b<n; b++)
	{
	  buf[0]='\0';
	  for (c=0; c<a; c++)
	    {
	      strcat (buf, result_array[b][c]);
	    }
	  falp[falp_size[0]]=(char*)vcalloc (strlen (buf)+1, sizeof (char));
	  sprintf ( falp[falp_size[0]++], "%s", buf);
	  vfree ( result_array[b]);
	}
      vfree (result_array);

    }

  falp[falp_size[0]]=(char*)vcalloc (2, sizeof (char));
  sprintf ( falp[falp_size[0]++], "*");
  free_char (alp, -1);
  return falp;
}

char** make_group_aa (int *ngroup, char *mode)
	{
/*mode:         indicates which matrix will be used for the grouping*/
/*n_group:      pointer to the number of groups                     */
/*return value: an array of strings containing the AA of each group */


	int **matrix;
	int a, b,c,is_in;
	char buf[28];
	char **group_list;
	char *matrix_name;
	int extend=0;
	matrix_name=(char*)vcalloc ( 100, sizeof (char));

	if (ngroup[0]==-1)extend=1;

	ngroup[0]=0;
	group_list=declare_char ( 100, 27);

	if (extend)
	  {
	    sprintf ( group_list[ngroup[0]++], "gG");
	    sprintf ( group_list[ngroup[0]++], "pP");
	    sprintf ( group_list[ngroup[0]++], "aA");
	    sprintf ( group_list[ngroup[0]++], "cC");
	    sprintf ( group_list[ngroup[0]++], "dD");
	    sprintf ( group_list[ngroup[0]++], "eE");

	    sprintf ( group_list[ngroup[0]++], "fF");
	    sprintf ( group_list[ngroup[0]++], "hH");
	    sprintf ( group_list[ngroup[0]++], "iI");
	    sprintf ( group_list[ngroup[0]++], "kK");
	    sprintf ( group_list[ngroup[0]++], "lL");
	    sprintf ( group_list[ngroup[0]++], "mM");
	    sprintf ( group_list[ngroup[0]++], "nN");
	    sprintf ( group_list[ngroup[0]++], "qQ");
	    sprintf ( group_list[ngroup[0]++], "rR");

	    sprintf ( group_list[ngroup[0]++], "sS");
	    sprintf ( group_list[ngroup[0]++], "tT");
	    sprintf ( group_list[ngroup[0]++], "vV");
	    sprintf ( group_list[ngroup[0]++], "wW");
	    sprintf ( group_list[ngroup[0]++], "*");
	  }


	if ( mode && mode[0]=='_'){mode++;sprintf ( matrix_name, "%s", mode);}


	if (mode==NULL || mode[0]=='\0' || strstr (mode, "mat_"))
	  {
	    if (mode==NULL || mode[0]=='\0')sprintf ( matrix_name, "idmat");
	    else if (strstr (mode, "mat_"))sprintf ( matrix_name, "%s", strstr (mode, "mat_")+4);



	    matrix=read_matrice ( matrix_name);

	    for ( a=0;a< 26; a++)
	      {
		if ( matrix[a][a]>0)
		  {
		    for ( c=0,b=0;b< 26; b++)
		      {

			if ( matrix[a][b]>0 && matrix[b][b]>0)
			  {
			    buf[c++]=b+'A';
			    buf[c++]=b+'a';
			  }
		      }
		    buf[c]='\0';
		    for ( is_in=0,b=0; b< ngroup[0]; b++)if ( strcmp (buf, group_list[b])==0)is_in=1;
		    if (is_in==0)sprintf ( group_list[ngroup[0]++], "%s", buf);

		  }
	      }
	    free_int (matrix, -1);
	    vfree (matrix_name);
	  }

	else if ( strstr (mode, "sim") || strm (mode, "idmat") || mode==NULL)
	  {
	    sprintf ( group_list[ngroup[0]++], "aA");
	    sprintf ( group_list[ngroup[0]++], "bB");
	    sprintf ( group_list[ngroup[0]++], "cC");
	    sprintf ( group_list[ngroup[0]++], "dD");
	    sprintf ( group_list[ngroup[0]++], "eE");
	    sprintf ( group_list[ngroup[0]++], "fF");
	    sprintf ( group_list[ngroup[0]++], "gG");
	    sprintf ( group_list[ngroup[0]++], "hH");
	    sprintf ( group_list[ngroup[0]++], "iI");
	    sprintf ( group_list[ngroup[0]++], "jJ");
	    sprintf ( group_list[ngroup[0]++], "kK");
	    sprintf ( group_list[ngroup[0]++], "lL");
	    sprintf ( group_list[ngroup[0]++], "mM");
	    sprintf ( group_list[ngroup[0]++], "nN");
	    sprintf ( group_list[ngroup[0]++], "oO");
	    sprintf ( group_list[ngroup[0]++], "pP");
	    sprintf ( group_list[ngroup[0]++], "qQ");
	    sprintf ( group_list[ngroup[0]++], "rR");
	    sprintf ( group_list[ngroup[0]++], "sS");
	    sprintf ( group_list[ngroup[0]++], "tT");
	    sprintf ( group_list[ngroup[0]++], "uU");
	    sprintf ( group_list[ngroup[0]++], "vV");
	    sprintf ( group_list[ngroup[0]++], "wW");
	    sprintf ( group_list[ngroup[0]++], "xX");
	    sprintf ( group_list[ngroup[0]++], "yY");
	    sprintf ( group_list[ngroup[0]++], "zZ");
	    vfree (matrix_name);
	  }
	else if ( strm (mode, "simple"))
	     {
	       sprintf ( group_list[ngroup[0]++], "avilmAVILM");
	       sprintf ( group_list[ngroup[0]++], "dekrDEKR");
	       sprintf ( group_list[ngroup[0]++], "stcnqhSTCNQH");
	       sprintf ( group_list[ngroup[0]++], "wfyWFY");
	       sprintf ( group_list[ngroup[0]++], "gG");
	       sprintf ( group_list[ngroup[0]++], "pP");
	       vfree (matrix_name);
	     }

	else if ( strm (mode, "mafft"))
	     {


	       sprintf ( group_list[ngroup[0]++],"agjopstAGJOPST");
	       sprintf ( group_list[ngroup[0]++],"ilmvILMV");
	       sprintf ( group_list[ngroup[0]++],"bdenqzBDENQZ");
	       sprintf ( group_list[ngroup[0]++],"hkrHKR");
	       sprintf ( group_list[ngroup[0]++],"fwyFWY");
	       sprintf ( group_list[ngroup[0]++],"cC");
	       vfree (matrix_name);
	     }
	else if ( strm (mode, "clustalw"))
	     {

		 sprintf ( group_list[ngroup[0]++],"astaASTA");
		 sprintf ( group_list[ngroup[0]++],"bneqkBNEQK");
		 sprintf ( group_list[ngroup[0]++],"cnhqkCNHQK");
		 sprintf ( group_list[ngroup[0]++],"dndeqDNDEQ");
		 sprintf ( group_list[ngroup[0]++],"eqhrkEQHRK");
		 sprintf ( group_list[ngroup[0]++],"fmilvFMILV");
		 sprintf ( group_list[ngroup[0]++],"gmilfGMILF");
		 sprintf ( group_list[ngroup[0]++],"hhyHHY");
		 sprintf ( group_list[ngroup[0]++],"ifywIFYW");
		 sprintf ( group_list[ngroup[0]++],"jcJC");
		 sprintf ( group_list[ngroup[0]++],"kpKP");
		 vfree (matrix_name);
	     }
	else if ( strm (mode, "polarity"))
	     {

	       sprintf ( group_list[ngroup[0]++],"eqrsdnkhtEQRSDNKHT");
	       sprintf ( group_list[ngroup[0]++],"pP");
	       sprintf ( group_list[ngroup[0]++],"gG");
	       sprintf ( group_list[ngroup[0]++],"cC");
	       sprintf ( group_list[ngroup[0]++],"fywFYW");
	       sprintf ( group_list[ngroup[0]++],"iavlmIAVLM");
	       vfree (matrix_name);
	     }
	else if ( strm (mode, "vasiliky"))
	     {
		 ngroup[0]=0;
		 sprintf ( group_list[ngroup[0]++], "rkRK");
		 sprintf ( group_list[ngroup[0]++], "deDE");
		 sprintf ( group_list[ngroup[0]++], "qhQH");
		 sprintf ( group_list[ngroup[0]++], "vilmVILM");
		 sprintf ( group_list[ngroup[0]++], "fyFY");
		 sprintf ( group_list[ngroup[0]++], "sS");
		 sprintf ( group_list[ngroup[0]++], "wW");
		 sprintf ( group_list[ngroup[0]++], "aA");
		 sprintf ( group_list[ngroup[0]++], "cC");
		 sprintf ( group_list[ngroup[0]++], "gG");
		 sprintf ( group_list[ngroup[0]++], "nN");
		 sprintf ( group_list[ngroup[0]++], "pP");
		 sprintf ( group_list[ngroup[0]++], "tT");
		 vfree (matrix_name);

	     }
	else if ( strm (mode, "clustalw_col"))
	     {
		 sprintf ( group_list[ngroup[0]++], "staSTA");
		 sprintf ( group_list[ngroup[0]++], "neqkNEQK");
		 sprintf ( group_list[ngroup[0]++], "nhqkNHQK");
		 sprintf ( group_list[ngroup[0]++], "ndeqNDEQ");
		 sprintf ( group_list[ngroup[0]++], "qhrkQHRK");
		 sprintf ( group_list[ngroup[0]++], "milvMILV");
		 sprintf ( group_list[ngroup[0]++], "milfMILF");
		 sprintf ( group_list[ngroup[0]++], "hyHY");
		 sprintf ( group_list[ngroup[0]++], "fywFYW");
		 sprintf ( group_list[ngroup[0]++], "gG");
		 sprintf ( group_list[ngroup[0]++], "pP");
		 sprintf ( group_list[ngroup[0]++], "cC");
		 vfree (matrix_name);

	     }
	else if ( strm (mode, "clustalw_dot"))
	     {
		 sprintf ( group_list[ngroup[0]++], "csaCSA");
		 sprintf ( group_list[ngroup[0]++], "atvATV");
		 sprintf ( group_list[ngroup[0]++], "sagSAG");
		 sprintf ( group_list[ngroup[0]++], "stnkSTNK");
		 sprintf ( group_list[ngroup[0]++], "stpaSTPA");
		 sprintf ( group_list[ngroup[0]++], "sgndSGND");
		 sprintf ( group_list[ngroup[0]++], "sndeqkSNDEQK");
		 sprintf ( group_list[ngroup[0]++], "ndeqhkNDEQHK");
		 sprintf ( group_list[ngroup[0]++], "neqhrkNEQHRK");
		 sprintf ( group_list[ngroup[0]++], "fvlimFVLIM");
		 sprintf ( group_list[ngroup[0]++], "hfyHFY");
		 vfree (matrix_name);
	     }
	else if ( strm (mode, "make_all"))
	     {
		 ngroup[0]=1;
		 sprintf ( group_list[0], "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
		 vfree (matrix_name);

	     }

	return group_list;
	}
char** make_group_aa_upgma (char*matrix, int max_n)
	{
	  char **group_list;
	  int **mat;
	  int *used;
	  int a, b, ba, bb, best, set, l, n;
	  l=26;

	  group_list=declare_char (l+1, l+1);
	  for (a=0; a<l; a++)group_list[a][0]='a'+a;
	  mat=read_matrice(matrix);
	  used=(int*)vcalloc ( l, sizeof (int));
	  n=l;

	  while (n>max_n)
	    {
	      for (set=0,a=0; a<l-1; a++)
		for (b=a+1; b<l; b++)
		  {
		    if (used[a]||used[b])continue;

		    if (set==0 || mat[a][b]>best)
		      {
			best=mat[a][b];
			ba=a;
			bb=b;
			set=1;
		      }
		  }

	      for (a=0; a<l; a++)
		{
		  mat[ba][a]=mat[a][ba]=(mat [ba][a]+mat[bb][a])/2;
		  used[bb]=1;
		}
	      strcat (group_list[ba], group_list[bb]);
	      vfree (group_list[bb]);
	      group_list[bb]=NULL;

	      n--;
	    }

	  for (n=0,a=0; a<l; a++)
	    {
	      if ( group_list[a])
		group_list[n++]=group_list[a];
	    }
	  vfree (used); free_int (mat, -1);
	  return group_list;
	}

int find_group_aa_distribution (char *col, int nseq,int n_group, char **gl,  int *distrib, char *mode )
	{
	static int *distribution;
	static char **lgl;
	static int ln_group;
	int a, b, c;
	int *d;
	char **gl2;
	int n_group2;



	if ( lgl==NULL)
		lgl=make_group_aa ( &ln_group, mode);

		if ( gl==NULL)
		{
		gl2=lgl;
		n_group2=ln_group;
		}
	else
		{
		gl2=gl;
		n_group2=n_group;
		}

	if ( distribution==NULL || ln_group<n_group)distribution=(int*)vcalloc ( n_group2, sizeof (int));
	if ( distrib==NULL)d=distribution;
	else d=distrib;


	for ( a=0; a< n_group2; a++)d[a]=0;

	for ( a=0; a< nseq; a++)
		{
		for ( b=0; b< n_group2; b++)
			d[b]+=is_in_set (col[a], gl2[b]);
		}
	c=d[0];
	for ( a=0; a< n_group2; a++)
		c=(d[a]>c)?d[a]:c;
	return c;
	}



int is_in_same_group_aa ( char r1, char r2, int n_group, char **gl, char *mode)
	{
	int a;
	static char **lgl;
	static int ln_group;

	char **gl2;
	int n_group2;

	/*use mode=idmat for similarity based on id*/

	r1=toupper(r1);
	r2=toupper(r2);
	if (mode==NULL)return (r1==r2)?1:0;

	if ( strm (mode, "clean"))
	     {
	     free_char (lgl, -1);
	     lgl=NULL;
	     ln_group=0;
	     return 0;
	     }
	else if ( strstr (mode, "cov"))
	  {
	    return 1;
	  }

	if ( lgl==NULL)
	  {
	    lgl=make_group_aa ( &ln_group, mode);
	  }

	if ( gl==NULL)
	  {
	    gl2=lgl;
	    n_group2=ln_group;
	  }
	else
	  {
	    gl2=gl;
	    n_group2=n_group;
	  }

	for ( a=0; a< n_group2; a++)
	  {
	    if ( is_in_set ( r1, gl2[a]) && is_in_set ( r2, gl2[a]))
	      {
		return 1;
	      }
	  }
	return 0;
	}


Alignment * gene2prot (Alignment *A){return A; }
char * test_gene2prot (Constraint_list *CL, int s1)
       {
	   int a, b,q, nal;
	   int F=-10000000; /*FORBIDEN STATE*/
	   int AL=0;       /*ALLOWED STATE*/
	   int SPLICE_PENALTY=1000;
	   int FRAME_PENALTY=1000;


	   int START,  ORF1,  ORF2, ORF3, s5NC;
	   int s3NC,ORF3_G1, ORF3_T2, ORF3_NC, ORF3_A3, ORF3_T4;
	   int U1_G1,   U1_T2,   U1_NC,   U1_A3,   U1_T4;
	   int U2_G1,   U2_T2,   U2_NC,   U2_A3,   U2_T4;
	   int U1,     U2, U3,  U4, U5, END;

	   int nstate=0;
	   int **transitions;
	   int **v_tab;
	   int **v_tab_p;
	   int **last_coding;
	   int **last_t4;
	   int *potential;
	   int v;

	   int orf1, orf2, orf3, ncp, p, state, pstate, e, best_state_p=0, best_state_v=0, best_pstate_p=0, best_pstate_v;
	   char *seq, *seq2, *seq3;
	   int l;
	   int  *is_coding;
	   int *is_t4;
	   char *codon;
	   int s, r, s2, r2, w2;

	   static int *entry;
	   int tot=0;

	   seq=(char*)vcalloc ( strlen ((CL->S)->seq[s1])+1, sizeof (char));
	   seq2=(char*)vcalloc ( strlen ((CL->S)->seq[s1])+1, sizeof (char));
	   seq3=(char*)vcalloc ( strlen ((CL->S)->seq[s1])+1, sizeof (char));
	   sprintf ( seq, "%s", (CL->S)->seq[s1]);
	   ungap (seq);

	   l=strlen (seq);
	   for ( a=0; a< l; a++) seq[a]=tolower ( seq[a]);
	   for ( a=0; a< l; a++) seq[a]=(seq[a]=='t')?'u': seq[a];


	   potential=(int*)vcalloc (l+1, sizeof (int));

	   for (nal=0, s=0; s<(CL->S)->nseq; s++)
	     {
	       for ( r=1; r<=(CL->S)->len[s]; r++)
		 {
		   for ( b=1; b<CL->residue_index[s1][r][0]; b++)
		     {

		       s2=CL->residue_index[s][r][b+SEQ2];
		       r2=CL->residue_index[s][r][b+R2];
		       w2=CL->residue_index[s][r][b+WE];
		       if (s==s1)potential[r-1]+=w2;
		       else if ( s2==s1)potential[r2-1]+=w2;
		       tot+=w2;
		       nal++;
		     }
		 }
	     }


	   SPLICE_PENALTY=10000;
	   FRAME_PENALTY=1000;


	   nstate=0;
	   START=nstate++;  ORF1=nstate++;  ORF2=nstate++; ORF3=nstate++; s5NC=nstate++;
	   s3NC=nstate++;
	   ORF3_G1=nstate++;U1_G1=nstate++;U2_G1=nstate++;
	   ORF3_T2=nstate++;U1_T2=nstate++;U2_T2=nstate++;
	   ORF3_NC=nstate++;U1_NC=nstate++;U2_NC=nstate++;
	   ORF3_A3=nstate++;U1_A3=nstate++;U2_A3=nstate++;
	   ORF3_T4=nstate++;U1_T4=nstate++;U2_T4=nstate++;


	   U1=nstate++;     U2=nstate++; U3=nstate++;  U4=nstate++; U5=nstate++;
	   END=nstate++;

	   is_coding=(int*)vcalloc ( nstate, sizeof (int));
	   is_coding[ORF1]=is_coding[ORF2]=is_coding[ORF3]=is_coding[U1]=is_coding[U2]=1;
	   is_coding[U3]=is_coding[U4]=is_coding[U5]=1;

	   is_t4=(int*)vcalloc ( nstate, sizeof (int));
	   is_t4[ORF3_T4]=is_t4[U1_T4]=is_t4[U2_T4]=1;
	   transitions=declare_int ( nstate, nstate);
	   for (a=0; a< nstate; a++)
		   for ( b=0; b< nstate; b++)transitions[a][b]=F;

	   transitions[START][ORF1]=AL;
	   transitions[START][s5NC]=AL-FRAME_PENALTY;
	   transitions[s5NC][s5NC]=AL;

	   transitions[s5NC][ORF1]=AL-FRAME_PENALTY;

	   transitions[ORF1][ORF2]=AL;
	   transitions[ORF2][ORF3]=AL;
	   transitions[ORF3][U1]=AL;
	   transitions[ORF3][ORF1]=AL;
	   transitions[ORF3][ORF3_G1]=AL-SPLICE_PENALTY;


	   transitions[ORF3_G1][ORF3_T2]=AL;
	   transitions[ORF3_T2][ORF3_NC]=AL;
	   transitions[ORF3_NC][ORF3_NC]=AL;
	   transitions[ORF3_NC][ORF3_A3]=AL;
	   transitions[ORF3_A3][ORF3_T4]=AL;
	   transitions[ORF3_T4][ORF1]=AL-SPLICE_PENALTY;

	   transitions[U1][U2]=AL;
	   transitions[U1][U1_G1]=AL-SPLICE_PENALTY;
	   transitions[U1_G1][U1_T2]=AL;
	   transitions[U1_T2][U1_NC]=AL;
	   transitions[U1_NC][U1_NC]=AL;
	   transitions[U1_NC][U1_A3]=AL;
	   transitions[U1_A3][U1_T4]=AL;
	   transitions[U1_T4][U3]=AL-SPLICE_PENALTY;
	   transitions[U3][U4]=AL;
	   transitions[U4][ORF1]=AL;

	   transitions[U2][U2_G1]=AL-SPLICE_PENALTY;
	   transitions[U2_G1][U2_T2]=AL;
	   transitions[U2_T2][U2_NC]=AL;
	   transitions[U2_NC][U2_NC]=AL;
	   transitions[U2_NC][U2_A3]=AL;
	   transitions[U2_A3][U2_T4]=AL;
	   transitions[U2_T4][U5]=AL-SPLICE_PENALTY;
	   transitions[U5][ORF1]=AL;

	   transitions[ORF3][s3NC]=AL-FRAME_PENALTY;
	   transitions[ORF3][END]=AL;
	   transitions[s3NC][END]=AL;


	   v_tab=declare_int ( l+1,nstate);
	   v_tab_p=declare_int ( l+1,nstate);
	   last_coding=declare_int ( l+1,nstate);
	   last_t4=declare_int ( l+1,nstate);

	   for (a=0; a< l; a++) potential[a]-=200;

	   codon=(char*)vcalloc ( 4, sizeof (char));
	   best_pstate_p=START;
	   best_pstate_v=0;
	   nal=0;
	   for ( p=1; p<=l; p++)
	       {
	       if  (translate_dna_codon (seq+(p-1), 'x')=='x' || p>(l-2))orf1=F;
	       else orf1=potential[p-1];

	       if  (p<2 || translate_dna_codon (seq+(p-2), 'x')=='x' || p>(l-1))orf2=F;
	       else orf2=potential[p-1];


	       if  (p<3 || translate_dna_codon (seq+(p-3), 'x')=='x' || p>l)orf3=F;
	       else orf3=potential[p-1];

	       if ( best_int (3, 1, &a, orf1, orf2, orf3)!=F)ncp=-best_int (3, 1, &a, orf1, orf2, orf3);
	       else ncp=1000;

	       for ( state=0; state< nstate; state++)
	           {

		       if      ( state==ORF1)e=orf1;
		       else if ( state==ORF2)e=orf2;
		       else if ( state==ORF3)e=orf3;
		       else if ( state>=U1 && state<=U3)
			   {
			   e=0;
			   }
		       else if ( state==U4)
		          {
			      codon[2]=seq[p-1];
			      codon[1]=seq[last_coding[p-1][U3]-1];
			      codon[0]=seq[last_coding[p-2][U1_T4]-1];
			      if ( translate_dna_codon (codon, 'x')=='x')e=F;
			      else e=0;
			  }
		       else if ( state==U5)
		          {
			      codon[2]=seq[p-1];
			      codon[1]=seq[last_coding[p-1][U2_T4]-1];
			      q=seq[last_coding[p-1][U2_T4]];
			      codon[0]=seq[last_coding[q-1][U1]-1];
			      if ( translate_dna_codon (codon, 'x')=='x')e=F;
			      else e=0;
			  }

		       else if (state>=ORF3_G1 && state<=U2_G1)e=(p<l-1 && seq[p-1]=='g' && seq[p]=='u')?ncp:F;
		       else if ( state>=ORF3_T2 && state<=U2_T2)
			   {
			   e=(p>1 && seq[p-2]=='g' && seq[p-1]=='u')?ncp:F;
			   }
		       else if ( state>=ORF3_A3 && state<=U2_A3)e=(seq[p-1]=='a')?ncp:F;
		       else if ( state>=ORF3_T4 && state<=U2_T4)e=(seq[p-1]=='u')?ncp:F;
		       else e=ncp;

		       for ( pstate=0; pstate<nstate; pstate++)
		           {
			       if (e==F ||  transitions[pstate][state]==F || v_tab[p-1][pstate]==F)v=F;
			       else v=e+transitions[pstate][state]+v_tab[p-1][pstate];

			       if ( pstate==0 || v>best_pstate_v)
			          {best_pstate_v=v;best_pstate_p=pstate;}
			   }
		      v_tab[p][state]=best_pstate_v;
		      v_tab_p[p][state]=best_pstate_p;

		      if (!is_coding[state])last_coding[p][state]=last_coding[p-1][best_pstate_p];
		      else if (is_coding[state])last_coding[p][state]=p;

		      if (!is_t4[state])
		         {
			     if (is_coding[state] && last_t4[p-1][best_pstate_p]==0)last_t4[p][state]=p;
			     else last_t4[p][state]=last_t4[p-1][best_pstate_p];
			 }
		      else if (is_t4[state])last_t4[p][state]=p;

		      if (state==0 ||best_pstate_v>best_state_v ){best_state_p=state; best_state_v=best_pstate_v;}
		   }
	       }
	   tot=0;
	   for ( p=l; p>0; p--)
	           {
		       if ( best_state_p>=ORF1 &&  best_state_p<=ORF3){seq2[tot++]=tolower (seq[p-1]);}
		       else if ( best_state_p>=U1 && best_state_p<=U5){seq2[tot++]=tolower (seq[p-1]);}
		       if (best_state_p==ORF1)seq[p-1]=toupper (seq[p-1]);
		       else if (best_state_p==ORF2 || best_state_p==ORF3)seq[p-1]=tolower (seq[p-1]);
		       else if ( best_state_p==ORF3_NC || best_state_p==U1_NC ||  best_state_p==U2_NC) seq[p-1]='.';
		       else if ( best_state_p==U1 || best_state_p==U2 || best_state_p==U3 || best_state_p==U4 || best_state_p==U5) seq[p-1]=best_state_p-U1+'1';
		       else seq[p-1]=toupper (seq[p-1]);
		       best_state_p=v_tab_p[p][best_state_p];
		   }

	   for ( a=0, b=tot-1; b>=0; b--, a++)
	       seq3[a]=seq2[b];

	   fprintf ( stderr, "\n%s\n", seq);
	   fprintf ( stderr, "\nN coding=%d\n", tot);
	   for ( a=0; a< tot; a+=3)
	        {
		b=translate_dna_codon (seq3+a, 'x');
		fprintf ( stderr, "%c",b);
		if ( b=='x'){fprintf ( stderr, "\n");myexit (EXIT_SUCCESS);}
		}

	    fprintf ( stderr, "\n");
	    myexit (EXIT_SUCCESS);
	    return 0;



       }
Alignment * dna_aln2_3frame_cdna_aln(Alignment *A,int *ns,int **l_s)
{
  Alignment *B;
  int a;
  B=realloc_aln2 (NULL,6,strlen(A->seq_al[l_s[0][0]])+strlen(A->seq_al[l_s[1][0]]));
  for ( a=0; a< 3; a++)
    {
      B->seq_al[a]=translate_dna_seq (A->seq_al[l_s[0][0]]+a, 0, 'o',B->seq_al[a]);
      B->seq_al[a+3]=translate_dna_seq (A->seq_al[l_s[1][0]]+a, 0, 'o',B->seq_al[a+3]);
    }
  for ( a=1; a<3; a++)
    {
      if ( strlen(B->seq_al[a])<strlen(B->seq_al[0])) B->seq_al[a]=strcat ( B->seq_al[a], "x");
      if ( strlen(B->seq_al[a+3])<strlen(B->seq_al[3])) B->seq_al[a+3]=strcat ( B->seq_al[a+3], "x");
    }

  B->nseq=6;
  B->len_aln=strlen (B->seq_al[0]);
  return B;
}

//JM_ADD
//For normal distribution scan
#ifndef PI
#define PI 3.141592653589793238462643
#endif

double normal(double x, double mean, double std)
{
	return (1/(std*sqrt(2.0*PI)))*exp((-0.5*(x-mean)*(x-mean))/(std*std));
}

int ** get_sim_aln_array_normal_distribution ( Alignment *A, char *mode, int *STD, int *CENTER)
	{
	int **w;
	int a, b;


	w=declare_int ( A->nseq, A->nseq);

	for ( a=0; a< A->nseq-1; a++)
	  {
	    for ( b=a+1; b< A->nseq; b++)
	      {

		w[a][b]=w[b][a]=generic_get_seq_sim_normal_distribution ( A->seq_al[a], A->seq_al[b], (A->cdna_cache)?A->cdna_cache[0]:NULL, mode, STD, CENTER);
	      }
	  }
	return w;
	}
int generic_get_seq_sim_normal_distribution ( char *seq1, char *seq2, int*cache, char *mode, int *STD, int *CENTER)
{
  return get_seq_sim_distribution ( seq1,seq2,GAP_LIST, mode, STD, CENTER);
}

int get_seq_sim_distribution ( char *string1, char *string2, char *ignore, char *in_mode, int *STD, int *CENTER)
	{
	int len1;
	int a;
	int pos0, gap=0;
	int p1, p2;
	int r=0,r1=0,r2=0;
	char *p;
	char mode[1000];

	double sim;


	sprintf ( mode, "%s", in_mode);

	/*mode: <mat>__<sim_mode>
	  mat: idscore to get the alignment done
	       any legal cw matrix
	  sim_mode: sim1->identities/matches
                    sim2->identities/min len
	*/


	if ( (p=strstr (mode, "_"))!=NULL)
	  {
	    p[0]='\0';
	    p++;
	  }


	if (strstr (mode, "idscore"))
	  {
	    static int **mat;
	    if (!mat) mat=read_matrice ("blosum62mt");
	    return idscore_pairseq (string1, string2, -12, -1, mat,mode);
	  }

	len1=strlen (string1);
	for ( sim=pos0=0,a=0; a< len1; a++)
		{
		  r1=string1[a];
		  r2=string2[a];
		  p1=1-is_in_set (r1, ignore);
		  p2=1-is_in_set (r2, ignore);
		  if (p1 && p2)
			{
			    pos0++;
			    if (is_in_same_group_aa(r1,r2,0, NULL, mode))
			    {
			      sim += normal(a, *CENTER, *STD);
			    }
			}
		  else if (p1+p2==1)
		    {
		      gap++;
		    }
		}

	if ( p==NULL || strm (p, "sim1") || strm (p, "sim"))
	  {
	    r=(pos0==0)?0:(sim*MAXID);
	  }
/*	else if ( strm (p, "sim2"))
	  {
	    r=(pos1==0 || pos2==0)?0:(sim*MAXID)/MIN(pos1,pos2);
	  }
	else if ( strm (p, "sim3"))
	  {
	    r=(pos1==0 || pos2==0)?0:(sim*MAXID)/MAX(pos1,pos2);
	  }
	else if ( strm (p, "gap1"))
	  {
	    r=(len1==0)?MAXID:(gap*MAXID)/len1;
	    r=MAXID-r;
	  }
	else if ( strm (p, "logid"))
	  {
	    r=logid_score (pos0, sim);
	  }*/
	return r;

	}

Alignment *aln2clean_pw_aln (Alignment *A, OveralnP *F)// char *mode, int t, int f, int p1,int p2, int p3, char *fsa_mode)
{
  int **C, **T;
  int a, b, c;
  Alignment *B;


  if (F->t==0)F->t=2;

  C=declare_int ( A->nseq, A->len_aln);
  T=declare_int ( A->nseq, A->len_aln);
  B=copy_aln (A, NULL);

  for (a=0; a< A->nseq;a++)
    {
      for (b=0; b<A->nseq; b++)
	{
	  int *w;
	  w=pw_aln2clean_aln_weight (A->seq_al[a], A->seq_al[b], 1,F);//f,p1, p2, p3, fsa_mode);
	  for (c=0; c<A->len_aln; c++)
	    {
	      if (A->seq_al[a][c]=='-')continue;
	      C[a][c]+=w[c];
	      T[a][c]++;
	    }
	  vfree (w);
	}
    }



  for (a=0; a<A->nseq; a++)
    {
      for (b=0; b<A->len_aln; b++)
	{
	  int c;
	  c=A->seq_al[a][b];
	  if ( c=='-');
	  else if (T[a][b]==0);
	  else
	    {
	      int r;
	      r=(C[a][b]*10)/T[a][b];
	      r=(r==10)?9:r;
	      if (!F->mode || strm (F->mode, "number"))
		B->seq_al[a][b]='0'+r;
	      else if ( F->mode && (strm (F->mode, "unalign") ||strm (F->mode, "unalign2")))
		B->seq_al[a][b]='0'+r;
	      else if ( F->mode && strm (F->mode, "lower") )
		{
		  if (r<=F->t)B->seq_al[a][b]=tolower (B->seq_al[a][b]);
		  else B->seq_al[a][b]=toupper (B->seq_al[a][b]);
		}
	    }
	}
    }

  if (F->mode && strm (F->mode, "unalign"))
    {
      A=unalign_aln (A, B, F->t);
      free_aln (B);
      B=copy_aln (A, NULL);
    }
  else if (F->mode && strm (F->mode, "unalign2"))
    {
      A=unalign_aln_2 (A, B, F->t);
      free_aln (B);
      B=copy_aln (A, NULL);
    }



  free_int (C, -1);
  free_int (T, -1);

  return B;
}

char **pw_aln2clean_pw_aln_fsa1 (char ** aln, OveralnP *F);
char **pw_aln2clean_pw_aln_fsa2 (char ** aln, OveralnP *F);

int  * pw_aln2clean_aln_weight ( char *seq1, char *seq2, int w, OveralnP *F)
{
  char **aln;
  int *weight;
  int l, a;

  if ( (l=strlen (seq1)) !=strlen (seq2))
    {
      HERE ("\n%s\n%s\n", seq1, seq2);
      printf_exit ( EXIT_FAILURE, stderr, "\nERROR: Comparing unaligned sequences [FATAL:%s]", PROGRAM);

    }

  aln=declare_char (2, l+1);
  sprintf ( aln[0], "%s", seq1);
  sprintf ( aln[1], "%s", seq2);


  aln=pw_aln2clean_pw_aln (aln, F);

  weight=(int*)vcalloc (l+1, sizeof (int));
  for (a=0; a<l; a++)
    {
      if ( aln[0][a] || seq1[a]=='x' || seq1[a]=='X' || seq2[a]=='x' || seq2[a]=='X')weight[a]=w;
    }
  free_char (aln, -1);

  return weight;
}


char **pw_aln2clean_pw_aln (char ** aln, OveralnP *F)
{

  if ( strm (F->model, "fsa2"))return pw_aln2clean_pw_aln_fsa2 (aln,F);
  else if ( strm (F->model, "fsa1"))return pw_aln2clean_pw_aln_fsa1 (aln,F);
  else return pw_aln2clean_pw_aln_fsa1 (aln,F);
}

char **pw_aln2clean_pw_aln_fsa2 (char ** aln, OveralnP *FO)
{
  int a, b, c, d, l, id;
  int c1, c2, e0, e1,tb, obs;
  int T0, T1,T2;
  int **mat, **tran, **p, **t, *s, *ids;
  int ns, ps, cs;
  int S, M1, M2, m1, m2,B1, B2,G1,G2, K;
  int F=-9999999;
  int MID_EXON_FACTOR=50;
  int best;
  static int **smat;
  int model_type=1;
  int *translate;

  if ( getenv ("MID_EXON_FACTOR"))MID_EXON_FACTOR=atoi (getenv ("MID_EXON_FACTOR"));



  if (!smat)smat=read_matrice ( "blosum62mt");

  l=strlen (aln[0]);

  if ( l!=strlen (aln[1]))
    {
      printf_exit ( EXIT_FAILURE, stderr, "\nERROR: unaligned strings");
    }



  s=(int*)vcalloc (l, sizeof (int));
  ids=(int*)vcalloc (l, sizeof (int));

  //record the id level of each posotion
  for (b=0; b<l; b++)
    {
      c1=tolower(aln[0][b]);c2=tolower(c2=aln[1][b]);

      if (c1=='-' || c2=='-' || c1=='X' || c2=='X' || c1!=c2)ids[b]=0;
      else ids[b]=1;
    }

  //record the state of each position: M, m, T, gap
  for (id=0,b=0,a=0;a<l; a++)
    {
      c1=aln[0][a];c2=aln[1][a];
      if (islower (c1))s[a]=3;
      else if (c1=='-' || c2=='-' || c1=='X' || c2=='X')s[a]=2;
      else
	{
	  int sc;
	  sc=smat[c1-'A'][c2-'A'];
	  if (sc>=2){id++; s[a]=1;}
	  else {s[a]=0;}
	  b++;
	}
    }

  if (b==0)
    {
      vfree(s);vfree (ids);
      return aln;
    }



  FO->p1=(FO->p1==0)?5:FO->p1;
  FO->p2=(FO->p2==0)?15:FO->p2;
  FO->p3=(FO->p3==0)?0:FO->p3;
  FO->p4=(FO->p4==0)?100:FO->p4;


  T1=100*(float)id/(float)b;
  T2=(FO->f==0)?30:T1*(float)((float)FO->f/(float)100);
  T2=MAX(T2,20);

  //0: unaligned
  //1: aligned
  //2: gap
  //3: exon boundary

  ns=0;
  S=ns++;
  M1=ns++;//1 matched  aligned
  m1=ns++;//2 mmatched aligned
  M2=ns++;//3 matched  unaligned
  m2=ns++;//4 mmatched unaligned
  B1=ns++;//5 transition aligned
  B2=ns++;//6 transition unaligned

  mat=declare_int (ns, 4);
  tran=declare_int (ns, ns);
  p=declare_int (l+1, ns);
  t=declare_int (l+1, ns);

  //emission Values
  mat[M1][0]=F; //non id
  mat[M1][1]=T1;//id
  mat[M1][2]=0; //gap
  mat[M1][3]=F; //transition

  mat[M2][0]=F;
  mat[M2][1]=T2;
  mat[M2][2]=0;
  mat[M2][3]=F;

  mat[m1][0]=100-T1;
  mat[m1][1]=F;
  mat[m1][2]=0;
  mat[m1][3]=F;

  mat[m2][0]=100-T2;
  mat[m2][1]=F;
  mat[m2][2]=0;
  mat[m1][3]=F;

  mat[B1][0]=F;
  mat[B1][1]=F;
  mat[B1][2]=F;
  mat[B1][3]=0;

  mat[B2][0]=F;
  mat[B2][1]=F;
  mat[B2][2]=F;
  mat[B2][3]=0;

  //transition values
  tran[S][m1]=0;
  tran[S][m2]=0;
  tran[S][M1]=0;
  tran[S][M2]=0;
  tran[S][B1]=0;
  tran[S][B2]=0;


  tran[M1][m1]= 0;
  tran[M1][m2]=-FO->p4;
  tran[M1][M1]=+FO->p2;
  tran[M1][M2]= F;
  tran[M1][S ]= F;
  tran[M1][B1]= 0;
  tran[M1][B2]=-FO->p1;

  tran[M2][m1]= F;
  tran[M2][m2]=+FO->p3;
  tran[M2][M1]= F;
  tran[M2][M2]= 0;
  tran[M2][S] = F;
  tran[M2][B1]= F;
  tran[M2][B2]= 0;


  tran[m1][m1]= 0;
  tran[m1][m2]= F;
  tran[m1][M1]= 0;
  tran[m1][M2]= F;
  tran[m1][S] = F;
  tran[m1][B1]= 0;
  tran[m1][B2]=-FO->p1;

  tran[m2][m1]=  F;
  tran[m2][m2]=  0;
  tran[m2][M1]= -FO->p4;
  tran[m2][M2]= +FO->p3;
  tran[m2][S] =  F;
  tran[m2][B1]=  F;
  tran[m2][B2]=  0;

  tran[B1][m1]=  0;
  tran[B1][m2]=  F;
  tran[B1][M1]=  0;
  tran[B1][M2]=  F;
  tran[B1][S]=   F;
  tran[B1][B1]=  F;
  tran[B1][B2]=  F;

  tran[B2][m1]= -FO->p1;
  tran[B2][m2]=  0;
  tran[B2][M1]= -FO->p1;
  tran[B2][M2]=  0;
  tran[B2][S]=   F;
  tran[B2][B1]=  F;
  tran[B2][B2]=  F;

  translate=(int*)vcalloc (ns, sizeof (int));
  translate[M1]=1;
  translate[m1]=1;
  translate[M2]=0;
  translate[m2]=0;
  translate[B1]=1;
  translate[B2]=0;

  for (a=1;a<=l; a++)
    {
      obs=s[a-1];

      for (cs=0; cs<ns; cs++)
	{
	  for (ps=0; ps<ns; ps++)
	    {
	      c=p[a-1][ps]+mat[cs][obs]+tran[ps][cs];
	      if (ps==0 || c>=best){t[a][cs]=ps;best=p[a][cs]=c;}
	    }

	}
    }


  for (a=0; a<ns; a++)
    {
      if (a==0 || p[l][a]>=best){tb=a;best=p[l][a];}
    }

  for (a=l; a>0; a--)
    {
      int v;
      int p2;

      p2=a-1;
      aln[0][p2]=aln[1][p2]=translate[tb];
      tb=t[a][tb];

    }

  free_int (p, -1);
  vfree(s);
  free_int (t, -1);
  free_int (mat, -1);
  free_int (tran, -1);
  vfree (translate);
  return aln;
}
char **pw_aln2clean_pw_aln_fsa1 (char ** aln, OveralnP *FO)
{
  int a, b, c, d, l, id;
  int c1, c2, e0, e1,tb, obs;
  int T0, T1,T2;
  int **mat, **tran, **p, **t, **s;
  int ns, ps, cs;
  int S, M1, M2, m1, m2, K;
  int F=-9999999;
  int best;
  static int **smat;
  int *translate;


  if (!smat)smat=read_matrice ( "blosum62mt");

  l=strlen (aln[0]);

  if ( l!=strlen (aln[1]))
    {
      printf_exit ( EXIT_FAILURE, stderr, "\nERROR: unaligned strings");
    }


  s=declare_int (l+1, 2);
  for (id=0,b=0,a=0;a<l; a++)
    {
      c1=aln[0][a];c2=aln[1][a];

      if ( c1=='-' || c2=='-' || c1=='x' || c1=='X' || c2=='x' || c2=='X')continue;
      else
	{
	  int sc;
	  sc=smat[c1-'A'][c2-'A'];
	  if (sc>=2){id++; s[b][0]=1;}
	  else {s[b][0]=0;}
	  s[b][1]=a;
	  b++;

	}
    }
  if (b==0)
    {
      free_int (s, -1);
      return aln;
    }
  FO->f=(FO->f==0)?30:FO->f;
  FO->p1=(FO->p1==0)?90:FO->p1;
  FO->p2=(FO->p2==0)?15:FO->p2;
  FO->p3=(FO->p3==0)?0:FO->p3;

  l=b;//length of the ungapped aln
  T1=100*(float)id/(float)b;
  T2=FO->f;//T1*f;



  //0: unaligned
  //1: aligned


  ns=0;
  S=ns++;
  M1=ns++;//1 matched  aligned
  m1=ns++;//2 mmatched aligned
  M2=ns++;//3 matched  unaligned
  m2=ns++;//4 mmatched unaligned

  mat=declare_int (ns, 2);
  tran=declare_int (ns, ns);
  p=declare_int (l+1, ns);
  t=declare_int (l+1, ns);


  mat[M1][0]=F;
  mat[M1][1]=T1;

  mat[M2][0]=F;
  mat[M2][1]=T2;

  mat[m1][0]=100-T1;
  mat[m1][1]=F;

  mat[m2][0]=100-T2;
  mat[m2][1]=F;


  tran[S][m1]=0;
  tran[S][m2]=0;
  tran[S][M1]=0;
  tran[S][M2]=0;


  tran[M1][m1]= 0;
  tran[M1][m2]=-FO->p1;// -P;
  tran[M1][M1]=+FO->p2;
  tran[M1][M2]= F;
  tran[M1][S] = F;

  tran[M2][m1]= F;
  tran[M2][m2]=+FO->p3;
  tran[M2][M1]= F;
  tran[M2][M2]= 0;
  tran[M2][S]=  F;

  tran[m1][m1]= 0;
  tran[m1][m2]= F;
  tran[m1][M1]= 0;
  tran[m1][M2]= F;
  tran[m1][S]=  F;

  tran[m2][m1]= F;
  tran[m2][m2]= 0;
  tran[m2][M1]=-FO->p1;
  tran[m2][M2]=+FO->p3;
  tran[m2][S]=  F;

  translate=(int*)vcalloc (ns, sizeof (int));
  translate[M1]=1;
  translate[m1]=1;
  translate[M2]=0;
  translate[m2]=0;
  translate[S]=1;


  for (a=1;a<=l; a++)
    {
      obs=s[a-1][0];

      for (cs=0; cs<ns; cs++)
	{
	  for (ps=0; ps<ns; ps++)
	    {
	      c=p[a-1][ps]+mat[cs][obs]+tran[ps][cs];
	      if (ps==0 || c>=best){t[a][cs]=ps;best=p[a][cs]=c;}
	    }

	}
    }


  for (a=0; a<ns; a++)
    {
      if (a==0 || p[l][a]>=best){tb=a;best=p[l][a];}
    }
  for (a=l; a>0; a--)
    {
      int p2=s[a-1][1];
      aln[0][p2]=aln[1][p2]=translate[tb];

      tb=t[a][tb];
    }


  free_int (p, -1);
  free_int (s, -1);
  free_int (t, -1);
  free_int (mat, -1);
  free_int (tran, -1);
  vfree (translate);
  return aln;
}
float* analyze_overaln ( Alignment *iA, Alignment *iB, char *mode, int filter, int f, int p1,int p2, int p3)
{
  Alignment *C, *D;
  Alignment *A, *B;
  OveralnP *F;

  F=(OveralnP*)vcalloc (1, sizeof (OveralnP));
  F->p1=p1;
  F->p2=p2;
  F->p3=p3;
  F->f=f;
  F->t=filter;
  sprintf (F->mode, "%s", mode);


  float *r;
  A=copy_aln (iA, NULL);
  B=copy_aln (iB, NULL);

  C=aln2gap_cache (A,0);
  A=filter_aln_upper_lower (A, C, 0, 0);
  D=aln2clean_pw_aln (B, F);
  r=aln2pred (A,D,mode);
  free_aln (C);
  free_aln (D);
  free_aln (A);
  free_aln (B);
  return r;
}
float* aln2pred ( Alignment *A, Alignment*B, char *mode)
{
  int a, b, c, d, i, l, salp, s, n;
  static char **list, *buf1, *buf2, *alp, *alp_lu;
  static int ***r;
  int T, N;
  int fp, fn, tn, tp;
  int tfp, tfn, ttn, ttp;
  float sp, sn, sen2, best, result;
  int print=1;
  float *fresult;

  fresult=(float*)vcalloc ( 3, sizeof (float));

  if ( mode && strstr (mode, "case"))
    {
      A=aln2case_aln (A,"u","l");
      B=aln2case_aln (B,"u","l");
    }

  if (mode && strstr (mode, "printaln"))
    {
      Sequence *S;
      Alignment *C;
      S=aln2seq (A);
      C=copy_aln (B, NULL);
      for (a=0; a<B->nseq; a++)
	{
	  i=name_is_in_list (C->name[a], S->name, S->nseq, 100);
	  if ( i==-1)
	    for (b=0; b<C->len_aln; b++) C->seq_al[a][b]='-';
	  else
	    for (d=0,b=0; b<C->len_aln; b++)
	      {
		if ( !is_gap (C->seq_al[a][b]))
		  {
		    if (C->seq_al[a][b]==S->seq[i][d])C->seq_al[a][b]=toupper(C->seq_al[a][b]);
		    d++;
		  }
	      }
	}
      print_aln (C);
    }

  vfree (alp);vfree (alp_lu);
  alp=(char*)vcalloc ( 256, sizeof (char));
  alp_lu=(char*)vcalloc ( 256, sizeof (char));

  for (c=0; c<2; c++)
    {
      Alignment *AL;
      AL=(c==0)?A:B;
      for (salp=0,a=0; a<AL->nseq; a++)
	{
	  for (b=0; b<AL->len_aln; b++)
	    {
	      c=AL->seq_al[a][b];
	      if (!is_gap(c) && !alp[c])
		{
		  salp++;
		  alp_lu[salp]=c;
		  alp[c]=salp;
		}
	    }
	}
    }

  vfree (buf1); vfree(buf2);
  buf1=(char*)vcalloc ( A->len_aln+1, sizeof (char));
  buf2=(char*)vcalloc ( B->len_aln+1, sizeof (char));

  free_arrayN ((void **)r, 3);
  r=(int***)declare_arrayN(3, sizeof (int),A->nseq,salp+1,salp+1);
  free_char ( list, -1);
  list=declare_char ( A->nseq, 100);
  for (n=0,a=0; a< A->nseq; a++)
    {
      for ( b=0; b<B->nseq; b++)
	{
	  if ( strm (A->name[a], B->name[b]))
	    {
	      sprintf ( buf1, "%s", A->seq_al[a]);
	      sprintf ( buf2, "%s", B->seq_al[b]);
	      ungap (buf1); ungap (buf2);
	      if ((l=strlen (buf1))!=strlen (buf2))continue;
	      else
		{
		  sprintf ( list[n], "%s", A->name[a]);
		  for (c=0; c<l; c++)
		    {
		      int c1, c2;
		      c1=buf1[c];
		      c2=buf2[c];
		      r[n][alp[c1]][alp[c2]]++;
		    }
		  n++;
		}
	    }
	}
    }



  for ( s=1; s<=salp; s++)
    {
      char type[4];
      sprintf (type, "_%c_", alp_lu[s]);
      ttp=ttn=tfp=tfn=0;
      for (a=0; a<n; a++)
	{
	  tp=tn=fp=fn=0;
	  for (b=1; b<=salp; b++)
	    {
	      for (c=1; c<=salp; c++)
		    {
		      if ( b==s && c==s)     tp+=r[a][b][c];
		      else if ( b==s && c!=s)fn+=r[a][b][c];
		      else if ( b!=s && c==s)fp+=r[a][b][c];
		      else if ( b!=s && b!=s)tn+=r[a][b][c];
		    }

	    }

	  ttp+=tp;
	  ttn+=tn;
	  tfp+=fp;
	  tfn+=fn;
	  rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);
	  if ( mode && strstr (mode, "printstat"))fprintf ( stdout, ">%s S=%c sp=%6.2f sn=%6.2f sen2=%6.2f best=%6.2f\n", list[a],alp_lu[s],sp, sn, sen2, best);
	}

      rates2sensitivity (ttp, ttn, tfp, tfn, &sp, &sn, &sen2, &best);
      if (mode && strstr (mode, "printstat"))fprintf ( stdout, ">TOT S=%c sp=%6.2f sn=%6.2f re=%6.2f best=%6.2f\n", alp_lu[s],sp, sn, sen2, best);

      if ( mode && strstr (mode, type))
	{
	  fresult[0]=sn;
	  fresult[1]=sp;
	  fresult[2]=sen2;
	}
    }
  return fresult;
}

Alignment * mark_exon_boundaries  (Alignment *A, Alignment *E)
{
  char *buf, *buf2;
  int a, b, c, i, l;

  buf2=(char*)vcalloc ( E->len_aln+1, sizeof (char));
  buf =(char*)vcalloc ( E->len_aln+1, sizeof (char));

  for (a=0; a< A->nseq; a++)
    {
      i=name_is_in_list (A->name[a], E->name, E->nseq, 100);
      if ( i==-1) continue;
      sprintf (buf, "%s", E->seq_al[i]);
      ungap (buf);
      l=strlen (buf);
      //clean buf2
      for (c=0, b=0; b<l; b++)if (buf[b]!='o' && buf[b]!='b' && buf[b]!='j')buf2[c++]=toupper(buf[b]);
      buf2[c]='\0';

      //lowercase the boundaries of buf2;
      for ( c=0,b=0; b<l; b++)
	{
	  //ENSEMBL: o: 0, b:1 j:2
	  if (buf[b]=='b' || buf[b]=='o' && c>=1)buf2[c-1]=tolower(buf2[c-1]);
	  else if (buf[b]=='j' &&c<l)buf2[c+1]=tolower(buf2[c+1]);
	  else c++;
	}

      for (c=0,b=0; b<A->len_aln; b++)
	{
	  if (!is_gap(A->seq_al[a][b]))
	    {
	      A->seq_al[a][b]=buf2[c++];
	    }
	}
    }
  vfree (buf);
  vfree (buf2);
  return A;
}
//simple_trimseq2
//Creates Clusters
//In each cluser there is a path between every pair of sequence
//A path is made of edges connecting tow nodes with w>min_sim
static int trimseq2_getnext (int *gl,int g,int n,int *used, int ***sim,int minsim);
int ** simple_trimseq2 (int n, int **sim, int minsim)
{
  //used[se_index]=group_index (Groug Index: 1..N)
  //gl[GroupIndex][0]=group size
  //gl[GroupIndex][1..GrouSize]: sequences (indexed 0..N-1)
  //gl[0][0]: ngroups;


  int ***ssim;
  int *  used;
  int tot, a, b,c, s1,s2, ng, g;
  int ** gl;
  int ** gsim;

  used=(int*)vcalloc (n, sizeof (int));
  ssim=(int***)vcalloc (n, sizeof (int**));
  for (s1=0; s1<n; s1++)
    {
      ssim[s1]=declare_int (n,2);
      for (s2=0; s2<n; s2++)
	{
	  ssim[s1][s2][0]=s2;
	  ssim[s1][s2][1]=sim[s1][s2];
	}
      sort_int_inv(ssim[s1], 2,1,0,n-1);
    }

  ng=0;
  gl=declare_int (n+1,1);
  gsim=declare_int (n,2);

  tot=0;
  while (tot<n)
    {
      ng++;
      vfree(gl[ng]);gl[ng]=(int*)vcalloc (n+1, sizeof (int));
      for (s1=0; s1<n; s1++)
	if (!used[s1])
	  {
	    used[s1]=ng;
	    gl[ng][++gl[ng][0]]=s1;
	    break;
	  }

      while ((s1=trimseq2_getnext(gl[ng],ng,n,used,ssim,minsim)))
	{
	  used[s1-1]=ng;
	  gl[ng][++gl[ng][0]]=s1-1;
	}
      tot+=gl[ng][0];
    }

  gl[0][0]=ng;

  for (g=1; g<=ng;g++)
    {

      for (b=1; b<=gl[g][0]; b++){gsim[b-1][0]=gl[g][b];gsim[b-1][1]=0;}
      for (b=1; b<=gl[g][0]; b++)
	for (c=1; c<=gl[g][0]; c++)
	  {
	    int s1=gl[g][b];
	    int s2=gl[g][c];
	    if (s1!=s2)
	      {
		gsim[b-1][1]+=sim[s1][s2];
		gsim[c-1][1]+=sim[s1][s2];
	      }
	  }
      sort_int_inv (gsim, 2, 1, 0, gl[g][0]-1);

      for (b=1; b<=gl[g][0]; b++)gl[g][b]=gsim[b-1][0];
    }
  free_int (gsim, 2);
  free_arrayN((void*)ssim,3);


  return gl;
}

static int trimseq2_getnext (int *gl,int g,int n,int *used, int ***sim,int minsim)
{
  int s1, s2,a,b;
  int bseq=-1;
  int bscore=-1,cscore;

  for (a=1; a<=gl[0]; a++)
    {
      s1=gl[a];
      for (b=0; b<n; b++)
	{
	  s2     =sim[s1][b][0];
	  cscore =sim[s1][b][1];
	  if (!used[s2])
	    {
	      if (cscore>minsim && cscore>bscore)
		{

		  bseq  =s2;
		  bscore=cscore;
		}
	      break;
	    }
	}
    }
  return bseq+1;
}
void add_msa (Alignment *A,int seq, char *lu, char *file, char *mode);
char *msa2master_seq (Alignment *A, int seq, char *master);

char *msa2master_seq (Alignment *A, int seq, char *master)
{
  
  static int *aa=(int*)vcalloc(256, sizeof (int));
  int c, s,a;
  master=(char *) vreallocg (master,(A->len_aln+1)*sizeof (char),NOMEMSET,NORESIZE);
  sprintf (master, "%s", A->seq_al[seq]);
  
  for ( c=0; c<A->len_aln; c++)
    {
      if ( is_gap (A->seq_al[seq][c]))
	{
	  int best_naa=0;
	  char best_aa=0;
	  for (a=0; a<256; a++)aa[a]=0;
	  for (s=0; s<A->nseq; s++)
	    {
	      
	      char r=A->seq_al[s][c];
	      if (r!='1')//allow gaps in the consensus -> they get turned into X
		{
		  aa[r]++;
		  if (aa[r]>best_naa){best_naa=aa[r]; best_aa=r;}
		}
	    }

	  master[c]=(best_aa=='-')?'x':best_aa;
	}
    }
  
  return master;
}


void add_msa (Alignment *A,int seq, char *lu, char *file, char *mode)
{
  FILE *fp=vfopen(file,mode);
  int s, c,col;
  int len=strlen (lu);
  
  for (s=0; s<A->nseq; s++)
    {
      if ( s!=seq)
	{
	  fprintf ( fp, ">%s\n", A->name [s]);
	  
	  for (col=0,c=0; c<len; c++)
	    {
	      char r=lu[c];
	      if (r!='-')fprintf (fp, "%c",A->seq_al[s][col++]);
	      else fprintf (fp, "-");
	    }
	  fprintf ( fp, "\n");
	}
    }
  vfclose (fp);
  return;
}



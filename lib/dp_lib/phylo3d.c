
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <search.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

int compare4phylo3d(const void *a, const void *b);
// A structure to hold labeled data
typedef struct LabeledData{
    double score; // Model's score or probability for being positive
    int label;    // Actual label (0 or 1)
};
typedef struct LabeledData LabeleddData;


//// AUC for phylo3d
int splits2npp(int n, float **ns);
float splits2auc (int n, float **ns);
double calculateAUC(LabeledData *data, int n);
void clean_pp (float ***** ns, int ***nn);
int set_pp (char *family, int nseq,float ***** sl, char ***** splits, int ***nn, char *mode);

Alignment * multistrap  (Alignment *inA,Constraint_list *CL,int nbsF,char **inbsF) 
{
  //will only keep the sequences whose structure is known
  int a,b;
  char *meF,*mlF,*imdF, **bsF;
  Alignment *A=NULL;
  p3D *D;
  NT_node *TREE;
  int nrep;
  char *multistrap_mode=NULL;
  
  meF=mlF=imdF=NULL;
  bsF =(char**)vcalloc   (nbsF, sizeof (char*));
  TREE=(NT_node*)vcalloc (nbsF, sizeof (NT_node));
  bsF[0]=inbsF[0];
		 

   if (getenv ("REPLICATES_4_TCOFFEE"))
    {
      if ( strm (getenv ("REPLICATES_4_TCOFFEE"), "columns"))nrep=inA->len_aln;
      else nrep=atoigetenv ("REPLICATES_4_TCOFFEE");
    }
   else
     nrep=100;

   if (getenv ("MULTISTRAP_MODE"))multistrap_mode=csprintf (NULL,"%s",getenv("MULTISTRAP_MODE"));
   else multistrap_mode=csprintf (NULL, "average");
   
  //1 - generate the required replicates and possibly the reference trees
  for (a=1; a<nbsF; a++)
    {
      if (check_file_exists(bsF[a]))
	bsF[a]=inbsF[a];
      else if (strm (inbsF[a], "me"))
	{
	  bsF[a]=vtmpnam  (NULL);
	  meF=vtmpnam (NULL);
	  
	  aln2fastme_treeF (inA, nrep,meF,bsF[a]);
	}
      else if (strm (inbsF[a], "ml"))
	{
	  bsF[a]=vtmpnam  (NULL);
	  mlF=vtmpnam (NULL);
	  aln2iqtree_treeF (inA, nrep,mlF,bsF[a]);
	}
      else if (strm (inbsF[a], "imd"))
	{
	  bsF[a]=vtmpnam  (NULL);
	  imdF=vtmpnam (NULL);
	  
	  A=aln2trim3d(inA, CL);
	  
	  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");  
	  D=fill_p3D(A, CL);
	  D->N=filter_columns_with_gap  (D->col, A, D->max_gap);
	  if (D->maxd<MY_EPSILON)D->maxd=scan_maxd(D);
	  D->N=filter_columns_with_dist (A,D->pos, D->col,D->dm3d,D->maxd);
	  makerep (D,0);
	  aln2dm(D,A);  
	  fprintf ( stderr, "---- generate imd tree [%d bootstrap cycles]", nrep);
	  dist2fastme_treeF (D->dm, A->name, A->nseq, imdF);
	  
	  for (b=0; b<nrep; b++)
	    {
	      static char *treeF=vtmpnam(NULL);
	      makerep (D,2);
	      aln2dm(D,A);
	      dist2fastme_treeF (D->dm, A->name, A->nseq, treeF);
	      string2file(bsF[a],"a" "%s", file2string (treeF));
	    }
	}
    }
    
 
  //2 - collect the reference tree on which the bs will be estimated
  if (check_file_exists (bsF[0]));
  else if (strm (bsF[0], "ml"))
    {
      if   (!mlF)
	{
	  bsF[0]=vtmpnam(NULL);
	  aln2fastme_treeF (inA,0,bsF[0],NULL);
	}
      else bsF[0]=mlF;
    }
  else if ( strm (bsF[0],"me"))
    {
      if   (!meF)
	{
	  bsF[0]=vtmpnam(NULL);
	  aln2iqtree_treeF (inA,0,bsF[0],NULL);
	}
      else bsF[0]=meF;
    }
  else if ( strm (bsF[0],"imd"))
    {
      if   (!imdF)
	{
	  bsF[0]=vtmpnam(NULL);
	  A=aln2trim3d(inA, CL);
	  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");  
	  D=fill_p3D(A, CL);
	  D->N=filter_columns_with_gap  (D->col, A, D->max_gap);
	  if (D->maxd<MY_EPSILON)D->maxd=scan_maxd(D);
	  D->N=filter_columns_with_dist (A,D->pos, D->col,D->dm3d,D->maxd);
	  makerep (D,0);
	  aln2dm(D,A);  
	  dist2fastme_treeF (D->dm, A->name, A->nseq,bsF[0]);
	}
      else bsF[0]=imdF; 
    }
  //This is the reference tree in which the new BS will be added
  TREE[0]=main_read_tree(bsF[0]);

  for ( a=1; a<nbsF; a++)
    {
      Sequence  *BTS=NULL;
      Alignment *BTA=NULL;
      static char *treeF=vtmpnam (NULL);
      static char *tmpF =vtmpnam (NULL);
      NT_node t=main_read_tree (bsF[0]);
      reset_bs2(t);
      string2file (treeF,"w", "%s",tree2string(t));
      
      printf_system ("cat %s %s > %s",treeF, bsF[a], tmpF); 
      BTS=get_treelist(tmpF);
      BTA=seq2aln(BTS,NULL,NO_PAD);
      treelist2node_support(BTA);
      string2file (treeF, "w", "%s", BTA->seq_al[0]);
      TREE[a]=main_read_tree (treeF);
    }
  
  //The combined BS gets added to TREE[0]
  combine_bsN(nbsF,TREE, multistrap_mode);

  fprintf     ( stdout, "\n");
  //Print the other way round so that the combined BS is at the bottom
  for (a=nbsF-1; a>=0; a--)
    {
      if (a==0)fprintf ( stdout, "# Reference Tree [%s] with bootsraps added by combining the others, using the %s multistrap mode\n", inbsF[a], multistrap_mode);
      else fprintf     ( stdout, "# [%s] tree\n", inbsF[a]);
      fprintf ( stdout, "%s", tree2string (TREE[a]));
    }
  exit (EXIT_SUCCESS);
}

Alignment * multistrap_old  (Alignment *inA,char *RTF,Constraint_list *CL)
{
  //will only keep the sequences whose structure is known
  static char *treeF=vtmpnam (NULL);
  static char *tree_listF=vtmpnam (NULL);
  NT_node RT_INITIAL;
  NT_node RT_COMBINED;
  NT_node RT_PHYLO3D;
  Alignment *A=NULL;
  p3D *D;
  Sequence  *BTS=NULL;
  Alignment *BTA=NULL;
  FILE *fp;
  int a;

  A=aln2trim3d(inA, CL);
  D=fill_p3D(A, CL);
  treeF=vtmpnam (NULL);
  if (RTF &&  (strm (RTF, "iqtree") && !check_file_exists("iqtree")))
    {
      aln2iqtree_treeF (inA, D->replicates,treeF,NULL);
    }
  else if (!RTF || (strm (RTF, "fastme") && !check_file_exists("fastme")))
    {
      aln2fastme_treeF (inA, D->replicates,treeF,NULL);
    }
  else
    {
      printf_system ("cp %s %s", RTF, treeF);
    }
 

  RT_INITIAL=main_read_tree (treeF);
  RT_COMBINED=main_read_tree (treeF);
  RT_PHYLO3D=main_read_tree(treeF);
  reset_bs2(RT_PHYLO3D);
  
  
 

  
  
  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");  
  //Filter columns containing too many gaps
  D->N=filter_columns_with_gap  (D->col, A, D->max_gap);
  
  //Filter columns with distances
  if (D->maxd<MY_EPSILON)D->maxd=scan_maxd(D);
  D->N=filter_columns_with_dist (A,D->pos, D->col,D->dm3d,D->maxd);
  makerep (D,0);
  aln2dm(D,A);
  
  dist2fastme_treeF (D->dm, A->name, A->nseq, treeF);
  string2file(tree_listF,"w", "%s", tree2string (RT_PHYLO3D));
  
  fprintf ( stderr, "---- generate fastme/phylo3D bootstrap replicates [%d bootstrap cycles]", D->replicates);
  for (a=1; a<=D->replicates; a++)
    {
 
      makerep (D,2);
      aln2dm(D,A);
      dist2fastme_treeF (D->dm, A->name, A->nseq, treeF);
      string2file(tree_listF,"a" "%s", file2string (treeF));
    }
  fprintf ( stderr, "[DONE]\n");
  
  BTS=get_treelist(tree_listF);
  BTA=seq2aln(BTS,NULL,NO_PAD);
  BTA=treelist2node_support (BTA);

  string2file (treeF, "w", "%s", BTA->seq_al[0]);
  RT_PHYLO3D=main_read_tree (treeF);

  combine_bs(RT_COMBINED, RT_PHYLO3D, D->multistrap_mode);

  fprintf (stderr, "Line 1: %s tree with Combined bootstraped [%s] - Line 2 %s tree with phylo3D bootstraps Line 3 %s tree with original boostrap\n", RTF, (D->multistrap_mode)?D->multistrap_mode:"geometric", RTF, RTF);

  
  fprintf (stdout,"%s%s%s", tree2string (RT_INITIAL),tree2string (RT_PHYLO3D),tree2string (RT_COMBINED));
  exit (EXIT_SUCCESS);
}
/*
Old function, now moved in util_make_tree
int aln2fastme_treeF  (Alignment *A, int bs, char *treeF)
{
  static char *alnF=vtmpnam(NULL);
  static char *odm=vtmpnam (NULL);
  fprintf ( stderr, "---- generate fastme tree [%d bootstrap cycles]", bs);
  if (!check_program_is_installed ("fastme",NULL,NULL,"http://www.atgc-montpellier.fr/fastme",NO_REPORT))printf_exit ( EXIT_FAILURE,stderr, "\nERROR: fastme must be installed [FATAL]");;
  output_phylip_aln (alnF, A, "w");
  
  printf_system ("fastme -i %s  -o %s -m BioNJ -p LG -g 1.0 -s -n -z 5 -b %d -B bst -O %s >/dev/null 2>/dev/null", alnF, treeF,bs,odm);
  fprintf ( stderr, "[DONE]\n");
  return 1;
}

int dist2fastme_treeF (double **dm, char **name,int nseq, char *treeF)
{
  static char *dmF=vtmpnam(NULL);
  FILE *fp;
  int s1, s2;

  if (!check_program_is_installed ("fastme",NULL,NULL,"http://www.atgc-montpellier.fr/fastme",NO_REPORT))printf_exit ( EXIT_FAILURE,stderr, "\nERROR: fastme must be installed [FATAL]");
  
  fp=vfopen (dmF, "w");
  fprintf ( fp, "%d \n", nseq);
  for ( s1=0; s1<nseq;s1++)
    {
      fprintf (fp, "%s ",name[s1]);
      for (s2=0; s2<nseq; s2++)
	fprintf (fp, "%6.3f ", (float)((double)dm[s1][s2])/(float)100);
      fprintf (fp, "\n");
    }
  fprintf (fp, "\n");
  vfclose (fp);
  
  printf_system ("fastme -i %s -g 1.0 -s -n -z 5 -o %s >/dev/null 2>/dev/null", dmF,treeF);
  
  return 1;
}
*/

  
Alignment *phylo3d (Alignment *inA, Constraint_list *CL)
{
  //will only keep the sequences whose structure is known
  Alignment *A=aln2trim3d(inA, CL);
  p3D *D=fill_p3D(A, CL);
  int a;
  
  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");  
  
  
  
  //Filter columns containing too many gaps
  D->N=filter_columns_with_gap  (D->col, A, D->max_gap);
  
  //Filter columns with distances
  if (D->maxd<MY_EPSILON)D->maxd=scan_maxd(D);

  D->N=filter_columns_with_dist (A,D->pos, D->col,D->dm3d,D->maxd);
      
  //Compute the first tree or DM
  makerep(D,0);//The first replicates is always the original column list;
  aln2dm(D,A);
  addtree(D,A);
  
  for (a=1; a<=D->replicates;a++)
    {
      //0-> re-use all the columns
      //1-> draw with resampling among ALL the pairs
      //2-> draw with resampling among ALL the columns
      makerep(D,2);
      aln2dm(D,A);
      addtree(D,A);
    }
  inA->Tree=A->Tree;
  return inA;    
}

double scan_maxd (p3D *D)
{
  NT_node RT, T;
  float rf, brf;
  float bmaxd=1500;
  Sequence *S=aln2seq(D->A);
  Alignment *A=D->A;
  int strict;
  int bstrict=0;
  int bnsites=0;
  int scan3D_min=(getenv("min_maxd_4_TCOFFEE"))?atof(getenv("min_maxd_4_TCOFFEE")):5;
  int scan3D_max=(getenv("max_maxd_4_TCOFFEE"))?atof(getenv("max_maxd_4_TCOFFEE")):D->extremed;
  int start, end;
  int a;

  if (!getenv ("REFERENCE_TREE"))
    {
      RT= compute_cw_tree(D->A);
      if (verbose())fprintf ( stderr, "\n!#ref_tree= [cw] -- %s",tree2string (RT));
    }
  else if ( getenv ("REFERENCE_TREE") && isfile(getenv ("REFERENCE_TREE")))
    {
      RT=main_read_tree (getenv ("REFERENCE_TREE"));
      if (verbose())fprintf ( stderr, "\n!#ref_tree=%s", getenv ("REFERENCE_TREE"));
    }
  else
    {
      RT= compute_cw_tree(D->A);
      vfclose (tree2file (RT, S, "newick",vfopen(getenv ("REFERENCE_TREE"), "w")));
      if (verbose())fprintf ( stderr, "\n!# ref_tree=%s [cw]", getenv ("REFERENCE_TREE"));
    }
  
  if (verbose())fprintf ( stderr, "\n!# scan3D_max=%d", scan3D_max);
  
  brf=0;
  bmaxd=D->extremed;
 
   
  if (getenv ("soft_maxd_4_TCOFFEE")){start=0; end=0;}
  else if ( getenv ("strict_maxd_4_TCOFFEE")){start=1; end=1;}
  else{start=0; end=1;}
  for (strict=start; strict<=end; strict++)
    {
      for (a=scan3D_min; a<scan3D_max; a++)
	{
	  static char *treeF=vtmpnam (NULL);
	  D->maxd=(double)a*100;
	  cputenv ("strict_maxd_4_TCOFFEE=%d", strict);
	  makerep(D,0);
	  filter_columns_with_dist (A,D->pos, D->colrep, D->dm3d, D->maxd);
	  
	  if (aln2dm(D,A))
	    {
	      dist2nj_tree (D->dm, A->name, A->nseq, treeF);
	      T=main_read_tree(treeF);
	      rf=simple_tree_cmp(RT,T, S, 1);
	      
	      if (verbose())fprintf ( stderr, "\n!# scan        : +maxd %3d %-12s ==> RF vs reftree %5.2f Nsites: %5d", a, (strict)?"+strict_maxd":"+soft_maxd", (float) 100-rf,D->nsites);
			      
	      if ((!getenv ("first_maxd_4_TCOFFEE") && (rf>brf || (rf==brf && D->nsites>bnsites) || (rf==brf && (double)a>=bmaxd))) || ((rf>brf)))
		{
		  brf=rf;
		  bmaxd=(double)a;
		  bstrict=strict;
		  bnsites=D->nsites;
		  fprintf (stderr, "***");
		}
	    }
	  else
	    {
	      if (verbose())fprintf ( stderr, "\n!# Threshold = %3d Angstrom ==> Missing Values in the distance matrix", a) ;
	    }
	}
    }
  if (brf>0)
    {
      if (verbose())fprintf ( stderr, "\n!# scan result : +maxd %3d %-12s ==> RF vs reftree %5.2f Nsites: %5d\n", (int)bmaxd, (strict)?"+strict_maxd":"+soft_maxd", (float) 100-brf,bnsites);
    }
  else
    {
      if (verbose())fprintf ( stderr, "\n!# WARNING -- Missing Values -- Could not find any suitable threshold - Use max value +maxd %d Angstrom", (int)bmaxd);
    }
  free_sequence (S, -1);
  cputenv ("strict_maxd_4_TCOFFEE=%d", bstrict);
  return bmaxd*100;//in picometers
}

Alignment *phylo3d_gt (Alignment *inA, Constraint_list *CL)// phylo3D Guide tree: pw distances will be estimated from unaligned sequences
{
  Alignment *A=aln2trim3d(inA, CL);
  p3D *D=fill_p3D(A, CL);
  Sequence *S;
  Alignment *B;
  
  int n, s1, s2, tot, l;
  int **pos;
  char *align_method=getenv ("align_method_4_TCOFFEE");
  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");
  
  //set the parameters
  
  l=2*A->len_aln;
  D->colrep=declare_int ((l*l-l)/2+1,3);
  vfree(D->used_site);D->used_site=NULL;
  free_int(D->used_site_pair, -1);D->used_site_pair=NULL;
  
  S=aln2seq(A);

  if (D->maxd<MY_EPSILON)D->maxd=D->extremed*100+1;


  tot=((A->nseq*A->nseq)-A->nseq)/2;
  for (n=0,s1=0; s1<A->nseq-1; s1++)
    {
      int rs1, rs2;
      int cls1, cls2;
      cls1=name_is_in_hlist (A->name[s1], (CL->S)->name, (CL->S)->nseq);
      for ( s2=s1+1; s2<A->nseq; s2++, n++)
	{
	  int  *buf_pos_s1, *buf_pos_s2;
	  char *buf_s1, *buf_s2;
	  
	  cls2=name_is_in_hlist (A->name[s2], (CL->S)->name, (CL->S)->nseq); 
	  
	  if (!align_method)B=align_two_sequences (S->seq[s1], S->seq[s2], "blosum62mt", -8, -1, "myers_miller_pair_wise");
	  else B=align_two_structures  (CL->S, cls1, cls2,align_method); 
	  
	  pos=aln2pos_simple(B,2);
	  msa2column_list (B,D->colrep);
	  
	  //buffer data
	  buf_s1=A->seq_al[s1];
	  buf_s2=A->seq_al[s2];
	  buf_pos_s1=D->pos[s1];
	  buf_pos_s2=D->pos[s2];
	  
	  //swapp data
	  A->seq_al[s1]=B->seq_al[0];
	  A->seq_al[s2]=B->seq_al[1];
	  D->pos[s1]=pos[0];
	  D->pos[s2]=pos[1];
	  
	  
	  
	  
	  D->dm[s1][s2]=D->dm[s2][s1]=pair2dist(D,s1,s2);
	  if (verbose())output_completion (stderr,n,tot,1, "Guide Tree Computation");
	  
	  //restaure data
	  A->seq_al[s1]=buf_s1;
	  A->seq_al[s2]=buf_s2;
	  D->pos[s1]=buf_pos_s1;
	  D->pos[s2]=buf_pos_s2;
	  free_int (pos, -1);
	  free_aln (B);
	}
    }
  
  addtree (D, A);
  inA->Tree=A->Tree;
  return inA;
}

int free_3dD (p3D*D)
{
  vfree (D->tree_mode);
  free_int (D->pos, -1);
  free_arrayN((void***)D->dm3d,3);
  free_arrayN((void**)D->dm,2);
  free_arrayN ((void**)D->col,2);
  vfree    (D->used_site);
  free_int(D->used_site_pair, -1);
  vfree (D);
  return 1;
}

p3D* fill_p3D (Alignment *A, Constraint_list *CL)
{
  int a,b;
  p3D *D=(p3D*)vcalloc ( sizeof (p3D), 1);
  D->A=A;
  D->CL=CL;
  D->S =CL->S;
  
  if (getenv ("max_gap_4_TCOFFEE")){D->max_gap  =atofgetenv("max_gap_4_TCOFFEE");}
  else D->max_gap=0.5;
  if (getenv ("MULTISTRAP_MODE"))D->multistrap_mode=csprintf (D->multistrap_mode, "%s",getenv("MULTISTRAP_MODE"));
  
  if (getenv ("TREE_MODE_4_TCOFFEE"))D->tree_mode=csprintf (D->tree_mode, "%s",getenv("TREE_MODE_4_TCOFFEE"));
  else D->tree_mode=csprintf (D->tree_mode, "nj");
  if (getenv ("REPLICATES_4_TCOFFEE"))
    {
      if ( strm (getenv ("REPLICATES_4_TCOFFEE"), "columns"))D->replicates=A->len_aln;
      else D->replicates=atoigetenv ("REPLICATES_4_TCOFFEE");

    }
  else D->replicates=0;

  

  if (getenv ("enb_4_TCOFFEE")){D->enb  =atoi("enb_4_TCOFFEE");}
  else D->enb=3;

    
  if (getenv ("columns4tree_4_TCOFFEE"))D->col=file2column_list (getenv ("columns4tree_4_TCOFFEE"), NULL);
  else D->col=msa2column_list (A, NULL);
  

  

  D->pos=aln2pos_simple (A, A->nseq);

  
  D->dm3d=aln2dm3d (A, CL, &D->extremed);
  if (getenv ("maxd_4_TCOFFEE"))
    {
      if (strm (getenv ("maxd_4_TCOFFEE"), "scan"))D->maxd=-1;
      else D->maxd  =(double)atof(getenv("maxd_4_TCOFFEE"))*100;
    }
  else D->maxd=D->extremed*100+1;
  


  if (D->maxd>=0 && D->maxd<MY_EPSILON)D->maxd=D->extremed*100;
  
  D->dm =declare_double(A->nseq, A->nseq);
  D->nsites=-1;
  D->used_site=(int*)vcalloc (A->len_aln, sizeof (int));
  
  D->nsitepairs=-1;
  D->used_site_pair=declare_int (A->len_aln, A->len_aln);
  
  if (A->Tree) free_aln (A->Tree); A->Tree=NULL;
  A->Tree=declare_aln2(D->replicates+1, 0);
  (A->Tree)->dmF_list=(char**)vcalloc (D->replicates+1, sizeof (char*));
  for (a=0; a<=D->replicates; a++)(A->Tree)->dmF_list[a]=vtmpnam (NULL);
  (A->Tree)->nseq=0;
  return D;
}


p3D * makerep (p3D *D, int mode)
{
  if      ( mode==0)D->colrep=col2rep    (D->col, D->colrep, D->N);
  else if ( mode==1)D->colrep=col2bsrep1 (D->col, D->colrep, D->N);
  else if ( mode==2)D->colrep=col2bsrep2 (D->col, D->colrep, D->N);
  
  return D;
}


int **col2rep (int **colin,int **colout, int ni)
{
  int i;
  if (!colout)colout=declare_int (ni+1,3);
  for (i=0; i<ni; i++)
    {
      colout[i][0]=colin[i][0];
      colout[i][1]=colin[i][1];
    }
  colout[i][0]=-1;
  return colout;
}
int **col2bsrep1 (int **colin,int **colout, int ni)
{
  int i;
  if (!colout)colout=declare_int (ni+1,3);
  for (i=0; i<ni; i++)
    {
      int ri=rand()%ni;
      colout[i][0]=colin[ri][0];
      colout[i][1]=colin[ri][1];
    }
  colout[i][0]=-1;
  return colout;
}
int **col2bsrep2 (int **colin,int **colout, int ni)
{
  /*select maxp sites with re-sampling among the maxp sites in colin*/
  /*keep the all against all that are declared in colin*/
  /*The total number is guarranteed, but not the total number of pairs*/
  
    
  int i, j, k;
  int p1,p2, pmax;
  static int  ns;
  static int *l1;
  static int *ls;
  static int **pairs;
  static int *bs;
  static int rni;
  if (!colout)colout=declare_int (ni+1,3);

  if (rni!=ni)//flush out everything
    {
      if (l1)vfree(l1);
      if (ls)vfree (ls);
      ns=0;
      if (pairs)free_int(pairs, -1);
      l1=NULL;
      ls=NULL;
      pairs=NULL;
      rni=ni;
    }
  if (rni==0)return NULL;
  
  if (!l1)
    {
      //estimate the max #columns
      
      int mni=0;//index of the right-most column
      
      for (mni=0,i=0; i<ni; i++)
	{
	  if (colin[i][0]>mni)mni=colin[i][0];
	  if (colin[i][1]>mni)mni=colin[i][1];
	}
      l1=(int*)vcalloc (mni+1, sizeof (int));//cache for the columns
      ls=(int*)vcalloc (mni+1, sizeof (int));//ls: list of sites => list of used columns, collected from tghe cache
      for (i=0; i<ni; i++)
	{
	  l1[colin[i][0]]=1;
	  l1[colin[i][1]]=1;
	}
      for (ns=0,i=0; i<=mni; i++)if (l1[i])ls[ns++]=i;
      
      pairs=declare_int (mni+1,mni+1);//Look up 2D array of all used columns - sparse because of the distance filtering
      for (i=0; i<ni; i++)
	{
	  int p1=colin[i][0];
	  int p2=colin[i][1];
	  pairs[p1][p2]=pairs[p2][p1]=1;
	}
      bs=(int*)vcalloc (ns, sizeof (int));
    }
  
  for (i=0; i< ns; i++)bs[i]=ls[rand()%ns];
  
  for (k=0,i=0; i<ns-1; i++)
    {
      for (j=i+1; j<ns; j++)
	{
	  p1=bs[i];
	  p2=bs[j];
	  if (pairs[p1][p2])
	    {
	      colout[k][0]=p1;
	      colout[k][1]=p2;
	      k++;
	    }
	}
    }
  
  colout[k][0]=-1;
  return colout;
}



int **col2scramble_col (int **col, int ni)
{
  int i;
  for (i=0; i<ni; i++)col[i][2]=rand()%(ni*2);
  sort_int (col, 3, 2, 0, ni-1);
  return col;
}
int col2n (int **col)
{
  int i=0;
  if (!col) return 0;
  else if (!col[0])return 0;
  
  while (col[i++][0]!=-1);
  return i-1;
}
int   col2max(int **col)
{
  int i;
  int max=0;
  while (col[i][0]!=-1)
    {
      int c1=col[i][0];
      int c2=col[i][1];
      if ( c1>max)max=c1;
      if ( c2>max)max=c2;
      i++;
    }
  return max;
}
Alignment * addtree (p3D *D,Alignment *A)
{
  static char *treeF=vtmpnam (NULL);
  static int ml=get_longest_string (A->name, A->nseq, NULL, NULL)+1;
  FILE *fp;
  int s1, s2;
  Alignment *TREEA=A->Tree;
  int dotree=1;
  if (!getenv ("THREED_TREE_DM"))dotree=0;
   

  if (A->nseq>2 && dotree)
    {
      dist2nj_tree (D->dm, A->name, A->nseq, treeF);
      TREEA->seq_al[TREEA->nseq]=file2string(treeF);
      sprintf (TREEA->name[TREEA->nseq], "%d",TREEA->nseq+1); 	    
     
    }
  else
    {
      sprintf (TREEA->name[TREEA->nseq], "%d",TREEA->nseq+1); 
      TREEA->seq_al[TREEA->nseq]=csprintf (TREEA->seq_al[TREEA->nseq], "uncomputedtree");
    }
  fp=vfopen (TREEA->dmF_list[TREEA->nseq], "w");
  
  fprintf ( fp, "%d \n", A->nseq);
  if (getenv ("PRINT_NSITES"))
    {
      float p1=(float)(D->nsites*100)/(float)A->len_aln;
      float p2=(float)(D->nsitepairs*100)/(float)(A->len_aln*A->len_aln-A->len_aln);
     
      
      fprintf (fp, "!# MAXD: %.2f  Angstrom",(float)D->maxd/100);
      if (D->nsites >-1)   fprintf (fp, " -- NSITES: %d  %.2f%%",D->nsites, p1);
      if (D->nsitepairs>-1)fprintf (fp, " -- NSITEPAIRS: %d %.2f%%", D->nsitepairs,p2);
      fprintf (fp, "\n");
    }
  for ( s1=0; s1<A->nseq;s1++)
    {
      fprintf (fp, "%-*.*s ", ml,ml,A->name[s1]);
      for (s2=0; s2<A->nseq; s2++)
	fprintf (fp, "%6.3f ", (float)((double)D->dm[s1][s2])/(float)100);
      fprintf (fp, "\n");
    }
  fprintf (fp, "\n");
  
  vfclose (fp);
  TREEA->nseq++;
  return TREEA;
}
int filter_columns_with_dist(Alignment *B, int **pos,int **col, int***dm, double maxd)
{
  if ( getenv ("strict_maxd_4_TCOFFEE") && atoi(getenv ("strict_maxd_4_TCOFFEE"))==1)return filter_columns_with_dist_strict(B,pos,col,dm,maxd);
  else return filter_columns_with_dist_relaxed(B,pos,col,dm,maxd);
} 
  
int filter_columns_with_dist_strict(Alignment *B, int **pos,int **col, int***dm, double maxd)
{
  //Keep every pair of columns with ALL distances fulfilling the condition*/
  int i,ni,s;
  i=ni=0;
  if (maxd<MY_EPSILON)return col2n(col);
 
  while (col[i][0]!=-1)
    {
      int count;
      int max;
      int c1=col[i][0];
      int c2=col[i][1];
      
      for (max=0,count=0,s=0; s<B->nseq; s++)
	{
	  int r1=pos[s][c1]-1;
	  int r2=pos[s][c2]-1;
	  if (r1<r2){int rb=r1; r1=r2;r2=rb;}
	  
	  if (r1>=0 && r2>=0)
	    {
	      max++;
	      if (dm[s][r1][r2]<=maxd)count++;
	    }
	  
	}
      if ( count==max)
	{
	  col[ni][0]=c1;
	  col[ni][1]=c2;
	  ni++;
	}
      i++;
    }
  col[ni][0]=-1;
  
  return ni;
}
int filter_columns_with_dist_relaxed(Alignment *B, int **pos,int **col, int***dm, double maxd)
{
  //Keep every pair of columns with at least one distance below the threshold
  int i,ni,s;
  i=ni=0;
  if (maxd<MY_EPSILON)return col2n(col);
  
  while (col[i][0]!=-1)
    {
      int c1=col[i][0];
      int c2=col[i][1];
      
      for (s=0; s<B->nseq; s++)
	{
	  int r1=pos[s][c1]-1;
	  int r2=pos[s][c2]-1;
	  if (r1<r2){int rb=r1; r1=r2;r2=rb;}
	  
	  if (r1>=0 && r2>=0 && dm[s][r1][r2]<=maxd)
	    {
	      col[ni][0]=c1;
	      col[ni][1]=c2;
	      ni++;
	      break;
	    }
	}
      i++;
    }
  col[ni][0]=-1;
  return ni;
}
int filter_columns_with_gap (int **col, Alignment *B, float max_gap)
{
  float *gap=(float*)vcalloc (B->len_aln, sizeof (float));
  int c,s, i, ni;

  for (c=0; c<B->len_aln; c++)
    {
      for (s=0; s<B->nseq; s++)if (B->seq_al[s][c]=='-')gap[c]+=1;
      gap[c]/=(float)B->nseq;
    }
  i=ni=0;
  while (col[i][0]!=-1)
    {
      int c1=col[i][0];
      int c2=col[i][1];
      if (gap[c1]<=max_gap && gap[c2]<=max_gap)
	{
	  col[ni][0]=c1;
	  col[ni][1]=c2;
	  ni++;
	}
      i++;
    }
  vfree (gap);
  col[ni][0]=-1;
  return ni;
}
	  
	  

 int** msa2column_list (Alignment *B, int **col)
{
  int i, j, n;
  if (!col)col=declare_int(((B->len_aln*B->len_aln)-B->len_aln)/2+1,3);
  for (n=0,i=1; i<B->len_aln-1; i++)
    for (j=i+1; j<B->len_aln; j++, n++)
      {
	col[n][0]=i;
	col[n][1]=j;
      }
  col[n][0]=-1;
  return col;
}
 int   **  file2column_list (char *file, int **list)
{
  int i,j,n;
  char ***l;
  
  if (!(l=file2list (file, " \t")))return NULL;
  
  
  i=j=0;
  while (l[i++]);
  if (!list)list=declare_int (i+1, 3);
  i=j=0;
  while (l[i])
    {
      
      int n=atoi(l[i][0]);
  
      if (n>=3 && l[i][1][0]!='#')
	{
	  
	  list[j][0]=atoi (l[i][1])-1;
	  list[j][1]=atoi (l[i][2])-1;
	  j++;
	}
      i++;
    }  
  list[j][0]=-1;
  free_arrayN((void***)l, 3);
  
  return list;
}
   
int*** aln2dm3d (Alignment *A, Constraint_list*CL, double *max)
{
  int s, r, sA, i;
  Sequence *S=CL->S;
  int ***dm=(int***)vcalloc (A->nseq, sizeof (int**));
  max[0]=0;
  for (sA=0; sA<A->nseq; sA++)
    {
      s=name_is_in_hlist (A->name[sA], S->name, S->nseq);
      dm[sA]=(int**)vcalloc(S->len[s]+1,sizeof (int*));
      for (r=1; r<=S->len[s]; r++)dm[sA][r]=(int*)vcalloc (r, sizeof(int));
      for (r=1; r<=S->len[s]; r++)
	{
	  for (i=1; i<CL->residue_index[s][r][0]; i+=ICHUNK)
	    {
	      int nr =CL->residue_index[s][r][i+R2];
	      int we =CL->residue_index[s][r][i+WE];
	      int ns =CL->residue_index[s][r][i+SEQ2];
	      
	      if (s!=ns){HERE ("Warning: Contact library contains inter-sequence data contacts: %d %d",i,ns);}
	      if (r>nr)dm[sA][r-1][nr-1]=we;
	      else dm[sA][nr-1][r-1]=we;
	      if ((double)we>max[0])max[0]=(double)we;
	    }
	}
    }
  max[0]/=100;//put max back into Angstrom
  return dm;
}


Alignment *aln2trim3d (Alignment *A, Constraint_list *CL)
{
  //Keep only sequences with contact information
  Alignment *B=NULL;
  static char *tmpF=vtmpnam (NULL);
  int ns, s,r,rs;
  FILE*fp;
  if (!CL) return B;
  
  fp=vfopen (tmpF, "w");
  for (ns=0,s=0; s<A->nseq; s++)
    {
      if ((rs=name_is_in_hlist (A->name[s], (CL->S)->name, (CL->S)->nseq)!=-1))
	{
	  for (r=1; r<=(CL->S)->len[s];r++)
	    {
	      if (CL->residue_index[s][r][0])
		{
		  fprintf (fp, ">%s\n%s\n",A->name[s], A->seq_al[s]);
		  ns++;
		  break;
		}
	    }
	}
    }
  vfclose (fp);
  if (ns)B=quick_read_fasta_aln (B,tmpF);
  return B;
}

int aln2dm (p3D *D, Alignment *A)
{
  // Turns A into either a distance matrix or a similarity matrix
  int rv=1;
  int s1, s2;
  int c1,c2, t;
    
  if (D->used_site)
    {
      for (c1=0; c1<A->len_aln; c1++)
	{
	  D->used_site[c1]=0;
	  for ( c2=0; c2<A->len_aln; c2++)
	    D->used_site_pair[c1][c2]=0;
	}
    }
  
  
  for (s1=0; s1<A->nseq-1; s1++)
    for (s2=s1+1; s2<A->nseq; s2++)
      {
	D->dm[s1][s2]=D->dm[s2][s1]=pair2dist (D,s1, s2);
	if (D->dm[s1][s2]<-99)rv=0;
      }
  if ( D->used_site)
    {
      for (D->nsites=0,D->nsitepairs=0,c1=0; c1<A->len_aln; c1++)
	{
	  D->nsites+=D->used_site[c1];
	  for (c2=0; c2<A->len_aln; c2++)
	    {
	      D->nsitepairs+=D->used_site_pair[c1][c2];
	    }
	}
    }
  
  return rv;
}
double pair2dist(p3D *D, int s1, int s2)
{ 
  //s1 and s2 match the D->dm3d 
  //check the offseet of the residue index - starts at 1 or 0
  double score=0;
  double max=0;
  double rscore, rmax;
  double w1, w2;
  int i;
  int ns=0;
  if (s1>s2){int sb=s1;s1=s2;s2=sb;}

  i=0;
  while (D->colrep[i][0]!=-1)
    {
      int c1=D->colrep[i][0];
      int c2=D->colrep[i][1];
      i++;
      if (c1<c2){int cb=c1;c1=c2;c2=cb;}
      

      int r11=D->pos[s1][c1]-1;
      int r12=D->pos[s1][c2]-1;
      int r21=D->pos[s2][c1]-1;
      int r22=D->pos[s2][c2]-1;
      
      if (r11<0 || r12<0)continue;
      if (r21<0 || r22<0)continue;
      
      if (abs((r12-r11))<D->enb)continue;
      if (abs((r22-r21))<D->enb)continue;
      
      
      
      w1=(double)D->dm3d[s1][r11][r12];
      w2=(double)D->dm3d[s2][r21][r22];
      
      if (w1<MY_EPSILON || w2<MY_EPSILON)continue;
      if (w1>D->maxd || w2>D->maxd)continue;

      ns++;
      phylo3d2score (w1, w2, &rscore, &rmax);
      
      score+=rscore;
      max+=rmax;
      if (D->used_site)D->used_site[c1]=D->used_site[c2]=1;
      if (D->used_site_pair)D->used_site_pair[c1][c2]=D->used_site_pair[c2][c1]=1;

    }
  if (max<MY_EPSILON)score=-100;
  else if (getenv ("THREED_TREE_DM"))
    {
      if(atoigetenv ("THREED_TREE_DM")==1)
	score=score/max;
      else
	score=sqrt(score/max);
      
    }
  else score=(double)100*((double)1 - (score/max));

  
  return score;
}
double phylo3d2score (double w1, double w2, double *rscore, double *rmax)
{
  double we=0;
  double sc=0;
  static int setparam;
  static int   distance_mode;
  static double distance_modeE;
  static int no_weights;
  //4 - initial measure (1-rlativeDelta)²=> concave function, different from Kimura
  //7 - more Kimura like: (1-rlativeDelta²)=>convex
  //8 - linear rdelta (1-RelativeDelta)
  if (!setparam)
    {
      setparam=1;
      if (getenv ("THREED_TREE_MODE"))distance_mode=atoigetenv ("THREED_TREE_MODE");
      else distance_mode=4;
      
      if (getenv ("THREED_TREE_MODE_EXP"))distance_modeE=atofgetenv ("THREED_TREE_MODE_EXP");
      else distance_modeE=2;

      if (getenv ("THREED_TREE_NO_WEIGHTS"))no_weights=atoigetenv ("THREED_TREE_NO_WEIGHTS");
      else no_weights=1;
    }

 
			   
   //first attempt-- Major issue because non symetrical and therefore not a distance
   if (!distance_mode)
     {
       static int warn;
       we=w1;
       sc=((MIN((w1/w2),(w2/w1))));
       if (warn==0)
	 {
	   add_warning ( stderr, "distance_mode==0 should not be used");
	   warn=1;
	 }
       
     }
   //same as before but symetrical: the distance ratio weighted by the average distance
   else if (distance_mode==1)
     {
       we=(w1+w2)/2;
       sc=((MIN((w1/w2),(w2/w1))));
     }
   //the absolute difference of distance normalized by the average distance
   else if (distance_mode==2)
     {
       
       we=(w1+w2)/2;
       sc=(double) 1-(FABS((w1-w2))/we);
       
     }
   //the absolute difference of distance normalized by the average distance and weighted by the average distance
   else if (distance_mode==3)
     {
       we=(w1+w2)/2;
       sc=(double)1-(FABS((w1-w2))/we);
     }
   else if (distance_mode ==4)
     {
       we=((w1>w2)?w1:w2);
       sc=(double) 1-(FABS((w1-w2))/we);
       
     }
   else if (distance_mode ==5)
     {
       we=(w1+w2);
       sc=(double)1-(FABS((w1-w2))/we);
       
     }
   else if (distance_mode ==6)
     {
       we=1;
       sc=(double)(FABS((w1-w2)));
     }
   else if ( distance_mode ==7)//Convexb -- Kimura like: 1-D^2 -- assiumes that D is over-estimated
     {
       double rdelta;
       we=((w1>w2)?w1:w2);
       rdelta=FABS((w1-w2))/we;
              
       rdelta*=rdelta;
       rscore[0]=1-rdelta;
       rmax[0]=1;
       return rscore[0];
     }
   else if ( distance_mode ==8)//Linear -- 1-D -- No assumption on D
     {
       double rdelta;
       we=((w1>w2)?w1:w2);
       rdelta=FABS((w1-w2))/we;
       rscore[0]=1-rdelta;
       rmax[0]=1;
       return rscore[0];
     }	 
   else if  (distance_mode ==9)//concave -- (1-D)^2 -- Asumes D is iunder- estimated
     {
       double rdelta;
       we=((w1>w2)?w1:w2);
       rdelta=FABS((w1-w2))/we;
       rscore[0]=1-rdelta;
       rscore[0]*=rscore[0];
       rmax[0]=1;
       return rscore[0];
     }
   else if (distance_mode==10)
     {
       static int setv=0;
       if (!setv){cputenv ("THREED_TREE_DM=1");setv=1;}
       rscore[0]=FABS((w1-w2));
       rmax[0]=1;
       
       return rscore[0];
     }
   else if (distance_mode==11)
     {
       static int setv=0;
       if (!setv){cputenv ("THREED_TREE_DM=2");setv=1;}
       rscore[0]=(w1-w2)*(w1-w2);
       rmax[0]=1;
       return rscore[0];
     }
	 
   
   if (no_weights)
     we=1;
   
   //Compute the score and its normalization
   //If 
   rmax[0]=we;
   rscore[0]=pow(sc,distance_modeE)*we;
   return rscore[0]/rmax[0];   
}


//////////////////////////////////////////////////////////////////////
///                                                              /////
///                                                              /////
///       Combined Boostrap Analysis                            /////
///                                                              /////
///                                                              /////
//////////////////////////////////////////////////////////////////////

Sequence * get_phylo3d_seq  (char*family);
//NT_node get_phylo3d_bm_tree (char*family,char *type, int nbs, Sequence *S);
NT_node get_phylo3d_bm_tree (char*family,int type, int nbs, int bstype, Sequence *S);
int tree2splits4phylo3d_bm (NT_node T, int ns,float **sl,char**splits, int* n);


int phylo3d_bm ( char *name)
{
  NT_node ***T;
  
  Sequence *S;
  float *****sl;
  char  *****splits;
  int   ***nn;
  char **testlist=(char**)vcalloc(100, sizeof (char*));
  int n_ppformula, n_bsformula, n_cmode,npp,pp, ref, bs, split, cmode;
  char **ppformula, **bsformula, **cmodelist;
  char *refS;
  float *auc;
  
  S=get_phylo3d_seq(name);

  
  T        =(NT_node***  )declare_arrayN(3,sizeof (NT_node),3,3,3);
  nn       =(int    ***  )declare_arrayN(3,sizeof (int)    ,3,3,3);
  sl       =(float  *****)declare_arrayN(5,sizeof (float)  ,3,3,3,S->nseq*3,20);
  splits   =(char   *****)declare_arrayN(5,sizeof (char)   ,3,3,3,S->nseq*3,S->nseq+1);
  ppformula=(char**)vcalloc (100, sizeof (char*));
  bsformula=(char**)vcalloc (100, sizeof (char*));
  cmodelist=(char**)vcalloc (100, sizeof (char*));
  auc      =(float*)vcalloc (100, sizeof (float));


  //Collect all the splits of all the provided trees
  //assuming: <family><reference trees><boostrap method>
  //assuming: <family>_<IMD|ME|ML>splits_<IMD|ME|ML>bs<25|100|[empty=200]>.trees
  //assuming: <family>_<IMD|ME|ML>.trees for the original trees based on 200 columns, and boostrapped on their MSA
  //In each file, the top tree is the reference <IMD|ME|ML> tree based on 200 columns, followed with a 100 replicates dine according to the boostrap method

  //T[Ref=I,M,L][ncol=25,100,200][bs_support_trees=I,E,L]
  
  if (1==1)// Read all the trees
    {
      int a, b, c;    
      for (a=0; a<3; a++)
	for ( b=0; b<3; b++)
	  for (c=0; c<3; c++)
	    T[a][b][c]=get_phylo3d_bm_tree(name,a,b,c, S);
    }
  // tree T[ref][ncol4bs][method4bs] contains nn[a][b][c] splits
  // each split is a 010001 string in splits[a][b][c]
  // the various characteristics of the splits are in sl[a][b][c][0-nn[a][b][c]][0-..]
  
  if (1==1)//collect all the splits
    {
      int a, b, c;
      for (a=0; a<3; a++)//ref
	for ( b=0; b<3; b++)//ncol for support
	  for (c=0; c<3; c++)//method for support
	    {
	      tree2splits4phylo3d_bm(T[a][b][c],S->nseq,sl[a][b][c],splits[a][b][c], &nn[a][b][c]);
	    }
    }
  
  /*Set up formula to determine reference branches*/
  // Methods used to decide on the pp + BS Threshold - as measured on 200 columns
  // bs defines the threshold for a brnach to be set as PP, ppformula defines the number of methods tha must feature the same branch with a support > BS
  n_ppformula=0;
  for (bs=0; bs<3; bs++)
    {
      char *bss;
      if      (bs==0)bss=csprintf (NULL, "000");
      else if (bs==1)bss=csprintf (NULL, "080");
      else if (bs==2)bss=csprintf (NULL, "100");
      
      ppformula[n_ppformula++]=csprintf (NULL, "I__%s",bss);
      ppformula[n_ppformula++]=csprintf (NULL, "_E_%s",bss);
      ppformula[n_ppformula++]=csprintf (NULL, "__L%s",bss);
      ppformula[n_ppformula++]=csprintf (NULL, "IE_%s",bss);
      ppformula[n_ppformula++]=csprintf (NULL, "I_L%s",bss);
      ppformula[n_ppformula++]=csprintf (NULL, "IEL%s",bss);
      ppformula[n_ppformula++]=csprintf (NULL, "_EL%s",bss);
      vfree(bss);
    }
  
  //Set up methods to estimate BS: combined methods + number BS
  //ppformula defines the replicates used estimate the combined support and the number of columns on which these replicates are based
  n_bsformula=0;
  for(bs=0; bs<3; bs++)
    {
      char *bss;
      if      (bs==0)bss=csprintf (NULL, "025");
      else if (bs==1)bss=csprintf (NULL, "100");
      else if (bs==2)bss=csprintf (NULL, "200");
      
      bsformula[n_bsformula++]=csprintf(NULL,"I__%s",bss);
      bsformula[n_bsformula++]=csprintf(NULL,"_E_%s",bss);
      bsformula[n_bsformula++]=csprintf(NULL,"__L%s",bss);
      bsformula[n_bsformula++]=csprintf(NULL,"IE_%s",bss);
      bsformula[n_bsformula++]=csprintf(NULL,"I_L%s",bss);
      bsformula[n_bsformula++]=csprintf(NULL,"IEL%s",bss);
      bsformula[n_bsformula++]=csprintf(NULL,"_EL%s",bss);
      vfree(bss);
    }

  n_cmode=0;
  cmodelist[n_cmode++]=csprintf (NULL,"avg");
  cmodelist[n_cmode++]=csprintf (NULL,"geo");
  cmodelist[n_cmode++]=csprintf (NULL,"max");
  cmodelist[n_cmode++]=csprintf (NULL,"min");
  
  

  refS= (char*)vcalloc (3, sizeof (char*));
  for (cmode=0; cmode<n_cmode; cmode++)
    {
      
      for (pp=0; pp<n_ppformula; pp++)
	{
	  // pp have been set according the ppformula rules
	  // Now each pp branched is marked as sl[reftree][nbs][bstype][split][3]=1
	  
	  set_pp(name,S->nseq,sl, splits,nn,ppformula[pp]);
	  
	  
	  for (ref=0; ref<3; ref++)// Reference Tree - Estimated on 200 sites with I/L/E
	    {
	      refS=csprintf (refS, "___");
	      if      (ref==0) refS[0]='I'; //IMD
	      else if (ref==1) refS[1]='E'; //ME
	      else if (ref==2) refS[2]='L'; //ML
	      
	      //Depending on how pp were collected, the ref tree may have less npp (e.g. if reftree is I and PP=L inter E)
	      //Topology (sl[a][0][a]==topology(sl[a][1][a]==topology(sl[a][2][a]
	      //work on sl[a][2][a] by convention as it is the true reference tree with its native bs
	      
	      npp =splits2npp (nn[ref][2][ref], sl[ref][2][ref]);
	      if ( npp>0 && npp<S->nseq)
		{
		  for (bs=0; bs<n_bsformula; bs++)
		    {
		      int bsbin=bsformula[bs][3]-'0';
		      
		      for (split=0; split<nn[ref][2][ref]; split++)
			{
			  //float newbs=0;
			  //int nbs=0;
			  
			  //if (bsformula[bs][0]=='I'){newbs+=sl[ref][bsbin][0][split][2];nbs++;}
			  //if (bsformula[bs][1]=='E'){newbs+=sl[ref][bsbin][1][split][2];nbs++;}
			  //if (bsformula[bs][2]=='L'){newbs+=sl[ref][bsbin][2][split][2];nbs++;}
			  //sl[ref][bsbin][ref][split][4]=newbs/(float)nbs;
			  int ncbs=0;
			  float cbs[3];
			  
			  if (bsformula[bs][0]=='I')cbs[ncbs++]=sl[ref][bsbin][0][split][2];
			  if (bsformula[bs][1]=='E')cbs[ncbs++]=sl[ref][bsbin][1][split][2];
			  if (bsformula[bs][2]=='L')cbs[ncbs++]=sl[ref][bsbin][2][split][2];
			  
			  sl[ref][bsbin][ref][split][4]=bs2combo(cbs, ncbs, cmodelist[cmode]);
			}
		      auc[bs]=splits2auc(nn[ref][bsbin][ref], sl[ref][bsbin][ref]);
		    }
		  //Display AUC - the float values of the AUCs
		  fprintf ( stdout,"AUC family: %-20s cmode: %s ref: %s ppmode: %s npp: %d bs: ", name, cmodelist[cmode],refS, ppformula[pp],npp); 
		  for (bs=0; bs<n_bsformula; bs++)
		    {
		      fprintf (stdout,"%s %.2f %d %d ",bsformula[bs],auc[bs], (auc[bs]>0.99999)?1:0, 4*(bs+1)-2+10);
		    }
		  fprintf (stdout, "\n");
		}
	    }
	}
    }
  exit (0);

}

float bs2combo(float *cbs, int n, char *mode)
{
  int a;
  float newbs;
  
  newbs=0;
  if (n==0)return 0;
  else if (strm (mode, "average"))
    {
      for (a=0; a<n; a++)newbs+=cbs[a];
      newbs/=(float)n;
    }
  else if ( strm (mode, "min"))
    {
      newbs=cbs[0];
      for ( a=0; a<n; a++)
	if (cbs[a]<newbs)newbs=cbs[a];
    }
  else if ( strm (mode, "max"))
    {
      newbs=cbs[0];
      for (a=0; a<n; a++)
	if (cbs[a]>newbs)newbs=cbs[a];
    }
  else if ( strm (mode, "geometric"))
    {
      for (a=0; a<n; a++)
	{
	  if (cbs[a]<0.0001) return 0;
	  else newbs+= log((double)cbs[a]);
	}
      newbs=(float)exp(double(newbs/(float)n));
    }
  return newbs;
  
}  
	    
	

Sequence * get_phylo3d_seq  (char*family)
{
  Sequence *S;
  char *treelist=csprintf (NULL, "%s_%s.trees",family, "IMD");
  S=get_treelist (treelist);
  return tree2seq(newick_string2tree(S->seq[0]), NULL);
}




NT_node get_phylo3d_bm_tree (char*family,int type, int ncol, int bstype, Sequence *S)
{
  static char *treelist;
  static char *typeS;
  static char *ncolS;
  static char *bstypeS;
  NT_node T;

  if (type == bstype && ncol==2)
    {
      if      (type==0)treelist=csprintf (treelist, "%s_IMD.trees",family);
      else if (type==1)treelist=csprintf (treelist, "%s_ME.trees",family);
      else if (type==2)treelist=csprintf (treelist, "%s_ML.trees",family);
    }
  else
    {  
      if      (type==0)typeS=csprintf (typeS, "IMDsplits");
      else if (type==1)typeS=csprintf (typeS, "MEsplits");
      else if (type==2)typeS=csprintf (typeS, "MLsplits");
      
      
      if      (ncol==0)ncolS=csprintf (ncolS, "25");
      else if (ncol==1)ncolS=csprintf (ncolS, "100");
      else if (ncol==2)ncolS=csprintf (ncolS, "");

      

      
      if      (bstype==0)bstypeS=csprintf (bstypeS, "IMDbs");
      else if (bstype==1)bstypeS=csprintf (bstypeS, "MEbs");
      else if (bstype==2)bstypeS=csprintf (bstypeS, "MLbs");

      treelist=csprintf (treelist,"%s_%s_%s%s.trees", family,typeS,bstypeS,ncolS);
    }

  if ( check_file_exists (treelist)==0)
    {
      printf_exit ( EXIT_FAILURE,stderr, "\nERROR: %s does not exist [FATAL]", treelist);
    }

  T=treelistF2node_support (treelist);
  T=prune_tree(T,S);
  T=recode_tree(T,S);
  return T;

}
 


int tree2splits4phylo3d_bm (NT_node T, int ns,float **sl,char**splits, int* n)
{
  if (!T)  return 0;
  if (!sl) return 0;

  tree2splits4phylo3d_bm(T->right, ns, sl, splits,n);
  tree2splits4phylo3d_bm(T->left , ns, sl, splits,n);

  if (!T->right) return 1;
  else if (T->parent && !(T->parent)->parent)return 1;
  else if ( T->dist==0)return 1;
  else
    {
      int flip=T->lseq2[0];
      int a;
      int depth=0;
      for (a=0; a< ns; a++)
	{
	  splits[n[0]][a]=(flip==1)?T->lseq2[a]+'0':1-T->lseq2[a]+'0';
	  depth+=T->lseq2[a];
	}
      
      splits[n[0]][a]='\0';
      sl[n[0]][0]=T->dist;
      sl[n[0]][1]=MIN((ns-depth), depth);
      sl[n[0]][2]=T->bootstrap;
      n[0]++;

    }
  return 1;
}






void clean_pp (float *****sl, int ***nn)
{
  int a, b, c,d;
  for ( a=0; a<3; a++)
    for (b=0; b<3; b++)
      for (c=0; c<3; c++)
	for ( d=0; d<nn[a][b][c]; d++)
	  sl[a][b][c][d][3]=0;
  return;
}
int set_pp (char *family, int nseq,float ***** sl, char ***** splits, int ***nn, char *teststring) 
{
  //imd is 0 me is 1, ml is 2
  // 25 columns is 0, 100 columns is 2, 200 columns is 2
  //set the pp according to the agreement of imd, me, ml on the 200 trees (ns[imd/me/ml][2][imd/me/ml]
  // sl..[0] branch lentgh
  // sl..[1] depth
  // sl..[2] native bs
  // sl..[3] PP mark (1 if PP)
  // sl..[4] will contained the combined bs
  
  int imd, me, ml, bs;
  int a, b, c, d, e,f,g,t;
  int passed;
  char *rb;
  int   depth;
  int npp=0;

  imd=me=ml=0;
  if (teststring[0]=='I')imd=1;
  if (teststring[1]=='E')me=1;
  if (teststring[2]=='L')ml=1;
  bs=atoi(teststring+3);

  
  clean_pp(sl, nn);
  for (a=0; a<nn[0][2][0]; a++)
    for ( b=0; b<nn[1][2][1]; b++)
      for ( c=0; c<nn[2][2][2]; c++)
	{
	  char  *b1=splits[0][2][0][a];
	  float bs1   = sl[0][2][0][a][2];
	  
	  char  *b2=splits[1][2][1][b];
	  float bs2=    sl[1][2][1][b][2];
	  
	  char  *b3=splits[2][2][2][c];
	  float bs3=    sl[2][2][2][c][2];
	  
	  passed=0;
	  if ( imd && me && ml)
	    {
	      if (bs1>=bs && bs2>=bs && bs3>=bs && strm (b1,b2) && strm (b1, b3))
		{
		  passed=1;
		  rb=b1;
		}
	    }
	  else if ( imd && me)
	    {
	       if (bs1>=bs && bs2>=bs  && strm (b1,b2))
		 {
		   passed=1;
		   rb=b1;
		 }
	    }
	  else if ( imd && ml)
	    {
	      if (bs1>=bs && bs3>=bs  && strm (b1,b3))
		 {
		   passed=1;
		   rb=b1;
		 }
	    }
	  else if ( me && ml)
	    {
	      if (bs2>=bs && bs3>=bs  && strm (b2,b3))
		 {
		   passed=1;
		   rb=b2;
		 }
	    }
	  else if (imd)
	    {
	      if (bs1>=bs)
		 {
		   passed=1;
		   rb=b1;
		 }
	    }
	  else if (me)
	    {
	      if (bs2>=bs)
		 {
		   passed=1;
		   rb=b2;
		 }
	    }
	  else if (ml)
	    {
	      if (bs3>=bs)
		 {
		   passed=1;
		   rb=b3;
		 }
	    }
	  
	  if (passed) 
	    {
	      for (d=0; d<3; d++)
		for (e=0; e<3; e++)
		  for (f=0; f<3; f++)
		    for (g=0; g<nn[d][e][f]; g++)
		      {
			if (strm (rb, splits[d][e][f][g]))sl[d][e][f][g][3]=1;
		      }
	    }
	}
  // All the original trees may OR may not contain the split that has just been identified
  // This reporting makes it possible to collect the PP associated with a ppmode that collects specific splits across the references
  // for instance I+E may reprt 200 branches, but if L is the refernce tree, it may only contains 150 of these branches thus leading to 150 trees being reported as PP for ref: __L ppmode: IE_080
  for (a=0; a<3; a++)
    {
      for (b=0; b<nn[a][0][0]; b++)
	{
	  if (sl[a][0][0][b][3]==1)
	    {
	      char ref[10];
	      depth=sl[a][0][0][b][1];
	      if (a==0)      sprintf (ref, "I__");
	      else if ( a==1)sprintf (ref, "_E_");
	      else if ( a==2)sprintf (ref, "__L");
	      
	      fprintf ( stdout, "PPLIST family: %-10s nseq: %d ref: %s ppmode: %s split: %s depth: %3d rdepth: %.2f\n",family, nseq,ref,teststring,splits[a][0][0][b],depth, (float)((float)(2*depth)/(float)(nseq)));
	    }
	}
    }
  
  return 1;
}

int splits2npp  (int nnodes, float **nodes)
{
  int a;
  int npp=0;
  for (a=0; a<nnodes; a++)
    npp+=(nodes[a][3]>0.01)?1:0;
  return npp;
}
  




///////////////////////////////// Test //////////////////////////////////////////////
int phylo3d_bm_test ( char *name)
{
  NT_node ***T;
  
  Sequence *S;
  float *****sl;
  char  *****splits;
  int   ***nn;
  char **testlist=(char**)vcalloc(100, sizeof (char*));
  int n_ppformula, n_bsformula, npp,pp, ref, bs, split;
  char **ppformula, **bsformula;
  char *refS;
  float *auc;
  
  S=get_phylo3d_seq(name);

  
  T        =(NT_node***  )declare_arrayN(3,sizeof (NT_node),3,3,3);
  nn       =(int    ***  )declare_arrayN(3,sizeof (int)    ,3,3,3);
  sl       =(float  *****)declare_arrayN(5,sizeof (float)  ,3,3,3,S->nseq*3,20);
  splits   =(char   *****)declare_arrayN(5,sizeof (char)   ,3,3,3,S->nseq*3,S->nseq+1);
  ppformula=(char**)vcalloc (100, sizeof (char*));
  bsformula=(char**)vcalloc (100, sizeof (char*));
  auc      =(float*)vcalloc (100, sizeof (float));

  if (1==1)// Read all the trees
    {
      int a, b, c;    
      for (a=0; a<3; a++)
	for ( b=0; b<3; b++)
	  for (c=0; c<3; c++)
	    T[a][b][c]=get_phylo3d_bm_tree(name,a,b,c, S);
    }
  
  if (1==1)//collect all the splits
    {
      int a, b, c;
      for (a=0; a<3; a++)
	for ( b=0; b<3; b++)
	  for (c=0; c<3; c++)
	    {
	      tree2splits4phylo3d_bm(T[a][b][c],S->nseq,sl[a][b][c],splits[a][b][c], &nn[a][b][c]);
	    }
    }
  
  /*Set up formula to determine reference branches*/
  // Methods used to decide on the pp + BS Threshold - as measured on 200 columns
  n_ppformula=0;
  ppformula[n_ppformula++]=csprintf (NULL, "IEL080");
  set_pp (name, S->nseq,sl, splits,nn,ppformula[pp]);
  n_bsformula=0;
  bsformula[n_bsformula++]=csprintf(NULL,"I__025");
  bsformula[n_bsformula++]=csprintf(NULL,"_E_025");
  bsformula[n_bsformula++]=csprintf(NULL,"__L025");
  bsformula[n_bsformula++]=csprintf(NULL,"IE_025");
  bsformula[n_bsformula++]=csprintf(NULL,"I_L025");
  bsformula[n_bsformula++]=csprintf(NULL,"_EL025");
  bsformula[n_bsformula++]=csprintf(NULL,"IEL025");
  
  ref=0;
  npp =splits2npp (nn[0][0][0], sl[0][0][0]);
  
  
  if (npp>0 && npp<S->nseq)
    {
      fprintf ( stdout,"%-20s -- %2d pp %d seq -- %s%d%d%d ---", name,npp, S->nseq,"IEL",0, 0, 0);
      for (bs=0; bs<n_bsformula; bs++)
	{
	  int bsbin=bsformula[bs][3]-'0';
	  
	  float newbs=0;
	  int  nbs=0;
	  for (split=0; split<nn[ref][2][ref]; split++)
	    {
	      nbs=0;
	      newbs=0;
	      if (bsformula[bs][0]=='I'){newbs+=sl[ref][bsbin][0][split][2];nbs++;}
	      if (bsformula[bs][1]=='E'){newbs+=sl[ref][bsbin][1][split][2];nbs++;}
	      if (bsformula[bs][2]=='L'){newbs+=sl[ref][bsbin][2][split][2];nbs++;}
	      sl[ref][bsbin][ref][split][4]=newbs/nbs;
	      //if ( strm (bsformula[bs], "IML025"))HERE ("%2d %2d --- %3d %s\n",split, (int)sl[ref][bsbin][ref][split][2],(int)sl[ref][bsbin][ref][split][4],splits[ref][bsbin][ref][split]);
	    }
	  
	  fprintf (stdout, " %.2f",splits2auc(nn[ref][bsbin][ref], sl[ref][bsbin][ref]));
	}
      fprintf (stdout, "\n");
    }
  
  exit (0);

}

//////////////////////  AUC


float splits2auc(int nnodes, float **nodes)
{
  static LabeledData *data=(LabeledData *) vcalloc (nnodes*2, sizeof (LabeledData));
  int a;
  
  for (a=0; a<nnodes; a++)
    {
      //fprintf (stderr, "\t %d %f\n", (int)nodes[a][3],nodes[a][4]); 
      data[a].label=(nodes[a][3]>0.01)?1:0;
      data[a].score=(double)nodes[a][4];
    }
  int n = nnodes;
    
    // Sort the data by scores in descending order
    qsort(data, n, sizeof(LabeledData), compare4phylo3d);
    //for (a=0; a<n; a++)
    //  HERE ("Sample %d %f", data[a].label,data[a].score);
    
   
    // Calculate and print the AUC
    double auc = calculateAUC(data, n);
    return (float)auc;
}
int compare4phylo3d(const void *a, const void *b) {
    LabeledData *dataA = (LabeledData *)a;
    LabeledData *dataB = (LabeledData *)b;
    if (dataB->score > dataA->score) return 1;
    if (dataB->score < dataA->score) return -1;
    return 0;
}
double calculateAUC_1(LabeledData *data, int n);
double calculateAUC_2(LabeledData *data, int n);
double calculateAUC(LabeledData *data, int n) 
  {
    double auc1,auc2,auc;
    //auc1=calculateAUC_1(data,n);
    auc2=calculateAUC_2(data,n);
    //if ( auc!=calculateAUC_2(data,n))
    //  exit (0);
    //HERE ("%.2f %.2f", (float) auc1, (float) auc2);
    return auc2;
  }

double calculateAUC_1(LabeledData *data, int n)
{
  
    double auc = 0.0, prev_fpr = 0.0, prev_tpr = 0.0, tpr = 0.0, fpr = 0.0;
    int tp = 0, fp = 0, positives = 0, negatives = 0;
    
    for (int i = 0; i < n; i++) {
        if (data[i].label == 1) positives++;
        else negatives++;
    }
    
    for (int i = 0; i < n; i++) {
        if (data[i].label == 1) tp++;
        else fp++;

        tpr = (double)tp / positives;
        fpr = (double)fp / negatives;

        if (fp == 0 && tp > 0) { // Start of ROC curve (0,0) to (0,1)
            prev_tpr = tpr;
            continue;
        }

        if (i > 0 || (fp > 0 && tp > 0)) { // Avoid division by zero and ensure movement on ROC
            auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2.0;
            prev_fpr = fpr;
            prev_tpr = tpr;
        }
    }

    return auc;
}

double calculateAUC_3(LabeledData *data, int n)
{
  int i, item,positives, negatives;
  int tt, tf, ft, ff;
  int total_true_0, total_true_1;
  double sens, spec, tpf, fpf, tpf_prev, fpf_prev, roc_area;

  i=item=positives=negatives=tt=tf=ff=ft;
  tpf=fpf=tpf_prev=fpf_prev=sens=spec=roc_area=0;
  
  for (i = 0; i < n; i++) {
        if (data[i].label == 1) positives++;
        else negatives++;
  }

  for (i=0; i<n; i++)
    {
      HERE ( "------- Labels: %d", data[i].label);
    }
  
    tt=0;
    tf=positives;
    ft=0;
    ff=negatives;

    sens = ((double) tt) / ((double) (tt+tf));
    spec = ((double) ff) / ((double) (ft+ff));
    tpf = sens;
    fpf = 1.0 - spec;
    
    roc_area = 0.0;
    tpf_prev = tpf;
    fpf_prev = fpf;

  for (item=0; item<n; item++)
    {
      tt+= data[item].label;
      tf-= data[item].label;
      ft+= 1 - data[item].label;
      ff-= 1 - data[item].label;
      
      sens = ((double) tt) / ((double) (tt+tf));
      spec = ((double) ff) / ((double) (ft+ff));

     
      
      tpf  = sens;
      fpf  = 1.0 - spec;

      if ( item < n-1 )
	if ( data[item].score != data[item-1].score )
	  {
	    
	    roc_area+= 0.5*(tpf+tpf_prev)*(fpf-fpf_prev);
	    HERE ("         Sen: %.2f Spe: %.2f ---- %.2f --- %.2f  --- %d ",sens, spec, 0.5*(tpf+tpf_prev)*(fpf-fpf_prev),data[item].score, data[item].label );
	    tpf_prev = tpf;
	    fpf_prev = fpf;
	    
	  }

      if ( item == n-1 )
	{
	  HERE ("*************");
	  roc_area+= 0.5*(tpf+tpf_prev)*(fpf-fpf_prev);
	  HERE ("         Sen: %.2f Spe: %.2f ---- %.2f --- %.2f  --- %d ",sens, spec, 0.5*(tpf+tpf_prev)*(fpf-fpf_prev),data[item].score, data[item].label );
	 
	}
    }
   
  return roc_area;
}


//////7 AUC New
double calculateAUC_2(LabeledData *data, int n) {
    double auc = 0.0, prev_fpr = 0.0, prev_tpr = 0.0;
    int tp = 0, fp = 0, positives = 0, negatives = 0;
    
    for (int i = 0; i < n; i++) {
        if (data[i].label == 1) positives++;
        else negatives++;
    }
    
    for (int i = 0; i < n; )
      {
        int j = i;
        int tie_fp = 0, tie_tp = 0;
        
        // Handle ties by aggregating them
        while (j < n && data[j].score == data[i].score) {
            if (data[j].label == 1) tie_tp++;
            else tie_fp++;
            j++;
        }

        tp += tie_tp;
        fp += tie_fp;

	double tpr = tp / (double)positives;
	double fpr = fp / (double)negatives;
	     
	if (fp == 0 && tp > 0) // Start of ROC curve (0,0) to (0,1)
	  {
	    prev_tpr = tpr;
	  }
	else if (i>0 || (fp>0 && tp>0))
	  {
	    auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2.0;
	    prev_fpr = fpr;
	    prev_tpr = tpr;
	  }
	
        i = j; // Move to the next group of scores
    }
    
    return auc;
}

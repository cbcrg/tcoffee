#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
#include "util_declare.h"

void free_pair_wise()
{
  //Free static allocated memory
  free_proba_pair_wise();
}

/************************************************************************/
/*                                                                      */
/*            CONSTRAINT_LIST                                           */
/*                                                                      */
/*                                                                      */
/************************************************************************/
int *** duplicate_residue_index (int ***r)
{
  int a,b,c;
  int d1,d2,d3;

  int ***nr;

  d1=read_array_size_new(r);
  nr=(int***)vcalloc ( d1, sizeof (int**));
  for (a=0; a<d1; a++)
    {
      d2=read_array_size_new (r[a])-1;
      nr[a]=(int**)vcalloc ( d2+1, sizeof (int*));
      for (b=0; b<d2; b++)
	{
	  d3=read_array_size_new (r[a][b]);
	  for (c=0; c<d3; c++)nr[a][b][c]=r[a][b][c];
	}
    }
  return nr;
}
int *** declare_residue_index (Sequence *S)
{
  int ***r;
  int a,b,c;

  if ( !S)return NULL;
  r=(int***)vcalloc ( S->nseq, sizeof (int**));
  for ( a=0; a<S->nseq; a++)
    {
      r[a]=(int**)vcalloc (S->len[a]+2, sizeof (int*));//The empty terminator makes a scan possible without knowing len
      for ( b=0; b<=S->len[a]; b++)
	{
	  r[a][b]=(int*)vcalloc ( 1, sizeof (int));
	  r[a][b][0]=1;
	}
    }
  return r;
}

Constraint_list * declare_constraint_list_simple ( Sequence *S)
{
  return declare_constraint_list (S, NULL, NULL, 0, NULL, NULL);
}

Constraint_list * declare_constraint_list ( Sequence *S, char *name, int *L, int ne,FILE *fp, int **M)
    {
    Constraint_list *CL;

    CL=( Constraint_list*)vcalloc (1, sizeof ( Constraint_list));


    CL->S=S;
    CL->M=M;

    if ( name!=NULL)
	{
	sprintf ( CL->list_name, "%s", name);

	}
    CL->cpu=1;
    CL->fp=fp;
    if (L)
      {
	HERE ("The USE of L is now Deprecated with Constraint Lists");
	exit (0);
      }
    CL->ne=ne;
    CL->entry_len=LIST_N_FIELDS;
    CL->el_size=sizeof (CLIST_TYPE);
    CL->matrices_list=declare_char(20,20);


    CL->weight_field=WE;
    if ( S)CL->seq_for_quadruplet=(int*)vcalloc ( S->nseq, sizeof (int));
    CL->Prot_Blast=( Blast_param*)vcalloc ( 1, sizeof ( Blast_param));
    CL->DNA_Blast=( Blast_param*)vcalloc ( 1, sizeof ( Blast_param));
    CL->Pdb_Blast=( Blast_param*)vcalloc ( 1, sizeof ( Blast_param));
    CL->TC=(TC_param*)vcalloc (1, sizeof (TC_param));

    //New data structure
    CL->residue_index=declare_residue_index (S);



    return CL;
    }

Constraint_list *free_constraint_list4lib_computation (Constraint_list *CL)
{
  if (!CL)return NULL;

  free_arrayN(CL->residue_index, 3);
  free_int (CL->M, -1);

  vfree (CL);
  return CL;
}
Constraint_list *duplicate_constraint_list4lib_computation (Constraint_list *CL)
{
  Constraint_list *SCL;
  SCL=( Constraint_list*)vcalloc (1, sizeof ( Constraint_list));
  SCL[0]=CL[0];
  SCL->S=CL->S;
  SCL->RunName=CL->RunName;

  SCL->max_L_len=0;
  SCL->M=NULL;
  SCL->ne=0;

  return SCL;
}
Constraint_list *duplicate_constraint_list_soft (Constraint_list *CL)
{
  /*Duplication that does not copy the long lists*/
  return copy_constraint_list (CL,SOFT_COPY);
}
Constraint_list *duplicate_constraint_list (Constraint_list *CL)
     {
     /*Duplicate everything in the constraint_list*/
     return copy_constraint_list (CL,HARD_COPY);
     }
Constraint_list *copy_constraint_list (Constraint_list *CL, int mode)
    {
    Constraint_list *NCL;
    Sequence *S;
    int a, b;



    /*Sequences*/


      S=(mode==HARD_COPY)?duplicate_sequence (CL->S):CL->S;


      if (mode==HARD_COPY)
	NCL=declare_constraint_list (S, NULL, NULL,0, NULL, NULL);
      else
	{
	  NCL=( Constraint_list*)vcalloc ( 1, sizeof (Constraint_list));
	  NCL[0]=CL[0];
	}


      NCL->copy_mode=mode;
      if (mode==SOFT_COPY)NCL->pCL=CL;
      NCL->S=S;
      /*master*/
      if (mode==HARD_COPY && CL->master)
	{NCL->master=(int*)vcalloc ( S->nseq, sizeof(int));
	for ( a=0; a< S->nseq; a++)
	  NCL->master[a]=CL->master[a];
	}
      else if (mode==SOFT_COPY)
	{
	  NCL->seq_for_quadruplet=CL->seq_for_quadruplet;
	}
      NCL->o2a_byte=CL->o2a_byte;

      /*struc List*/
      NCL->STRUC_LIST=(mode==HARD_COPY)?duplicate_sequence (CL->STRUC_LIST):CL->STRUC_LIST;
      sprintf ( NCL->align_pdb_param_file, "%s", CL->align_pdb_param_file);
      sprintf ( NCL->align_pdb_hasch_mode, "%s", CL->align_pdb_hasch_mode);


      NCL->W=(mode==HARD_COPY)?duplicate_weights (CL->W):CL->W;
      NCL->DM=(mode==HARD_COPY)?duplicate_distance_matrix (CL->DM):CL->DM;
      NCL->ktupDM=(mode==HARD_COPY)?duplicate_distance_matrix (CL->ktupDM):CL->ktupDM;
      NCL->RunName=CL->RunName;

      if (  mode==HARD_COPY && CL->translation){NCL->translation=(int*)vcalloc ((CL->S)->nseq, sizeof (int)); for ( a=0; a< (CL->S)->nseq; a++)NCL->translation[a]=CL->translation[a];}
      else{NCL->translation=CL->translation;}

      NCL->out_aln_format=(mode==HARD_COPY)?duplicate_char (CL->out_aln_format, -1, -1):CL->out_aln_format;
      NCL->n_out_aln_format=CL->n_out_aln_format;

    /*Packing Sequence: To use with domain analysis*/
      NCL->packed_seq_lu=(mode==HARD_COPY)?duplicate_int (CL->packed_seq_lu, -1, -1):CL->packed_seq_lu;
    /*DATA*/
      if (CL->fp)(mode==HARD_COPY)?NCL->fp=vtmpfile():CL->fp;

      if ( mode==HARD_COPY)
	{
	  NCL->residue_index=duplicate_residue_index (NCL->residue_index);
	}
      else NCL->residue_index=CL->residue_index;


     if ( mode==HARD_COPY)
       {
	 NCL->M=copy_int ( CL->M,NCL->M);
       }
     else
       NCL->M=CL->M;


    /*List Information*/
      NCL->ne=CL->ne;
      sprintf ( NCL->list_name, "%s", CL->list_name);
      NCL->entry_len=CL->entry_len;
      NCL->el_size=CL->el_size;

    /*Normalisation information*/
      NCL->filter_lib=CL->filter_lib;
      NCL->normalise=CL->normalise;
      NCL->overweight=CL->overweight;
      NCL->max_ext_value=CL->max_ext_value;
      NCL->max_value=CL->max_value;

    /*Pair wise alignment method*/
      NCL->pw_parameters_set=CL->pw_parameters_set;
      NCL->gop=CL->gop;
      NCL->f_gop=CL->f_gop;
      NCL->gep=CL->gep;
      NCL->f_gep=CL->f_gep;

      NCL->nomatch=CL->nomatch;

      NCL->TG_MODE=CL->TG_MODE;
      NCL->F_TG_MODE=CL->F_TG_MODE;

      sprintf ( NCL->dp_mode, "%s", CL->dp_mode);
      NCL->maximise=CL->maximise;
      sprintf ( NCL->matrix_for_aa_group, "%s", CL->matrix_for_aa_group);
      sprintf ( NCL->method_matrix, "%s", CL->method_matrix);

      NCL->diagonal_threshold=CL->diagonal_threshold;
      NCL->ktup=CL->ktup;

      NCL->use_fragments=CL->use_fragments;
      NCL->fasta_step=CL->fasta_step;
      NCL->lalign_n_top=CL->lalign_n_top;
      NCL->sw_min_dist=CL->sw_min_dist;
      NCL->matrices_list=(mode==HARD_COPY)?duplicate_char (CL->matrices_list, -1, -1):CL->matrices_list;
      NCL->n_matrices=CL->n_matrices;

      sprintf (NCL->distance_matrix_mode, "%s", CL->distance_matrix_mode);
      sprintf (NCL->distance_matrix_sim_mode, "%s", CL->distance_matrix_sim_mode);

      sprintf (NCL->tree_mode, "%s", CL->tree_mode);
      NCL->tree_aln=(mode==HARD_COPY)?copy_aln (CL->tree_aln, NULL):CL->tree_aln;
    /*Functions used for dynamic programming and Evaluation*/
      NCL->no_overaln=CL->no_overaln;
      NCL->profile_mode=CL->profile_mode;
      sprintf ( NCL->profile_comparison, "%s",CL->profile_comparison);
      NCL->get_dp_cost=CL->get_dp_cost;
      NCL->evaluate_residue_pair=CL->evaluate_residue_pair;
      NCL->pair_wise=CL->pair_wise;

      NCL->weight_field=CL->weight_field;
      NCL->max_n_pair=CL->max_n_pair;

    /*threading parameters*/
      NCL->Prot_Blast=(mode==HARD_COPY)?duplicate_blast_param ( CL->Prot_Blast):CL->Prot_Blast;
      NCL->DNA_Blast =(mode==HARD_COPY)?duplicate_blast_param ( CL->DNA_Blast):CL->DNA_Blast;
      NCL->Pdb_Blast =(mode==HARD_COPY)?duplicate_blast_param ( CL->Pdb_Blast):CL->Pdb_Blast;
      NCL->TC =(mode==HARD_COPY)?duplicate_TC_param ( CL->TC):CL->TC;

    /*Split parameters*/
      NCL->split=CL->split;
      NCL->split_nseq_thres= CL->split_nseq_thres;
      NCL->split_score_thres= CL->split_score_thres;
    /*Structural status*/
      NCL->check_pdb_status=CL->check_pdb_status;
    /*log*/
      sprintf ( NCL->method_log, "%s",CL->method_log);
      sprintf ( NCL->evaluate_mode, "%s",CL->evaluate_mode);
    /* Gene Prediction*/
      sprintf ( NCL->genepred_score, "%s",CL->genepred_score);

    /*Parameters for domain extraction*/
      NCL->moca=(mode==HARD_COPY)?duplicate_moca ( CL->moca):CL->moca;



    /*Functions for hiding forbiden pairs of residues*/
      /* Copy only for soft_copy*/
      if (mode==SOFT_COPY)
	{
	  NCL->forbiden_pair_list=CL->forbiden_pair_list;
	}
   /*extention properties:*/
      NCL->nseq_for_quadruplet=CL->nseq_for_quadruplet;
      if (mode==HARD_COPY && CL->seq_for_quadruplet)
	{NCL->seq_for_quadruplet=(int*)vcalloc ( S->nseq, sizeof(int));
	for ( a=0; a< S->nseq; a++)
	  NCL->seq_for_quadruplet[a]=CL->seq_for_quadruplet[a];
	}
      else if (mode==SOFT_COPY)
	{
	  NCL->seq_for_quadruplet=CL->seq_for_quadruplet;
	}

   /*extention properties: Do only a soft copy*/
      /* Not To be copied yet */
      if ( mode==SOFT_COPY)
	{
	  NCL->extend_jit=CL->extend_jit;
	  NCL->extend_threshold=CL->extend_threshold;
	  sprintf ( NCL->extend_clean_mode, "%s", CL->extend_clean_mode);
	  sprintf ( NCL->extend_compact_mode, "%s", CL->extend_compact_mode);
	}

    /*Lookup table parameteres*/
      NCL->chunk= CL->chunk;
      /* Do NOT copy NCL->seq_indexed, NCL->start_index, NCL->max_L_len, NCL->chunk*/
      /*
	if ( mode==SOFT_COPY)
	{
	  NCL->seq_indexed=CL->seq_indexed;
	  NCL->start_index=CL->start_index;
	  NCL->end_index=CL->start_index;
	  NCL->max_L_len=CL->max_L_len;
	  }
      */
    /*PDB STRUCTURE ALIGNMENTS*/
      /* Do only a soft copy */
      if ( mode==SOFT_COPY)
	{
	  NCL->T=CL->T;
	}
    /*MISC*/
       NCL->cpu=CL->cpu;
       NCL->local_stderr=CL->local_stderr;
       sprintf (NCL->multi_thread, "%s", CL->multi_thread);
       if (mode==SOFT_COPY)
	 {
	   NCL->comment=CL->comment;
	 }
       else
	 {
	   if (CL->comment)NCL->comment=(char*)vcalloc(strlen (NCL->comment)+1, sizeof (char));
	   sprintf (NCL->comment, "%s", CL->comment);
	 }
       

    return NCL;
    }
Constraint_list *free_constraint_list_full (Constraint_list *CL)
{
  free_sequence (free_constraint_list (CL), -1);
  return NULL;
}
Sequence *free_constraint_list (Constraint_list *CL)
    {
    Sequence *S;
    int a, b;
    Constraint_list *pCL;


    /*Prepare the selective freeing of the CL data structure:
      If the CL has been obtained from copy, every pointer that is identical to the parent CL (CL->pCL)
      will not be saved.
    */


    if ( !CL)return NULL;
    else S=CL->S;

    if ( CL->copy_mode==SOFT_COPY && !CL->pCL)
      {vfree(CL); return S;}
    else if ( CL->copy_mode==SOFT_COPY)
      {

	pCL=CL->pCL;
	CL->residue_index=NULL;

	if ( CL->M                      ==pCL->M                       )CL->M=NULL;

	if (CL->start_index             ==pCL->start_index             )CL->start_index=NULL;
	if (CL->end_index             ==pCL->end_index                 )CL->end_index=NULL;

	if ( CL->fp                     ==pCL->fp                      )CL->fp=NULL;
	if ( CL->matrices_list          ==pCL->matrices_list           )CL->matrices_list=NULL;


	if ( CL->STRUC_LIST             ==pCL->STRUC_LIST              )CL->STRUC_LIST=NULL;
	if ( CL->W                      ==pCL->W                       )CL->W=NULL;
	if ( CL->DM                     ==pCL->DM                      )CL->DM=NULL;
	if ( CL->ktupDM                 ==pCL->ktupDM                      )CL->ktupDM=NULL;


	if ( CL->translation            ==pCL->translation             )CL->translation=NULL;
	if ( CL->moca                   ==pCL->moca                    )CL->moca=NULL;
	if ( CL->Prot_Blast             ==pCL->Prot_Blast              )CL->Prot_Blast=NULL;
	if ( CL->DNA_Blast              ==pCL->DNA_Blast               )CL->DNA_Blast=NULL;
	if ( CL->Pdb_Blast              ==pCL->Pdb_Blast               )CL->Pdb_Blast=NULL;
	if ( CL->seq_for_quadruplet     ==pCL->seq_for_quadruplet      )CL->seq_for_quadruplet=NULL;
	if ( CL->TC                      ==pCL->TC                       )CL->TC=NULL;

      }


    /*End of selective freeing of the CL data structure*/



    if ( CL->residue_index)free_arrayN(CL->residue_index, 3);

    if ( CL->M)free_int (CL->M, -1);
    if ( CL->fp)vfclose (CL->fp);
    if ( CL->matrices_list)free_char(CL->matrices_list,-1);


    if ( CL->start_index)free_int ( CL->start_index,-1);
    if ( CL->end_index)free_int ( CL->end_index,-1);




    if ( CL->STRUC_LIST)free_sequence ( CL->STRUC_LIST, (CL->STRUC_LIST)->nseq);
    if ( CL->W)free_weights (CL->W);

    CL->DM=free_distance_matrix (CL->DM);
    CL->ktupDM=free_distance_matrix (CL->ktupDM);

    if ( CL->translation)vfree(CL->translation);
    if ( CL->moca)free_moca (CL->moca);
    if ( CL->Prot_Blast)free_blast_param ( CL->Prot_Blast);
    if ( CL->DNA_Blast) free_blast_param ( CL->DNA_Blast);
    if ( CL->Pdb_Blast) free_blast_param ( CL->Pdb_Blast);
    if ( CL->TC) free_TC_param ( CL->TC);

    if (CL->seq_for_quadruplet)vfree (CL->seq_for_quadruplet);

    vfree(CL);
    return S;
    }

Distance_matrix * free_distance_matrix ( Distance_matrix *DM)
{
  if (!DM)return NULL;
  free_int ( DM->similarity_matrix,-1);
  free_int ( DM->distance_matrix,-1);
  free_int ( DM->score_similarity_matrix,-1);
  vfree (DM);
  return NULL;
}
Distance_matrix * duplicate_distance_matrix ( Distance_matrix *DMin)
{
  Distance_matrix *DM;
  if (!DMin) return NULL;

  DM=(Distance_matrix*)vcalloc ( 1, sizeof (Distance_matrix));
  DM->similarity_matrix=duplicate_int ( DMin->similarity_matrix, -1, -1);
  DM->distance_matrix=duplicate_int ( DMin->distance_matrix, -1, -1);
  DM->score_similarity_matrix=duplicate_int ( DMin->score_similarity_matrix, -1, -1);
  return DM;
}

/************************************************************************/
/*                                                                      */
/*            MOCA Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Moca * duplicate_moca ( Moca *m)
      {
	Moca *nm;

	if ( m==NULL)return m;

	nm=(Moca*)vcalloc ( 1, sizeof (Moca));

	nm->moca_scale=m->moca_scale;
	nm->evaluate_domain=m->evaluate_domain;
	nm->moca_threshold=m->moca_threshold;
	nm->cache_cl_with_domain=m->cache_cl_with_domain;
	if ( m->forbiden_residues)nm->forbiden_residues=copy_int  (m->forbiden_residues,nm->forbiden_residues);
	nm->make_nol_aln=m->make_nol_aln;


	return nm;
      }
Moca * free_moca ( Moca *m)
      {
	if ( m->forbiden_residues)free_int ( m->forbiden_residues, -1);
	vfree ( m);
	return NULL;
      }
/************************************************************************/
/*                                                                      */
/*            TC_param Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
TC_param * duplicate_TC_param ( TC_param*B)
{
  TC_param *N;
  N=(TC_param*)vcalloc (1, sizeof ( TC_param));
  memcpy(B, N, sizeof(TC_param));
  return N;
  }
TC_param * free_TC_param ( TC_param*B)
{
  vfree (B);
  return NULL;
}
/************************************************************************/
/*                                                                      */
/*            Blast_param Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Blast_param * duplicate_blast_param ( Blast_param*B)
{
  Blast_param *N;
  N=(Blast_param*)vcalloc (1, sizeof ( Blast_param));
  sprintf ( N->blast_server, "%s", B->blast_server);
  sprintf ( N->db, "%s", B->db);
  N->min_id=B->min_id;
  N->max_id=B->min_id;
  N->min_cov=B->min_cov;
  return N;
}
Blast_param * free_blast_param ( Blast_param*B)
{
  vfree (B);
  return NULL;
}

/************************************************************************/
/*                                                                      */
/*            PDB Functions                                             */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Structure* declare_structure ( int n, char **array)
    {
    Structure *S;
    int a;

    S=(Structure*)vcalloc (1, sizeof (Structure));
    S->n_fields=1;
    S->nseq=n;

    S->struc=(int***)vcalloc ( n, sizeof (int**));
    S->len=(int*)vcalloc ( n, sizeof (int));
    for ( a=0; a< n; a++)
        {
	S->len[a]=strlen(array[a]);
	S->struc[a]=declare_int ( strlen ( array[a])+2, 1);
	}
    return S;
    }

Structure *extend_structure ( Structure *S)
    {
    int a, b;


    for ( a=0; a< S->nseq; a++)
        {
	for ( b=0; b< S->len[a]; b++)
		S->struc[a][b]=(int*)vrealloc ( S->struc[a][b],( S->n_fields+1)*sizeof (int));
	}
    S->n_fields++;
    return S;
    }
Sequence * reset_sequence_len (Sequence *S)
{
  int min,max,a,l;

  if ( !S || !S->nseq)return 0;
  max=min=strlen (S->seq[0]);
  for (a=0; a<S->nseq; a++)
    {
      l=strlen (S->seq[a]);
      min=MIN(l,min);
      max=MAX(l,max);
      S->len[a]=l;
    }
  S->max_len=max;
  S->min_len=min;
  return S;
}
Sequence * declare_sequence ( int min, int max, int nseq)
    {
    Sequence *LS;



    LS=( Sequence*)vcalloc (1, sizeof ( Sequence));

    LS->seq_comment=declare_char ( nseq,COMMENT_SIZE);
    LS->aln_comment=declare_char ( nseq,COMMENT_SIZE);

    LS->file=declare_char( nseq,STRING+1);
    LS->seq=declare_char ( nseq, max+1);
    LS->name=declare_char( nseq,MAXNAMES+1);

    LS->len=(int*)vcalloc ( nseq, sizeof (int));
    LS->max_len=max;
    LS->min_len=min;
    LS->nseq=nseq;
    LS->max_nseq=nseq;
    LS->type=(char*)vcalloc(30, sizeof (char));
    LS->T=(Template**)declare_arrayN(2, sizeof (Template), nseq, 1);


    LS->dc=declare_int (nseq, 2);
    return LS;
    }


Sequence * realloc_sequence   (Sequence *OUT, int new_nseq, int max_len)
{


	if ( new_nseq<OUT->max_nseq)
		return OUT;

	OUT->min_len =MIN(OUT->min_len,max_len);
	OUT->max_len =MAX(OUT->max_len,max_len);

	OUT->seq_comment =new_realloc_char ( OUT->seq_comment, new_nseq,COMMENT_SIZE);
	OUT->aln_comment =new_realloc_char ( OUT->aln_comment, new_nseq,COMMENT_SIZE);
	OUT->name    =new_realloc_char ( OUT->name,    new_nseq,MAXNAMES+1);
	OUT->seq     =new_realloc_char ( OUT->seq,     new_nseq,OUT->max_len+1);


	if (OUT->genome_co != NULL)
		OUT->genome_co =(Genomic_info*)vrealloc(OUT->genome_co, new_nseq * sizeof(Genomic_info));

	OUT->file=new_realloc_char ( OUT->file,    new_nseq,STRING+1);
	OUT->len=(int*)vrealloc( OUT->len,     (new_nseq+1)*sizeof (int));

	OUT->T=(Template**)realloc_arrayN (2, (void **)OUT->T,sizeof (Template), new_nseq, 1);
	OUT->dc=(int **)realloc_arrayN (2, (void **)OUT->dc,sizeof (int), new_nseq, 2);

	OUT->max_nseq=new_nseq;
	return OUT;
}
Sequence * duplicate_sequence (Sequence *S )
{
	Sequence *LS;
	int a, b;
	char*tmp=vtmpnam(NULL);
	FILE *fp;
	
	fp=vfopen (tmp, "w");
	for (a=0; a<S->nseq; a++)
	  fprintf (fp, ">%s %s\n%s\n", S->name[a],S->seq_comment[a], S->seq[a]);
	vfclose (fp);
	
	LS=get_fasta_sequence(tmp, NULL);

	for (a=0; a<S->nseq; a++)
	  {
	    LS->file[a]=csprintf (LS->file[a],"%s", S->file[a]);
	    LS->aln_comment[a]=csprintf (LS->aln_comment[a],"%s", S->aln_comment[a]);
	    LS->seq_comment[a]=csprintf (LS->seq_comment[a],"%s", S->seq_comment[a]);
	    LS->dc[a][0]=S->dc[a][0];
	    LS->dc[a][1]=S->dc[a][1];
	    LS->len[a]=S->len[a];
	    LS->T[a][0]=S->T[a][0];
	  }
	sprintf (LS->type,"%s", S->type);
	sprintf (LS->template_file,"%s", S->template_file);
	if (S->W)LS->W=duplicate_weights (S->W);
	
	LS->max_len=S->max_len;
	LS->min_len=S->min_len;

	LS->blastdb=S->blastdb;
	LS->max_nseq=S->nseq;
	LS->blastdbS=S->blastdbS;
	LS->MasterS=S->MasterS;

	return LS;
}
Sequence * duplicate_sequence_old (Sequence *S )
{
	Sequence *LS;
	int a, b;
// 	printf("DUP\n");
	if (S==NULL)return S;
	LS=declare_sequence (S->min_len, S->max_len, S->nseq);

	if (S->genome_co != NULL)
		LS->genome_co =(Genomic_info*) vcalloc(S->nseq, sizeof(Genomic_info));
	for (b=0, a=0; a<S->nseq; ++a)
	  {
	    if (S->seq && S->seq[a])
	      {

		sprintf ( LS->file[b], "%s", S->file[a]);

		if ( S->seq_comment && S->seq_comment[a])
		  sprintf ( LS->seq_comment[b], "%s", S->seq_comment[a]);
		if ( S->aln_comment && S->aln_comment[a])
		  sprintf ( LS->aln_comment[b], "%s", S->aln_comment[a]);

		if ( S->seq && S->seq[a])
		  sprintf  ( LS->seq[b],  "%s", S->seq[a]);

		if ( S->name&& S->name[a])
		  sprintf ( LS->name[b], "%s", S->name[a]);

		if (S->genome_co != NULL)
			{
				unsigned int len;

				Genomic_info *tmp, *tmp_ori;
				unsigned int i;
// 				for ( i=0; i< S->nseq; a++)
// 				{
				tmp = &(LS->genome_co[a]);
				tmp_ori = &(S->genome_co[a]);
				tmp->start = tmp_ori->start;
				tmp->strand = tmp_ori->strand;
				tmp->end = tmp_ori->end;
				tmp->seg_len = tmp_ori->seg_len;
// 				printf("A: %i\n", a);
				tmp->seg_name =(char*) vcalloc(strlen(tmp_ori->seg_name)+1, sizeof(char));
				strcpy(tmp->seg_name, tmp_ori->seg_name);
// 				printf("NAME %s\n", tmp->seg_name);
// 				printf("-%s %i-\n", tmp_ori->seg_name, tmp_ori->start);
// 				}
			}


			LS->dc[b][0]=S->dc[a][0];
			LS->dc[b][1]=S->dc[a][1];
			LS->len[b]=S->len[a];
			LS->T[b][0]=S->T[a][0];
			b++;

		}
	}

	LS->max_len=S->max_len;
	LS->min_len=S->min_len;
	LS->nseq=b;

	if (S->W)LS->W=duplicate_weights (S->W);
	LS->blastdb=S->blastdb;
	sprintf ( LS->type, "%s", S->type);
	sprintf ( LS->template_file, "%s", S->template_file);
	LS->max_nseq=S->nseq;
	LS->blastdbS=S->blastdbS;
	LS->MasterS=S->MasterS;
	return LS;
}

void free_sequence ( Sequence *LS, int nseq)
	{


	if ( !LS) return;
	
	free_char ( LS->file, -1);
	
	free_char ( LS->seq_comment, -1);
	free_char ( LS->aln_comment, -1);

	free_char ( LS->seq, -1);
	free_char ( LS->name,-1);
	
	free_int  (LS->dc, -1);

	free_arrayN((void*)LS->T, 2);
	vfree (LS->type);
	vfree (LS->len);
	free_weights (LS->W);
	//Don't free LS->blastdb;
	vfree (LS);

	}
/************************************************************************/
/*                                                                      */
/*            Weights Functions                                         */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Weights* declare_weights ( int nseq)
	{
	Weights *W;

	W=(Weights*)vcalloc ( 1, sizeof ( Weights));
	W->comments=(char*)vcalloc ( 1000, sizeof (char));
	W->nseq=nseq;
	W->mode=(char*)vcalloc (FILENAMELEN, sizeof (char));
	W->seq_name= declare_char ( W->nseq*2, 200);
	W->PW_SD=declare_float ( W->nseq, W->nseq);
	W->PW_ID=declare_float ( W->nseq, W->nseq);
	W->SEQ_W=(float*)vcalloc ( W->nseq, sizeof ( float));
	return W;
	}
Weights* duplicate_weights (Weights *W)
     {
       Weights *NW;
       int a, b, c;

       NW=declare_weights (W->nseq);
       sprintf ( NW->comments, "%s", W->comments);
       sprintf ( NW->mode, "%s", W->mode);
       for (a=0, c=0; a< W->nseq; a++)
	 {
	   if ( W->seq_name[a])
	     {
	       sprintf ( NW->seq_name[c], "%s", W->seq_name[a]);
	       NW->SEQ_W[c]=W->SEQ_W[a];
	       for(b=0; b< W->nseq; b++)
		 {
		   NW->PW_SD[c][b]=W->PW_SD[a][b];
		   NW->PW_ID[c][b]=W->PW_ID[a][b];
		 }
	       c++;
	     }
	 }
       return NW;
     }
Weights* free_weights ( Weights* W)
	{

	if ( !W)return NULL;

	vfree(W->comments);


	vfree(W->mode);
	free_char(W->seq_name, -1);
	free_float(W->PW_SD,-1);
	free_float(W->PW_ID, -1);
	vfree(W->SEQ_W);
	vfree(W);
	return NULL;
	}


Alignment* copy_aln ( Alignment *A, Alignment *B)
        {
	  int a, b;
	  int nnseq;
	  int nlen, mlen;
	  /*	  c[100]=10;*/



	  if ( A==NULL){free_aln(B); return NULL;}
	  	  
	  nnseq=A->nseq;
	  nlen=A->len_aln;
	  for (a=0; a<A->nseq; a++)
	    {
	      int cl=(A->seq_al[a])?strlen(A->seq_al[a]):0;
	      if (cl>nlen)nlen=cl;
	    }
	  nlen++;
	  
	  if (B)
	    B=realloc_alignment2 (B, nnseq+1, nlen);
	  else
	    B=declare_aln2 (nnseq+1, nlen);
	  B->S=A->S;


	  /*SIZES*/
	  B->max_len=A->max_len;
	  B->min_len=A->min_len;
	  B->declared_len=nlen;
	  B->max_n_seq=nnseq+1;

	  B->nseq=A->nseq;
	  B->len_aln=A->len_aln;


/*sequence Information*/
	    if ( A->generic_comment)
	      {
		vfree(B->generic_comment);
		B->generic_comment=(char*)vcalloc (strlen(A->generic_comment)+1, sizeof (char));
		sprintf ( B->generic_comment, "%s", A->generic_comment);
	      }
	    if ( (A->S)==NULL){vfree (B->len); B->len=(int*)vcalloc ( A->max_n_seq, sizeof (int));}
	    ga_memcpy_int ( A->len, B->len, B->nseq);

	    
	    B->seq_comment=copy_char ( A->seq_comment,  B->seq_comment);
	    B->aln_comment=copy_char ( A->aln_comment,  B->aln_comment);

	    B->name=copy_char ( A->name,     B->name);

	    B->file=copy_char ( A->file,     B->file);
	    B->tree_order=copy_char ( A->tree_order,     B->tree_order);
	    B->expanded_order=A->expanded_order;
	    
	    //Copy seq_al
	    
	    for (a=0; a<A->nseq; a++)
	      {
		if (A->seq_al[a])B->seq_al[a]=csprintf (B->seq_al[a], "%s", A->seq_al[a]);
	      }
	    
	    //B->seq_al[A->nseq]=csprintf (B->seq_al[A->nseq], "%s", A->seq_al[A->nseq-1]);
	   
	    //This is meant to take care of the consensus
	    //
	    
	    
	    B->order=copy_int  ( A->order,    B->order);
	    B->S=A->S;
	    if (A->seq_cache)
	        {
		  B->seq_cache=copy_int  ( A->seq_cache,    B->seq_cache);
		}

	    if (A->cdna_cache)
	        {
		  B->cdna_cache=copy_int  ( A->cdna_cache,    B->cdna_cache);
		}

	    B->P=copy_profile (A->P);

	    B->Dp_result=A->Dp_result;

/*Score*/

	    if ( (A->S)==NULL){vfree (B->score_seq); B->score_seq=(int*)vcalloc ( A->max_n_seq, sizeof (int));}
	    ga_memcpy_int(  A->score_seq,B->score_seq,B->nseq);
	    B->score_res=A->score_res;

	    B->score_aln=A->score_aln;
	    B->score=A->score;
	    B->ibit=A->ibit;
	    B->cpu=A->cpu;
	    B->finished=A->finished;
	    if (A->dm)
	      {
		int i,j;
		B->dm=declare_double (A->nseq, A->nseq);
		for (i=0; i<A->nseq; i++)
		  for (j=0; j<A->nseq; j++)
		    B->dm[i][j]=A->dm[i][j];
		B->dm=A->dm;
	      }

/*Output Options*/
	    B->output_res_num=A->output_res_num;
	    B->residue_case=A->residue_case;
	    B->expand=A->expand;

	    B->CL=A->CL;
	    B->random_tag=A->random_tag;

/*Make the function Recursive */
	    if ( A->A)
	      {
		B->A=copy_aln (A->A, NULL);
	      }
	    else B->A=NULL;
	    /*Deal with Trees*/
	    //if (A->Tree)
	    //B->Tree=copy_aln (A->Tree, NULL);
	    if (A->tname){B->tname=(char*)vcalloc (strlen (A->tname)+1, sizeof (char)); sprintf (B->tname, "%s", A->tname);}
	    return B;
	}

Alignment* shrink_aln ( Alignment *A, int nseq, int *list)
        {
	Alignment *B=NULL;
	int a,seq;

	B=copy_aln (A, B);
	for ( a=0; a< nseq; a++)
	    {
	    seq=list[a];
	    sprintf ( A->seq_comment[a], "%s",B->seq_comment[seq]);
	    sprintf ( A->aln_comment[a], "%s",B->aln_comment[seq]);

	    sprintf ( A->seq_al [a], "%s",B->seq_al [seq]);
	    A->order[a][0]=B->order[seq][0];
	    A->order[a][1]=B->order[seq][1];
	    A->order[a][2]=B->order[seq][2];
	    A->order[a][3]=B->order[seq][3];
	    A->order[a][4]=B->order[seq][4];

	    A->score_seq[a]=B->score_seq[seq];
	    A->len[a]=B->len[seq];
	    }
	A->nseq=nseq;
	A->len_aln=strlen (A->seq_al[0]);
	free_aln (B);
	return A;
	}
Alignment* extract_sub_aln2 ( Alignment *B, int ns, char **ls)
        {
	  int *list;
	  Alignment *A;

	  list=name_array2index_array(ls, ns, B->name, B->nseq);
	  A=extract_sub_aln ( B,ns, list);
	  vfree (list);
	  return A;
	}
Alignment* extract_sub_aln ( Alignment *B, int nseq, int *list)
        {
	Alignment *A=NULL;
	int a,b,n,seq;

	A=declare_aln2(nseq, B->len_aln+1);
	for ( n=0,a=0; a< nseq; a++)
	    {
	    seq=list[a];
	    if ( seq==-1)continue;
	    else n++;
	    sprintf ( A->seq_comment[a], "%s",B->seq_comment[seq]);
	    sprintf ( A->aln_comment[a], "%s",B->aln_comment[seq]);
	    sprintf ( A->name[a], "%s",B->name[seq]);


	    for (b=0; b<=B->len_aln; b++)A->seq_al [a][b]=B->seq_al [seq][b];
	    A->order[a][0]=B->order[seq][0];
	    A->order[a][1]=B->order[seq][1];
	    A->order[a][2]=B->order[seq][2];
	    A->order[a][3]=B->order[seq][3];
	    A->order[a][4]=B->order[seq][4];

	    A->score_seq[a]=B->score_seq[seq];
	    A->len[a]=B->len[seq];
	    }
	A->nseq=n;
	A->len_aln=B->len_aln;
	return A;
	}

Alignment *declare_aln2 ( int nseq, int len)
        {
	  Sequence *S;
	  Alignment *A;

	  S=(Sequence*)vcalloc ( 1, sizeof ( Sequence));
	  S->nseq=nseq;
	  S->max_len=len;

	  A=declare_aln (S);
	  A->S=NULL;
	  vfree(S);
	  return A;
	}



Alignment *declare_aln ( Sequence *S){return declare_Alignment(S);}

Alignment *declare_Alignment ( Sequence *S)
	{
	Alignment *LA;
	int a;

	/*ordre:
	  [x][0]= which is the xth seq of aln
	  [x][1]= how many deleted residues before the first one
	*/


	LA=( Alignment*)vcalloc (1, sizeof ( Alignment));
	aln_stack (LA, DECLARE_ALN);
	if ( S==NULL)
	    {
	      LA->declared_len=MAX_LEN_ALN;
	      LA->max_n_seq=MAX_N_SEQ;
	    }
	else
	  {
	    LA->declared_len=2*S->max_len+1;
	    LA->max_n_seq=S->nseq+1;
	  }
	LA->S=S;


	LA->seq_comment=declare_char (LA->max_n_seq, COMMENT_SIZE);
	LA->aln_comment=declare_char (LA->max_n_seq, COMMENT_SIZE);


	LA->seq_al=declare_char ( LA->max_n_seq,LA->declared_len );
	LA->name=declare_char (LA->max_n_seq, MAXNAMES+1);


	LA->file=declare_char (LA->max_n_seq, STRING);
	LA->tree_order=declare_char (LA->max_n_seq, STRING);
	LA->order= declare_int (LA->max_n_seq , 5);
	
	//order[a][0]: sequence index in S
	//order[a][1]: offset of the sequence
	//order[a][2]: used by sw_gotoh_pair_wise
	//order[a][3]: used by sw_gotoh_pair_wise
	//order[a][4]: weight, -1
	
	LA->score_seq= (int*)vcalloc (LA->max_n_seq, sizeof (int));

	for ( a=0; a< LA->max_n_seq; a++)LA->order[a][0]=a;
	LA->len_aln=0;
	LA->score_aln=0;
	LA->len=(int*)vcalloc (LA->max_n_seq, sizeof (int));

	if (S && S->name)for ( a=0; a<S->nseq; a++)
	  {
	    sprintf ( LA->name[a], "%s", S->name[a]);

	  }

	return LA;

	}

Alignment *Realloc_Alignment4nseq ( Alignment *A, int nseq)
{
  int a;
  nseq++;
  A->seq_comment=(char **)vrealloc(A->seq_comment, nseq*sizeof (char*)); 
  A->aln_comment=(char **)vrealloc(A->aln_comment, nseq*sizeof (char*)); 
 
  A->seq_al=(char **)vrealloc(A->seq_al, nseq*sizeof (char*)); 
  A->name=(char **)vrealloc(A->name, nseq*sizeof (char*)); 

  A->file=(char **)vrealloc(A->file, nseq*sizeof (char*)); 
  A->tree_order=(char **)vrealloc(A->tree_order, nseq*sizeof (char*)); 
  
  
  A->score_seq= (int*)vrealloc (A->score_seq,sizeof (int)*nseq);
  A->len= (int*)vrealloc (A->len,sizeof (int)*nseq);
  A->order=(int **)vrealloc(A->order, nseq*sizeof (int*)); 
  
  for (a=0; a<nseq; a++)
    {
      if (!A->aln_comment[a])A->aln_comment[a] =(char*)vcalloc (COMMENT_SIZE,sizeof (char));
      if (!A->seq_comment[a])A->seq_comment[a] =(char*)vcalloc (COMMENT_SIZE,sizeof (char));
      if (!A->file[a])       A->file[a]        =(char*)vcalloc (STRING,      sizeof (char));
      if (!A->tree_order[a]) A->tree_order[a]  =(char*)vcalloc (STRING,      sizeof (char));
      if (!A->order[a])      A->order[a]       =(int* )vcalloc (5,           sizeof (int));
      if (!A->name[a])       A->name[a]        =(char*)vcalloc (MAXNAMES+1,  sizeof (char));
    
    }

  A->max_n_seq=nseq;
  return A;
}
Alignment * realloc_aln ( Alignment *A, int new_len){return realloc_alignment(A, new_len);}
Alignment * realloc_alignment ( Alignment *A, int new_len)
	{
	if (A==NULL)A=declare_Alignment (NULL);

	return realloc_alignment2( A, A->max_n_seq,new_len);
	}

Alignment * realloc_aln2 ( Alignment *A, int n_nseq, int n_len)
{
  if (n_len==0 && n_nseq==0)return A;
  else if (n_nseq && !n_len)return Realloc_Alignment4nseq(A,n_nseq);
  else   return realloc_alignment2(A, n_nseq, n_len);
}




Alignment * realloc_alignment2 ( Alignment *A, int n_nseq, int n_len)
	{
	int a;
	int len, nseq;
	int delta_len, delta_nseq;

        if ( A==NULL) A=declare_Alignment(NULL);

	n_len++;
	n_nseq++;

	len=A->declared_len;
	nseq=A->max_n_seq;

	n_len=MAX(len, n_len);
	n_nseq=MAX(nseq,n_nseq);
	delta_nseq=MAX(0,n_nseq-nseq);
	delta_len =MAX(0,n_len-len);

	if ( delta_nseq<=0 && delta_len<=0)return A;


	else
	    {
	        A->len=(int*)vrealloc( A->len, sizeof (int)*n_nseq);
		for (a=nseq; a< n_nseq; a++)A->len[a]=0;

		A->declared_len =n_len;
		A->max_n_seq    =n_nseq;


		A->seq_comment=new_realloc_char ( A->seq_comment, n_nseq, -1);
		A->aln_comment=new_realloc_char ( A->aln_comment, n_nseq, -1);

		A->name   =new_realloc_char ( A->name, n_nseq, -1);


		A->file   =new_realloc_char ( A->file, n_nseq, -1);

		A->tree_order   =new_realloc_char ( A->tree_order, n_nseq, -1);
		A->seq_al =new_realloc_char ( A->seq_al, n_nseq, n_len);
		A->order  =new_realloc_int  ( A->order, n_nseq, -1);

		if ( A->seq_cache) A->seq_cache=new_realloc_int  ( A->seq_cache, n_nseq,n_len);
		if ( A->cdna_cache)A->cdna_cache=new_realloc_int  ( A->cdna_cache, n_nseq,n_len);


		A->score_seq=(int*)vrealloc( A->score_seq, sizeof (int)*(n_nseq));
		for ( a=nseq; a< n_nseq; a++)A->score_seq[a]=0;



	    }
	return A;
	}


long aln_stack (Alignment *A, int mode)
{
  static long *list;
  static int size;
  static int max_size;


  if (A==NULL) return 0;
  else if ( mode==DECLARE_ALN)
    {
      if ( size==max_size)
	{
	  max_size+=1000;
	  list=(long int*)vrealloc (list, max_size*sizeof (long));
	}
      list[size++]=(long)A;
      return 0;
    }
  else if (mode==FREE_ALN)
    {
      int a, b;
      for (a=0; a<size; a++)
	{
	  if (list[a]==(long)A)
	    {
	      for (b=a+1; b<size;b++)
		{
		  list[b-1]=list[b];
		}
	      list[b-1]=0;
	      size--;
	      return 1;
	    }
	}
      return 0;
    }
  else if ( mode==EXTRACT_ALN)
    {
      return list[size--];
    }
  else
    {
      printf_exit (EXIT_FAILURE, stderr, "ERROR: Unknown mode for aln_stack");
      return 0;
    }
}

Sequence* free_aln ( Alignment *LA){ return free_Alignment(LA);}
Alignment *free_data_in_aln (Alignment *A)
{
  //Frees only the sequence data (keeps profile information)
  A->seq_al=free_char (A->seq_al, -1);
  A->seq_comment=free_char (A->seq_comment, -1);
  A->aln_comment=free_char (A->aln_comment, -1);
  A->name=free_char (A->name, -1);
  A->expanded_order=free_char (A->expanded_order, -1);

  A->order=free_int (A->order, -1);
  A->seq_cache=free_int (A->seq_cache, -1);
  A->cdna_cache=free_int (A->cdna_cache, -1);
  A->score_res=free_int (A->score_res, -1);
  free_sequence (A->S, -1);
  if (A->Tree)free_Alignment (A->Tree);
  if (A->RepColList)free_int (A->RepColList, -1);
  if (A->dm)free_double (A->dm,-1);
  A->S=NULL;
  return A;

}
Sequence* free_Alignment ( Alignment *LA)
	{
	  /* Does not free the A->S field (sequences of A)*/

	  if (!LA) return NULL;
	  Sequence *S;
	  //aln_stack checks the alignment has not already been freed
	  if ( LA==NULL || !aln_stack(LA,FREE_ALN)){return NULL;}

	  S=LA->S;
	  free_char ( LA->file, -1);
	  free_char ( LA->seq_al, -1);
	  free_int  ( LA->seq_cache, -1);
	  free_int  ( LA->cdna_cache, -1);
	  free_char ( LA->name,-1);

	  free_char ( LA->tree_order,-1);
	  vfree ( LA->generic_comment);
	  free_char ( LA->seq_comment, -1);
	  free_char ( LA->aln_comment, -1);

	  free_int  ( LA->order, -1);

	  free_double (LA->dm, -1),
	  vfree ( LA->score_seq);
	  vfree ( LA->len);

	  free_profile (LA->P);
	  if ( LA->A){free_Alignment (LA->A);LA->A=NULL;}
	  if ( LA->Tree){free_Alignment (LA->Tree);LA->Tree=NULL;}
	  if ( LA->RepColList)free_int (LA->RepColList, -1);
	  if (LA->tname)vfree (LA->tname);
	  vfree ( LA);
	  return S;
	}

Alignment * update_aln_random_tag ( Alignment *A)
{
  static int tag;
  if ( !A) return A;

  A->random_tag=++tag;
  return A;
}

Profile   *copy_profile   (Profile *P1)
{

  Profile *P;

  if ( !P1) return NULL;
  P=declare_profile ( P1->alphabet, P1->max_len);
  P->count=copy_int (P1->count, P->count);
  P->count2=copy_int (P1->count2, P->count2);
  P->count3=copy_int (P1->count3, P->count3);

  return P;

}


Profile   *declare_profile(char *alphabet, int len)
{
  Profile *P;
  P=( Profile*)vcalloc ( 1, sizeof ( Profile));
  P->alp_size=strlen(alphabet);
  P->max_len=len;
  P->alphabet=(char*)vcalloc ( strlen (alphabet)+2, sizeof (char));
  sprintf ( P->alphabet, "%s", alphabet);

  P->count=declare_int( P->alp_size+2, len);
  P->count2=declare_int(100, len);
  P->count3=declare_int(100, len);

  return P;
}
Profile * free_profile ( Profile *P)
{
  if (!P) return NULL;
  else
    {
      vfree (P->alphabet);
      free_int ( P->count, -1);
      free_int ( P->count2, -1);
      vfree (P);
    }
  return NULL;
}


/************************************************************************/
/*                                                                      */
/*             ALLOCATION                                               */
/*                                                                      */
/*                                                                      */
/************************************************************************/


double alloc_mem;
double max_mem;
double tot_mem;
Memcontrol *memlast;


void mem_profile (char *msg)
{
 
  if (getenv ("PROFILE_MEM_4_TCOFFEE"))print_mem_usage(stdout,msg);
}
FILE* print_mem_usage (FILE *fp, char *comment)
{
  fprintf ( fp, "# %s Memory Usage: Current= %.3f Mb, Max= %.3f Mb\n", (comment)?comment:"",(float)((float)alloc_mem/(1024*1024)),(float)((float)tot_mem/(1024*1024)) );
  return fp;
}
void set_max_mem (int m)
{
  max_mem=m*1024*1024;
}

int verify_memory (int s)
{
  static int flushed;
    
  alloc_mem+=s;
  
  tot_mem=(alloc_mem>tot_mem)?alloc_mem:tot_mem;

  if (max_mem && alloc_mem>max_mem && !flushed)
    {
      flushed=1;
      fprintf (stderr, "\nERROR: Current Memory usage of %d Mb exceeds allowed maximum (%d Mb). [FATAL:%s]\n",(int)(alloc_mem/(1024*1024)),(int)(max_mem/(1024*1024)), PROGRAM);
      myexit (EXIT_FAILURE);
    }
  else
    return 1;
  return 0;

}

int my_assert ( void *p, int index)
{
  static int warning;

  if (!warning)
    {
      fprintf ( stderr, "\n****************************************************************\n");
      fprintf ( stderr, "\n          DEBUG MODE [Rebuild For Better Performances]          \n");
      fprintf ( stderr, "\n*****************************************************************\n");
      warning=1;
    }

  if ( !is_dynamic_memory(p)) return 1;
  else if ( read_array_size_new (p)<=index)
    {
      fprintf ( stderr, "\nFaulty Allocation: Size=%d Access=%d\n", read_array_size (p,0),index);
      return 0;
    }
  else
    {
      return 1;
    }
}



void * vmalloc ( size_t size)
	{
	void * x;
	Memcontrol *M;

	verify_memory (size+2*sizeof (Memcontrol));
	
	if ( size==0)
	    return NULL; /*crash ("\n0 bytes in vmalloc\n");*/
	else
	    {

	      x= malloc (size + 2*sizeof (Memcontrol));
	      //x=dlmalloc (size + 2*sizeof (Memcontrol));
	    if ( x==NULL)
		{
		  printf_exit (EXIT_FAILURE,stderr, "\nFAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)\n");

		}
	    else
	      {
		M=(Memcontrol*)x;
		M[0].size=size;
		M[0].size_element=0;
		sprintf ( M[0].check, "dy");
		M+=2;
		x=M;
		return x;
	      }
	    }
	return NULL;
	}



void *vcalloc (size_t nobj, size_t size)
{
  return sub_vcalloc (nobj,size, MEMSET0);
}
void *vcalloc_nomemset ( size_t nobj, size_t size)
{
  return sub_vcalloc (nobj, size, NO_MEMSET0);
}
void *sub_vcalloc ( size_t nobj, size_t size, int MODE)
	{
	void *x;
	Memcontrol *M;

	if ( nobj<=0 || size<=0)return NULL;/*crash ("\n0 bytes in vmalloc\n");*/
	else x=vmalloc (nobj*size);


	M=(Memcontrol*)x;
	M-=2;M[0].size_element=size;
	M+=2;
	x=M;

	if ( x==NULL)
		{
		crash ( "\nFAILED TO ALLOCATE REQUIRED MEMORY (vcalloc)\n");
		return NULL;
		}
	else
	  {
	    if ( MODE==MEMSET0)
	      {
		x=memset (x,0, nobj*size);
	      }
	    else
	      {
		if (nobj)x=memset (x, 0, size);
	      }
	    return x;
	  }
	}
void *vreallocg ( void *p, size_t size, int set, int resize)
	{
	  //makes it possibl�e to realloc w/o reintilaizing and w/o freeing left over: MEMSET or NOMEMSET and RESIZE or NORESIZE
	void *x;
	Memcontrol *M;
	size_t i_size;
	int a;
	void *ip=p;
	
	if ( p==NULL)
	  {
	    x=vmalloc (size);
	    if (set==MEMSET)memset (x, 0, size);
	    return x;
	  }
	else
	  {
	    M=(Memcontrol*)p;
	    M-=2;
	    i_size=M[0].size;
	    p=M;


	    if ( size<=0){return NULL;vfree (p);return NULL;}
	    else if (resize==NORESIZE && size<=i_size)return ip;
	    else
	      {
		verify_memory (size - i_size);
		x=realloc ( p, size+2*sizeof(Memcontrol));

		if ( x==NULL){crash ( "\nFAILED TO ALLOCATE REQUIRED MEMORY (realloc)\n");return NULL;}
		M=(Memcontrol*)x;
		M[0].size=size;
		M+=2;
		x=M;
		if (set==MEMSET)for ( a=i_size; a< size; a++)((char*)x)[a]=0;
		return x;
	      }
	  }
	return NULL;
	}

void *vrealloc ( void *p, size_t size)
	{
	void *x;
	Memcontrol *M;
	size_t i_size;
	int a;


	if ( p==NULL)
	  {
	    return vcalloc (size/sizeof(char), sizeof (char));
	  }
	else
	  {
	    M=(Memcontrol*)p;
	    M-=2;
	    i_size=M[0].size;
	    p=M;


	    if ( size<=0){return NULL;vfree (p);return NULL;}
	    else
	      {
		verify_memory (size - i_size);
		x=realloc ( p, size+2*sizeof(Memcontrol));

		if ( x==NULL){crash ( "\nFAILED TO ALLOCATE REQUIRED MEMORY (realloc)\n");return NULL;}
		M=(Memcontrol*)x;
		M[0].size=size;
		M+=2;
		x=M;
		for ( a=i_size; a< size; a++)((char*)x)[a]=0;
		return x;
	      }
	  }
	return NULL;
	}

void vfree ( void *p)
     {
       Memcontrol *M;
       size_t size;

       if ( !p)return;
       else
	 {
	   M=(Memcontrol*)p;
	   M-=2;
	   size=M[0].size;

	   p=M;
	   free(p);

	   verify_memory (-(size+2*sizeof(Memcontrol)));
	 }
     }
void vfree_all (void *p)
{
  Memcontrol *n;
  while (memlast)
    {
      n=memlast->p;
      vfree (memlast+2);
      memlast=n;
    }
}
/*********************************************************************/
/*                                                                   */
/*                          SIZES                                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
#define WRITE_SIZE(type,function)\
void function ( int x, type *array, int os)\
     {\
     fprintf(stderr, "\nwrite_size is a deprecated function [Warning:%s]\n", PROGRAM);return;\
     }
WRITE_SIZE(short,write_size_short)
WRITE_SIZE(char,write_size_char)
WRITE_SIZE(int,write_size_int)
WRITE_SIZE(float,write_size_float)
WRITE_SIZE(double,write_size_double)

#define READ_ARRAY_SIZE(type, function)\
int function (void *array, size_t size)\
    {\
      return read_array_size (array, size);\
    }
READ_ARRAY_SIZE(short,read_size_short)
READ_ARRAY_SIZE(char,read_size_char)
READ_ARRAY_SIZE(int,read_size_int)
READ_ARRAY_SIZE(float,read_size_float)
READ_ARRAY_SIZE(double,read_size_double)


int arrlen (void *array)
    {
    Memcontrol *p;
    if (array==NULL)return 0;
    p=(Memcontrol *)array;
    p-=2;
    if ( p[0].size_element ==0)
      {
	int *x;
	fprintf ( stderr, "\nERROR in read_array_size: trying to read the size of a malloc-ed block [FATAL]");
	x[1000]=1;//THis helps triggering a trace in Valgring
	exit (0);
      }
    else return (int)p[0].size/p[0].size_element;
    }

int read_array_size_new (void *array)
{
  if (!array)return 0;
  return read_array_size ( array, 0);
}
int read_array_size (void *array, size_t size)
    {
    Memcontrol *p;
    if (array==NULL)return 0;
    p=(Memcontrol *)array;
    p-=2;
    if ( p[0].size_element ==0 && size==0)
      {
	int *x;
	fprintf ( stderr, "\nERROR in read_array_size: trying to read the size of a malloc-ed block [FATAL]");
	x[1000]=1;
	exit (0);
      }
    else if ( size ==0) return (int)p[0].size/p[0].size_element;

    return (int)p[0].size/size;

    }
int is_dynamic_memory ( void *array)
{
  Memcontrol *p;
  if (array==NULL)return 0;
  p=(Memcontrol *)array;
  p-=2;
  if ( strm (p[0].check, "dy"))return 1;
  return 0;
}

/************************************************************************/
/*                                                                      */
/*             DECLARE 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/

void * free_arrayN(void *p, int n)
{
  int a, s;
  void **i;


  if ( p==NULL) return NULL;
  else if ( n==1)vfree ((void *)p);
  else
    {
      i=(void**)p;
      s=read_array_size ( (void *)p, sizeof ( void *));
      for ( a=0; a< s; a++) free_arrayN ((void *)i[a], n-1);
      vfree (p);
    }
  return NULL;
}

void * declare_arrayNnomemset (int ndim, size_t size, ...)
{
  va_list ap;
  int *array;
  void **p;
  int a;

  va_start (ap, size);

  array=(int*)vcalloc (ndim, sizeof (int));
  for ( a=0; a< ndim; a++)
    {
      array[a]=va_arg (ap,int);
      if ( array[a]<0){va_end(ap);return NULL;}

    }
  va_end (ap);

  if ( ndim==2)
    {

      p=(void**)vcalloc_nomemset (array[0], sizeof ( void*));
      for (a=0; a< array[0]; a++)
	{
	p[a]=vcalloc_nomemset (array[1], size);
	}
    }
  else
    {
      p=(void**)declare_arrayN2nomemset (ndim, array, size);
    }
  vfree (array);
  return p;
}

void *declare_arrayN2nomemset ( int ndim, int *A, size_t size)
{
  int a;
  void **p;

  if ( ndim>1)
    {
      p=(void**)vcalloc_nomemset (A[0], sizeof (void*));
      for ( a=0; a<A[0]; a++)
    	p[a]=declare_arrayN2(ndim-1, A+1, size);
    }
  else
    {
      p=(void**)vcalloc_nomemset (A[0], size);
    }
  return p;
}

void * declare_arrayN (int ndim, size_t size, ...)
{
  va_list ap;
  int *array;
  void **p;
  int a;


  va_start (ap, size);

  array=(int*)vcalloc (ndim, sizeof (int));
  for ( a=0; a< ndim; a++)
    {
      array[a]=va_arg (ap,int);
      if ( array[a]<0){va_end(ap);return NULL;}

    }
  va_end (ap);

  if ( ndim==2)
    {

      p=(void**)vcalloc_nomemset (array[0], sizeof ( void*));
      for (a=0; a< array[0]; a++)
	p[a]=vcalloc (array[1], size);
    }
  else
    {
      p=(void**)declare_arrayN2 (ndim, array, size);
    }
  vfree (array);
  return p;
}

void *declare_arrayN2 ( int ndim, int *A, size_t size)
{
  int a;
  void **p;

  if ( ndim>1)
    {
      p=(void**)vcalloc_nomemset (A[0], sizeof (void*));
      for ( a=0; a<A[0]; a++)
    	p[a]=declare_arrayN2(ndim-1, A+1, size);
    }
  else
    {
      p=(void**)vcalloc (A[0], size);
    }
  return p;
}

void **declare_array (int first, int second, size_t size)
{
  return (void **)declare_arrayN (2,size,first, second);
}


#define DECLARE_ARRAY(type,wf,rf,function)\
type**  function (int first, int second)\
  {\
    return (type **)declare_arrayN (2,sizeof(type), first, second);\
   }

int **declare_int2 (int f, int *s, int d)
{
  int **r;
  int a;
  r=(int**)vcalloc ( f, sizeof (int*));
  for (a=0; a<f; a++)
    r[a]=(int*)vcalloc (s[a]+d, sizeof (int));
  return r;
}



char *resize_string (char *buf)
{
  char *nbuf;
  int l;

  l=strlen (buf);
  nbuf=(char*)vcalloc (l+1, sizeof (char));
  sprintf (nbuf, "%s",buf);
  vfree (buf);
  return nbuf;
}


DECLARE_ARRAY(short,write_size_short,read_size_short,declare_short)
DECLARE_ARRAY(char,write_size_char,read_size_char,declare_char)
DECLARE_ARRAY(int,write_size_int,read_size_int,declare_int)
DECLARE_ARRAY(float,write_size_float,read_size_float,declare_float)
DECLARE_ARRAY(double,write_size_double,read_size_double,declare_double)

void **declare_array_nomemset (int first, int second, size_t size)
{
  return (void **)declare_arrayNnomemset (2,size,first, second);
}
#define DECLARE_ARRAY_NMS(type,wf,rf,function)\
type**  function (int first, int second)\
  {\
    return (type **)declare_arrayNnomemset (2,sizeof(type), first, second);\
   }
DECLARE_ARRAY_NMS(short,write_size_short,read_size_short,declare_short_nomemset)
DECLARE_ARRAY_NMS(char,write_size_char,read_size_char,declare_char_nomemset)
DECLARE_ARRAY_NMS(int,write_size_int,read_size_int,declare_int_nomemset)
DECLARE_ARRAY_NMS(float,write_size_float,read_size_float,declare_float_nomemset)
DECLARE_ARRAY_NMS(double,write_size_double,read_size_double,declare_double_nomemset)


Alignment ** declare_aln_array ( int first)
    {
    Alignment ** array;
    int a;


    array=(Alignment**)vcalloc (first, sizeof (Alignment*));
    for ( a=0; a< first; a++)
      array[a]=declare_Alignment (NULL);
    return array;
    }




/************************************************************************/
/*                                                                      */
/*             Realloc 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/
void ** realloc_arrayN(int ndim,void **main_array,size_t size, ...)
{
  va_list ap;
  int *array;
  void **p;
  int a;

  /*dim size==-1: keep current size*/
  /*Dim sizes are the absolute size (not the extension*/
  /*If array getting shorter, memory is Not Claimed back*/

  array=(int*)vcalloc (ndim, sizeof (int));
  va_start (ap, size);
  for ( a=0; a< ndim; a++)
    {
      array[a]=va_arg (ap,int);
      if ( array[a]<-1){va_end(ap);return NULL;}
    }
  va_end (ap);

  p=realloc_arrayN2 (ndim, main_array, array,size);
  vfree (array);
  return p;
}

void **realloc_arrayN2 ( int ndim, void ** p, int *A, size_t size)
{
  int a;
  int o, n;


  if (ndim==0)return NULL;

  if ( ndim>1)
    {
      o=read_array_size (p,sizeof (void*));
      if (A[0]>o)p=(void**)vrealloc (p, sizeof (void*)*A[0]);
      n=(A[0]==-1)?o:A[0];
      for ( a=0; a<n; a++)
    	p[a]=(void*)realloc_arrayN2(ndim-1, (void**)p[a], A+1, size);
    }
  else
    {
      o=read_array_size (p, size);
      if (A[0]>o)p=(void**)vrealloc (p, size*A[0]);
    }
  return p;
}



void ** realloc_array (void **array,size_t size, int first, int second, int ext1, int ext2)
{
  int a;
  int d1, d2;
  if ( array==NULL)return declare_array (((first==-1)?0:first)+ext1, ((second==-1)?0:second)+ext2, size);
  else if ( first==-1)
    {
      first=read_array_size (array, sizeof (void *));
    }
  if (second==-1)second=read_array_size(array[0], size);

  d1=first+ext1;
  d2=second+ext2;

  for ( a=d1; a<first; a++)vfree (array[a]);
  array=(void**)vrealloc (array, sizeof (void*)*d1);

  if ( d2!=second)
    {
      for (a=0; a<d1 && a<first; a++)
	array[a]=vrealloc ( array[a], size*d2);
    }
  for ( a=first; a< d1; a++)
    array[a]=vrealloc ( array[a], size*d2);
  return array;
}





#define REALLOC_ARRAY(type,wf,rf,function1,function2,function3)\
type ** function1 ( type **array, int first, int second, int ext1, int ext2)\
    {\
return (type **)realloc_array ((void **)array, sizeof (type),first, second, ext1, ext2);\
\
    }
REALLOC_ARRAY(short,write_size_short,read_size_short,realloc_short,declare_short,free_short)
REALLOC_ARRAY(char,write_size_char,read_size_char,realloc_char,declare_char,free_char)
REALLOC_ARRAY(int,write_size_int,read_size_int,realloc_int,declare_int,free_int)
REALLOC_ARRAY(float,write_size_float,read_size_float,realloc_float,declare_float,free_float)
REALLOC_ARRAY(double,write_size_double,read_size_double,realloc_double,declare_double,free_double)



#define NEW_REALLOC_ARRAY(type,wf,rf,function1,function2,function3)\
type ** function1 ( type **array, int ext1, int ext2)\
{\
	int first, l1;\
	int second, l2;\
	first=rf(array,sizeof (type*));\
	second=rf(array[0],sizeof (type));\
	if (ext1==-1) ext1=first;\
	if (ext2==-1) ext2=second;\
		\
	int i,j;\
	if (ext1>first)\
	{\
		array=(type**)vrealloc(array, sizeof(type*)*ext1);\
		for (i=first; i<ext1; ++i)\
			array[i]=(type*)vcalloc(ext2,sizeof(type));\
	}\
	else if (ext1<first)\
	{\
		for (i=ext1; i<first; ++i)\
			vfree(array[i]);\
		array=(type**)vrealloc(array, sizeof(type*)*ext1);\
	}\
	if (ext2!=second)\
	{\
		if (ext1>first)\
			ext1=first;\
		for (i=0; i<ext1; ++i)\
			array[i]=(type*)vrealloc(array[i], sizeof(type)*ext2);\
	}\
	return  array;\
}




NEW_REALLOC_ARRAY(short,write_size_short,read_size_short,new_realloc_short,declare_short,free_short)
NEW_REALLOC_ARRAY(char,write_size_char,read_size_char,new_realloc_char,declare_char,free_char)
NEW_REALLOC_ARRAY(int,write_size_int,read_size_int,new_realloc_int,declare_int,free_int)
NEW_REALLOC_ARRAY(float,write_size_float,read_size_float,new_realloc_float,declare_float,free_float)
NEW_REALLOC_ARRAY(double,write_size_double,read_size_double,new_realloc_double,declare_double,free_double)
Alignment ** realloc_aln_array ( Alignment **array, int ext1)
    {
    int a;
    int first;

    if ( array==NULL)
	 {
	 array=declare_aln_array(ext1);
	 return array;
	 }
     first=read_array_size ( array, sizeof ( Alignment *));
     if ( ext1>0)
	{
	array=(Alignment**)vrealloc ( array, (sizeof (Alignment*))*(first+ext1));
        for ( a=first; a<first+ext1; a++)array[a]=declare_Alignment (NULL);
	}
    else if ( ext1==0);
    else if ( ext1<0)
         {
	 for ( a=first-1; a>=(first+ext1);a--)free_Alignment (array[a]);
	 array=(Alignment**)vrealloc ( array, (sizeof (Alignment*))*(first+ext1));
	 }
     return array;
    }

/************************************************************************/
/*                                                                      */
/*            free 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/

#define FREE_ARRAY(type,wf,rf,function) \
type ** function (type **array, int first)\
    {\
      return (type **)free_arrayN((void*)array, 2);\
    }
FREE_ARRAY(short,write_size_short,read_size_short,free_short)
FREE_ARRAY(char,write_size_char,read_size_char,free_char)
FREE_ARRAY(int,write_size_int,read_size_int,free_int)
FREE_ARRAY(float,write_size_float,read_size_float,free_float)
FREE_ARRAY(double,write_size_double,read_size_double,free_double)


Alignment ** free_aln_array (Alignment **array)
   {
   int a;
   int len;


   if ( array==NULL)return NULL;
   len=read_array_size ( array, sizeof (Alignment *));
   for ( a=1; a< len; a++)free_Alignment(array[a]);
   vfree ( array);
   return  NULL;
   }

Fname *declare_fname (int size)
   {
   Fname *F;

   size+=strlen (get_home_4_tcoffee())+FILENAMELEN+1;

   F=(Fname*)vcalloc ( 1, sizeof (Fname));
   F->name  =(char*)vcalloc ( size, sizeof (char));
   F->path  =(char*)vcalloc ( size, sizeof (char));
   F->suffix=(char*)vcalloc ( size, sizeof (char));
   F->full=(char*)vcalloc ( size, sizeof (char));
   return F;
   }

Fname *free_fname ( Fname *F)
   {
   vfree (F->name);
   vfree (F->path);
   vfree (F->suffix);
   return NULL;
   }

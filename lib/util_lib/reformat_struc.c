#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h" 
#include "define_header.h"



#define FATAL "fatal:reformat_struc"

char * process_repeat (char *aln, char *seq, char *pdb)
{
  char *tf, *file, *name;
  Alignment *A, *A2;
  Sequence *S, *P;
  int r1, r2, is_r1, is_r2, l1=0, l2=0, a;
  int *pos;
  FILE *fp;
  
  A=main_read_aln (aln, NULL);
  for (a=0; a<A->nseq; a++)ungap(A->seq_al[a]);
  
  S=main_read_seq (seq);
  P=get_pdb_sequence (pdb);
  

  A2=align_two_sequences (S->seq[0], P->seq[0], "pam250mt", -10, -1, "myers_miller_pair_wise");
  
  pos=(int*)vcalloc ( A2->len_aln+1, sizeof (int));

  for (l1=0, l2=0,a=0; a< A2->len_aln; a++)
    {
      
      r1=A2->seq_al[0][a];
      r2=A2->seq_al[1][a];
      
      is_r1=1-is_gap(r1);
      is_r2=1-is_gap(r2);
      
      l1+=is_r1;
      l2+=is_r2;
      
      if (!is_r1);
      else
	{
	  pos[l1]=l2;
	}
    }
  tf=vtmpnam(NULL);
  fp=vfopen (tf, "w");
  for (a=0; a<A->nseq; a++)
    {
      int *coor, b, c;
      
      name=A->name[a];
      file=vtmpnam (NULL);
      coor=string2num_list2 (name, "-");

      //Check the compatibility between the guide sequence and the coordinates
      for ( c=0,b=coor[1]-1; b<coor[2]; b++, c++)
	{

	  if (tolower(A->seq_al[a][c])!=tolower(S->seq[0][b]))
	    printf_exit (EXIT_FAILURE, stderr, "Incompatibility between the repeat [%s] and the master Sequence [%s]\n%s",A->name[a], seq, A->seq_al[a]);
	}
   
      printf_system ( "extract_from_pdb %s -coor %d %d -seq_field SEQRES> %s", pdb,pos[coor[1]-1],pos[coor[2]-1], file);
      fprintf (fp, ">%s _P_ %s\n", name, file);
      vfree (coor);
    }
  vfclose (fp);
  return tf;
}

void bfactor2x_in_pdb (char *infile, char *outfile, float *v)
{
  char **ll=file2lines (infile);
  int nn=1;
  int i, ci;
  FILE *out;
  out=vfopen (outfile, "w");
  i=-1;
  while (ll[nn])
    {
      int b;
      char **lw=string2list (ll[nn]);
      if (strm (lw[1], "ATOM"))
	{
	  
	  ci=atoi (lw[6]);
	  //if (i==-1)i=ci;
	  //ci-=i;
	  fprintf (out, "ATOM  ", lw[1]);//ATOM
	  fprintf (out, "%5s ",lw[2]);//ATOM NUmber
	  fprintf (out, "%-4s", lw[3]);//ATOM Name
	  fprintf (out, " "   );//alt locator -17 -
	  fprintf (out, "%3s ", lw[4]);//Res Name
	  fprintf (out, "%s", lw[5]);//Chain
	  fprintf (out, "%4s", lw[6]);//Residue Sequence Number
	  fprintf (out, " "   );//code for res insertion
	  fprintf (out, "   "   );//space before XYZ
	  fprintf (out, "%8s", lw[7]);//X
	  fprintf (out, "%8s", lw[8]);//Y
	  fprintf (out, "%8s", lw[9]);//Z
	  fprintf (out, "%6s", lw[10]);//Occupancy
	  
	  if (v[ci]<100)fprintf (out,"%6.2f\n", v[ci]);
	  else
	    fprintf (out," 99.99\n");
	}
      else
	{
	  fprintf (out, "%s\n", ll[nn]);
	}
      free_char (lw, -1);
      nn++;
    }
  vfclose (out);
  free_char (ll, -2);
}  
char *     normalize_pdb_file  (char *name_in, char *seq, char *out_file)
{
  char command[1000];
  Sequence *S;
  Alignment *A;
  int a;
  int start, end, npdb, r1, r2;
  char name[100];

  if ( !name_in) return NULL;
  else
    {
      sprintf ( name, "%s", name_in);
    }
  
 if ( !is_pdb_file(name))
    {
      fprintf(stdout, "\nERROR[normalize_pdb_file]: %s is not a pdb file[FATAL:%s]\n", name, PROGRAM);
      myexit (EXIT_FAILURE);
    }

  S=get_pdb_sequence (name);
  A=align_two_sequences (S->seq[0],seq,"idmat",-3,0, "fasta_pair_wise");
 
 for (start=-1, end=-1,npdb=0,a=0; a< A->len_aln; a++)
    {
      r1=1-is_gap(A->seq_al[0][a]);
      r2=1-is_gap(A->seq_al[1][a]);
      
      npdb+=r1;

      if (r1 && r2 && start==-1)start=npdb;
      if (r1 && r2)end=npdb;
    }

  free_aln(A);  
  free_sequence (S, -1);
  
  sprintf ( command, "extract_from_pdb -infile %s -atom ALL -chain FIRST -coor %d %d -nodiagnostic > %s", check_file_exists(name), start, end, out_file);
  my_system ( command);
  return out_file;
  }

Ca_trace * trim_ca_trace (Ca_trace *T, char *seq )
{
  /*This function:
    -removes from Ca trace all the residues that are not in the sequence
    -add in the Ca trace residues unmatched in the structure (it gives them a NULL structure)
  */
  Alignment *ALN;
  Atom *A;  
  int a,l, s, r, is_r, is_s;
  int *seq_cache, *struc_cache;
  int tot_l=0;

  char buf1[10000];
  char buf2[10000];

  /*  lower_string (T->seq);
      lower_string (seq);
  */


  sprintf (buf1, "%s", T->seq);
  sprintf (buf2, "%s", seq);
  lower_string (buf1);
  lower_string (buf2);
  

  if ( strm (buf1,buf2))return T;
  else
    {
      ALN=align_two_sequences (T->seq,seq, "est_idmat",-1, 0,"fasta_pair_wise"); 
      struc_cache=(int*)vcalloc (ALN->len_aln+1, sizeof (int));
      seq_cache  =(int*)vcalloc (ALN->len_aln+1, sizeof (int));
      
     for ( r=0, s=0,a=0; a< ALN->len_aln; a++)
       {
	 is_r=!is_gap(ALN->seq_al[0][a]);
	 is_s=!is_gap(ALN->seq_al[1][a]);
	 
	 r+=is_r;
	 s+=is_s;
	 
	 if ( is_s && is_r)
	   {
	     struc_cache[r-1]=s-1;
	     seq_cache[s-1]=r-1;
	   }
	 else if ( is_s && !is_r)
	   {
	     seq_cache[s-1]=-1;
	     
	   }
	 else if ( !is_s && is_r)
	   {
	     struc_cache[r-1]=-1;
	   }
       }
     
      T->ca=(Atom**)vrealloc ( T->ca, sizeof (Atom*)*(ALN->len_aln+1));
      T->peptide_chain=(Amino_acid**)vrealloc ( T->peptide_chain, (sizeof (Amino_acid*))*(ALN->len_aln+1));
      T->seq=(char*)vrealloc ( T->seq,   ALN->len_aln+1);

      for ( a=T->len; a< ALN->len_aln; a++)
	{
	  T->peptide_chain[a]=(Amino_acid*)vcalloc (1, sizeof (Amino_acid));
	}
      
      
      /*Read every atom*/
      for ( a=0; a< T->n_atom; a++)
	{
	  
	  A=(T->structure[a]);
	  if ( struc_cache[A->res_num-1]==-1)continue;
	  else
	    {
	      /*set the struc residue to its sequence index*/
	      A->res_num=struc_cache[A->res_num-1]+1;
	      if (strm (A->type, "CA")) {T->ca[A->res_num-1]=A;tot_l++;}
	      if ( strm (A->type, "CA"))(T->peptide_chain[A->res_num-1])->CA=A;
	      if ( strm (A->type, "C"))(T->peptide_chain[A->res_num-1] )->C=A;
	      if ( strm (A->type, "CB"))(T->peptide_chain[A->res_num-1])->CB=A;
	      if ( strm (A->type, "N"))(T->peptide_chain[A->res_num-1] )->N=A;
	    }
	}
 
      l=strlen(seq);     
      for ( a=0;a< l; a++)
       {

	 if ( seq_cache[a]==-1)
	   {
	     tot_l++;
	     T->ca[a]=NULL;
	     
	     if (!T->peptide_chain[a])T->peptide_chain[a]=(Amino_acid*)vcalloc (1, sizeof (Amino_acid));
	     T->peptide_chain[a]->CA=NULL;
	     T->peptide_chain[a]->C =NULL;
	     T->peptide_chain[a]->CB=NULL;
	     T->peptide_chain[a]->N=NULL;
	     T->seq[a]='x';
	   }
	 else
	   {
	     T->seq[a]=seq[a];
	   }
       }
      T->len=ALN->len_aln;
     
      /*
	T->len=tot_l;
      */
      
     free_aln (ALN);
     vfree(seq_cache);
     vfree(struc_cache);
     
	  
    }
  return T;
}

Ca_trace * read_ca_trace (char *name, char *seq_field )
{
  char *tp_name=NULL;
  char command[10000];
  

  if ( !is_simple_pdb_file (name))
    {
      tp_name=vtmpnam (NULL);
      sprintf ( command, "extract_from_pdb -seq_field %s  -infile %s -atom ALL -chain FIRST -mode simple> %s",seq_field, check_file_exists(name), tp_name);
      if ( getenv4debug ("DEBUG_EXTRACT_FROM_PDB"))fprintf ( stderr, "\n[DEBUG_EXTRACT_FROM_PDB:read_ca_trace] %s\n", command);
      my_system ( command);
    }
  else
    tp_name=name;
  
  return simple_read_ca_trace (tp_name);
}

Ca_trace * simple_read_ca_trace (char *tp_name )
    {
	/*This function reads a pdb file into a Ca_trace structure*/
	 
	 

	int a, c, n;
	FILE *fp;
	Atom *A, *CurrentAtom, *ListAtom;
	char res;
	char *buf;
	Ca_trace *T=NULL;
	int res_num=0, last_res_num=0;
	int cres;
	
	
	buf=(char*)vcalloc ( VERY_LONG_STRING, sizeof (char));
	n=count_n_line_in_file (tp_name );
	
	if ( !T)
	    {
	    T=(Ca_trace*)vcalloc ( 1, sizeof ( Ca_trace));
	    declare_name (T->name);
	    }
	
	/*1 Get the complete sequence: replace missing residues with Xs*/
	for (a=0; a< VERY_LONG_STRING; a++)buf[a]='x';
	res=res_num=0;
	
	fp=vfopen (tp_name, "r");
        while ( (c=fgetc(fp))!='>');
        fscanf ( fp, "%*s" );
        while ( (c=fgetc(fp))!='\n');
        fscanf ( fp, "%*s" );
        while ( (c=fgetc(fp))!='\n');	
	while ((c=fgetc(fp))!=EOF)
	  {
	       ungetc(c, fp);

	       fscanf (fp, "%*s %*s %c %*c %d %*f %*f %*f\n",&res,&res_num);
	       if (res)
		 {
		  
		   res=tolower (res);
		   buf[res_num-1]=res;
		   last_res_num=res_num;
		 }
	       res=res_num=0;
	  }
	buf[last_res_num]='\0';
	vfclose (fp);
	/*Sequence Read*/
	

	T->len=strlen (buf);
	T->seq=(char*)vcalloc ( T->len+1, sizeof (char));
	buf=lower_string (buf);
	sprintf ( T->seq, "%s", buf);
	n+=T->len;
	T->structure=(Atom**)vcalloc ( n, sizeof (Atom*));
	for ( a=0; a< n; a++)T->structure[a]=(Atom*)vcalloc ( 1, sizeof (Atom));
	T->ca=(Atom**)vcalloc ( T->len+1, sizeof ( Atom*));
	a=0;

	fp=vfopen (tp_name, "r");
        while ( (c=fgetc(fp))!='>');
        fscanf ( fp, "%*s" );
        while ( (c=fgetc(fp))!='\n');
        fscanf ( fp, "%*s" );
        while ( (c=fgetc(fp))!='\n');	
	
	cres=-1;
	while ((c=fgetc(fp))!=EOF)
	   {
	       ungetc(c, fp);
	       A=T->structure[a];
	       A->num=a;
	       fscanf (fp, "%*s %s %s %*c %d %f %f %f\n",A->type, A->res,&A->res_num, &A->x, &A->y, &A->z);
	       res=A->res[0];
	       res=tolower (res);
	       T->seq[A->res_num-1]=res;
	       A->r=atom2radius(A->type);
	       
	       if (( strm ( A->type, "CA")) || ( strm ( A->type, "C3'")))
		 {
		   T->ca[A->res_num-1]=A;
		   A->L=ListAtom;
		 }
	       if (A->res_num!=cres)
		 {
		   cres=A->res_num;
		   ListAtom=CurrentAtom=A;
		 }
	       else
		 {
		   CurrentAtom->N=A;
		   CurrentAtom=CurrentAtom->N;
		   CurrentAtom->N=NULL;
		 }
	       a++;
	   }
	T->n_atom=a;

	T->peptide_chain=(Amino_acid**)vcalloc (T->len+1, sizeof (Amino_acid*));
	for ( a=0; a<=T->len; a++) T->peptide_chain[a]=(Amino_acid*)vcalloc (1, sizeof (Amino_acid));				   
	for ( a=0; a< T->n_atom; a++)
	  {
	    A=T->structure[a];

	    if ( strm (A->type, "CA"))(T->peptide_chain[A->res_num-1])->CA=A;
	    if ( strm (A->type, "C3'"))(T->peptide_chain[A->res_num-1])->CA=A;
	    if ( strm (A->type, "C"))(T->peptide_chain[A->res_num-1] )->C=A;
	    if ( strm (A->type, "CB"))(T->peptide_chain[A->res_num-1])->CB=A;
	    if ( strm (A->type, "N"))(T->peptide_chain[A->res_num-1] )->N=A;
	  }
	    

	vfclose (fp);	
	vfree (buf);
	return T;
    }
Ca_trace * hasch_ca_trace    ( Ca_trace *T)
    {
      
      T=hasch_ca_trace_nb (T);
      T=hasch_ca_trace_bubble (T);
      T=hasch_ca_trace_transversal (T);
      return T;
    }
Ca_trace * hasch_ca_trace_transversal ( Ca_trace *TRACE)
    {
	/*This function gets the Coordinates of a protein and computes the distance of each Ca to its 
	  
	  given a Ca,
	  Compute the distance between, CA-x and CA+x with x=[1-N_ca]
	  T->nb[a][0]-->Number of distances.
	  T->nb[a][1... T->nb[a][0]]-->ngb index with respect to the Ca chain
	  T->d_nb[a][1... T->d_nb[a][0]]-->ngb index with respect to the Ca chain
	*/
	
	int a, b, d;
	float dist;
	Atom *A, *B;
		
	Struct_nb *T;
	Pdb_param *PP;


	TRACE->Transversal=(Struct_nb*)vcalloc ( 1, sizeof (Struct_nb));
	
	T=TRACE->Transversal;
	PP=TRACE->pdb_param;

	if ( !T->nb)T->nb=declare_int (TRACE->len+1, 1);
	if ( !T->d_nb)T->d_nb=declare_float (TRACE->len+1, 1);

	for (d=0,a=0; a< TRACE->len; a++)
	    {

	        for ( b=1; b<=PP->N_ca; b++)
		    {
		    if ( (a-b)<0 || (a+b)>=TRACE->len)continue;
		    A=TRACE->ca[a-b];
		    B=TRACE->ca[a+b];
		    dist=get_atomic_distance ( A, B);
		    
		    T->nb[a]=(int*)vrealloc ( T->nb[a], (++T->nb[a][0]+1)*sizeof (int));
		    T->nb[a][T->nb[a][0]]=b;

		    T->d_nb[a]=(float*)vrealloc ( T->d_nb[a], (T->nb[a][0]+1)*sizeof (float));
		    T->d_nb[a][T->nb[a][0]]=dist;
		    
		    d++;
		    }
		T->max_nb=MAX (T->max_nb, T->nb[a][0]);
	    }
	return TRACE;

    }			       

Ca_trace * hasch_ca_trace_nb ( Ca_trace *TRACE)
    {
	/*This function gets the Coordinates of a protein and computes the distance of each Ca to its 
	  T->N_ca Ngb.
	  The first Ngb to the left and to the right are excluded
	  Ngd to the left get negative distances
	  Ngb to the right receive positive distances
	  T->nb[a][0]-->Number of ngb.
	  T->nb[a][1... T->nb[a][0]]-->ngb index with respect to the Ca chain
	  T->d_nb[a][1... T->d_nb[a][0]]-->ngb index with respect to the Ca chain
	*/
	


	int a, b, d;
	float dist;
	Atom *A, *B;
		
	Struct_nb *T;
	Pdb_param *PP;
	
	TRACE->Chain=(Struct_nb*)vcalloc ( 1, sizeof (Struct_nb));
	
	T=TRACE->Chain;
	PP=TRACE->pdb_param;

	if ( !T->nb)T->nb=declare_int (TRACE->len+1, 1);
	if ( !T->d_nb)T->d_nb=declare_float (TRACE->len+1, 1);

	for (d=0,a=0; a< TRACE->len; a++)
	    {
		for ( b=MAX(0,a-PP->N_ca); b< MIN( a+PP->N_ca, TRACE->len); b++)
		        {
			if (FABS(a-b)<2)continue;
			A=TRACE->ca[a];
			B=TRACE->ca[b];
			if ( !A || !B)continue;
			dist=get_atomic_distance ( A, B);
			if (b<a) dist=-dist;
		

			T->nb[a]=(int*)vrealloc ( T->nb[a], (++T->nb[a][0]+1)*sizeof (int));
			T->nb[a][T->nb[a][0]]=b;
			
			T->d_nb[a]=(float*)vrealloc ( T->d_nb[a], (T->nb[a][0]+1)*sizeof (float));
			T->d_nb[a][T->nb[a][0]]=dist;
			d++;
			}
		    
		T->max_nb=MAX (T->max_nb, T->nb[a][0]);
	    }
	return TRACE;

    }	
Ca_trace * hasch_ca_trace_bubble ( Ca_trace *TRACE)
    {
      
	int a, b;
	float dist;
	Atom *A, *B;
	float **list;
	Pdb_param *PP;
	Struct_nb *T;

	PP=TRACE->pdb_param;
	TRACE->Bubble=(Struct_nb*)vcalloc ( 1, sizeof (Struct_nb));
	T=TRACE->Bubble;
	
		
	
	if ( !T->nb)T->nb=declare_int (TRACE->len+1, 1);
	if ( !T->d_nb)T->d_nb=declare_float (TRACE->len+1, 1);
	list=declare_float ( TRACE->n_atom, 3);
		

	for (a=0; a< TRACE->len; a++)
	    {
		for ( b=0; b< TRACE->len; b++)
		    {
			A=TRACE->ca[a];
			B=TRACE->ca[b];
			if ( !A || !B)continue;
			dist=get_atomic_distance ( A, B);
			
			if ( dist<PP->maximum_distance && FABS((A->res_num-B->res_num))>2)
			   {
			     T->nb[a][0]++;			       			       
			     T->nb[a]=(int*)vrealloc ( T->nb[a], (T->nb[a][0]+1)*sizeof (int));
			     T->nb[a][T->nb[a][0]]=(TRACE->ca[b])->num;
			     
			     T->d_nb[a]=(float*)vrealloc ( T->d_nb[a], (T->nb[a][0]+1)*sizeof (float));
			     T->d_nb[a][T->nb[a][0]]= ((a<b)?dist:-dist);
			   }						
		    }
		T->max_nb=MAX (T->max_nb, T->nb[a][0]);
	
	    }
	
	for ( a=0; a< TRACE->len; a++)
	    {
	    for ( b=0; b< T->nb[a][0]; b++)
	        {
		    list[b][0]=T->nb[a][b+1];
		    list[b][1]=T->d_nb[a][b+1];
		    list[b][2]=(TRACE->structure[T->nb[a][b+1]])->res_num;
		}
	    
	    sort_float ( list, 3,2, 0, T->nb[a][0]-1);
	    for ( b=0; b< T->nb[a][0]; b++)
	        {
		    T->nb[a][b+1]=list[b][0];
		    T->d_nb[a][b+1]=list[b][1];
		}
	    }

	free_float ( list, -1);
	return TRACE;

    }

int get_best_contact (Ca_trace *T,int enb, int *used, int *x, int *y, int *s);
float   **traces2best_contacts (Ca_trace *T1,Ca_trace*T2, float probe)
{
  int a, b, ba, bb;
  float **dist, bd,d;
  Atom *A, *B;
  //returns a non overlapping set of the best contacting pairs
  dist=declare_float (T1->len, T2->len);
  

  bd=0;
  for (a=0; a<T1->len; a++)
    for (b=0; b<T2->len; b++)
      {
	 A=T1->ca[a];
	 B=T2->ca[b];
	 d=dist[a][b]=residues2contacts(A,B, probe);
	 if (d>bd){bd=d;ba=a;bb=b;}
      }
  
  if (bd==0);
  else 
    {
      int   * used =(int*)vcalloc (T1->len, sizeof (int));
      float** bdist=declare_float (T1->len, T2->len);
      
      
      bdist[ba][bb]=bdist[bb][ba]=bd;
      used[ba]=used[bb]=1;
      
      while (bd>0)
	{
	  bd=0;
	  for (a=0; a<T1->len;a++)
	    {
	      if (used[a])continue;
	      for (b=0; b<T2->len; b++)
		{
		  if (used[b])continue;
		  d=dist[a][b];
		  if (d>bd){bd=d; ba=a; bb=b;}
		}
	    }
	  if (bd>0){used[ba]=1; used[bb]=1; bdist[ba][bb]=dist[ba][bb];}
	}
      free_float (dist,-1);
      vfree (used);
      dist=bdist;
    }
  return dist;
}
float   **trace2best_contacts (Ca_trace *T, float probe)
{
  int a, b, ba, bb;
  int enb=3;
  float **dist, bd,d;
  Atom *A, *B;
  //returns a non overlapping set of the best contacting pairs
  dist=declare_float (T->len, T->len);
  

  bd=0;
  for (a=0; a<T->len-enb; a++)
    for (b=a+enb; b<T->len; b++)
      {
	 A=T->ca[a];
	 B=T->ca[b];
	 d=dist[a][b]=dist[b][a]=residues2contacts(A,B, probe);
	 if (d>bd){bd=d;ba=a;bb=b;}
      }
  if (bd==0);
  else 
    {
      int   * used =(int*)vcalloc (T->len, sizeof (int));
      float** bdist=declare_float (T->len, T->len);
      
      
      bdist[ba][bb]=bdist[bb][ba]=bd;
      used[ba]=used[bb]=1;
      
      while (bd>0)
	{
	  bd=0;
	  for (a=0; a<T->len-enb;a++)
	    {
	      if (used[a])continue;
	      for (b=a+enb; b<T->len; b++)
		{
		  if (used[b])continue;
		  d=dist[a][b];
		  if (d>bd){bd=d; ba=a; bb=b;}
		}
	    }
	  if (bd>0){used[ba]=1; used[bb]=1; bdist[ba][bb]=bdist[bb][ba]=dist[ba][bb];}
	}
      free_float (dist,-1);
      vfree (used);
      dist=bdist;
    }
  return dist;
}
float   **traces2closest_contacts (Ca_trace *T1,Ca_trace *T2, float probe)
{
  int a, b, ba, bb;
  float **dist, bd,d;
  Atom *A, *B;
  dist=declare_float (T1->len, T2->len);
  

  bd=99;
  for (a=0; a<T1->len; a++)
    for (b=0; b<T2->len; b++)
      {
	 A=T1->ca[a];
	 B=T2->ca[b];
	 d=dist[a][b]=dist[b][a]=get_closest_vdw_contact (A,B);
	 if (d<bd && d>UNDEFINED){bd=d;ba=a;bb=b;}
      }

  
  if (bd>=99);
  else 
    {
      int   * used =(int*)vcalloc (T1->len, sizeof (int));
      float** bdist=declare_float (T1->len, T2->len);
      
      bdist[ba][bb]=bdist[bb][ba]=bd;
      used[ba]=used[bb]=1;
      
      while (bd!=99)
	{
	  bd=99;
	  for (a=0; a<T1->len;a++)
	    {
	      if (used[a])continue;
	      for (b=0; b<T2->len; b++)
		{
		  if (used[b])continue;
		  d=dist[a][b];
		  if (d<bd && d>UNDEFINED){bd=d; ba=a; bb=b;}
		}
	    }
	  if (bd!=99)
	    {
	      used[ba]=1; used[bb]=1; bdist[ba][bb]=(dist[ba][bb]<0)?0.01:dist[ba][bb];
	    }
	}
      free_float (dist,-1);
      vfree (used);
      dist=bdist;
    }
  return dist;
}
float   **trace2closest_contacts (Ca_trace *T, float probe)
{
  int a, b, ba, bb;
  int enb=3;
  float **dist, bd,d;
  Atom *A, *B;
  dist=declare_float (T->len, T->len);
  

  bd=99;
  for (a=0; a<T->len-enb; a++)
    for (b=a+enb; b<T->len; b++)
      {
	 A=T->ca[a];
	 B=T->ca[b];
	 d=dist[a][b]=dist[b][a]=get_closest_vdw_contact (A,B);
	 if (d<bd && d>UNDEFINED){bd=d;ba=a;bb=b;}
      }

  
  if (bd>=99);
  else 
    {
      int   * used =(int*)vcalloc (T->len, sizeof (int));
      float** bdist=declare_float (T->len, T->len);
      
      bdist[ba][bb]=bdist[bb][ba]=bd;
      used[ba]=used[bb]=1;
      
      while (bd!=99)
	{
	  bd=99;
	  for (a=0; a<T->len-enb;a++)
	    {
	      if (used[a])continue;
	      for (b=a+enb; b<T->len; b++)
		{
		  if (used[b])continue;
		  d=dist[a][b];
		  if (d<bd && d>UNDEFINED){bd=d; ba=a; bb=b;}
		}
	    }
	  if (bd!=99)
	    {
	      used[ba]=1; used[bb]=1; bdist[ba][bb]=bdist[bb][ba]=(dist[ba][bb]<0)?0.01:dist[ba][bb];
	    }
	}
      free_float (dist,-1);
      vfree (used);
      dist=bdist;
    }
  return dist;
}
float   **traces2count_contacts (Ca_trace *T1, Ca_trace*T2, float probe)
{
  int a, b;
  Atom *A, *B;
  float **dist,d;

  
  dist=declare_float ( T1->len, T2->len);
  for (a=0; a< T1->len; a++)
    {
      for ( b=0; b< T2->len; b++)
	{
	  
	  A=T1->ca[a];
	  B=T2->ca[b];
	  dist[a][b]=residues2contacts (A,B,probe);
	}
    }
  return dist;
}

float   **trace2count_contacts (Ca_trace *T, float probe)
{
  int a, b;
  Atom *A, *B;
  float **dist,d;
  int enb=3;
  
  dist=declare_float ( T->len, T->len);
  for (a=0; a< T->len-enb; a++)
    {
      for ( b=a+enb; b< T->len; b++)
	{
	  
	  A=T->ca[a];
	  B=T->ca[b];
	  dist[a][b]=dist[b][a]=residues2contacts (A,B,probe);
	}
    }
  return dist;
}
float   **traces2contacts (Ca_trace *T1,Ca_trace *T2, float probe)
{
  //identifies contacts betwee AA through ANY atom
  int a, b;
  Atom *A, *B;
  float **dist,d;
  
  float MaxD=50;
  //Hard coded maximum that corrsponds to the largest possible distance between the Ca of two contacting AA
  //just a trck to make things run faster, should not impact results

  dist=declare_float ( T1->len, T2->len);
  for (a=0; a< T1->len; a++)
    {
      for ( b=0; b< T2->len; b++)
	{
	  
	  A=T1->ca[a];
	  B=T2->ca[b];
	  d=get_atomic_distance (A,B);
	  
	  if (d>=UNDEFINED && d<=MaxD)
	    {
	      while (A)
		{
		  Atom *CB=B;
		  while (CB)
		    {
		      float max=A->r+CB->r+probe;
		      float cd;
		      
		      d=get_atomic_distance (A,CB);
		      cd=d-(CB->r+A->r);
		      if (cd<0)cd=0.01;
		      if (d<max && d>=UNDEFINED){dist[a][b]=cd; A=CB=NULL;}
		      else CB=CB->N;
		    }
		  A=(A)?A->N:NULL;
		}
	    }
	}
    }
  return dist;
}
float   **trace2contacts (Ca_trace *T, float probe)
{
  //identifies contacts betwee AA through ANY atom
  int a, b;
  Atom *A, *B;
  float **dist,d;
  int enb=3;
  float MaxD=50;
  //Hard coded maximum that corrsponds to the largest possible distance between the Ca of two contacting AA
  //just a trck to make things run faster, should not impact results

  dist=declare_float ( T->len, T->len);
  for (a=0; a< T->len-enb; a++)
    {
      for ( b=a+enb; b< T->len; b++)
	{
	  
	  A=T->ca[a];
	  B=T->ca[b];
	  d=get_atomic_distance (A,B);
	  
	  if (d>=UNDEFINED && d<=MaxD)
	    {
	      while (A)
		{
		  Atom *CB=B;
		  while (CB)
		    {
		      float max=A->r+CB->r+probe;
		      float cd;
		      
		      d=get_atomic_distance (A,CB);
		      cd=d-(CB->r+A->r);
		      if (cd<0)cd=0.01;
		      if (d<max && d>=UNDEFINED){dist[a][b]=dist[b][a]=cd; A=CB=NULL;}
		      else CB=CB->N;
		    }
		  A=(A)?A->N:NULL;
		}
	    }
	}
    }
  return dist;
}

float ** traces2ca_distances(Ca_trace *T1,Ca_trace *T2, float max)
   {
     //Measures Ca distances
     int a, b;
     Atom *A, *B;
     float **dist;
     
     dist=declare_float ( T1->len, T2->len);
     
     for (a=0; a< T1->len; a++)
       {
	 for ( b=0; b< T2->len; b++)
	   {
	     float d;
	     dist[a][b]=UNDEFINED;
	     A=T1->ca[a];
	     B=T2->ca[b];
	     d=get_atomic_distance ( A,      B);
	     if (d<max || max<=0.000001)dist[a][b]=d;
	   }
       }
     return dist;
   }  

float ** trace2ca_distances(Ca_trace *T, float max)
   {
     //Measures Ca distances
       int a, b;
       Atom *A, *B;
       float **dist;
       
       dist=declare_float ( T->len, T->len);
    
       for (a=0; a< T->len-1; a++)
	    {
		for ( b=a+1; b< T->len; b++)
		    {
		      float d;
		      dist[a][b]=dist[b][a]=UNDEFINED;
		      A=T->ca[a];
		      B=T->ca[b];
		      d=get_atomic_distance ( A,      B);
		      if (d<max || max<=0.000001)dist[a][b]=dist[b][a]=d;
		    }
	    }
       return dist;
   }

int residues2contacts (Atom *iA, Atom*iB, float probe)
{
  Atom *A, *B;
  int x, y, nx, ny, bx,by, nc;
  int *luA, *luB;
  float **dm, d, bd;
  
  
  
  d=get_atomic_distance (iA,iB);
  if ( d<=UNDEFINED || d>30)return 0;

  A=iA;
  B=iB;

  nx=0;ny=0;
  while (A){A=A->N; nx++;}
  while (B){B=B->N, ny++;}
  
  dm=declare_float (nx, ny);
  luA=(int*)vcalloc (nx, sizeof (int));
  luB=(int*)vcalloc (ny, sizeof (int));
  
  x=y=0;
  A=iA;
  bd=99;
  while (A)
    {
      B=iB;
      y=0;
      while (B)
	{
	  d=dm[x][y]=get_atomic_distance (A,B)-(A->r+B->r);
	  B=B->N;
	  y++;
	  if (d<probe && d>UNDEFINED && d<bd){bd=d; bx=x; by=y;}
	}
      A=A->N;
      x++;
    }
  if (bd>=99)nc=0;
  else
    {
      nc=1;luA[bx]=1;luB[by]=1;
      while (bd<99)
	{
	  bd=99;
	  for (x=0; x<nx; x++)
	    {
	      if (luA[x]==1)continue;
	      for (y=0; y<ny; y++)
		{
		  if (luB[y]==1) continue;
		  d=dm[x][y];
		  if (d<probe && d>UNDEFINED && d<bd){bd=d; bx=x;by=y;}
		}
	    }
	  if (bd<99){nc++; luA[bx]=1; luB[by]=1;}
	}
    }
  //free_float (dm, -1);
  //vfree (luA); vfree (luB);
  return nc;
} 
float    get_closest_vdw_contact   ( Atom *iA, Atom*iB)
{
  float d, best_d;
  Atom *A, *B;
  if (!iA || !iB) return UNDEFINED;
  
  best_d=get_atomic_distance (iA, iB)-(iA->r+iB->r);
  if (best_d>30) return best_d;
  
  A=iA;
  best_d=99;
  while (A)
    {
      B=iB;
      while (B)
	{
	  d=get_atomic_distance (A, B);
	  d-=(A->r+B->r);
	  if (d<best_d && d>=UNDEFINED){best_d=d;}
	  B=B->N;
	}
      A=A->N;
    }
  if (best_d==99)return UNDEFINED;
  else return best_d;
}
  
float    get_atomic_distance ( Atom *A, Atom*B)
   {
       float dx, dy, dz, d;
       
       if ( !A || !B)
	 {
	   return UNDEFINED;
	 }

       
       dx=A->x - B->x;
       dy=A->y - B->y;
       dz=A->z - B->z;
       

       d=(float) sqrt ( (double) ( dx*dx +dy*dy +dz*dz));
       return d;
   }

float atom2radius (char *t)
{
  //radius taken from rasmol: http://www.rasmol.org/software/RasMol_2.7.3/src/abstree.h 
  //converted from rasmol units as detailed in https://www.umass.edu/microbio/rasmol/rasbonds.htm
  //07/12/16

  char a, b;
  static float **lu;
  if (!t)return 0;
  else if ( strm (t, "O1P"))return atom2radius ("O"); //NMR
  else if ( strm (t, "O2P"))return atom2radius ("O"); //NMR
  else if ( strm (t, "O3P"))return atom2radius ("O"); //NMR
  else if ( strm (t, "O4P"))return atom2radius ("O"); //NMR
  else if ( strm (t, "O5P"))return atom2radius ("O"); //NMR
  else if ( strm (t, "AS"))return atom2radius ("As"); 

  a=t[0]; b=t[1];
  if (!lu)
    {
      int x, y;
      lu=(float **)declare_float (257, 257);
      for (x=0; x<256; x++)
	for (y=0; y<256; y++)
	  lu[x][y]=-1;
     
      lu['A']['c']=2.12;
      lu['A']['g']=1.548;
      lu['A']['l']=1.5;
      lu['A']['m']=1.66;
      lu['A']['r']=2.768;
      lu['A']['s']=0.828;
      lu['A']['t']=1.12;
      lu['A']['u']=1.448;
      lu['B'][' ']=1.548;
      lu['B']['a']=2.408;
      lu['B']['e']=0.628;
      lu['B']['i']=1.728;
      lu['B']['k']=1.64;
      lu['B']['r']=1.748;
      lu['C'][' ']=1.548;
      lu['C']['a']=1.948;
      lu['C']['d']=1.748;
      lu['C']['e']=1.86;
      lu['C']['f']=1.628;
      lu['C']['l']=1.748;
      lu['C']['m']=1.648;
      lu['C']['o']=1.128;
      lu['C']['r']=1.128;
      lu['C']['s']=3.008;
      lu['C']['u']=1.148;
      lu['D']['y']=1.628;
      lu['E']['r']=1.588;
      lu['E']['s']=1.62;
      lu['E']['u']=1.96;
      lu['F'][' ']=1.3;
      lu['F']['e']=1.948;
      lu['F']['m']=1.608;
      lu['F']['r']=3.24;
      lu['G']['a']=1.548;
      lu['G']['d']=1.688;
      lu['G']['e']=3.996;
      lu['H'][' ']=1.1;
      lu['H']['e']=2.2;
      lu['H']['f']=1.4;
      lu['H']['g']=1.98;
      lu['H']['o']=1.608;
      lu['I'][' ']=1.748;
      lu['I']['n']=1.448;
      lu['I']['r']=1.22;
      lu['K'][' ']=2.388;
      lu['K']['r']=1.9;
      lu['L']['a']=1.828;
      lu['L']['i']=1.22;
      lu['L']['r']=1.58;
      lu['L']['u']=1.528;
      lu['M']['d']=1.6;
      lu['M']['g']=1.5;
      lu['M']['n']=1.188;
      lu['M']['o']=1.748;
      lu['N'][' ']=1.4;
      lu['N']['a']=2.2;
      lu['N']['b']=1.328;
      lu['N']['d']=1.788;
      lu['N']['e']=2.02;
      lu['N']['i']=1.24;
      lu['N']['o']=1.588;
      lu['N']['p']=1.708;
      lu['O'][' ']=1.348;
      

      lu['O']['s']=1.58;
      lu['P'][' ']=1.88;
      lu['P']['a']=1.6;
      lu['P']['b']=2.16;
      lu['P']['d']=1.44;
      lu['P']['m']=1.76;
      lu['P']['o']=1.208;
      lu['P']['r']=1.62;
      lu['P']['t']=1.548;
      lu['P']['u']=1.668;
      lu['R']['a']=2.568;
      lu['R']['b']=2.648;
      lu['R']['e']=1.3;
      lu['R']['h']=1.22;
      lu['R']['n']=2.3;
      lu['R']['u']=1.2;
      lu['S'][' ']=1.808;
      lu['S']['b']=1.12;
      lu['S']['c']=1.32;
      lu['S']['d']=1.808;
      lu['S']['e']=0.9;
      lu['S']['i']=2.2;
      lu['S']['m']=1.74;
      lu['S']['n']=1.668;
      lu['S']['r']=2.02;
      lu['T']['a']=1.22;
      lu['T']['b']=1.66;
      lu['T']['c']=1.8;
      lu['T']['e']=1.26;
      lu['T']['h']=1.84;
      lu['T']['i']=1.948;
      lu['T']['l']=1.708;
      lu['T']['m']=1.568;
      lu['U'][' ']=1.748;
      lu['V'][' ']=1.06;
      lu['W'][' ']=1.26;
      lu['X']['e']=2.1;
      lu['Y'][' ']=1.608;
      lu['Y']['b']=1.54;
      lu['Z']['n']=1.148;
      lu['Z']['r']=1.42;
      
    }
  if (!b)b=' ';
  else if (b>='A' && b<='Z')b=' ';
  else if (b>='0' && b<='9')b=' ';
  else if (b=='\'')b=' ';

  if (lu[a][b]<0)
    {
      if (strlen (t)>=2 && lu[a][tolower(t[1])]>0)return lu[a][tolower(t[1])];
      else
	add_warning ( stderr, "\nWARNING: Atom [%s] is unknown Defaulf Van Der Waals Value of 1.5 Angstrom will be used",t);
      return 1;
    }
  else 
    return lu[a][b];

}
    
float atom2radius_old (char *t)
{
  char a, b;
  
  if (!t)return 0;
  a=t[0]; b=t[1];
  //radius taken from rasmol: https://www.umass.edu/microbio/rasmol/rasbonds.htm
  if      (a=='C' && b=='a')return 1.940;
  else if (a=='C' && b=='d')return 1.748;
  else if (a=='C')return 1.548;
  else if (a=='N')return 1.400;
  else if (a=='O')return 1.348;
  else if (a=='P')return 1.880;
  else if (a=='S')return 1.808;
  else if (a=='F')return 1.948;
  else if (a=='Z')return 1.148;
  else if (a=='I')return 1.748;
  else if (a=='H')return 1.100;
  else return 2;
  
}
  
char * map_contacts ( char  *file1, char *file2, float T)
{
  
  Ca_trace *ST1, *ST2;
  int *contact_list;
  int a;
  

  ST1=read_ca_trace (file1, "ATOM");
  ST2=read_ca_trace (file2, "ATOM");

  contact_list=identify_contacts (ST1, ST2, T);
  for ( a=0; a<ST1->len; a++)
    {
      ST1->seq[a]=(contact_list[a]==1)?toupper(ST1->seq[a]):tolower(ST1->seq[a]);
    }

  return ST1->seq;
}

float ** print_contacts ( char  *file1, char *file2, float T)
{
  
  Ca_trace *ST1, *ST2;
  float **dist, d;
  int a, b;
  Atom *A, *B;
  char *list;
  int *list1=NULL, *list2=NULL;
  int *cache1, *cache2;
  
  
   
  if ((list=strstr (file1, "_R_"))!=NULL)
      {
	list[0]='\0';
	list+=3;
	list1=string2num_list2(list, "_");
      }
  
  if ((list=strstr (file2, "_R_"))!=NULL)
      {
	list[0]='\0';
	list+=3;
	list2=string2num_list2(list, "_");
      }	

  fprintf ( stdout, "\n#%s (%s) struc2contacts_01", PROGRAM, VERSION);
  fprintf ( stdout, "\nStructure %s vs %s", file1, file2);
  ST1=read_ca_trace (file1, "SEQRES");
  ST2=read_ca_trace (file2, "SEQRES");
  
  
  cache1=(int*)vcalloc (ST1->len+1, sizeof (int));
  cache2=(int*)vcalloc (ST2->len+1, sizeof (int));
  
  if (list1)for ( a=1; a<list1[0]; a++)cache1[list1[a]]=1;
  else for ( a=1; a<ST1->len; a++)cache1[a]=1;
  
  if (list2)for ( a=1; a<list2[0]; a++)cache2[list2[a]]=1;
  else for ( a=1; a<ST2->len; a++)cache2[a]=1;
  
  
  dist=declare_float (ST1->len+1,ST2->len+1);
  vfree (list1); vfree(list2);
  
  for ( a=0; a< ST1->n_atom; a++)
    {
      A=ST1->structure[a];
      if ( !cache1[A->res_num])continue;
      for ( b=0; b<ST2->n_atom; b++)
	    {
	      
	      B=ST2->structure[b];
	      if( !cache2[B->res_num])continue;
	    
	      d=get_atomic_distance (A,B);
	      
	      if (dist[A->res_num][B->res_num]==0 || dist[A->res_num][B->res_num]>d)dist[A->res_num][B->res_num]=d;
	    }
    }

  for ( a=1; a<=ST1->len; a++)
    {
      A=ST1->ca[a-1];
      if ( !A || !cache1[A->res_num])continue;
      for ( b=1; b<= ST2->len; b++)
	{
	  B=ST2->ca[b-1];
	  if( !B || !cache2[B->res_num])continue;
	  if(dist[a][b]!=0)fprintf ( stdout, "\nResidue %3d [%s] vs %3d [%s] %9.4f Angstrom",A->res_num,A->res,B->res_num,B->res,dist[a][b]);
	}
    }
  fprintf ( stdout, "\n");
  vfree (cache1);vfree(cache2);
  free_float (dist, -1);
  return NULL;
}
  

int * identify_contacts (Ca_trace *ST1,Ca_trace *ST2, float T) 
{
  int a, b;
  float d;
  int *result;
  
  
    
  result=(int*)vcalloc ( ST1->len+1, sizeof (int));

  
  for ( a=0; a< ST1->n_atom; a++)
    for ( b=0; b<ST2->n_atom; b++)
  
      {
	
	d=get_atomic_distance (ST1->structure[a], ST2->structure[b]);
	if (d<T)
	  {
	    /*Atom *A, *B;
	    A=ST1->structure[a]; B=ST2->structure[b];
	    fprintf ( stderr, "\n%d %s %s Vs %d %s %s: %f", A->res_num, A->res, A->type, B->res_num, B->res, B->type, d); */
	    result[(ST1->structure[a])->res_num-1]=1;
	  }
      }
  return result;
}

Sequence *seq2contacts4ligands ( Sequence *S, float T)
{
  int a;
  Sequence *NS;
  
 
  
  NS=duplicate_sequence (S);
  for ( a=0; a< S->nseq; a++)
    {
      NS->seq[a]=string2contacts ( S->seq[a], S->name[a], S->seq_comment[a], T);
    }
  
  return NS;
}

char *string2contacts (char *seq,char *name, char *comment, float T)
{
  char **nlist;
  char *r;
  char *result;
  int a, b, n;
  char *struc_name;
  static char *struc_file;
  static char *ligand_file;
  Alignment *A;
  char r0, r1;
  
  int l, ln;
  char command[1000];
  /*>seq__struc Ligand1 Chain1 Ligand 2 cahin2
    Chain: index or ANY if unknown
    Ligand: name of pdb file
  */
  
  if ( !struc_file)
    {
      struc_file=vtmpnam  (NULL);
      ligand_file=vtmpnam (NULL);
    }
  

  
  result=(char*)vcalloc ( strlen (seq)+1, sizeof (char));
  for ( a=0; a< strlen (seq); a++)result[a]='0';
  
  nlist=string2list (comment);
  if ( !nlist)return result;
  else
    n=atoi(nlist[0]);

  struc_name=strstr(name, "_S_");
  if (!struc_name && !is_pdb_struc (name))
    {
      
      return result;
    }
  else if ( !struc_name && is_pdb_struc (name))
    {
      struc_name=name;
    }
  else
    {
      struc_name+=3;
      if ( check_file_exists (struc_name) && is_simple_pdb_file(struc_name))
	{
	 sprintf (command, "cp %s %s", name, struc_file);
	}
      else
	{
	  sprintf ( command, "extract_from_pdb -infile %s -atom ALL -mode simple -force >%s",name, struc_file);
	}
      my_system (command);
    }
  
  
  
  for ( a=1, ln=1;a<n; ln++)
	{
	  fprintf ( stderr, "\n\tProcess: Structure %s Against Ligand %s T=%.2f Angstrom",struc_name, nlist[a], T);
	 
	  
	  if (check_file_exists ( nlist[a]) && is_simple_pdb_file(nlist[a]))
	      {		
		sprintf (command, "cp %s %s", nlist[a], ligand_file);
		a++;
	      }
	  else if ( check_file_exists ( nlist[a]))
	    {
	      sprintf (command, "extract_from_pdb -infile %s -ligand ALL -atom ALL -mode simple -force >%s", nlist[a], ligand_file);
	      a++;
	    }
	  else
	    {
	      sprintf ( command, "extract_from_pdb -infile %s -chain %s -ligand %s -ligand_only -atom ALL -mode simple -force >%s",struc_name, nlist[a+1],nlist[a], ligand_file);
	      a+=2;	      
	    }
	  my_system (command);
	  
	  if ( T>0)
	    {
	      r=map_contacts (struc_file,ligand_file,T);
	      
	      toggle_case_in_align_two_sequences (KEEP_CASE);
	      A=align_two_sequences (seq,r,"pam250mt", -10, -1, "myers_miller_pair_wise");     
	      toggle_case_in_align_two_sequences (CHANGE_CASE);
	      
	      
	      for ( l=0,b=0; b< A->len_aln; b++)
		{
		  r0=A->seq_al[0][b];r1=A->seq_al[1][b];
		  if (!is_gap(r0))
		    {
		      if (isupper(r1))result[l]=(result[l]!='0')?'9':'0'+ln;
		      l++;
		    }
		}
	      
	      free_aln (A);
	      fprintf ( stderr, " [DONE]");
	    }
	  else if ( T==0)
	    {
	      print_contacts( struc_file,ligand_file,T);
	    }
	  
	}
      fprintf ( stderr, "\n");
      
      return result;
}


char **struclist2nb (char *name,char *seq, char *comment, float Threshold, char *atom, char *output)
{
  char *list, **pdb;
  int a;
  char **R, *tmpf;
  
  tmpf=vtmpnam (NULL);
  
  list=strstr ( comment, "_P_")+3;
  if ( !strstr (comment, "_P_"))return NULL;
  else
    {
      pdb=string2list ( strstr ( comment, "_P_")+3);
    }
  
  for (a=1; a<atoi(pdb[0]); a++)
    {
      char com[1000];
      sprintf ( com, "_P_ %s", pdb[a]);
      R=struc2nb (name, seq, com, Threshold, atom, tmpf);
    }
  if (output==NULL|| strstr (output, "fasta"))
    {
      for (a=0; a< strlen (seq); a++)
	if (isupper (R[a][a]))fprintf (stdout, ">%s R=%d T=%.2f %s\n%s\n", name, a+1, Threshold, comment, R[a]);
    }
  else
    {
      FILE *fp;
      char c;
      
      fp=vfopen (tmpf, "r");
      while ( (c=fgetc(fp))!=EOF)fprintf (stdout, "%c", c);
      vfclose (fp);
    }
  return NULL;
}

char **struc2nb (char *name,char *seq, char *comment, float Threshold, char *atom, char *output)
{
  char *struc_file;
  char *struc_name;
  Ca_trace *T;
  Atom *A1, *A2;
  int a, b;
  short **hasch;
  FILE *fp;
  float d;
  char command[10000];
  static char **R;
  
  struc_file=vtmpnam (NULL);
  declare_name (struc_name);
  
  sscanf ( strstr(comment, "_P_"), "_P_ %s", struc_name);
  //struc_name=strstr(name, "_S_");
  
  if (!R)
    {
      int l;
      l=strlen (seq);
      R=declare_char (l+1, l+1);
      for ( a=0; a<l; a++)
	{
	  sprintf ( R[a], "%s", seq);
	  lower_string (R[a]);
	}
    }
  if ( atom) atom=translate_string (atom, "_", " ");
  if ( !struc_name) return NULL;
  
  else 
    {
      //struc_name+=3;
      if ( check_file_exists (struc_name) && is_simple_pdb_file(struc_name))
	{
	  sprintf (command, "cp %s %s", name, struc_file);
	}
      else
	{
	  sprintf ( command, "extract_from_pdb -infile %s -atom %s -mode simple -force >%s",struc_name,(atom==NULL)?"ALL":atom, struc_file);
	}
  
      my_system (command);
    }
  T=read_ca_trace (struc_file, "ATOM");
  hasch=declare_short (T->len, T->len);
  
  if (!R)
    {
      int l;
      l=strlen (seq);
      R=declare_char (l+1, l+1);
      for ( a=0; a<l; a++)
	{
	  sprintf ( R[a], "%s", seq);
	  lower_string (R[a]);
	}
    }
  for ( a=0; a< T->n_atom; a++)
    for ( b=0; b< T->n_atom; b++)
      {
	A1=T->structure[a];A2=T->structure[b];
	d=get_atomic_distance (A1, A2);

	if ( d<Threshold)
	  {
	    hasch[A1->res_num-1][A2->res_num-1]=1;
	  }
      }
  fp=vfopen (output, "a");
  fprintf ( fp, "#Distance_map_format_01\n#Sequence %s with T= %.2f", struc_name, Threshold);
  for ( a=0; a<T->len; a++)
    {
      int c;
      c=change_residue_coordinate ( T->seq,seq,a);
      
      if ( c!=-1 && seq)
	{
	  char r1, r2;
	  r1=(T->seq)[a];
	  r2=seq[c];
	  r1=tolower(r1);r2=tolower(r2);
	  if ( r1!=r2) continue;
	  R[c][c]=toupper (R[c][c]);
	  fprintf (fp, "\n%s Residue %d ",struc_name,c+1);
	  for ( b=0; b<T->len; b++)
	    {
	      int d;
	      char r3, r4;
	      r3=(T->seq)[a];r4=seq[c];
	      r3=tolower(r3);r4=tolower(r4);
	      if ( r3!=r4) continue;
	      
	      d=change_residue_coordinate (T->seq,seq,b);
	      
	      if ( hasch[a][b] && d!=-1)
		{
		  fprintf (fp, "%d ",d+1);
		  R[c][d]=toupper (R[c][d]);
		}
	    }
	  fprintf (fp, ";");
	}
    }
  free_short (hasch, -1);
  fprintf (fp, "\n");
  return R;
}

short **seq2nb (char *seq, char *pdb, float Threshold, char *atom)
{
  char *struc_file;
 
  Ca_trace *T;
  Atom *A1, *A2;
  int a, b;
  short **hasch;
  short **result; //Contains the result for every residue of seq
  float d;
  char command[10000];

  // get a clean pdb file
  struc_file=vtmpnam (NULL);
  if ( check_file_exists (struc_file) && is_simple_pdb_file(struc_file))
    {
      sprintf (command, "cp %s %s", pdb, struc_file);
    }
  else
    {
      sprintf ( command, "extract_from_pdb -infile %s -atom %s -mode simple -force >%s",pdb,(atom==NULL)?"ALL":atom, struc_file);
    }
  my_system (command);
  
  //Read and hasch the PDB file
  T=read_ca_trace (struc_file, "ATOM");
  hasch=declare_short (T->len, T->len);
  result=declare_short (strlen (seq)+1, strlen (seq)+1);
  for ( a=0; a< T->n_atom; a++)
    for ( b=0; b< T->n_atom; b++)
      {
	A1=T->structure[a];A2=T->structure[b];
	d=get_atomic_distance (A1, A2);

	if ( d<Threshold)
	  {
	    hasch[A1->res_num-1][A2->res_num-1]=1;
	  }
      }
  
  for ( a=0; a<T->len; a++)
    {
      int c;
      c=change_residue_coordinate ( T->seq,seq,a);
      if ( c!=-1)
	{
	  for ( b=0; b<T->len; b++)
	    {
	      int d;
	      d=change_residue_coordinate (T->seq,seq,b);
	      
	      if ( hasch[a][b] && d!=-1)
		{
		  result[c][result[c][0]++]=d;
		}
	    }
	}
    }
  free_short (hasch, -1);
  return result;
}

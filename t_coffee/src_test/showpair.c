#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

static void make_p_ptrs(int *tptr, int *pl, int naseq, int l);
static void make_n_ptrs(int *tptr, int *pl, int naseq, int len);
static void put_frag(int fs, int v1, int v2, int flen);
static int frag_rel_pos(int a1, int b1, int a2, int b2);
static void des_quick_sort(int *array1, int *array2, int array_size);
static void pair_align(int seq_no, int l1, int l2);


/*
*	Prototypes
*/

/*
*	 Global variables
*/
/*extern int *seqlen_array;
  extern char **seq_array;
  extern int  dna_ktup, dna_window, dna_wind_gap, dna_signif; params for DNA 
 extern int prot_ktup,prot_window,prot_wind_gap,prot_signif; params for prots 
  extern int 	nseqs;
  extern Boolean 	dnaflag;
  extern double 	**tmat;
  extern int 	max_aa;
  extern int  max_aln_length;
*/

static int *seqlen_array;
static char **seq_array;
 
static int 	nseqs;
static int 	dnaflag;
static int      max_aln_length;
static int      max_aa;

static int 	next;
static int 	curr_frag,maxsf;
static int 	**accum;
static int 	*diag_index;
static char 	*slopes;

int ktup,window,wind_gap,signif;    		      /* Pairwise aln. params */
int *displ;
int *zza, *zzb, *zzc, *zzd;

static Boolean percent=1;


static void make_p_ptrs(int *tptr,int *pl,int naseq,int l)
{
	static int a[10];
	int i,j,limit,code,flag;
	int residue;
	
	/*tptr--> pointer to the last occurence of the same residue or ktuple:

	  abcdeabef
	  
	  tptr: 0 0 0 0 0 1 2 5 0
	  pl[a]=6
	  pl[b]=7
	*/

	
	for (i=1;i<=ktup;i++)
           a[i] = (int) pow((double)(max_aa+1),(double)(i-1));

	limit = (int)   pow((double)(max_aa+1),(double)ktup);
	for(i=1;i<=limit;++i)
		pl[i]=0;
	for(i=1;i<=l;++i)
		tptr[i]=0;
	
	

	for(i=1;i<=(l-ktup+1);++i) {
		code=0;
		flag=FALSE;
		for(j=1;j<=ktup;++j) {
			residue = seq_array[naseq][i+j-1];
			if((residue<0) || (residue > max_aa)){
				flag=TRUE;
				break;
			}
			code += ((residue) * a[j]);
		}
		if(flag)
			continue;
		++code;
		if(pl[code]!=0)tptr[i]=pl[code];
		pl[code]=i;
	}
}


static void make_n_ptrs(int *tptr,int *pl,int naseq,int len)
{
	static int pot[]={ 0, 1, 4, 16, 64, 256, 1024, 4096 };
	int i,j,limit,code,flag;
	int residue;
	
	limit = (int) pow((double)4,(double)ktup);
	
	for(i=1;i<=limit;++i)
		pl[i]=0;
	for(i=1;i<=len;++i)
		tptr[i]=0;
	
	for(i=1;i<=len-ktup+1;++i) {
		code=0;
		flag=FALSE;
		for(j=1;j<=ktup;++j) {
			residue = seq_array[naseq][i+j-1];
			if((residue<0) || (residue>4)){
				flag=TRUE;
				break;
			}
			code += ((residue) * pot[j]);  /* DES */
		}
		if(flag)
			continue;
		++code;
		if(pl[code]!=0)
			tptr[i]=pl[code];
		pl[code]=i;
	}
}


static void put_frag(int fs,int v1,int v2,int flen)
{
	int end;
	accum[0][curr_frag]=fs;
	accum[1][curr_frag]=v1;
	accum[2][curr_frag]=v2;
	accum[3][curr_frag]=flen;
	
	if(!maxsf) {
		maxsf=1;
		accum[4][curr_frag]=0;
		return;
	}
	
        if(fs >= accum[0][maxsf]) {
		accum[4][curr_frag]=maxsf;
		maxsf=curr_frag;
		return;
	}
	else {
		next=maxsf;
		while(TRUE) {
			end=next;
			next=accum[4][next];
			if(fs>=accum[0][next])
				break;
		}
		accum[4][curr_frag]=next;
		accum[4][end]=curr_frag;
	}
}


static int frag_rel_pos(int a1,int b1,int a2,int b2)
{
	int ret;
	
	ret=FALSE;
	if(a1-b1==a2-b2) {
		if(a2<a1)
			ret=TRUE;
	}
	else {
		if(a2+ktup-1<a1 && b2+ktup-1<b1)
			ret=TRUE;
	}
	return ret;
}


static void des_quick_sort(int *array1, int *array2, int array_size)
/*  */
/* Quicksort routine, adapted from chapter 4, page 115 of software tools */
/* by Kernighan and Plauger, (1986) */
/* Sort the elements of array1 and sort the */
/* elements of array2 accordingly */
/*  */
{
	int temp1, temp2;
	int p, pivlin;
	int i, j;
	int lst[50], ust[50];       /* the maximum no. of elements must be*/
								/* < log(base2) of 50 */

	lst[1] = 1;
	ust[1] = array_size;
	p = 1;

	while(p > 0) {
		if(lst[p] >= ust[p])
			p--;
		else {
			i = lst[p] - 1;
			j = ust[p];
			pivlin = array1[j];
			while(i < j) {
				for(i=i+1; array1[i] < pivlin; i++)
					;
				for(j=j-1; j > i; j--)
					if(array1[j] <= pivlin) break;
				if(i < j) {
					temp1     = array1[i];
					array1[i] = array1[j];
					array1[j] = temp1;
					
					temp2     = array2[i];
					array2[i] = array2[j];
					array2[j] = temp2;
				}
			}
			
			j = ust[p];

			temp1     = array1[i];
			array1[i] = array1[j];
			array1[j] = temp1;

			temp2     = array2[i];
			array2[i] = array2[j];
			array2[j] = temp2;

			if(i-lst[p] < ust[p] - i) {
				lst[p+1] = lst[p];
				ust[p+1] = i - 1;
				lst[p]   = i + 1;
			}
			else {
				lst[p+1] = i + 1;
				ust[p+1] = ust[p];
				ust[p]   = i - 1;
			}
			p = p + 1;
		}
	}
	return;

}





static void pair_align(int seq_no,int l1,int l2)
{
	int pot[8],i,j,l,m,flag,limit,pos,tl1,vn1,vn2,flen,osptr,fs;
	int tv1,tv2,encrypt,subt1,subt2,rmndr;
	int residue;
	
	if(dnaflag) {
		for(i=1;i<=ktup;++i)
			pot[i] = (int) pow((double)4,(double)(i-1));
		limit = (int) pow((double)4,(double)ktup);
	}
	else {
		for (i=1;i<=ktup;i++)
           		pot[i] = (int) pow((double)(max_aa+1),(double)(i-1));
		limit = (int) pow((double)(max_aa+1),(double)ktup);
	}
	
	tl1 = (l1+l2)-1;
	
	for(i=1;i<=tl1;++i) {
		slopes[i]=displ[i]=0;
		diag_index[i] = i;
	}
	

/* increment diagonal score for each k_tuple match */
/* Attempt at guessing the best band by looking at identities*/

	for(i=1;i<=limit;++i) 
	        {
		vn1=zzc[i];
		while(TRUE) 
		        {
			if(!vn1) break;
			vn2=zzd[i];
			while(vn2 != 0) 
			        {
				osptr=vn1-vn2+l2;
				++displ[osptr]; /*PLUG THE Pos Dependant Scheme Here!!!! (For Id only)*/
				vn2=zzb[vn2];
				}
			vn1=zza[vn1];
			}
		}

/* choose the top SIGNIF diagonals */

	des_quick_sort(displ, diag_index, tl1);

	j = tl1 - signif + 1;
	if(j < 1) j = 1;
 
/* flag all diagonals within WINDOW of a top diagonal */

	for(i=tl1; i>=j; i--) 
		if(displ[i] > 0) {
			pos = diag_index[i];
			l = (1  >pos-window) ? 1   : pos-window;
			m = (tl1<pos+window) ? tl1 : pos+window;
			for(; l <= m; l++) 
				slopes[l] = 1;
		}

	for(i=1; i<=tl1; i++)  displ[i] = 0; /*reset the diagonals score*/

	
	next=curr_frag=maxsf=0;	
	for(i=1;i<=(l1-ktup+1);++i) 
	        {
		encrypt=flag=0;
		for(j=1;j<=ktup;++j)
		{residue = seq_array[seq_no][i+j-1];if((residue<0) || (residue>max_aa)){flag=TRUE; break;}encrypt += ((residue)*pot[j]);}
		if(flag) continue;
		else flag=FALSE;
		
		++encrypt;	
		vn2=zzd[encrypt];
      
		/*now trying to match i-ktup and vn2-ktup*/
		while(TRUE) 
		        {
			if(!vn2) 
			        {
				flag=TRUE;
				break;
				}
			osptr=i-vn2+l2;      /*osptr=Diagonal under investigation*/
			if(slopes[osptr]!=1) /*Get the next diagonal if that one is not flagged*/
			        {
				vn2=zzb[vn2];
				continue;
			        }
			flen=0;
			fs=ktup;
			next=maxsf;	
			
			/* A-loop*/
			while(TRUE) 
			    {
			    if(!next)
				  {
				  ++curr_frag;
				  if(curr_frag>=2*max_aln_length) 
					  {

					  return;
					  }
				  displ[osptr]=curr_frag;
				  put_frag(fs,i,vn2,flen); /*sets the coordinates of the fragments*/
				  }
			    else 
				  {
				  tv1=accum[1][next];
				  tv2=accum[2][next];
				  if(frag_rel_pos(i,vn2,tv1,tv2)) 
				        {
					if(i-vn2==accum[1][next]-accum[2][next]) 
					    {
					    if(i>accum[1][next]+(ktup-1))
						fs=accum[0][next]+ktup;
					    else  
					          {
						  rmndr=i-accum[1][next];
						  fs=accum[0][next]+rmndr;
					          }
					    flen=next;
					    next=0;
					    continue;
					    }
					else 
					    {
					    if(displ[osptr]==0)
					    subt1=ktup;
					    else 
					       {
					       if(i>accum[1][displ[osptr]]+(ktup-1))
					       subt1=accum[0][displ[osptr]]+ktup;
					       else 
					          {
						  rmndr=i-accum[1][displ[osptr]];
						  subt1=accum[0][displ[osptr]]+rmndr;
						  }
					       }
					    subt2=accum[0][next]-wind_gap+ktup;
					    if(subt2>subt1) 
					        {
						flen=next;
						fs=subt2;
						}
					    else 
					        {
						flen=displ[osptr];
						fs=subt1;
						}
					    next=0;
					    continue;
					    }
				           
					}
				  else 
				        {
					next=accum[4][next];
					continue;
					}
				  }
			    break;
			}
		/*
		* End of Aloop
		*/
		
			vn2=zzb[vn2];
		}
	}

}		 

int ** show_pair(int istart, int iend, int jstart, int jend, int *in_seqlen_array, char **in_seq_array, int dna_ktup, int dna_window, int dna_wind_gap, int dna_signif,int prot_ktup, int prot_window,int prot_wind_gap,int prot_signif, int in_nseqs,int in_dnaflag, int in_max_aa, int in_max_aln_length  )
{
	int i,j,dsr;
	double calc_score;
	int **tmat;
	
	seqlen_array=vcalloc ( in_nseqs+1, sizeof(int));
	for ( i=0; i< in_nseqs; i++)seqlen_array[i+1]=in_seqlen_array[i];
	
	
	seq_array=declare_char ( in_nseqs+1, in_max_aln_length);
	for ( i=0; i< in_nseqs; i++)sprintf (seq_array[i+1], "%s",in_seq_array[i]); 

	
	nseqs=in_nseqs;
	dnaflag=in_dnaflag;
	max_aa=in_max_aa;
	max_aln_length=in_max_aln_length;

	
	tmat=declare_int ( nseqs+1, nseqs+1);
	accum=declare_int( 5, 2*max_aln_length+1);

	displ      = (int *) vcalloc( (2*max_aln_length +1), sizeof (int) );
	slopes     = (char *)vcalloc( (2*max_aln_length +1) , sizeof (char));
	diag_index = (int *) vcalloc( (2*max_aln_length +1) , sizeof (int) );

	zza = (int *)vcalloc( (max_aln_length+1),sizeof (int) );
	zzb = (int *)vcalloc( (max_aln_length+1),sizeof (int) );

	zzc = (int *)vcalloc( (max_aln_length+1), sizeof (int) );
	zzd = (int *)vcalloc( (max_aln_length+1), sizeof (int) );

        if(dnaflag) {
                ktup     = dna_ktup;
                window   = dna_window;
                signif   = dna_signif;
                wind_gap = dna_wind_gap;
        }
        else {
                ktup     = prot_ktup;
                window   = prot_window;
                signif   = prot_signif;
                wind_gap = prot_wind_gap;
        }

	for(i=istart+1;i<=iend;++i) 
	        {
		if(dnaflag)
			make_n_ptrs(zza,zzc,i,seqlen_array[i]);
		else
			make_p_ptrs(zza,zzc,i,seqlen_array[i]);
		for(j=MAX(jstart+1, i+1);j<=jend;++j) 
		        {
			    if (i!=j)
			       {
				   if(dnaflag)
				       make_n_ptrs(zzb,zzd,j,seqlen_array[j]);
				   else
				       make_p_ptrs(zzb,zzd,j,seqlen_array[j]);

				   pair_align(i,seqlen_array[i],seqlen_array[j]);

				   if(!maxsf)
				       calc_score=0.0;
				   else {
				       calc_score=(double)accum[0][maxsf];
				       if(percent) {
					   dsr=(seqlen_array[i]<seqlen_array[j]) ?
					       seqlen_array[i] : seqlen_array[j];
					   calc_score = (calc_score/(double)dsr) * 100.0;
				       }
				   }

				   tmat[i-1][j-1] = (int)(((100.0 - calc_score)/100.0)*1000);
				   tmat[j-1][i-1] = (int)(((100.0 - calc_score)/100.0)*1000);
				   fprintf ( stderr, "\r[%d %d]=> %.2f",i, j, (float)calc_score );
			       }
			}
		}

	free_int ( accum, -1);

	vfree(displ);
	vfree(slopes);
        vfree(diag_index);

	vfree(zza);
	vfree(zzb);
	vfree(zzc);
	vfree(zzd);
	return tmat;
}


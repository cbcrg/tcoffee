#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <sys/times.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "motif_lib_header.h"

Motif*	build_motif_tree ( Motif *motif_tree, int len, int n, char *al,int ns, int n_motifs)
        {
	int a;
	    
	if ( n==len)
	   {
	   motif_tree->motif[n]='\0';
	   
	   motif_tree->motif_offset=motif_tree->n_motifs;
	   motif_tree->leaf_list[motif_tree->n_motifs++]=motif_tree;
	   return motif_tree;
	   }
	 
	if ( motif_tree==NULL)
	    {
	    motif_tree=declare_motif(len, al, ns, n_motifs);
	    motif_tree->leaf_list=calloc ( n_motifs, sizeof (Motif*));
	    }  
	
	for ( a=0; a< ns; a++)
	    {
	    (motif_tree->child[a])=declare_motif(len, al, ns, n_motifs);
	    (motif_tree->child[a])->leaf_list=motif_tree->leaf_list;
	    
	    sprintf ( (motif_tree->child[a])->motif, "%s", motif_tree->motif);
	   
	    (motif_tree->child[a])->motif[n]=al[a];
	    (motif_tree->child[a])->parent=motif_tree;
	    (motif_tree->child[a])->n_motifs=motif_tree->n_motifs;
	    
	    motif_tree->child[a]=build_motif_tree ( motif_tree->child[a], len,n+1,al,ns, n_motifs);
	    motif_tree->n_motifs=(motif_tree->child[a])->n_motifs;
	    }
	return motif_tree;
	}	

Motif* find_motif_in_tree ( char *motif,int p, int len, Motif *mtree)
        {
	int a, b;
	
	for ( a=0; a< mtree->n_children; a++)
	    {
	     if ((mtree->child[a])->motif[p]==motif[p])
	        {

		if ( p==len-1)return mtree->child[a];
		else  return find_motif_in_tree ( motif,p+1,len,mtree->child[a]);
		}
	    }
	return NULL;
	}
	

int bin2int ( int *in, int len)
    {
    int tot=0, a;
    static *array;

    if ( array==NULL)
       {
       array=calloc ( len, sizeof (int));
       for ( a=0; a< len; a++)
           {
	   array[a]=(int)pow (2, (double) a);
	   }
       }
    
    for ( a=0; a< len; a++)
        {
	tot+=in[a]*array[a];
	}
    return tot;
    }
Motif * declare_motif (int len, char *al, int ns, int n_motif)
    {
    Motif * m;

    m=calloc ( 1, sizeof (Motif));
    m->motif=calloc ( len+1, sizeof (char));   
    m->child=calloc (ns, sizeof ( Motif*));
    m->n_children=ns;
    return m;
    }

int motif_is_in_seq ( char *seq, char *motif)
    {
    int a, b, c, l,w;
    static char *buf;

    if ( buf==NULL)buf=vcalloc (100, sizeof (char));

    l=strlen ( seq);
    w=strlen ( motif);
    for ( a=0; a< l-(w-1); a++)
        {
	sprintf ( buf, "%s",seq);
	buf[a+w]='\0';
	if ( same_motif ( buf+a, motif))return 1;

	}
    return 0;
    }

int same_motif ( char *m1, char *m2)
    {
    int l, x;
    static int **lu;
    int a, g, c, t,u, i, y, s;
    
    if (lu==NULL)
       {
       a='A'-'*';
       g='G'-'*';
       c='C'-'*';
       t='T'-'*';
       u='U'-'*';
       i='I'-'*';
       y='Y'-'*';
       s=0;

       
       lu=calloc ( 100, sizeof (int*));
	   for ( x=0; x< 100; x++)lu[x]=calloc ( 100, sizeof (int));
	   
       lu[a][a]=1;
       lu[g][g]=1;
       lu[c][c]=1;
       lu[t][t]=1;
       lu[u][u]=1;
       lu[u][t]=lu[t][u]=1;
       lu[i][i]=1;
       lu[y][y]=1;
       lu[i][g]=lu[g][i]=1;
       lu[i][a]=lu[a][i]=1;
       lu[y][c]=lu[c][y]=1;
       lu[y][t]=lu[t][y]=1;
       lu[y][u]=lu[u][y]=1;
       lu[s][a]=lu[a][s]=1;
       lu[s][g]=lu[g][s]=1;
       lu[s][c]=lu[c][s]=1;
       lu[s][t]=lu[t][s]=1;
       lu[s][u]=lu[u][s]=1;
       }
       

    if ( (l=strlen (m1))!=(strlen (m2)))return 0;
    else
       {
       for ( a=0; a< l; a++)
           {
	   if ( lu[m1[a]-'*'][m2[a]-'*']);
	   else return 0;
	   }
       }

    return 1;
    }
Oligo * read_oligo_list ( char *fname)
    {
    Oligo *O;
    FILE *fp;
    int a, b;
    

    
    O=vcalloc (1, sizeof (Oligo));
    O->ALPHABET=vcalloc ( 100, sizeof (char));
    O->EALPHABET=vcalloc ( 100, sizeof (char));
    O->AMBIGUITIES=vcalloc ( 100, sizeof (char));
    
    fp=vfopen ( fname, "r");
    fscanf ( fp, "ALPHABET %s\n", O->ALPHABET);
    fscanf ( fp, "AMBIG_ALPHABET %s\n", O->AMBIGUITIES);
    if ( O->AMBIGUITIES[0]=='@')O->AMBIGUITIES[0]='\0';
    fscanf ( fp, "WORD_SIZE %d\n", &O->WSIZE);
    fscanf ( fp, "NSEQ %d\n", &O->NSEQ);
    fscanf ( fp, "LEN %d\n", &O->LEN);
    fscanf ( fp, "SCORE %d", &a);
    sprintf ( O->EALPHABET, "%s%s", O->ALPHABET, O->AMBIGUITIES);
   
    O->seq=declare_char ( O->NSEQ, O->LEN+1);
    for ( a=0; a< O->NSEQ; a++)
	{
	fscanf ( fp, "%*s\n%s\n",O->seq[a]);
	}
    vfclose (fp);
    return O;
    }
int weight_motif ( char *motif)
   {
   int l,a;
   static int *lu;
   float tot=1;
   
   if (!lu)
      {
      lu=vcalloc (100, sizeof (int));
      lu['A'-'*']=1;
      lu['G'-'*']=1;
      lu['C'-'*']=1;
      lu['T'-'*']=1;
      lu['U'-'*']=1;
      lu['I'-'*']=2;
      lu['Y'-'*']=2;
      lu['*'-'*']=4;
      }

   l=strlen (motif);
   for ( a=0; a< l; a++)
       {
       tot=lu[motif[a]-'*']*tot;
       }
   tot=(1/tot)*100;
   
   return ((int)tot==0)?1:(int)tot;
   }

FILE *print_motif_and_bin ( FILE *fp, char *p, Motif *M, int size)
      {
      int c;
      fprintf ( fp, "%s%s ", p,M->motif);
      for ( c=0; c<size; c++)fprintf (fp,"%d", M->bin_val[c]);
      fprintf (fp, "\n");
      return fp;
      }
int get_bin_edit_distance ( int *b1, int *b2, int l)
    {
    int a,d=0;

    for (d=0, a=0; a< l; a++)
	d+=(b1[a]!=b2[a]);
    return d;
    }
int get_bin_edit_distance_fn ( int *b1, int *b2, int l)
    {
    int a,d=0;

    for (d=0, a=0; a< l; a++)
        {
	if ( b1[a] && !b2[a])return -1;
	d+=(!b1[a] && b2[a]);
	}
    return d;
    }
int get_bin_edit_distance_fp ( int *b1, int *b2, int l)
    {
    int a,d=0;

    for (d=0, a=0; a< l; a++)
	{
	if (!b1[a] && b2[a])return -1;    
	d+=(b1[a] && !b2[a]);	
	}
    return d;
    }
int bin_is_included_in_bin ( int *b1, int *b2, int l)
   {
   int a;
   
   for ( a=0; a< l; a++)if ( b1[a] && !b2[a])return 0;
   return 1;
   }
int count_bits ( int *b, int l)
   {
   int a, t;
   
   for (t=0, a=0; a< l; a++)t+=b[a];
   return t;
   }
/**************HASHING UTILITIES**************************/
Bin_node * find_val_in_bin_tree ( Bin_node ** Node_list,Bin_node *N,int index, int *bin_val, int n, int len, int *used)
   {
   Bin_node *C;

   if (used[0]==0){used[0]=1;N=Node_list[0];N->left_child=NULL; N->right_child=NULL;N->index=0;}
  
   if ( n==len)
      {
      N->n++;
      if (N->n==0)N->n++;
      if (N->n>N->max_list)
         {
	 N->list=vrealloc ( N->list, N->n*sizeof (int));
	 N->max_list=N->n;	 
	 }
      N->list[N->n-1]=index;
      return N;
      }
   else 
      {
      
      if      ( N->left_child &&  !bin_val[n]){C=N->left_child ;}
      else if ( N->right_child &&  bin_val[n]){C=N->right_child;}
      else
         {
	 if (    !bin_val[n]){C=N->left_child =Node_list[used[0]++];}
	 else if (bin_val[n]){C=N->right_child=Node_list[used[0]++];}
	      	      
	 C->left_child=NULL;
	 C->right_child=NULL;
	 C->n=0;
	 C->score=0;
	 C->used=0;
	 C->index=used[0]-1;
	 }
     
     
      return find_val_in_bin_tree ( Node_list,C,index, bin_val, n+1,    len, used);
 
      }
   }

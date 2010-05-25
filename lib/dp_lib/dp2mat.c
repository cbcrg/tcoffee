#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "seq2mat.h"

double p_add(int n, ...);
double get_max_double ( int n, double *list);

void print_matrice ( double **m, int l1, int l2, char *s);
double is_defined ( double a, int g);
double is_defined2 ( double a);
/******DP MATRICES*******/
int I= 0;
int D =1;
int S =2;
int C =3;
/******DP MATRICES*******/
int maximise=1;

main ( int argc, char *argv[])
	{
	double *****dp_matrices;
	Sequence *S;
	int gep, gop, T;
	FILE *OUT;
	char seq_file[1000];
	char matrix[100];
	int a;
	
	sprintf ( matrix, "blosum62mt");
	gop=220;
	gep=100;
	OUT=stdout;
	T=90;
	
	for ( a=1; a< argc; a++)
		{
		
		if ( strcmp ( argv[a], "-matrix")==0)
			{
			sprintf ( matrix, "%s", argv[a+1]);
			a++;
			}
		else if (strcmp ( argv[a], "-seq")==0)
			{
			sprintf ( seq_file, "%s", argv[a+1]);
			a++; 
			}
		else if ( strcmp ( argv[a], "-gop")==0)
			{
			gop=atoi(argv[a+1]);
			a++;
			}
		else if ( strcmp ( argv[a], "-gep")==0)
			{
			gep=atoi(argv[a+1]);
			a++;
			}
		else if ( strcmp ( argv[a], "-T")==0)
			{
			T=atoi(argv[a+1]);
			a++;
			}
		else if ( strcmp ( argv[a], "-out")==0)	
			{
			if ( strcmp ( argv[a+1], "stdout")==0)
				OUT=stdout;
			else
				OUT=fopen ( argv[a+1], "w");
			a++;
			}
		}
	S=read_sequences ( seq_file);
	fprintf ( OUT, "%d\n", S->nseq);
	for ( a=0; a<S->nseq; a++)
		fprintf ( OUT, "%s %d %s\n", S->name[a], S->len[a], S->seq[a]);
		
	OUT=make_matrices (S, gep, gop, matrix, T, OUT);
	fclose (OUT);
	}

FILE *make_matrices ( Sequence *S, int gep, int gop, char *matrix_name, int T, FILE *OUT)
	{
	int a, b, c, d, e, f;
	double *** dpm;
	double ***tm;
	int **matrix;
	Alignment *A;
	int score=0;
	
	
	dpm=calloc (4, sizeof (double**));
	for ( c=0; c<4; c++)dpm[c]=declare_double ( S->max_len+2,S->max_len+2);	
	
	A=declare_Alignment(S);	
	
		 
	matrix=read_matrice (matrix_name);	 		
	for ( b=0; b<S->nseq-1; b++)
		for ( c=b+1; c<S->nseq; c++)
			{
			make_dp (S->seq[b], S->seq[c], dpm,matrix,-gop,-gep);
			trace_back (S->seq[b], S->seq[c], dpm, A);
			OUT=output_matrix ( S, b, c, dpm, T, OUT);
			
			fprintf ( stderr, "\n%d %d Score=%d",b, c, (int)A->score);
			print_aln (A);
			
			}
	return  OUT;
	}

FILE* output_matrix ( Sequence *SEQ, int s1, int s2, double ***dpm, int T, FILE *OUT)
	{
	int l1, l2;
	double delta, min, max;
	int a, b, c;
	
	l1=strlen ( SEQ->seq[s1]);
	l2=strlen (SEQ->seq[s2]);
	
	min=max=dpm[S][1][1];
	for ( a=1; a<=l1; a++)
		for ( b=1; b<=l2; b++)
			{
			max=(max<dpm[S][a][b])?dpm[S][a][b]:max;
			min=(min>dpm[S][a][b])?dpm[S][a][b]:min;
			}
	delta=(max-min)-((max-min)*(100-T))/(100);
	
	fprintf ( OUT, "#%d %d\n", s1, s2);
	
	
	for ( a=1; a<=l1; a++)
		for ( b=1; b<=l2; b++)
			{
			if ((dpm[S][a][b] -min)>=delta)
				fprintf ( OUT, "%d %d %d\n", a, b, (int)dpm[S][a][b]);
				
			}
	return OUT;
	}
	
void make_dp ( char *seq1, char *seq2, double ***M, int **matrix,int gop, int gep)
 	{
	int i, j, offset, t, a, b,d;
	int l1, l2, max_len;
	static double ***m_d; 
	static double ***m_r;
	static double ***m_d2;
	static len_array;
	double cost;
	int t_gop;
	int i_gop,d_gop;
	int t_gep, i_gep, d_gep;
	
	
	
	maximise=1;
	l1=strlen ( seq1);
	l2=strlen ( seq2);
	max_len =((l1>l2)?l1:l2)+2;
	
	if (max_len>len_array)
		{
		if ( m_d==NULL)
			{
			m_d=calloc ( 4, sizeof (double**));
			m_d2=calloc ( 4, sizeof (double**));
			m_r=calloc ( 4, sizeof (double**));
			}
		for ( a=0; a<4; a++)
			{
			
			free_double ( m_d[a], len_array);
			m_d[a]=declare_double ( max_len, max_len);
			free_double ( m_d2[a], len_array);
			m_d2[a]=declare_double ( max_len, max_len);
			free_double ( m_r[a], len_array);
			m_r[a]=declare_double ( max_len, max_len);
			}
		len_array=max_len;
		}
	
			

	t_gop=gop;		
	t_gep=gep;
	
	m_d[C][0][0]=UNDEFINED;
	m_d[D][0][0]=UNDEFINED;
	m_d[I][0][0]=UNDEFINED;
	
	m_d2[C][0][0]=UNDEFINED;
	m_d2[D][0][0]=UNDEFINED;
	m_d2[I][0][0]=UNDEFINED;
	
	t=t_gop;
	for ( j=1; j<=l2; j++)
		{
		t+=t_gep;
		
		m_d [C][0][j]=t;
		m_d [D][0][j]=t+gop;
		m_d [I][0][j]=UNDEFINED;
		m_d [S][0][j]=UNDEFINED;
		
		}
	t=t_gop;
	for ( i=1; i<=l1; i++)
		{
		t+=t_gep;
		m_d [C][i][0]=t;
		m_d [D][i][0]=UNDEFINED;
		m_d [I][i][0]=t+gop;
		m_d [S][i][0]=UNDEFINED;
		for ( j=1; j<=l2; j++)
			{
			
			
			i_gop=(i==l1)?t_gop:gop;
			i_gep=(i==l1)?t_gep:gep;
			m_d[I][i][j]=best (2, maximise,&offset,m_d[I][i][j-1], m_d[C][i][j-1]+(double)i_gop)+(double)i_gep;
			m_d2[I][i][j]=((offset==0)?m_d[I][i][j-1]:(m_d[C][i][j-1]+i_gop))-i_gop;
			
				
			d_gop=((j==l2)?t_gop:gop);
			d_gep=((j==l2)?t_gep:gep);
			m_d[D][i][j]=best (2, maximise,&offset,m_d[D][i-1][j], m_d[C][i-1][j]+(double)d_gop)+(double)d_gep;
			m_d2[D][i][j]=((offset==0)?m_d[D][i-1][j]:(m_d[C][i-1][j]+d_gop))-d_gop;
			
			
			
			m_d [S][i][j]=((m_d[C][i-1][j-1]==UNDEFINED)?0:m_d[C][i-1][j-1])+(double)matrix[seq1[i-1]-'a'][seq2[j-1]-'a'];
			m_d2[S][i][j]=(m_d[C][i-1][j-1]==UNDEFINED)?0:m_d[C][i-1][j-1];
			
				
			m_d[C][i][j]=best	(3, maximise,&offset, m_d[D][i][j], m_d[I][i][j], m_d[S][i][j]);
			}
			
		}
	
	
	t=t_gop;
	for ( j=1; j<=l2; j++)
		{
		t+=t_gep;
		
		m_d [C][0][j]=t;
		m_d [D][0][j]=t;
		m_d [I][0][j]=t;
		m_d [S][0][j]=t;
		
		m_d2[C][0][j]=m_d [C][0][j-1]-t_gop;
		m_d2[D][0][j]=m_d [D][0][j-1]-t_gop;
		m_d2[I][0][j]=m_d [I][0][j-1]-t_gop;
		m_d2[S][0][j]=m_d [S][0][j-1]-t_gop;
		
		}
	t=t_gop;
	for ( i=1; i<=l1; i++)
		{
		t+=t_gep;
		m_d [C][i][0]=t;
		m_d [D][i][0]=t;
		m_d [I][i][0]=t;
		m_d [S][i][0]=t;
		
		m_d2[C][i][0]=m_d[C][i-1][0]-t_gop;
		m_d2[D][i][0]=m_d[D][i-1][0]-t_gop;
		m_d2[I][i][0]=m_d[I][i-1][0]-t_gop;
		m_d2[S][i][0]=m_d[S][i-1][0]-t_gop;
		}
	
	
	
	
	
	
	m_r[C][l1+1][l2+1]=UNDEFINED;
	m_r[D][l1+1][l2+1]=UNDEFINED;
	m_r[I][l1+1][l2+1]=UNDEFINED;
	m_r[S][l1+1][l2+1]=UNDEFINED;	
	
	
	t=t_gop;
	for ( j=l2; j>0; j--)
		{
		t+=t_gep;
		m_r[C][l1+1][j]=t;
		m_r[I][l1+1][j]=UNDEFINED;
		m_r[D][l1+1][j]=t+gop;
		m_r[S][l1+1][j]=UNDEFINED;
		}
	t=t_gop;
	for ( i=l1; i>0; i--)
		{
		t+=t_gep;
		m_r[C][i][l2+1]=t;
		m_r[I][i][l2+1]=t+gop;
		m_r[D][i][l2+1]=UNDEFINED;
		m_r[S][i][l2+1]=UNDEFINED;
		
		
		for ( j=l2; j>0; j--)
			{
			
			i_gop=(i==1)?t_gop:gop;
			i_gep=(i==1)?t_gep:gep;
			m_r[I][i][j]=best (2, maximise,&offset,m_r[I][i][j+1], m_r[C][i][j+1]+(double)i_gop)+(double)i_gep;
				
			d_gop=(j==1)?t_gop:gop;
			d_gep=(j==1)?t_gep:gep;
			m_r[D][i][j]=best (2, maximise,&offset,m_r[D][i+1][j], m_r[C][i+1][j]+(double)d_gop)+(double)d_gep;
			
			
			
			m_r[S][i][j]=((m_r[C][i+1][j+1]==UNDEFINED)?0:m_r[C][i+1][j+1])+(double)matrix[seq1[i-1]-'a'][seq2[j-1]-'a'];
			m_r[C][i][j]=best	(3, maximise,&offset, m_r[D][i][j], m_r[I][i][j], m_r[S][i][j]);
			
			}
		}
	
	t=t_gop;
	for ( j=l2; j>0; j--)
		{
		t+=t_gep;
		m_r[C][l1+1][j]=t;
		m_r[I][l1+1][j]=t;
		m_r[D][l1+1][j]=t;
		m_r[S][l1+1][j]=t;
		}
	t=t_gop;
	for ( i=l1; i>0; i--)
		{
		t+=t_gep;
		m_r[C][i][l2+1]=t;
		m_r[I][i][l2+1]=t;
		m_r[D][i][l2+1]=t;
		m_r[S][i][l2+1]=t;
		}
	cost=m_r[C][1][1];
	
	
	for ( i=0; i<=l1+1; i++)
		for (j=0; j<=l2+1; j++)
			for ( d=0; d<4; d++)M[d][i][j]=0;
			
	
	for ( i=0; i<=l1+1; i++)M[S][i][l2+1]=M[I][i][l2+1]=UNDEFINED;
	for ( j=0; j<=l2+1; j++)M[S][l1+1][j]=M[D][l1+1][j]=UNDEFINED;
		
	
	
	if ( m_d[C][l1][l2]!=cost)
		{
		fprintf ( stderr, "\nPB: COST down=%d COST up=%d", (int)m_d[C][l1][l2], (int)cost);	
		exit (0);
		}
			
	for ( i=1; i<=l1+1; i++)
		for (j=1; j<=l2+1; j++)
			{
			
			if ( M[I][i][j]!=UNDEFINED)
				{
				M[I][i][j]=(is_defined2(m_d2[I][i-1][j])+is_defined2(m_r[I][i][j]));
				if ( M[I][i][j]>cost)fprintf ( stderr, "\nPB in I %d %d: %d = %d+%d COST=%d", i, j, (int)M[I][i][j],(int)m_d2[I][i-1][j],(int)m_r[I][i][j],(int)cost);
				}
			if ( M[D][i][j]!=UNDEFINED)
				{
				M[D][i][j]=(is_defined2(m_d2[D][i][j-1])+is_defined2(m_r[D][i][j]));
				if ( M[D][i][j]>cost)fprintf ( stderr, "\nPB in D %d %d: %d = %d+%d COST=%d", i, j,(int)M[D][i][j],(int) m_d2[D][i][j-1],(int)m_r[D][i][j],(int)cost);
				}
			if ( M[S][i][j]!=UNDEFINED)
				{
				M[S][i][j]=(is_defined2(m_d2[S][i][j])  +is_defined2(m_r[S][i][j]));
				if ( M[S][i][j]>cost)fprintf ( stderr, "\nPB in S %d %d: %d = %d+%d COST=%d", i, j,(int)M[S][i][j],(int)m_d2[S][i][j] ,(int)m_r[S][i][j],(int)cost);
				}
				
			M[C][i][j]=best(3, maximise, &offset, M[D][i][j], M[I][i][j], M[S][i][j]);
			}
	
			
	
	
	/*
	print_matrice ( M[I], l1+2, l2+2, "I");
	print_matrice ( M[D], l1+2, l2+2, "D");
	print_matrice ( M[S], l1+2, l2+2, "S");
	*/
	}


void trace_back (char *seq1, char *seq2, double ***m, Alignment *A)
	{
	int a, b, c;
	int l,l1, l2, mode;
	int array[10000];
	int i, j, c1, c2;
	double cost;
	l1=strlen ( seq1);
	l2=strlen ( seq2);
	
	/************************/
	/* I nsertion in Seq2   */
	/* D eletion  in Seq1	*/
	/* S ubstitution	*/
	
	/* 0=> x/x              */
	/* 1=> x/-		*/
	/* 2=> -/x              */
	/************************/
	
	i=l1+1;
	j=l2+1;
	l=0;
	
	cost=best ( 3, maximise, &mode,  m[S][i-1][j-1], m[D][i-1][j], m[I][i][j-1]);	
	while (!((i==1) || (j==1)))
		{
		
		if ( cost!=best ( 3, maximise, &mode,  m[S][i-1][j-1], m[D][i-1][j], m[I][i][j-1]))
			{
			fprintf ( stderr, "\ncost1=%d cost2=%d i=%d j=%d", cost, (int)best ( 3, maximise, &mode,  m[S][i-1][j-1], m[D][i-1][j], m[I][i][j-1]),i, j);
			fprintf ( stderr, "\n>A\n%s\n>B\n%s\n", seq1, seq2);
			
			}
			
		array[l++]=mode;
		
		if (mode==0 || mode==1)i--;
		if (mode==0 || mode==2)j--;
		}
	A->score=cost;
	for ( ; i>1; i--)array[l++]=1;
	for ( ; j>1; j--)array[l++]=2;	
	l--;
	
	for ( c1=0, c2=0,a=l; a>=0; a--)
		{
		if ( array[a]==0 || array[a]==1)A->seq_al[0][l-a]=seq1[c1++];
		if ( array[a]==0 || array[a]==2)A->seq_al[1][l-a]=seq2[c2++];		 
		if ( array[a]==1)A->seq_al[1][l-a]='-';
		if ( array[a]==2)A->seq_al[0][l-a]='-';
		}
	if ( c1!=l1)
		{
		fprintf ( stderr, "\nRESIDUES MISSED IN seq1");
		exit (0);
		}
	
	if ( c2!=l2)
		{
		fprintf ( stderr, "\nRESIDUES MISSED IN seq1");
		exit (0);
		}
			
	A->seq_al[0][l+1]=A->seq_al[1][l+1]='\0';
	
	A->nseq=2;
	A->len_aln=l+1;
	}
	

	
		

/**********************************DEBUG ROUTINES******************************/

double get_cost ( Sequence *SEQ,int M, int r1, int r2,int s1, int s2, double ***** m)
	{
	int a, b, c, e,s;
	double v1, v2, v, tot_score, score;
	
	score=0;
	for ( s=0; s<SEQ->nseq; s++)
		{
		if ( s==s1|| s==s2)
			{
			if ( s1>s2)tot_score=is_defined2 ( m[s2][s1][M][r2][r1]);
			else
				tot_score=is_defined2 ( m[s1][s2][M][r1][r2]);
			}	
			
		else
		    {
		    
		    for ( tot_score=UNDEFINED,a=1; a<=SEQ->len[s]; a++)
		    	{
		    	if (s>s1)v1=is_defined2(m[s1][s][M][r1][a]);
		    	else v1=is_defined2(m[s][s1][M][a][r1]);
		    	
		    	if (s>s2)v2=is_defined2(m[s2][s][M][r2][a]);
		    	else v2=is_defined2(m[s][s2][M][a][r2]);
		    	
		    	v=v1+v2;
		    	
		    	if ( tot_score==UNDEFINED)tot_score=v;
		    	else tot_score=p_add (2, tot_score, v);
		    	
		    	}
		    
		    }
		 if ( score==UNDEFINED)score=tot_score;
		 else score=p_add ( 2, score, tot_score);
		    
		 }
	
	
	return score;
	}	
				
void print_matrice ( double **m, int l1, int l2, char *s)
	{
	int a, b;
	
	fprintf ( stderr, "\n%s\n",s);
	for ( a=0; a<l1; a++)
		{
		fprintf ( stderr, "\n");
		for (b=0; b<l2; b++)
			if ( (int)m[a][b]==UNDEFINED)fprintf ( stderr, "%4c ",'*');
			else fprintf ( stderr, "%4d ", (int)m[a][b]);
		}
	fprintf ( stderr, "\n");
	}		 	
double is_defined ( double a, int g)
	{
	if ( a==UNDEFINED)return UNDEFINED;
	else return a+(double)g;
	}
double is_defined2 ( double a)
	{
	if ( a==UNDEFINED)return 0;
	else return a;
	}
	
double p_add(int n, ...)
	{
	va_list ap;
	static double *list1;
	int a;
	double max=0, tot=0;
	
	if (list1==NULL)list1=(double*)calloc (30, sizeof (double));
	
	if ( n>30)
		{
		fprintf ( stderr, "\nFAILED IN p_ADD");
		exit(0);
		}
	va_start(ap, n);
	for ( a=0; a< n; a++)
		{
		list1[a]=va_arg (ap, double);
		}
	va_end (ap);
	
	max=get_max_double ( n, list1);
	
	for ( a=0; a<n;a++)
		if ((list1[a]-=max)>=-100)tot=tot+exp (list1[a]);
	
	
	return max+log (tot);
	}

double get_max_double ( int n, double *list)
	{
	double max=0;
	int a;
	
	max=list[0];
	
	for (a=0; a< n; a++)
		{
		max=(max>=list[a])?max:list[a];
		}
	return max;
	} 
	

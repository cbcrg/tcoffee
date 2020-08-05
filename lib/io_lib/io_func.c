#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "matrices.h"

#define DEFAULT_COLOR -1
#define GAP_COLOR     -2
#define INK_COLOR     -3
	
Sequence * cw_read_sequences ( char *seq_name)
	{
	Sequence *S;


	char **seq=NULL;
	char **name=NULL;
	int *len=NULL;
	int nseq=0;
	int min_len=0;
	int max_len=0;
	int a;

	get_sequence ( seq_name, &nseq, &seq, &name, &len, &min_len, &max_len); 
	
	S=declare_sequence ( min_len, max_len, nseq);
	for ( a=0; a< nseq; a++)sprintf ( S->file[a], "%s", seq_name);
	
	for ( a=0; a<S->nseq; a++)
	    {
	    S->len[a]=len[a];
	    sprintf ( S->name[a],"%s", name[a]);
	    vfree ( name[a]);
	    sprintf ( S->seq[a], "%s", seq[a]);
	    vfree ( seq[a]);
	    }
	vfree (seq);
	vfree (name);
	vfree (len);
	S=get_sequence_type ( S);
	return S;
	}
char *    get_string_type   (char *S)
        {
	int a, l;
	int protein=0,  dna=0,rna=0, tot=0;
	char *type;
	static char *ltype;
	static int warning;
	
	if ( !ltype)
	  declare_name(ltype);
	
	declare_name(type);
	l=(int)strlen (S);
	
	if (l==0)
	  {
	    sprintf ( type, "UNKNOWN");
	    return type;
	  }
	
	for ( a=0; a<l; a++)     
	        {
		    if ( !is_gap(S[a]))
			{
			protein+=( is_aa (S[a]) && !is_dna(S[a]));
			dna+=    ( is_dna(S[a]));
			rna+=    ( is_rna(S[a]));
			tot++;
			}
		}

	protein=(protein*100)/tot;
	dna=(dna*100)/tot;
	rna=(rna*100)/tot;
	
	if ( l<20 && warning==0)
	  {
	    /*add_warning ( stderr, "WARNING: short sequences, use -type=DNA or PROTEIN");*/
	    warning=1;
	  }
	if ( l<20 && ltype && ltype[0])
	  {
	    if (dna==100)
	      {
		sprintf ( type, "%s", ltype);
		ltype[0]='\0';
	      }
	    else
	      sprintf ( type, "PROTEIN");
	    
	  }
	else if ( dna>98 && rna>0)sprintf ( type, "RNA");
	else if ( dna>98)sprintf ( type, "DNA");
	else sprintf ( type, "PROTEIN");
	
	sprintf ( ltype, "%s", type);
	return type;
	}

Alignment* get_aln_type (Alignment *A)
        {
	  if ( !A) return A;
	  
	  if ( A->S && !(A->S)->type)(A->S)->type=(char*)vcalloc (30, sizeof (char));
	  
	  if ( A->S && (A->S)->type[0]!='\0')
	    {
	    ;
	    }
	  else if (A->S!=NULL && (A->S)->type[0]=='\0') 
	    {
	      A->S=get_sequence_type (A->S);
	    }
	  else if (A->S==NULL) 
	      {
		A->S=aln2seq (A);
		A->S=get_sequence_type(A->S);
	      }
	  return A;
	}



char *unset_mode ()
{
  return set_mode (UNSET, NULL);
}
char *store_mode (char *val)
{
  return set_mode (SET, val);
}
char *retrieve_mode ()
{
  return set_mode (GET,NULL);
}
char *set_mode (int mode, char *val)
{
  static char type[100];
  if (mode==SET)
    {
      if (!val)printf_exit (EXIT_FAILURE, stderr, "Error:  programme mode unset in io_func.c:set_seq_type");
      sprintf ( type,"%s", val);
    }
  else if ( mode==GET)
    {
      ;
    }
  else if ( mode==UNSET)
    {
      type[0]='\0';
    }
  else
    {
      printf_exit (EXIT_FAILURE, stderr, "Error: unknown mode in function io_func.c:set_seq_type, use SET, GET or UNSET");
    }
  return type;
}
/************************************************************/

  
	
char *unset_seq_type ()
{
  return set_seq_type (UNSET, NULL);
}
char *store_seq_type (char *val)
{
  return set_seq_type (SET, val);
}
char *retrieve_seq_type ()
{
  return set_seq_type (GET,NULL);
}
char *set_seq_type (int mode, char *val)
{
  static char type[100];
  if (mode==SET)
    {
      if (!val)printf_exit (EXIT_FAILURE, stderr, "Error: sequence type unset in io_func.c:set_seq_type");
      sprintf ( type,"%s", val);
    }
  else if ( mode==GET)
    {
      ;
    }
  else if ( mode==UNSET)
    {
      type[0]='\0';
    }
  else
    {
      printf_exit (EXIT_FAILURE, stderr, "Error: unknown mode in function io_func.c:set_seq_type, use SET, GET or UNSET");
    }
  return type;
}
char * get_array_type (int n, char **seq)
{
  char *buf, *buf2;
  int a, tot=0;
  buf2=(char*)vcalloc ( 100, sizeof (char));
		 
  

  for ( tot=0,a=0; a<n; a++)tot+=(seq[a])?strlen (seq[a]):0;
  buf=(char*)vcalloc (tot+1, sizeof (char));
  for ( a=0; a<n; a++)strcat (buf, (seq[a])?(seq[a]):"");
  sprintf ( buf2, "%s", get_string_type(buf));
  vfree (buf);
  return buf2;
}
Sequence *get_sequence_type (Sequence *S)
{
  if ( !S) return NULL;
  else sprintf ( S->type, "%s", get_array_type (S->nseq, S->seq));
  return S;
}
Sequence *fast_get_sequence_type (Sequence *S)
{
  if ( !S) return NULL;
  else if (S->type && S->type[0])return S;
  else return get_sequence_type(S);
}


void get_sequence (char *seq_file,int *NSEQ, char ***SEQ, char ***SN, int **sl, int *min, int *max)
	{
	int a,b;
	int min_len;
	int max_len;
	int nseq;

	int **SL;
		
	nseq=NSEQ[0]= readseqs ( seq_file,  SEQ, SN, &SL);
	sl[0]=(int*)vcalloc ( nseq, sizeof (int));
	
	 
	min_len= max_len= (SL)[0][0];
	for ( a=0; a<NSEQ[0]; a++)
		{
		sl[0][a]=SL[a][0];
		for ( b=0; b<(SL)[a][0]; b++) 
		 	(SEQ[0])[a][b]=tolower ((SEQ[0])[a][b]);
		 } 
	for ( a=1; a<NSEQ[0]; a++)
		{
		min_len= ( min_len > (SL)[a][0])?(SL)[a][0]:min_len;
		max_len= ( max_len < (SL)[a][0])?(SL)[a][0]:max_len;
		}
	min[0]=min_len;
	max[0]=max_len;
	}

int ** get_matrix   ( char *name, char *format)
       {

       if ( strm ( "blast", format))return read_blast_matrix ( name);
       else if ( strm ( "clustalw", format))return read_matrice(name);
       else
           {
	   fprintf ( stderr, "\nError:\nUnknowm Format %s for Matrix %s[FATAL]", format, name);
	   myexit (EXIT_FAILURE);
	   }
       return NULL;
       }
void display_matrix (int **mat);
int ** read_matrice (char *mat_name_in)
	{
	int a,b,c, l;

	char *AA;
        FILE *fp;
	int **matrice;
	int **matrix2;
	char mat_name[200];
	int *vector=NULL;
	
	AA=(char*)vcalloc (256, sizeof (char));
	sprintf (AA, "abcdefghiklmnpqrstvwxyz");
	l=strlen(AA);

	
	if ( strcmp (mat_name_in, "list")==0)
	  {
	    fprintf ( stderr, "****List of available matrices\n");
	    fprintf ( stderr, "\tidmat\n");
	    fprintf ( stderr, "\tmd_40mt\n");
	    fprintf ( stderr, "\tmd_350mt\n");
	    fprintf ( stderr, "\tmd_250mt\n");
	    fprintf ( stderr, "\tpam120mt\n");
	    fprintf ( stderr, "\tpam160mt\n");
	    fprintf ( stderr, "\tpam350mt\n");
	    
	    fprintf ( stderr, "\tblosum30mt\n");
	    fprintf ( stderr, "\tblosum40mt\n");
	    fprintf ( stderr, "\tblosum45mt\n");
	    fprintf ( stderr, "\tblosum50mt\n");
	    fprintf ( stderr, "\tblosum55mt\n");
	    fprintf ( stderr, "\tblosum62mt\n");
	    fprintf ( stderr, "\tblosum80mt\n");
	    myexit (EXIT_SUCCESS);
	    
	  }
	if ( strm2 (mat_name_in, "pam", "PAM"))sprintf ( mat_name, "pam250mt");
	else if (strm2 (mat_name_in, "blosum", "BLOSUM"))sprintf ( mat_name, "blosum62mt");
	else if (strm3 (mat_name_in, "id", "ID", "idmat"))sprintf ( mat_name, "idmat");
	else sprintf ( mat_name, "%s", mat_name_in);
	
	/*Read Header Matrices*/
	if (strm(mat_name, "pam250mt"))vector=pam250mt;
	else if (strm(mat_name, "idmat"))vector=idmat;
	else if (strm(mat_name, "dna_idmat"))vector=idmat;
	else if (strm(mat_name, "est_idmat"))vector=est_idmat;
	else if (strm(mat_name, "md_350mt"))vector=md_350mt;
	else if (strm(mat_name, "md_250mt"))vector=md_250mt;
	else if (strm(mat_name, "md_120mt"))vector=md_120mt;
	else if (strm(mat_name, "md_40mt" ))vector= md_40mt;
	else if (strm(mat_name, "pam350mt" ))vector=pam350mt;
	else if (strm(mat_name, "pam160mt" ))vector=pam160mt;
	else if (strm(mat_name, "pam120mt" ))vector=pam120mt;
	
	else if (strm(mat_name, "blosum80mt" ))vector=blosum80mt;
	else if (strm(mat_name, "blosum62mt" ))vector=blosum62mt;
	else if (strm(mat_name, "exon2mt" ))vector=blosum62mt;
	else if (strm(mat_name, "blosum62mt3" ))vector=blosum62mt3;
	
	else if (strm(mat_name, "blosum62mt2" ))vector=blosum62mt2;
	else if (strm(mat_name, "blosum55mt" ))vector=blosum55mt;
	else if (strm(mat_name, "blosum50mt" ))vector=blosum50mt;
	else if (strm(mat_name, "blosum45mt" ))vector=blosum45mt;
	
	else if (strm(mat_name, "blosum40mt" ))vector=blosum40mt;
	else if (strm(mat_name, "blosum30mt" ))vector=blosum30mt;
	else if (strm(mat_name, "beta_mat" ))vector=beta_mat;
	else if (strm(mat_name, "alpha_mat" ))vector=alpha_mat;
	else if (strm(mat_name, "coil_mat" ))vector=coil_mat;
	
	else if (strm(mat_name, "rblosum80mt" ))vector=rblosum80mt;
	else if (strm(mat_name, "rblosum62mt" ))vector=rblosum62mt;
	else if (strm(mat_name, "rblosum30mt" ))vector=rblosum30mt;
	
	else if (strm(mat_name, "rpam250mt" ))vector=rpam250mt;
	else if (strm(mat_name, "rpam350mt" ))vector=rpam350mt;
	else if (strm(mat_name, "rpam160mt" ))vector=rpam160mt;
	else if (strm(mat_name, "rpam120mt" ))vector=rpam120mt;

	else if (strm(mat_name, "tmpam250mt" ))vector=tmpam250mt;
	else if (strm(mat_name, "rtmpam250mt" ))vector=rtmpam250mt;

	else if (strm(mat_name, "rbeta_mat" ))vector=rbeta_mat;
	else if (strm(mat_name, "ralpha_mat" ))vector=ralpha_mat;
	else if (strm(mat_name, "rcoil_mat" ))vector=rcoil_mat;
	else if (strm (mat_name, "jtttm250mt"))vector=jtttm250mt;
	
	else if (strm (mat_name, "promoter_tf1"))vector=promoter_tf1;
	else if (strm (mat_name, "blosumR"))vector=blosumR;

	/*Header Matrices*/
	if(vector)
	  {
	    matrice=declare_int ( 256, 256);
	    for (a=0; a<l; ++a)
		{
		 for (b=0;b<=a;++b)
		   {
		     matrice[a][b]=matrice[b][a]=(vector[(a*a+a)/2+b]);
		   }
		}
	  }

	  /*Hard coded Matrices*/
	
	else if ( strm (mat_name, "exon_mt"))
	  {
	    matrice=declare_int ( 256, 256);
	    sprintf (AA, "bojx");
	    l=4;
	    for ( a=0; a<l-1; a++)
	      {
		matrice[a][a]=10;
		for ( b=a+1; b< l-2; b++)
		  {
		    matrice[a][b]=matrice[b][a]=5;
		  }
	      }
	    
	  }
	else if ( strm ( mat_name, "tdamat"))
	  {
	    int x, y;
	    matrice=declare_int ( 256, 256);
	    sprintf ( AA, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	    l=strlen (AA);
	    for ( x=0; x<l; x++)
	      for ( y=0; y<l; y++)matrice[x][y]=-10;
	    for ( x=0; x<l; x++) {matrice[x][x]=0; matrice[x][GAP_CODE]=-3;}
	    return matrice;
	  }
	/*Blast Matrices*/
	else if (strm (mat_name, "strikeP"))
	  {
	    double min=100000;
	    
	    matrice=declare_int ( 256, 256);
	     for (a=0; a<26; ++a)
		{
		 for (b=0;b<26;++b)
		   {
		     if (min>strikeP_mat[a][b])min=strikeP_mat[a][b];
		   }
		}
	     for (a=0; a<26; ++a)
	       {
		 for (b=0;b<26;++b)
		   {
		     
		     matrice[a+'a'][b+'a']=(int)((double)10*(strikeP_mat[a][b]-min));
		     matrice[a+'A'][b+'A']=(int)((double)10*(strikeP_mat[a][b]-min));
		   }
	       }
	    return matrice;
	  }
	else if (strm (mat_name, "strikeR"))
	  {
	    matrice=declare_int ( 256, 256);
	    matrice['G']['A']=matrice['g']['a']=0;
	    matrice['G']['G']=matrice['g']['g']=0;
	    matrice['G']['C']=matrice['g']['c']=6;
	    matrice['G']['T']=matrice['g']['t']=2;
	    matrice['G']['U']=matrice['g']['u']=2;
	    matrice['C']['A']=matrice['g']['a']=0;
	    matrice['C']['G']=matrice['c']['g']=6;
	    matrice['C']['C']=matrice['c']['c']=0;
	    matrice['C']['T']=matrice['c']['t']=0;
	    matrice['C']['U']=matrice['c']['u']=0;
	    matrice['T']['A']=matrice['t']['a']=5;
	    matrice['T']['G']=matrice['t']['g']=0;
	    matrice['T']['C']=matrice['t']['c']=0;
	    matrice['T']['T']=matrice['t']['t']=0;
	    matrice['T']['U']=matrice['t']['u']=0;
	    matrice['U']['A']=matrice['u']['a']=5;
	    matrice['U']['G']=matrice['u']['g']=2;
	    matrice['U']['C']=matrice['u']['c']=0;
	    matrice['U']['T']=matrice['u']['t']=0;
	    matrice['U']['U']=matrice['u']['u']=0;
	    matrice['A']['A']=matrice['a']['a']=0;
	    matrice['A']['G']=matrice['a']['g']=0;
	    matrice['A']['C']=matrice['a']['c']=0;
	    matrice['A']['T']=matrice['a']['t']=5;
	    matrice['A']['U']=matrice['a']['u']=5;
	    return matrice;
	  }
	else if ( check_file_exists(mat_name) && is_blast_matrix (mat_name))
	  {
	    matrice=read_blast_matrix ( mat_name);
	  }

	
	else if ( check_file_exists(mat_name) && is_pavie_matrix (mat_name))
	  {
	    
	    matrice=read_pavie_matrix ( mat_name);
	    return matrice;
	  }
	
	else if ( check_file_exists(mat_name) && is_clustalw_matrix (mat_name))
		{
		  
		  fp=vfopen ( mat_name, "r");
		  while ( (c=fgetc (fp))!='$' && c!=EOF);
		  if ( c==EOF){vfclose (fp);return NULL;};
	
		  matrice=declare_int ( 256, 256);
		  fgetc(fp);
		  for ( a=0; a<l; a++)
		    {
		      for ( b=0; b<=a; b++)
			if (fscanf ( fp, "%d,", &matrice[a][b])==0){vfclose (fp);free_int (matrice, -1);return NULL;};
		      fscanf ( fp, "\n");
		    }		
		  fclose ( fp);
		}
	else {return NULL;}
	

	matrix2=declare_int ( 256, 256);
	for ( b=0; b<l; b++)
		for ( c=0; c<l; c++)
			{
			
			matrix2[toupper(AA[b])-'A'][toupper(AA[c])-'A']=matrice[b][c];
			matrix2[toupper(AA[b])-'A'][tolower(AA[c])-'A']=matrice[b][c];
			matrix2[tolower(AA[b])-'A'][toupper(AA[c])-'A']=matrice[b][c];
			matrix2[tolower(AA[b])-'A'][tolower(AA[c])-'A']=matrice[b][c];
			}
	if ( strm ( "exon2mt", mat_name))
	  {
	    char EE[4];
	    int v=(getenv ("EXONVALUE_4_TCOFFEE"))?atoi("EXONVALUE_4_TCOFFEE"):100;
	    sprintf (EE, "obj");
	    for (a=0; a<strlen (EE); a++)
	      {
		matrix2[toupper(EE[a])-'A'][toupper(EE[a])-'A']=v;
		matrix2[tolower(EE[a])-'A'][toupper(EE[a])-'A']=v;
		matrix2[toupper(EE[a])-'A'][tolower(EE[a])-'A']=v;
		matrix2[tolower(EE[a])-'A'][tolower(EE[a])-'A']=v;
	      }
	  }
	
	/*Correct for RNA: Add the U cost*/
	for (b=0; b<l; b++)
	  {
	    
	    matrix2['U'-'A'][toupper(AA[b])-'A']=matrix2['T'-'A'][toupper(AA[b])-'A'];
	    matrix2[toupper(AA[b])-'A']['U'-'A']=matrix2[toupper(AA[b])-'A']['T'-'A'];
	    matrix2['U'-'A'][toupper(AA[b])-'A']=matrix2['t'-'A'][tolower(AA[b])-'A'];
	    matrix2[tolower(AA[b])-'A']['u'-'A']=matrix2[tolower(AA[b])-'A']['t'-'A'];
					
	  }
	matrix2['U'-'A']['U'-'A']=matrix2['T'-'A']['T'-'A'];
	matrix2['u'-'A']['U'-'A']=matrix2['t'-'A']['T'-'A'];
	matrix2['U'-'A']['u'-'A']=matrix2['T'-'A']['t'-'A'];
	matrix2['u'-'A']['u'-'A']=matrix2['t'-'A']['t'-'A'];
	
	free_int (matrice, -1);
	return (matrix2);

	}

void display_matrix (int **mat)
{
  int a, b;
  
  for ( a=0; a< 26; a++)
    {
      fprintf ( stderr, "\n%c ", a+'a');
      for ( b=0; b< 26; b++)
	fprintf ( stderr, " %2d", mat[a][b]);
    }
}

int **neg_matrix2pos_matrix ( int **matrix)
    {
      int b,c,l, min, max;
      char AA[]="abcdefghiklmnpqrstvwxyzABCDEFGHIKLMNPQRSTVWXYZ";
      l=strlen(AA);
      min=max=matrix[AA[0]-'A'][AA[0]-'A'];
      for ( b=0; b<l; b++)
	for ( c=0; c<l; c++)
	  {
	    min=(matrix[AA[b]-'A'][AA[c]-'A']<min)?matrix[AA[b]-'A'][AA[c]-'A']:min;
	    max=(matrix[AA[b]-'A'][AA[c]-'A']<max)?matrix[AA[b]-'A'][AA[c]-'A']:max;	    	   
	  }
      if (min>0)return matrix;
      else
	 {
	   for ( b=0; b<l; b++)
	     for ( c=0; c<l; c++)
	       {
		 matrix[b][c]=matrix[b][c]-min;
	       }
	 }
      return matrix;
    }

      

/*****************************************************************/
void get_rgb_values ( int val, Color *C)
     {

       /*Colors can be modified, see the definition of 
	 COLOR_FILE in prgogrammes_define.h
       */
       
     static char **html_code;
     static float **ps_code;
     float *r, *g, *b;
     if (!html_code)
       {
	 int rgb255=0;
	 int a, b;
	 int new_scheme=3;
	 int color_file=0;
	 FILE *fp;
	 int classC, n, c;
	 char *cs;
	 
	 if ((cs=getenv ("COLOR_4_TCOFFEE"))) 
	    {
	      if ( strm (cs, "new")|| strm (cs, "NEW")||strm (cs, "monet"))new_scheme=3;
	      else if (strm (cs, "old")|| strm (cs, "OLD")||strm (cs, "delaunay"))new_scheme=0;
	      else if (strm (cs, "2"))new_scheme=2;
	      else if (strm (cs, "1"))new_scheme=1;
	      else 
		{
		  printf_exit ( EXIT_FAILURE,stderr, "\nERROR: color scheme %s is unknown",cs);
		}
	    }

	 color_file=(check_file_exists(COLOR_FILE))?1:0;
	 html_code=declare_char(20, 10);
	 ps_code=declare_float  (20, 3);
	 
	 if (!color_file && new_scheme==3 )
	   {
	     n=0;
	      /*0*/
	     sprintf (html_code[n], "#9B92FF");//blue
	     ps_code[n][0]=155;
	     ps_code[n][1]=146;
	     ps_code[n][2]=255;
	     n++;
	    
	     
	     
	     /*1*/
	     sprintf (html_code[n], "#B4FFB4");//green1
	     ps_code[n][0]=180;
	     ps_code[n][1]=255;
	     ps_code[n][2]=180;
	     n++;
	    

	     /*2*/
	     sprintf (html_code[n], "#BEFFBE");//green2
	     ps_code[n][0]=190;
	     ps_code[n][1]=255;
	     ps_code[n][2]=190;
	     n++;
	     
	     /*3*/
	     sprintf (html_code[n], "#C8FFC8");//green3
	     ps_code[n][0]=200;
	     ps_code[n][1]=255;
	     ps_code[n][2]=200;
	     n++;
	     
	     /*4*/
	     sprintf (html_code[n], "#FFFFB7");//orange 1
	     ps_code[n][0]=255;
	     ps_code[n][1]=255;
	     ps_code[n][2]=183;
	     n++;
	     
	     /*5*/
	     sprintf (html_code[n], "#FFFFAD");//orange 1
	     ps_code[n][0]=255;
	     ps_code[n][1]=255;
	     ps_code[n][2]=173;
	     n++;
	   
	     /*6*/
	     sprintf (html_code[n], "#FFFFAD");//orange 2
	     ps_code[n][0]=255;
	     ps_code[n][1]=255;
	     ps_code[n][2]=173;
	     n++;
	     
	     /*7*/
	     sprintf (html_code[n], "#FFE6E6");//red1
	     ps_code[n][0]=255;
	     ps_code[n][1]=230;
	     ps_code[n][2]=230;
	     n++;
	     
	     /*8*/
	     sprintf (html_code[n], "#FFDCDC");//red2
	     ps_code[n][0]=255;
	     ps_code[n][1]=220;
	     ps_code[n][2]=220;
	     n++;

	     /*9*/
	     sprintf (html_code[n], "#FFC8C8");//red3
	     ps_code[n][0]=255;
	     ps_code[n][1]=200;
	     ps_code[n][2]=200;
	     n++;

	   
	     

	   }
	 else if (!color_file && new_scheme==2 )
	   {
	     n=0;
	      /*0*/
	     sprintf (html_code[n], "#008080");//blue
	     ps_code[n][0]=0.4;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=1;
	     n++;
	     /*0
	     sprintf (html_code[n], "#BC4400");//blue
	     ps_code[n][0]=0.4;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=1;
	     n++;
	     */
	     
	     
	     /*2*/
	     sprintf (html_code[n], "#D05800");//blue green
	     ps_code[n][0]=0.6;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     /*2
	     sprintf (html_code[n], "#E46C03");//Green pea
	     ps_code[n][0]=0.6;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     /*

	     /*3*/
	     sprintf (html_code[n], "#FF8A21");//lighter pink
	     ps_code[n][0]=0.8;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*4*/
	     sprintf (html_code[n], "#FF9E35");//light pink
	     ps_code[n][0]=1.0;
	     ps_code[n][1]=1.0;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*5*/
	     sprintf (html_code[n], "#FFB249");//light orange
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.85;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*6*/
	     sprintf (html_code[n], "#FFC65D");//light orange/bown
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.7;
	     ps_code[n][2]=0;
	     n++;
	   
	     /*7*/
	     sprintf (html_code[n], "#FFDA71");//darker orange/brown
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.6;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*8*/
	     sprintf (html_code[n], "#FFEE85");//mysti rose
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*9*/
	     sprintf (html_code[n], "#FFFF99");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.2;
	     ps_code[n][2]=0;
	     n++;

	     /*9*/
	     sprintf (html_code[n], "#FFFFAD");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.2;
	     ps_code[n][2]=0;
	     n++;

	     /*9*/
	     sprintf (html_code[n], "#FFFFC1");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.2;
	     ps_code[n][2]=0;
	     n++;
	     

	   }
	 else if (!color_file && new_scheme==1 )
	   {
	     n=0;
	     /*0*/
	     sprintf (html_code[n], "#7AD1FF");//blue
	     ps_code[n][0]=0.4;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=1;
	     n++;
	     /*2*/
	     sprintf (html_code[n], "#91FF49");//blue green
	     ps_code[n][0]=0.6;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     /*2*/
	     sprintf (html_code[n], "#A5FF5D");//Green pea
	     ps_code[n][0]=0.6;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*3*/
	     sprintf (html_code[n], "#B9FF71");//lighter pink
	     ps_code[n][0]=0.8;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*4*/
	     sprintf (html_code[n], "#FFC5C3");//light pink
	     ps_code[n][0]=1.0;
	     ps_code[n][1]=1.0;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*5*/
	     sprintf (html_code[n], "#FAD1B9");//light orange
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.85;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*6*/
	     sprintf (html_code[n], "#FAE6B9");//light orange/bown
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.7;
	     ps_code[n][2]=0;
	     n++;
	   
	     /*7*/
	     sprintf (html_code[n], "#FAF1B9");//darker orange/brown
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.6;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*8*/
	     sprintf (html_code[n], "#F8FAB9");//mysti rose
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*9*/
	     sprintf (html_code[n], "#EDFAB9");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.2;
	     ps_code[n][2]=0;
	     n++;

	   }
	 else if (!color_file && new_scheme==0)//old code
	   {
	     n=0;
	     /*0*/
	     sprintf (html_code[n], "#6666FF");
	     ps_code[n][0]=0.4;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=1;
	     n++;
	     /*1*/
	     sprintf (html_code[n], "#00FF00");
	     ps_code[n][0]=0.6;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*2*/
	     sprintf (html_code[n], "#66FF00");
	     ps_code[n][0]=0.8;
	     ps_code[n][1]=1;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*3*/
	     sprintf (html_code[n], "#CCFF00");
	     ps_code[n][0]=1.0;
	     ps_code[n][1]=1.0;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*4*/
	     sprintf (html_code[n], "#FFFF00");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.85;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*5*/
	     sprintf (html_code[n], "#FFCC00");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.7;
	     ps_code[n][2]=0;
	     n++;
	   
	     /*6*/
	     sprintf (html_code[n], "#FF9900");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.6;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*7*/
	     sprintf (html_code[n], "#FF6600");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.4;
	     ps_code[n][2]=0;
	     n++;
	     
	     /*8*/
	     sprintf (html_code[n], "#FF3300");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0.2;
	     ps_code[n][2]=0;
	     n++;
	     
	     
	     /*9*/
	     sprintf (html_code[n], "#FF2000");
	     ps_code[n][0]=1;
	     ps_code[n][1]=0;
	     ps_code[n][2]=0;
	     n++;
	   }
	 else if (color_file)
	   {
	     fp=vfopen( COLOR_FILE, "r");
	     
	     while ((c=fgetc(fp))!='*');
	     while ((c=fgetc(fp))!='\n');
	     
		 c=0;
		 while ((c=fgetc(fp))!=EOF)
		   {
		     ungetc(c, fp);
		     if ( fscanf (fp, "%d", &classC)==0)break;
		     fscanf (fp, "%s %f %f %f", html_code[classC], &ps_code[classC][0], &ps_code[classC][1],&ps_code[classC][2]);
		     while ((c=fgetc(fp))!='\n' && c!=EOF);
		     if ( c==EOF)ungetc(c, fp);
		   }
		 vfclose(fp);
	   }
	 
	  for (a=0; a<n; a++)
	    {
	      for (b=0; b<3; b++)
		if (ps_code[a][b]>1){rgb255=1;b=3;a=n;}
	    }
	  if (rgb255)
	    for (a=0; a<n; a++)
	      {
		for (b=0; b<3; b++)
		  ps_code[a][b]/=255;
	      }
       }
     
    
     

     /*Conversions*/
     if ( val==10)val--;
       
     r=&C->r;
     g=&C->g;
     b=&C->b;

     if ( val==10)val--;
     sprintf ( C->html_color_class, "value%d",val); 
     
     
     if (val<=9 && val>=0)
       {

	 sprintf ( C->html_color, "%s", html_code[val]);
	 r[0]=ps_code[val][0];
	 g[0]=ps_code[val][1];
	 b[0]=ps_code[val][2];
       }

     else if (val==DEFAULT_COLOR || val==NO_COLOR_RESIDUE || val==NO_COLOR_GAP || (val>'A' && val<'z'))
       {
	C->html_color[0]='\0';
	sprintf ( C->html_color_class, "valuedefault");
	r[0]=1.;
	g[0]=1;
	b[0]=1;
	
	} 
     else if (val==GAP_COLOR)
       {
	C->html_color[0]='\0';
	sprintf ( C->html_color_class, "valuegap");
	r[0]=1.;
	g[0]=1;
	b[0]=1; 
	} 
     else if (val==INK_COLOR )
        {
	sprintf ( C->html_color, "000000");
	sprintf ( C->html_color_class, "valueink");
	r[0]=0.;
	g[0]=0;
	b[0]=0;
	}
     return;

    
     }

// void output_tm_mark(FILE_format *fps)
// {
//    Color *box_c=vcalloc ( 1, sizeof (Color));
//    Color *ink;
//    get_rgb_values_format (INK_COLOR,     (ink  =vcalloc ( 1, sizeof (Color))));
// 
//    get_rgb_values_format ( 5, box_c);
//    fps=print_format_char ( " IN ", box_c,ink,fps);
// 
//    get_rgb_values_format ( 9, box_c);
//    fps=print_format_char ( " HEL ", box_c,ink,fps);
// 
//    get_rgb_values_format ( 0, box_c);
//    fps=print_format_char ( "OUT", box_c,ink,fps);
// }

int output_color_format ( Alignment *B,Alignment *Sin,char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *))
    {
    int a, b, c;
    int max_name_len=15;
    int max_len=0;
    char *buf2, *buf3;
    
    static char *buf;
    int s;
    int *n_residues;
    static FILE_format *fps;
    Color *ink;
    Color *box_c;
    Color *white;
    Alignment *S;

    
    S=copy_aln (B, NULL);
    
    buf2=(char*)vcalloc (Sin->len_aln+1, sizeof (char));
    buf3=(char*)vcalloc (  B->len_aln+1, sizeof (char));
    for ( a=0; a<B->nseq; a++)
      {
	int i,n, r;
	
	i=name_is_in_list ( B->name[a], Sin->name, Sin->nseq, -1);
	if (i==-1)continue;
	sprintf (buf2, "%s", Sin->seq_al[i]);ungap(buf2);
	sprintf (buf3, "%s", S->seq_al[a]);ungap(buf3);
	
	
	
	if (strlen (buf2) !=strlen(buf3))
	  {
	    for (b=0; b<strlen (buf2); b++)HERE ("%d", buf2[b]);
	    fprintf ( stderr, "\nERROR: Incompatible cache ON sEQ: %s\n", S->name[a]);
	    fprintf ( stderr, "\n%s\n%s", buf2, buf3); 
	    fprintf ( stderr, "\n\n%s\n%s", Sin->seq_al[i],S->seq_al[a]); exit (EXIT_FAILURE);
	  }
	
	for (n=0,b=0;b<B->len_aln; b++)
	  {
	    r=S->seq_al[a][b];
	    if (!is_gap(r))
	      {
		S->seq_al[a][b]=buf2[n++];
	      }
	  }
      }

    S=aln2number(S);
    vfree (buf2);
    
    box_c=(Color*)vcalloc ( 1, sizeof (Color));    
    get_rgb_values_format (DEFAULT_COLOR, (white=(Color*)vcalloc ( 1, sizeof (Color))));	
    get_rgb_values_format (INK_COLOR,     (ink  =(Color*)vcalloc ( 1, sizeof (Color))));

    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=B->order[a][1];

    fps=vfopen_format( name);      
    if ( buf==NULL)
	{
	buf=(char*)vcalloc (10000, sizeof (int));
	}

    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	}	
    if ( max_len>max_name_len)max_len=max_name_len;
    
   sprintf (buf, "\n%s, %s (%s)\n%s\n",PROGRAM,VERSION,BUILD_INFO, AUTHOR);
   fps=print_format_string ( buf,white, ink, fps);

   if(B->output_tm == 1)
   {
	fps=print_format_char ( '\n', white, ink, fps);
	get_rgb_values_format ( 5, box_c);
	fps=print_format_string ( " IN ", box_c,ink,fps);
	
	get_rgb_values_format ( 9, box_c);
	fps=print_format_string ( " HEL ", box_c,ink,fps);
	
	get_rgb_values_format ( 0, box_c);
	fps=print_format_string ( " OUT ", box_c,ink,fps);
   }

   fps=print_format_string ( "\n\n",white,ink, fps);
   
   fps->line-=max_len;
   fps->line=fps->line-fps->line%3;
   
  
   
   
   for (a=0; a<B->len_aln; a+=fps->line)
	   {	   
	   
	   if ( (fps->n_line+(B->nseq+4))>fps->max_line_ppage && !((B->nseq+4)>fps->max_line_ppage))
		 {
		 fps=print_format_char ( fps->eop,white, ink, fps);		
		 }
	  
	   for (b=0; b<=B->nseq; b++)
	     {
	       sprintf (buf,"%-*.*s ",max_len+2, max_len,(b==B->nseq)?"":S->name[b]);
	       fps=print_format_string ( buf,white, ink, fps);
	       if(B->output_res_num)
		 {
		   sprintf (buf, " %4d ", n_residues[b]+1);
		   fps=print_format_string ( buf,white, ink, fps);
		 }
	       
	       for (fps->in_seq=1,c=a;c<a+fps->line && c<B->len_aln;c++)
		 {
		   if (b==B->nseq)
		     {
		       n_residues[b]++;
		       get_rgb_values_format (DEFAULT_COLOR,box_c);
		       s=analyse_aln_column ( B, c);
		     }
		   else
		     {
		       n_residues[b]+=!is_gap(B->seq_al[b][c]);
		       s=B->seq_al[b][c]; 
		       if (!is_gap(s) && S->seq_al[b][c]!=NO_COLOR_RESIDUE )
			 {
			   get_rgb_values_format ( S->seq_al[b][c], box_c);			    			
			 }
		       else
			 {
			   get_rgb_values_format (GAP_COLOR, box_c);		   	
			 }
		     }
		fps=print_format_char ( s,box_c, ink,fps);
		}
	      fps->in_seq=0;

	       if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ", n_residues[b]);
		 fps=print_format_string ( buf,white, ink, fps);
		 }

	      fps=print_format_char ( '\n', white, ink, fps);	      
	     
	     }
	    fps=print_format_string ( "\n\n",white, ink, fps);
	    }
    fps=print_format_string ( "\n\n\n",white, ink,fps);
    
    
    vfclose_format( fps);
    free_aln (S);
    vfree (n_residues);
    return 1;

    }

int output_raw_score (Alignment *A, Alignment *B, char *name)
{
  FILE *fp;
  int a, b,pos;
  int *ngap;
  int *len;
  
  if (!A)return NULL;
  ngap=(int*)vcalloc (A->len_aln, sizeof (int));
  len=(int*)vcalloc  (A->len_aln,sizeof (int)); 
  for (b=0; b<A->len_aln; b++)
    for (a=0; a<A->nseq; a++)
      {
	int g=is_gap(A->seq_al[a][b]);
	ngap[b]+=g;
	len[a]+=1-g;
      }
  
  fp=vfopen (name, "w");
  fprintf ( fp, "#ALN nseq: %d len: %d score: %d\n", A->nseq, A->len_aln, (B)?B->score:0);
  for (b=0; b<A->len_aln; b++)
    fprintf (fp, "#COLUMN pos: %d ngap: %d score: %d\n", b+1,ngap[b],(B)?B->seq_al[B->nseq][b]-'0':0);
  for (a=0; a<A->nseq; a++)
    fprintf ( fp, "#SEQUENCE name: %s len: %d score: %d\n", A->name[a], len[a], (B)?B->score_seq[a]:0);
  for (a=0; a<A->nseq; a++)
    for (pos=0,b=0; b<A->len_aln; b++)
      {
	int r=A->seq_al[a][b];
	int s=(B)?B->seq_al[a][b]:0;
	int g=is_gap(r);
	pos+=1-g;
	if (!g && s!=NO_COLOR_RESIDUE)
	  {
	    int s=B->seq_al[a][b];
	    if (s>='0')s-='0';
	    fprintf (fp, "#RESIDUE seq: %s pos: %d res: %c score: %d\n", A->name[a],pos, r,s);
	  }
	  }
  vfclose (fp);
  vfree (ngap);
  vfree(len);
  return 1;
}



int output_reliability_format ( Alignment *B,Alignment *S,char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *))
    {
    int a, b, c,l;
    int max_name_len=15;
    int max_len=0;
    static char *buf,*buf2;
    int s;
    static FILE_format *fps;
    Color *ink;
    Color *box_c;
    Color *white;
    int *n_residues;
   

    box_c=(Color*)vcalloc ( 1, sizeof (Color));    
    get_rgb_values_format (DEFAULT_COLOR, (white=(Color*)vcalloc ( 1, sizeof (Color))));	
    get_rgb_values_format (INK_COLOR,     (ink  =(Color*)vcalloc ( 1, sizeof (Color))));
    
    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=B->order[a][1];


    fps=vfopen_format( name);      
    if ( buf==NULL)
	{
	buf=(char*)vcalloc (10000, sizeof (int));
	buf2=(char*)vcalloc (10000, sizeof (int));
	}

    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	}

    if ( vfopen_format==vfopen_ascii)
      {
	fps->line+= max_len;
      }
    else if ( max_len>max_name_len)max_len=max_name_len;
    
    
    
    sprintf (buf, "%s, %s (%s)\n%s\nCPU TIME:%d sec.\n%s",PROGRAM,VERSION,BUILD_INFO, AUTHOR,  (B->cpu+get_time())/1000, (S->generic_comment)?S->generic_comment:"");
    fps=print_format_string ( buf,white, ink, fps);
    sprintf (buf, "SCORE=%d\n*\n", S->score_aln);
    fps=print_format_string ( buf,white, ink, fps);
    
   sprintf ( buf2, " BAD AVG GOOD");
   l=strlen(buf2);
   get_rgb_values_format ( DEFAULT_COLOR, box_c);
   fps=print_format_char ( buf2[0],box_c, ink, fps);
   for ( a=1; a<l-1; a++)
        {
	get_rgb_values_format ( MIN(9,a-1), box_c);
	fps=print_format_char ( buf2[a],box_c,ink,fps);
	}
   fps=print_format_char ( buf2[a], box_c, ink, fps);
   
   
   fps=print_format_string ( "\n*\n",white,ink, fps);
   for ( a=0;S->score_seq && a< B->nseq; a++)
       {
	 get_rgb_values_format (S->score_seq[a]/10, box_c);
       sprintf ( buf, "%-*.*s ", max_len+2,max_len, S->name[a]);
       fps=print_format_string ( buf,box_c, ink,fps);
       sprintf ( buf, ": %3d\n", S->score_seq[a]);
       fps=print_format_string ( buf,white, ink,fps);
       }
   //Print the Consensus score
   get_rgb_values_format (S->score_aln/10, box_c);
   sprintf ( buf, "%-*.*s ", max_len+2,max_len, S->name[S->nseq]);
   fps=print_format_string ( buf,box_c, ink,fps);
   sprintf ( buf, ": %3d\n", S->score_aln/10);
   fps=print_format_string ( buf,white, ink,fps);
   
   fps=print_format_string ( "\n",white, ink,fps);


   
   fps->line-=max_len;
   fps->line=fps->line-(fps->line%3);
   
   for (a=0; a<B->len_aln; a+=fps->line)
	   {	   
	   
	   if ( (fps->n_line+(B->nseq+4))>fps->max_line_ppage && !((B->nseq+4)>fps->max_line_ppage))
		 {
		 fps=print_format_char ( fps->eop,white, ink, fps);		
		 }
	  
	   for (b=0; b<=S->nseq; b++)
	     {
	     if ( b==S->nseq && print_format_string !=print_ascii_string) fps=print_format_string ( "\n",white, ink, fps);
	     sprintf (buf,"%-*.*s ",max_len+2,max_len,S->name[b]);
	     fps=print_format_string ( buf,white, ink, fps);
	     if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ", n_residues[b]+1);
		 fps=print_format_string ( buf,white, ink, fps);
		 }
	     
	     for (fps->in_seq=1,c=a;c<a+fps->line && c<B->len_aln;c++)
		{
		if (b==S->nseq)
		   {
		     
		   if (S->score_seq)
		     {
		       int s;
		       s=S->seq_al[b][c];
		       if ( s>='0' && s<='9')s-='0';
		       get_rgb_values_format (s,box_c);
		     }
		   else get_rgb_values_format (DEFAULT_COLOR,box_c);
		   n_residues[b]++;
		   s=analyse_aln_column ( B, c);
		   }
		else
		   {
		   n_residues[b]+=!is_gap(B->seq_al[b][c]);
		   //s=toupper(B->seq_al[b][c]);
		   s=GET_CASE(B->residue_case, B->seq_al[b][c]);
		   
		   if (!is_gap(s) && S->seq_al[b][c]!=NO_COLOR_RESIDUE )
			{
			get_rgb_values_format ( S->seq_al[b][c], box_c);			    			
			
			}
		   else
		        {
			get_rgb_values_format (GAP_COLOR, box_c);		   	
			
			}
		   
		   }
		fps=print_format_char ( s,box_c, ink,fps);
		}
	      fps->in_seq=0;

	      if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ",n_residues[b]);
		 fps=print_format_string ( buf,white, ink, fps);
		 }

	      fps=print_format_char ( '\n', white, ink, fps);	      
	     
	     }
	    fps=print_format_string ( "\n\n",white, ink, fps);
	    }
    fps=print_format_string ( "\n\n\n",white, ink,fps);
    vfclose_format( fps);
    return 1;

    }

int output_reliability_format_fasta ( Alignment *B,Alignment *S,char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *))
    {
    int a, b, c,l;
    int max_name_len=15;
    int max_len=0;
    static char *buf,*buf2;
    int s;
    static FILE_format *fps;
    Color *ink;
    Color *box_c;
    Color *white;
    int *n_residues;
   

    box_c=(Color*)vcalloc ( 1, sizeof (Color));    
    get_rgb_values_format (DEFAULT_COLOR, (white=(Color*)vcalloc ( 1, sizeof (Color))));	
    get_rgb_values_format (INK_COLOR,     (ink  =(Color*)vcalloc ( 1, sizeof (Color))));
    
    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=B->order[a][1];


    fps=vfopen_format( name);      
    if ( buf==NULL)
	{
	buf=(char*)vcalloc (10000, sizeof (int));
	buf2=(char*)vcalloc (10000, sizeof (int));
	}

    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	}

    if ( vfopen_format==vfopen_ascii)
      {
	fps->line+= max_len;
      }
    else if ( max_len>max_name_len)max_len=max_name_len;
   
   fps->line-=max_len;
   fps->line=fps->line-(fps->line%3);
 
	   for (b=0; b<S->nseq; b++)
	     {
	     sprintf (buf,">%-*.*s \n",max_len+2,max_len,S->name[b]);
	     fps=print_format_string ( buf,white, ink, fps);
	     for (fps->in_seq=1,c=0;c<B->len_aln;c++)
		{
		   n_residues[b]+=!is_gap(B->seq_al[b][c]);
		   s=B->seq_al[b][c];
		   if (!is_gap(s) && S->seq_al[b][c]!=NO_COLOR_RESIDUE )
			{
			get_rgb_values_format ( S->seq_al[b][c], box_c);			    	
			}
		   else
		        {
			get_rgb_values_format (GAP_COLOR, box_c);
			}
		   
		fps=print_format_char ( s,box_c, ink,fps);
		}
	      fps->in_seq=0;
	      fps=print_format_char ( '\n', white, ink, fps);	      
	     }
    vfclose_format( fps);
    return 1;
    }
    
/*****************************************************************************/
/*                       PDF         FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_color_pdf    ( Alignment *B,Alignment *S, char *name)
      {
      char *tmp_name;
      char command[LONG_STRING];

      
#ifndef PS2PDF 
      fprintf (stderr, "\nPDF FORMAT IS NOT SUPPORTED: INSTALL THE PROGRAM PS2PDF\n");
      myexit (EXIT_FAILURE);
#else
      tmp_name=vtmpnam(NULL);
      
      output_color_ps (B, S, tmp_name);
      sprintf ( command, "%s %s %s", PS2PDF, tmp_name, name);
      my_system  ( command); 
      vremove  ( tmp_name); 
#endif      
      
      
      return 1;
      }
int       output_reliability_pdf    ( Alignment *B,Alignment *S, char *name)
      {
      char *tmp_name;
      char command[LONG_STRING];

      

#ifndef PS2PDF 
      fprintf (stderr, "\nPDF FORMAT IS NOT SUPPORTED: INSTALL THE PROGRAM PS2PDF\n");
      myexit (EXIT_FAILURE);
#else
      tmp_name=vtmpnam(NULL);
      
      output_reliability_ps (B, S, tmp_name);
      sprintf ( command, "%s %s %s", PS2PDF, tmp_name, name);
      my_system  ( command); 
      vremove  ( tmp_name); 
#endif      
      
      
      return 1;
      }
/*****************************************************************************/
/*                       POST SCRIPT FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_color_ps     ( Alignment *B,Alignment *S, char *name)
      {
      output_color_format (B, S, name, vfopen_ps,print_ps_string,print_ps_char,get_rgb_values_ps, vfclose_ps);
      return 1;
      }
int       output_reliability_ps     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format (B, S, name, vfopen_ps,print_ps_string,print_ps_char,get_rgb_values_ps, vfclose_ps);
      return 1;
      }
FILE_format *print_ps_string( char *s, Color *box, Color *ink, FILE_format *fps)
      {
      int l;
      int a;
      
      l=strlen (s);
      
      for ( a=0; a< l; a++)
          {
	  fps=print_ps_char (s[a], box, ink, fps);
	  }
      return fps;
      }
      

FILE_format * print_ps_char ( int c, Color *box, Color *ink, FILE_format *f)
       {
       
       int ch;
       int cw;

       ch=f->font+3;
       cw=f->font-2;
        
       if ( c=='(' || c==')')return f;
       else if (c!='\n' && c!=f->eop)
          {
	  fprintf(f->fp,"%d %d moveto\n", f->x,f->y);
	  fprintf(f->fp,"0 %d rlineto\n%d 0 rlineto\n0 -%d rlineto\nclosepath\n",ch,cw,ch   );
	  fprintf(f->fp,"%3.1f %3.1f %3.1f setrgbcolor\nfill\n%3.1f %3.1f %3.1f setrgbcolor\n", box->r,box->g,box->b, ink->r, ink->g, ink->b);
	  fprintf(f->fp,"%d %d moveto\n(%c) show\n", f->x+1,f->y+3, c);
	  
	  f->x+=cw;
	  }
       else 
          {
	  f->n_line++;
	  if ( f->n_line==f->max_line_ppage || c==f->eop)
	     {
	     
	     f->n_line=0;
	     f->x=f->x0;
	     f->y=f->y0;
	     fprintf(f->fp,"showpage\n");
	     f->n_pages++;
	     fprintf ( f->fp, "%c%cPage:  %d %d\n",'%', '%', f->n_pages, f->n_pages);	    
	     }
	  else
	     {
	     f->x=f->x0;
	     f->y-=ch;	  
	     }
	  }
       return f;
       }
FILE_format* print_ps_line (int len, Color *c, FILE_format *f)
{
  static int x=0;
  static int y=0;
  int pix=1;
  int xs=0;
  static int width;

  if (!width)
    {
      width=1;
      fprintf (f->fp, "%d setlinewidth\n", pix);
    }
  
  if (len==-1)
    {
      y-=pix;
      x=xs;
    }
  else
    {
      fprintf (f->fp, "newpath\n");
      fprintf (f->fp, "%d %d moveto\n",x, y);
      fprintf(f->fp,"%3.1f %3.1f %3.1f setrgbcolor\n", c->r,c->g,c->b);
      x+=len*pix;
      fprintf (f->fp, "%d %d lineto\n", x, y);
      fprintf (f->fp, "stroke\n");
    }
  return f;
}

void get_rgb_values_ps ( int val, Color *C)
     {
     get_rgb_values ( val, C);
     }



FILE_format* vfopen_ps ( char *name)
      {
      FILE_format*fps;

      fps=(FILE_format*)vcalloc ( 1, sizeof ( FILE_format));
      fps->font=9;
      fps->max_line_ppage=60;
      fps->line=get_msa_line_length (0, 0);/*N char per line*/
      fps->x0=15;
      fps->y0=750;
      fps->eop='^';
      
      fps->fp=vfopen ( name, "w");
      fprintf(fps->fp,"%%!PS-Adobe-2.0\n/Courier findfont\n%d scalefont\nsetfont\n",fps->font);
      fprintf(fps->fp, "%%%%Pages: (atend)\n");
      fprintf(fps->fp,"newpath\n"); 
      ++(fps->n_pages);
      fprintf (fps->fp, "%%%%Page:  %d %d\n", fps->n_pages, fps->n_pages);
      fprintf (fps->fp,"%d %d translate\n",fps->x0, fps->y0);
      return fps;
      }

FILE_format* vfclose_ps ( FILE_format *fps)
      {
      
      fprintf(fps->fp,"showpage\n");	  
      fprintf ( fps->fp, "%%%%Pages:  %d\n", fps->n_pages);
      fprintf(fps->fp,"%%%%EOF"); 
      fprintf(fps->fp,"%%%%\n"); 
      vfclose ( fps->fp);
      vfree (fps);
      return NULL;
      }
/*****************************************************************************/
/*                       HTML FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
//JM_ADD
/*****************************************************************************/
void output_hit_matrix(char *fileName, float **ffpHitScoreMatrix, int nl)
{
	int i, j;
     	FILE *fp;

     	fp=vfopen(fileName, "w");
 	for(i = 0; i < nl; i++)
	{
		for(j = 0; j < i; j++)
			fprintf(fp, "%6.2f ", ffpHitScoreMatrix[j][i-j]);
		for(j = i; j < nl; j++)
			fprintf(fp, "%6.2f ", ffpHitScoreMatrix[i][j-i]);
		fprintf(fp, "\n");
	}
	vfclose(fp);
}
int output_hit_color_format (Alignment *B, float **ffPScoreTable, int nl, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));
int output_hit_color_html (Alignment *B, float **ffPScoreTable, int nl, char *name)
{
   output_hit_color_format (B, ffPScoreTable, nl, name, vfopen_html,print_html_string,print_html_char,get_rgb_values_html, vfclose_html);
   return 1;
}

int output_hit_color_format (Alignment *B, float **ffPScoreTable, int nl, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *))
{
	int a, b;
	int max_name_len=15;
	int max_len=0;
    
	static char *buf;
	static FILE_format *fps;
	Color *ink;
	Color *box_c;
	Color *white;

	box_c=(Color*)vcalloc ( 1, sizeof (Color));    
	get_rgb_values_format (DEFAULT_COLOR, (white=(Color*)vcalloc ( 1, sizeof (Color))));	
	get_rgb_values_format (INK_COLOR,     (ink  =(Color*)vcalloc ( 1, sizeof (Color))));

	if ( max_len==0)
	{
		for ( a=0; a< B->nseq; a++)
		{
			if ( strlen (B->name[a])>max_len)
				max_len= strlen ( (B->name[a]));
		}
	}
	if ( max_len>max_name_len)max_len=max_name_len;

     	if ( buf==NULL)
 		buf=(char*)vcalloc (10000, sizeof (int));
	int iEmptyChr = 32; //SPACE ASCIICODE
	int iColorValue;
	fps=vfopen_format(name);
   	for (a=0; a < nl; a++)
	{	   
	     	sprintf (buf,"%*d ", max_len+2, a);
	     	fps=print_format_string ( buf,white, ink, fps);
		for(b = 0; b < a; b++)
		{
		      iColorValue = (int)((ffPScoreTable[b][a-b]*9)/100);
		      get_rgb_values_format (iColorValue, box_c);
		      fps=print_format_char (iEmptyChr,box_c, ink,fps);
		}
		for(b = a; b < nl; b++)
		{
		      iColorValue = (int)((ffPScoreTable[a][b-a]*9)/100);
		      get_rgb_values_format (iColorValue, box_c);
		      fps=print_format_char (iEmptyChr,box_c, ink,fps);
		}
	      fps=print_format_char ('\n', white, ink, fps);
	}
	vfclose_format(fps);
	vfree(buf);
	vfree(box_c);
	return 1;
}

/*****************************************************************************/

int       output_color_html     ( Alignment *B,Alignment *S, char *name)
      {
      output_color_format (B, S, name, vfopen_html,print_html_string,print_html_char,get_rgb_values_html, vfclose_html);
      return 1;
      }
int       output_reliability_html     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format (B, S, name, vfopen_html,print_html_string,print_html_char,get_rgb_values_html, vfclose_html);
      return 1;
      }
FILE_format *print_html_string( char *s, Color *box, Color *ink, FILE_format *fhtml)
      {
      int l;
      int a;
      
      l=strlen (s);
      
      for ( a=0; a< l; a++)
          {
	  fhtml=print_html_char (s[a], box, ink, fhtml);
	  }
      fhtml=print_html_char (CLOSE_HTML_SPAN,NULL,NULL,fhtml);
      return fhtml;
      }
      

FILE_format * print_html_char ( int c, Color *box, Color *ink, FILE_format *f)
       {
       char html_color[100];
       int in_span, new_color;
       char string[1000];


        if (c==CLOSE_HTML_SPAN)
	 {
	   if (f->in_html_span)fprintf ( f->fp, "</span>");
	   f->in_html_span=0;
	   return f;
	 }
	
	
       in_span=f->in_html_span;
       new_color=1-(strm (box->html_color_class, f->previous_html_color));
     
       
 
       sprintf (f->previous_html_color, "%s", box->html_color_class);
       sprintf ( html_color, "class=%s", box->html_color_class);

       
       if ( c!=' ')sprintf ( string, "%c", c);
       else sprintf ( string, "&nbsp;");
       
       if ( !in_span &&                  c!='\n' && c!=f->eop)
          {
	  fprintf ( f->fp, "<span %s>%s",html_color,string );
	  f->in_html_span=1;
	  }
       else if (in_span && !new_color && c!='\n' && c!=f->eop)
	  {
	 
	  fprintf ( f->fp, "%s",string);
	  }
       else if (in_span &&  new_color && c!='\n' && c!=f->eop)
	  {
	  fprintf ( f->fp, "</span><span %s>%s",html_color,string);
	  }	   
       else if ( c=='\n')
          {
	  if ( f->in_html_span)fprintf ( f->fp, "</span>");
	  fprintf ( f->fp, "<br>");
	  sprintf ( f->previous_html_color, "no_color_set");
	  f->in_html_span=0;
	  f->n_line++;
	  }
       
       
       
       
     
       return f;
       }
      
void get_rgb_values_html ( int val, Color *C)
     {
     get_rgb_values ( val, C);
     }

FILE_format* vfopen_html ( char *name)
      {
      FILE_format*fhtml;
      Color *color;
      int a;

      color=(Color*)vcalloc ( 1, sizeof (Color));

      fhtml=(FILE_format*)vcalloc ( 1, sizeof ( FILE_format));
      fhtml->font=11;
      fhtml->max_line_ppage=100000;
      fhtml->line=get_msa_line_length (0, 0);/*N char per line*/
      fhtml->x0=15;
      fhtml->y0=800;
      fhtml->eop='^';
      sprintf ( fhtml->previous_html_color, "no_value_set");
      fhtml->fp=vfopen ( name, "w");

      fprintf(fhtml->fp,"<html>\n<style>\n");
     
      fprintf(fhtml->fp,"SPAN { font-family: courier new, courier-new, courier, monospace; font-weight: bold; font-size: %dpt;}\n", fhtml->font);
      fprintf(fhtml->fp,"SPAN { line-height:100%%}\n");
      fprintf(fhtml->fp,"SPAN {	white-space: pre}\n");
      
      for ( a=0; a< 10; a++)
          {
	  get_rgb_values_html ( a, color);    
	  if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
	  else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
	  }
      get_rgb_values_html (DEFAULT_COLOR, color);
      if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
      else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
      
      get_rgb_values_html (GAP_COLOR, color);
      if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
      else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
       
      get_rgb_values_html (INK_COLOR, color);
      if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
      else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
      
      
      
      fprintf(fhtml->fp,"</style>");
      fprintf(fhtml->fp,"<body>");
      
      return fhtml;
      }
FILE_format* vfclose_html ( FILE_format *fhtml)
      {
      if ( fhtml->in_html_span)fprintf(fhtml->fp,"</span>");
      fprintf(fhtml->fp,"</body></html>\n");	  
      vfclose ( fhtml->fp);
      vfree (fhtml);
      return NULL;
      }
/*****************************************************************************/
/*                       ascii FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_color_ascii     ( Alignment *B,Alignment *S, char *name)
      {

	output_color_format (B, S, name, vfopen_ascii,print_ascii_string,print_ascii_char,get_rgb_values_ascii, vfclose_ascii);
	return 1;
      }
int       output_reliability_ascii     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format (B, S, name, vfopen_ascii,print_ascii_string,print_ascii_char,get_rgb_values_ascii, vfclose_ascii);
      return 1;
      }

int       output_reliability_fasta     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format_fasta (B, S, name, vfopen_ascii,print_ascii_string,print_ascii_char,get_rgb_values_ascii, vfclose_ascii);
      return 1;
      }
      
FILE_format *print_ascii_string( char *s, Color *box, Color *ink, FILE_format *fascii)
      {
      int l;
      int a;
      
      l=strlen (s);
      
      for ( a=0; a< l; a++)
          {
	  fascii=print_ascii_char (s[a], box, ink, fascii);
	  }
      return fascii;
      }
      

FILE_format * print_ascii_char ( int c, Color *box, Color *ink, FILE_format *f)
       {             
       if (box->ascii_value>=0 && f->in_seq)fprintf ( f->fp, "%c", box->ascii_value);
       else fprintf ( f->fp, "%c",c);
       return f;
       }

      
void get_rgb_values_ascii ( int val, Color *C)
{

     
     if ( val==NO_COLOR_RESIDUE)C->ascii_value='-';
     else if ( val==NO_COLOR_GAP)C->ascii_value='*';
     else if ( val>9){C->ascii_value='#'; }
     else if ( val>=0 && val<=9) C->ascii_value=val+'0';
     else   C->ascii_value=val;
     }

FILE_format* vfopen_ascii ( char *name)
      {
      FILE_format*fascii;
      fascii=(FILE_format*)vcalloc ( 1, sizeof ( FILE_format));
      fascii->font=11;
      fascii->max_line_ppage=100000;
      fascii->line=get_msa_line_length (0,0);/*N char per line*/
      fascii->x0=15;
      fascii->y0=800;
      fascii->eop='^';
      fascii->fp=vfopen ( name, "w");
     
      
      return fascii;
      }
FILE_format* vfclose_ascii ( FILE_format *fascii)
      {
       vfclose ( fascii->fp);
      vfree (fascii);
      return NULL;
      }


/*****************************************************************************/
/*                       seq_score output                                    */
/*                                                                           */
/*****************************************************************************/

int       output_seq_reliability_ascii     ( Alignment *B,Alignment *S, char *name)
{
  FILE *fp;
  int a;
  int max_len=0;
  for ( a=0; a< B->nseq; a++)
    {if ( strlen (B->name[a])>max_len)
      max_len= strlen ( (B->name[a]));
    }
  

  fp=vfopen ( name, "w");
  fprintf ( fp, "ALN_SCORE %d\n", S->score_aln);
  for ( a=0; a< S->nseq; a++)fprintf (fp, "SEQ_SCORE %*.*s %3d\n",  max_len+2,max_len,S->name[a],S->score_seq[a]);
  vfclose (fp);
  
  return 1;
}
  
int aln2compressed_ps (Alignment *A,char *file)
{
  FILE_format *fps;
  int a, b, cs, cl, ns;
  
  static Color *C=(Color*)vcalloc ( 1, sizeof (Color));
  fps=vfopen_ps(file);
  for (a=0; a<A->nseq; a++)
    {
      for (cs=-1,cl=0,b=0; b<A->len_aln; b++)
	{
	  ns=A->seq_al[a][b]-'0';
	  if (cs==-1){cs=ns;cl=1;}
	  else if (ns==cs){cl++;}
	  else 
	    {
	      
	      get_rgb_values_ps(cs,C);
	      fps=print_ps_line (cl,C, fps);
	      cl=1;
	      cs=ns;
	    }
	}
      get_rgb_values_ps(cs,C);
      fps=print_ps_line (cl, C, fps);
      fps=print_ps_line (-1, NULL, fps);
    }
  vfclose_ps (fps);
}
int aln2compressed_pdf    ( Alignment *A,char *name)
      {
      char *tmp_name;
      char command[LONG_STRING];

      
#ifndef PS2PDF 
      fprintf (stderr, "\nPDF FORMAT IS NOT SUPPORTED: INSTALL THE PROGRAM PS2PDF\n");
      myexit (EXIT_FAILURE);
#else
      tmp_name=vtmpnam(NULL);
      
      aln2compressed_ps (A, tmp_name);
      sprintf ( command, "%s %s %s", PS2PDF, tmp_name, name);
      my_system  ( command); 
      vremove  ( tmp_name); 
#endif      
      
      
      return 1;
      }

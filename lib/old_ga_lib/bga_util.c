/*SAGA bga_util.c*/
#define MAIN 1
#include "ga.h"
/******************************************************************************************/			    	
Decoded_chromosome *decode_chromosome (Chromosome *C, Decoded_chromosome *D, Parameter *PARAM)
    {
    if ( C==NULL)crash ("Empty Chromosome [decode_chromosome\n");
    else
	{
	D=realloc_decoded_chromosome (D,C, NULL, PARAM);
	}
    return decode (C,D, PARAM);
    }

Chromosome *code_chromosome ( Decoded_chromosome *D, Chromosome *C, Parameter *PARAM)
    {
    if ( D==NULL)crash ("Empty Decoded_chromosome [code_chromosome]\n");
    else
	C=realloc_chromosome (C,NULL,D, PARAM);
    
    return code (D,C, PARAM);
    }
Decoded_chromosome *copy_decoded_chromosome ( Decoded_chromosome *D1, Decoded_chromosome *D2, Parameter *PARAM)
    {
    static Chromosome *C;

    C=code_chromosome ( D1, C, PARAM);
    return decode_chromosome ( C,D2, PARAM);
    }
Chromosome * copy_chromosome ( Chromosome *C1, Chromosome *C2, Parameter *PARAM)
    {
    static Decoded_chromosome *D;
    D=decode_chromosome ( C1, D, PARAM);
    return code_chromosome (D, C2, PARAM);
    }
Chromosome *extract_chrom ( Population *POP, int G, int i, Chromosome *C, Parameter *PARAM)
    {
    return copy_chromosome (((POP+G)->CHROMOSOME[i]),C, PARAM); 
    }

void insert_chrom  ( Chromosome *A, Population *POP, int G, int i, Parameter *PARAM)
    {
    copy_chromosome (A,((POP+G)->CHROMOSOME[i]), PARAM); 
    (POP+G)->raw_fitness[i][0]=A->score;
    }


void copy_individual ( Population *POP, int G0, int i0,int G1, int i1, Parameter *PARAM)
    {
    static Chromosome * C;
    
    if (C==NULL)
	C=declare_chromosome (PARAM);
	
    extract_chrom ( POP, G0,i0,C,PARAM);
    insert_chrom ( C,POP, G1, i1, PARAM);
    
    (POP+G1)->raw_fitness[i1][0]= (POP+G0)->raw_fitness[i0][0];
    (POP+G1)->raw_fitness[i1][1]= (POP+G0)->raw_fitness[i0][1];
    
    (POP+G1)->fitness[i1][0]= (POP+G0)->fitness[i0][0];
    (POP+G1)->fitness[i1][1]= (POP+G0)->fitness[i0][1];
    
    (POP+G1)->threshold[i1][0]= (POP+G0)->threshold[i0][0];
    (POP+G1)->threshold[i1][1]= (POP+G0)->threshold[i0][1];
    
    (POP+G1)->scaled_fitness[i1][0]= (POP+G0)->scaled_fitness[i0][0];
    (POP+G1)->scaled_fitness[i1][1]= (POP+G0)->scaled_fitness[i0][1];
    
    (POP+G1)->expected_os[i1][0]= (POP+G0)->expected_os[i0][0];
    (POP+G1)->expected_os[i1][1]= (POP+G0)->expected_os[i0][1];
    
    (POP+G1)->ranked[i1][0]= (POP+G0)->ranked[i0][0];
    (POP+G1)->ranked[i1][1]= (POP+G0)->ranked[i0][1];
    }
	
int is_dup (Chromosome *A, Population *L_POP, int **m_t, int tot_g1, Parameter *PARAM)
    {
    /* return 1 if the chrom is duplicated, return 0 if it is not*/
    int a;
    static int *check_array;
    
    if ( check_array==NULL)
	{
	 check_array=vmalloc ( sizeof ( int)*(PARAM->BGA)->MAXPOP);
	 }

    for ( a=0; a< (PARAM->BGA)->MAXPOP; a++)
	check_array[a]=-1;
    
    for ( a=0; a< tot_g1; a++)
	{
	 if ( is_same (A,(L_POP+1)->CHROMOSOME[m_t[a][2]], PARAM))
	    return 1;
	 else
	    check_array[m_t[a][2]]=0;
	} 	  	     
    for ( a=0; a< (PARAM->BGA)->MAXPOP; a++)
	{
	if ( check_array[a]== -1)
	    {
	    if ( is_same (A,(L_POP+0)->CHROMOSOME[a],PARAM))
		return 1;
	    }
	}
     return 0;
     }
     
static int sort_field;
void sort_float ( float **V,int N_F, int F, int left, int right)
	{
	sort_field=F;
	qsort ( V, right+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_int));
	}
int cmp_float ( const float **a, const float **b)
	{
	if ( a[0][sort_field]< b[0][sort_field])return-1;
	else if ( a[0][sort_field]==b[0][sort_field])return 0;
	else return 1;
	}

void sort_int_1D ( int *L, int n)
	{
	int **array;
	int a;
	
	array=declare_int ( n, 1);
	for ( a=0; a< n; a++)
		array[a][0]=L[a];
	sort_int ( array, 1, 0, 0, n-1);
	for ( a=0; a< n; a++)
		L[a]=array[a][0];
	free_int ( array, n);
	}

void sort_int ( int **V,int N_F, int F, int left, int right)
	{
	sort_field=F;
	qsort ( V, (right-left)+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_int));
	}

void sort_int_inv ( int **V,int N_F, int F, int left, int right)
	{
	int a,b,c;
	int **list;
	
	sort_field=F;
	qsort ( V, (right-left)+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_int));
	
	list=vcalloc ( (right-left)+1, sizeof (int*));
	for ( a=left; a< (right-left)+1; a++)
		{
		list[a-left]=vcalloc ( N_F, sizeof (int));
		for ( b=0; b< N_F; b++)
			{
			list[a-left][b]=V[a][b];
			}
		}
	for ( a=left; a< (right-left)+1; a++)
		{
		for ( b=0; b< N_F; b++)
			V[a][b]=list[(right-left)-a][b];
		}
	for ( a=0;a< (right-left)+1; a++)free (list[a]);
	free (list);
	}
	
	



int cmp_int ( const int**a, const int**b)
	{
	if ( a[0][sort_field]< b[0][sort_field])return-1;
	else if ( a[0][sort_field]==b[0][sort_field])return 0;
	else return 1;
	}
	

int ** duplicate_int ( int **array , int len, int field)
    {
    int **array2;
    int a,b;

    array2=declare_int ( len, field);
    
    for ( a=0; a< len; a++)
    	ga_memcpy_int( array[a],array2[a],field);
	
    return array2;
    }

float ** duplicate_float ( float **array , int len, int field)
    {
    float **array2;
    int a,b;

    array2=declare_float ( len, field);
    
    for ( a=0; a< len; a++)
	ga_memcpy_float( array[a],array2[a],field);
    return array2;
    }

void copy_int ( int **array1, int **array2, int len, int number_field)
    {
    int a,b;
    for ( a=0; a< len; a++)
	ga_memcpy_int( array1[a],array2[a],number_field);
    
    }

void copy_float ( float **array1, float **array2, int len, int number_field)
    {
    int a,b;
    for ( a=0; a< len; a++)
	ga_memcpy_float( array1[a],array2[a],number_field);
    
    }

double return_max_double ( double ** array, int len_array, int field)
    {
    double max;
    int a;
    
    max=array[0][field];
    for ( a=1; a< len_array; a++)
	max=( array[a][field]>max)?array[a][field]:max;

    return max;
    }

double return_min_double ( double ** array, int len_array, int field)
    {
    double min;
    int a;
    
    min=array[0][field];
    for ( a=1; a< len_array; a++)
	min=( array[a][field]<min)?array[a][field]:min;

    return min;
    }

float return_max_float ( float ** array, int len_array, int field)
    {
    float max;
    int a;
    
    max=array[0][field];
    for ( a=1; a< len_array; a++)
	max=( array[a][field]>max)?array[a][field]:max;

    return max;
    }

float return_2Dmax_float ( float ** array, int start, int len_array, int first_field, int number_field)
    {
    float max;
    int a,b;
    
    max=array[start][first_field];
    for ( a=start; a< start+len_array; a++)
	for (b=first_field; b< first_field+number_field; b++)     
	    max=( array[a][b]>max)?array[a][b]:max;

    return max;
    }

float return_min_float ( float ** array, int len_array, int field)
    {
    float min;
    int a;
    
    min=array[0][field];
    for ( a=1; a< len_array; a++)
	min=( array[a][field]<min)?array[a][field]:min;

    return min;
    }
    
float return_2Dmin_float ( float ** array, int start, int len_array, int first_field, int number_field)
    {
    float min;
    int a,b;
    
    min=array[start][first_field];
    for ( a=start; a< start+len_array; a++)
	for (b=first_field; b< first_field+number_field; b++)     
	    min=( array[a][b]<min)?array[a][b]:min;

    return min;
    }
int return_max_int ( int ** array, int len_array, int field)
    {
    int max;
    int a;
    
    max=array[0][field];
    for ( a=1; a< len_array; a++)
	max=( array[a][field]>max)?array[a][field]:max;

    return max;
    }

int return_max_int_hor ( short int ** array, int len_array, int field)
    {
    short int max;
    int a;
    
    max=array[field][0];
    for ( a=1; a< len_array; a++)
	max=( array[field][a]>max)?array[field][a]:max;

    return (int)max;
    }
int return_min_int ( int ** array, int len_array, int field)
    {
    int min;
    int a;
    
    min=array[0][field];
    for ( a=1; a< len_array; a++)
	min=( array[a][field]<min)?array[a][field]:min;

    return min;
    }

int return_min_int_hor (short int ** array, int len_array, int field)
    {
    short int min;
    int a;
    
    min=array[field][0];
    for ( a=1; a< len_array; a++)
	min=( array[field][a]<min)?array[field][a]:min;

    return (int)min;	
    }

int min_int ( int a, int b)
    {
    return ( a<b)?a:b;
    }
int max_int ( int a, int b)
    {
    return ( a>b)?a:b;
    }
int return_maxlen ( char ** array, int number)
    {
    int a;
    int max=0;
    for ( a=0; a< number; a++)
	max=( strlen ( array[a])>max)?strlen ( array[a]):max;

    return max;
    }

int return_minlen ( char ** array, int number)
    {
    int a;
    int min;

    min=strlen( array[0]);
    for ( a=1; a< number; a++)
	min=( strlen ( array[a])>min)?strlen ( array[a]):min;

    return min;
    }
 
float return_sum_float ( float **array, int len, int field)
    {
    int a;
    float b=0;
    
    for ( a=0; a< len; a++)
	b+=array[a][field];

    return b;
    }

float return_mean_float ( float **array, int len, int field)
    {
    return ( return_sum_float ( array, len, field)/len);
    }

float return_sd_float ( float **array, int len, int field,float mean)
    {
    int a;
    float b=0;
    
    for ( a=0; a< len; a++)
    	{
    	if ( (mean-array[a][field])!=0) 
	 	b+=(float) ( mean-array[a][field])*(mean-array[a][field]);
	}
	
    return ( (b=(float)sqrt(b/len))<1)?1:b;
    }
float return_mean_diff_float ( float **array, int len, int field,float mean)
    {
    int a;
    float b=0;
    
    for ( a=0; a< len; a++)
    	{
    	if ( (mean-array[a][field])!=0) 
	 	b+=sqrt((double)((float) ( mean-array[a][field])*(float)(mean-array[a][field])));
	}
	
    return ((float)b/(float)len);
    }
int return_sum_int ( int **array, int len, int field)
    {
    int a;
    int b=0;
    
    for ( a=0; a< len; a++)
	b+=array[a][field];

    return b;
    }

int return_mean_int ( int **array, int len, int field)
    {
    return ( return_sum_int ( array, len, field)/len);
    }

int return_sd_int ( int **array, int len, int field,int mean)
    {
    int a;
    int b=0;
    for ( a=0; a< len; a++)
	b+= pow (( mean-array[a][field]),2);

    return ((b=sqrt(b/len))<1)?1:b;
    }
    
void swap_int ( int *a, int *b, int n)
    {
    int t;
    int c;

    for ( c=0;c<n;c++) 
	{t=a[c];
	 a[c]=b[c];
	 b[c]=t;
	}
    }

void swap_float ( float *a, float *b, int n)
    {
    float t;
    int c;
    
    for ( c=0; c< n; c++)
	{
	t=a[c];
	a[c]=b[c];
	b[c]=t;
	}	
    }

void swap_double ( double *a, double *b, int n)
    {
    double t;
    int c;

    for (c=0; c< n; c++)
	{
        t=a[c];
	a[c]=b[c];
	b[c]=t;
	}
    }
	
int float_flip ( float proba)
    {
    int p;
    int a;
    
    
    proba= proba*100000;
    p=proba;
    a= addrand ( (unsigned long) 100000)+1;
    return ( a<=p)?1:0;
    }

char * extract_char ( char * array, int first, int len)
    {
    char *array2;
    int a;

    len= ( len<0)?0:len;    

    array2=vmalloc ( sizeof ( char)*(len+1));
    
    for ( a=0; a<len; a++)
	array2[a]=array[first++];

    array2[a]='\0';

    return array2;
    }    
    

void inverse_int ( int**array, int len, int field, int max, int min)
    {
    int a;
    for ( a=0; a< len; a++)
	array[a][field]=max-array[a][field]+min;
    }
void inverse_float ( float**array, int len, int field, int max, int min)
    {
    int a;
    for ( a=0; a< len; a++)
	array[a][field]=max-array[a][field]+min;
    }
void inverse_2D_float ( float **array, int start, int len, int start_field, int number_field, float max,float min)
    {
    int a, b;
    for ( a=0; a< start+len; a++)
	for ( b=start_field; b< start_field+ number_field; b++)
	    array[a][b]=max-array[a][b]+min;
    }

void reasses_local_bias (Population * POP , int G, Parameter *PARAM, int stab, int gen)
    {
    return;
    } 

char * generate_void ( int x)
    {
    int a;
    char *string;
    
    string = vcalloc ( x+1, sizeof ( char));

    for ( a=0; a< x; a++)
	string[a]=' ';

    string[a]='\0';
    return string;
    } 

FILE * vfopen ( char *name, char *mode)
    {
    FILE *fp;
    char new_name[1000];
    char c;
    int a=0;
    int n_check;
    
    if ( FILE_CHECK==1)
    	{
    	if ( (fp= fopen ( name, mode))==0)
    		{
    		if ( strcmp (mode, "r")==0)
    			{
    			printf ( "\nCOULD NOT READ %s\nENTER A NEW NAME FOR THE FILE: ", name);
    			a=0;	
    			while ( ( c=getchar())!='\n')
    				new_name[a++]=c;
    			new_name[a]='\0';
    			return vfopen ( new_name, mode);
    			}
    		else
    			{
    			printf ( "\nCANNOT WRITE %s\n", name);
    			return vfopen ( new_name, mode);	
    			}
    		}
     	else
     		{
     		return fp;
     		}
	}
     else
     	{
     	n_check=0;
     	while( ((fp= fopen ( name, mode))==0) && n_check<20)	
     		{
     	
		system ( "sleep 1");
     	
     		n_check++;
     		}
     	if ( fp==NULL)
     		{			 
     		printf ( "\nCOULD NOT %s %s\nFORCED EXIT (NON INTERACTIVE MODE)\n", (strcmp ( mode, "r")==0)?"READ":"WRITE", name);
     		crash ( "FORCED EXIT");
     		}
     		
     	else
     		return fp;
     	}					 	
    }

void input_parameter ( char *name1, char *name_type, int interactive, Parameter *PARAM)
	{
	char *string3;
	
	if ( interactive==0)
		return;
	else
		{
		printf ( "\nPRESENT NAME OF %s:\n##%s##", name_type, name1);
		printf ( "\nENTER A NEW NAME (RETURN TO KEEP THE PREVIOUS): ");
		if ( (string3=input_name())!=NULL)
    			{
    			sprintf ( name1, "%s", string3);
    			free ( string3);
    			}	
		}
	}		     
char *input_name ()
	{
	char *string;
	int a;
	char ch;
	
	string= vcalloc ( 500, sizeof ( char));
	
	a=0;
    	while ( ( ch=getchar())!='\n')
    			string[a++]=ch;
    	string[a]='\0'; 
    	
    	if ( string[0]=='\0')
    		{
    		free (string);
    		return NULL;
    		}
    	else
    		return string;
    	}			   
    
void set_name ( char *full_name, char *path, char *fname)
	{
	if ( is_full_name ( fname)==1)
		sprintf ( full_name, "%s", fname);
	else sprintf ( full_name, "%s/%s", path,fname);
	
	}
int is_full_name ( char *fname)
	{
	int len;
	int a=0;
	
	len=strlen( fname);
	for ( a=0; a< len; a++)
		if (fname[a]=='[' || fname[a]=='/')
			return 1;
			
	return 0;
	}		
							    
void delete_file ( char *fname)
	{
	char command[1000];
	FILE * fp;
	
	fp=fopen ( fname, "w");
	fprintf ( fp, "x");
    	fclose ( fp);
	

	sprintf ( command, "rm %s", fname);
	system ( command);
    	
    	}
Decoded_chromosome *get_mem_free_decoded_chrom ( Parameter *PARAM)
	{
	int a;
	
	if (PARAM->MEM==NULL)
		{
		PARAM->MEM= vcalloc ( 1, sizeof ( Mem));
		(PARAM->MEM)->C_BUF=vcalloc ( 100, sizeof (Decoded_chromosome *));
		}		 
	for ( a=0; a< (PARAM->MEM)->NC; a++)
		{
		if ( ((PARAM->MEM)->C_BUF[a])->used==0)
			{
			((PARAM->MEM)->C_BUF[a])->used=1;
			((PARAM->MEM)->C_BUF[a])->num=a;
			return ((PARAM->MEM)->C_BUF[a]);
			}
		}
		
	(PARAM->MEM)->C_BUF[(PARAM->MEM)->NC]=declare_decoded_chromosome (PARAM);
	((PARAM->MEM)->C_BUF[(PARAM->MEM)->NC])->used=1;
	((PARAM->MEM)->C_BUF[(PARAM->MEM)->NC])->num=(PARAM->MEM)->NC;
	(PARAM->MEM)->NC++;
	if ( (PARAM->MEM)->NC>=100)
		{
		printf ( "\nMEMORY LEAK: TOO MANY DEC CHROM ALLOCATED, FORCED EXIT");
		fprintf ((PARAM->BGA)->FP_ERROR, "\nMEMORY LEAK: TOO MANY DEC CHROM ALLOCATED, FORCED EXIT");
		crash ( "FORCED EXIT");
		}
		
	return ((PARAM->MEM)->C_BUF[(PARAM->MEM)->NC-1]);
	}
void free_mem_dc ( Decoded_chromosome *D, Parameter *PARAM)
	{
	int a;
	clean_decoded_chromosome ( D,PARAM);
	D->used=0;
	D->num=0;
	}
			
	
void ga_memcpy_int ( int *array1, int *array2, int n)
	{
	int a;
	
	for ( a=0; a< n; a++)
		array2[a]=array1[a];
	}			
			
void ga_memcpy_float ( float *array1, float *array2, int n)
	{
	int a;
	
	for ( a=0; a< n; a++)
		array2[a]=array1[a];
	}			
		
void crash(char *string)
	{
	int a;
	int b;

#ifndef CRASH_ON_EXIT
	fprintf ( stderr, "\n%s\n", string);
	exit (0);
#else	
	a=0;
	b=5/a;

	fprintf ( stderr, "\n%s", string);
	fprintf ( stderr, "%d", b);
#endif
	}

void crash2(int a, int b)
	{
	char string[100];
	
	if ( a>b)
		{
		sprintf ( string, "\n%d %d", a, b);
		crash ("  ");
		}
	}

void identify_operating_system ( Parameter *PARAM)
	{
	char command[1000];
	FILE *fp;
	char os[100];
	char fname[100];
	int a;
	
	a=addrand((unsigned long)10000);
	
	sprintf ( fname, "os_%d", a);
	sprintf ( command, "echo $OSTYPE>%s", fname);
	system ( command);
	fp=vfopen ( fname, "r");
	fscanf ( fp, "%s", os);
	fclose ( fp);
	sprintf ( command , "rm %s", fname);
	system ( command);
	
	if ( (strcmp (os, "irix")==0) || (strcmp (os, "IRIX")==0)) 
		sprintf ( PARAM->OS, "_sgi");
	else if ( (strcmp (os, "osf1")==0) || (strcmp (os, "OSF1")==0)) 
		sprintf ( PARAM->OS, "_osf");
	else 
		{
		crash("COULD NOT GET OS");
		}
	fprintf ( stderr, "\nOPERATING_SYSTEM: %s\n",PARAM->OS); 
	}	
		
void setenv_func ( char *string_name, char *string_value)
	{
	char command[1000];
	sprintf (command, "setenv %s %s", string_name,string_value);
	
	system ( command);
	}
void get_pwd ( char *name)
	{
	char *string;
	char command[1000];
	FILE *fp;
	 
	string=vcalloc ( L_tmpnam, sizeof (char));
	string=tmpnam(string);
	sprintf ( command, "pwd > %s", string);
	system (command);
	fp=fopen ( string, "r");
	fscanf ( fp, "%s",name);
	fclose (fp);
	sprintf ( command, "rm %s", string);
	system ( command);
	free (string);
	
	}
int int_fabs ( int k)
	{
	if (k>=0)	    
		return k;
 	else if ( k<0)
		return (k* -1);
 	}	 
int name_is_in_list ( char *name, char **name_list, int n_name, int len)
	{
	int a;
	for ( a=0; a< n_name; a++)
		if ( len!=-1)
			if ( strncmp ( name, name_list[a], len)==0)
				return a;
		else if ( strcmp ( name, name_list[a])==0)
				return a;
	return -1;
	}

FILE *get_number_list_in_file ( FILE *fp, int *list, int *n, int *max_len)
	{
	
	int c;
	
	while ( isspace((c=fgetc (fp))));
	ungetc(c, fp);
	while ( c!='\n')
		{
		while ( isspace((c=fgetc (fp))) && c!='\n');
		
		if ( c!='\n')
			{
			ungetc(c, fp);
			if ( n[0]>=max_len[0])
				list=realloc ( list, (n[0]+100)*sizeof (int));
				max_len[0]=(n[0]+100);
			
			fscanf ( fp, "%d",&list[n[0]++]);
			}
		}
	return fp;
	} 
FILE * find_token_in_file ( char *fname, FILE * fp, char *token)
	{
	
	static char *name;
	
	
	
	if ( name==NULL)name = vcalloc ( 1000, sizeof (char));
	if (fp==NULL)
		{
		fp=fopen ( fname, "r");
		
		}
	if ( fp==NULL) return NULL;
	
	while ( (fscanf ( fp, "%s", name))!=EOF)
		{
		
		if (strcmp ( name, token)==0)return fp;
		}
	fclose ( fp);
	return NULL;
	}
	
int check_file_exists ( char *fname)
	{
	FILE *fp;
	
	fp=fopen ( fname, "r");
	if ( fp==NULL) return 0;
	else
		fclose (fp);
	return 1;
	}	
	
	

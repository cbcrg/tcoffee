/*SAGA bga_util.c*/
#define MAIN 1
#include "ga.h"
/******************************************************************************************/			    	
Decoded_chromosome *decode_chromosome (Chromosome *C, Decoded_chromosome *D, Parameter *PARAM)
    {
    if ( C==NULL)crash ("Empty Chromosome [decode_chromosome\n");
    else
	{
	D=realloc_decoded_chromosome (D,C, PARAM);
	}
    return decode (C,D, PARAM);
    }

Chromosome *code_chromosome ( Decoded_chromosome *D, Chromosome *C, Parameter *PARAM)
    {
    if ( D==NULL)crash ("Empty Decoded_chromosome [code_chromosome]\n");
    else
	C=realloc_chromosome (C,D, PARAM);
    
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
     
 








int float_flip ( float proba)
    {
    int p;
    int a;
    
    
    proba= proba*100000;
    p=proba;
    a= addrand ( (unsigned long) 100000)+1;
    return ( a<=p)?1:0;
    }





void reasses_local_bias (Population * POP , int G, Parameter *PARAM, int stab, int gen)
    {
    return;
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
		


int int_fabs ( int k)
	{
	if (k>=0)	    
		return k;
 	else if ( k<0)
		return (k* -1);
 	}	 



	

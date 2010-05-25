/*SAGA bga_init*/ 

# include "ga.h"


int INTERACTIVE_PROMPT_SET;

main ( int argc, char **argv)
	{
	time_t r;
	int a, b, c, d, e;
	char time_string[100];
	
	Parameter *PARAM;    	    	    	
    	
	PARAM=declare_parameter();
	get_ref_time (PARAM);
	process_command_line ( argc,argv, PARAM);
	read_command_line    ( argc,argv, PARAM);
	
	fprintf ( stderr, "\nFINISHED READING PARAMETERS");

	declare_dos_2 ( PARAM);
	
	INTERACTIVE_PROMPT_SET=(PARAM->BGA)->INTERACTIVE;
	
	if ( (PARAM->BGA)->OUTPUT_SCREEN==0)
		{
		(PARAM->std0)=stderr;
		(PARAM->std1)=(PARAM->std2)=fopen ( "/dev/null", "w");
		}
	else if ( (PARAM->BGA)->OUTPUT_SCREEN==1)
		{
		(PARAM->std0)=(PARAM->std1)=stderr;
		(PARAM->std2)=fopen ( "/dev/null", "w");
		}
	else if ((PARAM->BGA)->OUTPUT_SCREEN==2)
		{
		(PARAM->std0)=(PARAM->std1)=(PARAM->std2)=stderr;
		}
		
	time(&r);
	sprintf ( time_string, "%d", (int) r);
	PARAM->START=r;
	
	if ( (PARAM->BGA)->RANDOM_SEED >0)
		{addrandinit ( (unsigned long) (PARAM->BGA)->RANDOM_SEED);
	 	if ( ((PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "\nRANDOM_SEED= %d",(PARAM->BGA)->RANDOM_SEED);
		} 
    	else if ( (PARAM->BGA)->RANDOM_SEED== -1)
		{
		
		b=0;
		for ( a= strlen ( time_string)-1; a>=strlen ( time_string)-3; a--)
	    	time_string[b++]= time_string[a];   
		time_string[b]='\0';
		(PARAM->BGA)->RANDOM_SEED= atoi ( time_string);
		addrandinit ( (unsigned long) (PARAM->BGA)->RANDOM_SEED);
		fprintf ((PARAM->std1), "\nRANDOM_SEED(undet)= %d",(PARAM->BGA)->RANDOM_SEED);
		}
		
	process_param (PARAM);
	input_data (PARAM);
	main_do_aln (PARAM);
	}
	 	
void main_do_aln(Parameter *PARAM)
    {
    Decoded_chromosome *R;
    Chromosome *C;
    char fname[200];
    
    
    if ( (PARAM->BGA)->EVALUATE_ONLY==1)
    	{
    	C=declare_chromosome (PARAM);
    	code_chromosome ((PARAM->BGA)->REFERENCE, C, PARAM);
    	save_individual ( C,NULL, 0, PARAM);
    	}
    else
    	{	 			
    	init_ga ( PARAM);
   	}
    }
    
void init_ga (Parameter *L_PARAM)
    {
    
    Population *L_POP;
    Statistic *L_STAT;

    fprintf ((L_PARAM->std1), "\nDECLARE_POPULATION");
    if ((L_PARAM->MEM)==NULL)
		{L_PARAM->MEM= vcalloc ( 1, sizeof ( Mem));
		(L_PARAM->MEM)->C_BUF=vcalloc ( 100, sizeof (Decoded_chromosome *));
		}
		
    (L_PARAM->MEM)->SD1=declare_decoded_chromosome ( L_PARAM);
    (L_PARAM->MEM)->SD2=declare_decoded_chromosome ( L_PARAM);
    (L_PARAM->MEM)->C1=declare_chromosome ( L_PARAM);
    (L_PARAM->MEM)->C2=declare_chromosome ( L_PARAM);
    
    L_STAT=declare_statistic(L_PARAM); 
    L_PARAM->POPULATION=L_POP=declare_population ((L_PARAM->BGA)->MAXPOP, L_PARAM);
    fprintf ((L_PARAM->std1), "*");    
    fprintf ((L_PARAM->std2), "\n\n###################################################\n\t\tSTART THE GA");
     
    do_ga ( L_POP,L_PARAM, L_STAT);
    }    
         

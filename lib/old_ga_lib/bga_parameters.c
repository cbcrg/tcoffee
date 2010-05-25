#define MAIN 1
#include "ga.h"

static char error [1000];
void process_command_line ( int nargc,char **list, Parameter *PARAM)
	{
	int a;
	int exit_flag=0;

	initialise_parameter_list (PARAM);
	for ( a=1; a< nargc; a++)
		{
		if ( is_cl_flag(list[a]))
			{
			
			if ( strcmp ( list[a], "-l")==0)
				{				
				input_main_parameter (list[++a], PARAM);
				}
			
			else if ( strcmp ( list[a], "-s")==0)
				{
				read_parameter_file (list[++a], PARAM);
				}
			else if ( (strcmp ( list[a], "-help")==0) || (strcmp ( list[a], "-?")==0) || (strcmp ( list[a], "-man")==0)||(strcmp ( list[a], "-help")==0))
				{
				output_help (PARAM, list[0]);
				exit (0);
				}
			
			else if ( strcmp (list[a], "-option")==0) 
				{
				option_help ( PARAM, list[a+1]);
				exit (0);
				}
			else
				{
				sprintf ( error, "\n%s is not a flag or a Parameter ( remove the semicolumn?)\n",list[a]);
				crash (	error);
				}
			}
		else
			{
			a++;
			}
		}
		
	
	
	for ( a=0;a< (PARAM->NP); a++)
		{
		if ( (PARAM->P_LIST+a)->read==0)
			{
			fprintf ( stderr, "\nPB: A PARAMETER VALUE COULD NOT BE READ: %s\n",(PARAM->P_LIST+a)->name); 
			exit_flag=1;
			}
		}
	if ( exit_flag)exit(0);
	}
	
void input_main_parameter  (char *fname, Parameter *PARAM)
	{
	FILE *fp;
	char name1[200], name2[200];
	
	
	
	fp=vfopen (fname, "r");
	while ( fscanf ( fp, "%s %s\n", name1, name2)==2)
		{
		
		read_parameter_file ( name2, PARAM);
		}
	fclose (fp);
	}

		

void initialise_parameter_list (Parameter *PARAM)
	{
	
	int a, b, c;
	
	PARAM->NP=0;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_ERROR_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_ERROR_FILE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_FINISH_FLAG");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_FINISH_FLAG);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_END_FLAG");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_END_FLAG);
	PARAM->NP++;

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_POP_READ_FLAG");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_POP_READ_FLAG);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_POP_DUMP_FLAG");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_POP_DUMP_FLAG);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_GENERATION_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_GENERATION_FILE);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_SCORE_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_SCORE_FILE);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_CPU_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_CPU_FILE);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_GEN_BEST_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_GEN_BEST_FILE);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_BILAN_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_BILAN_FILE);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_CURRENT_OPERATOR_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_CURRENT_OPERATOR_FILE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "RESCALING_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->RESCALING_MODE);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "RESCALING_C");
	(PARAM->P_LIST+PARAM->NP)->TYPE='F';
	(PARAM->P_LIST+PARAM->NP)->F=&((PARAM->BGA)->RESCALING_C);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "PRESELECT_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->PRESELECT_MODE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "SELECT_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->SELECT_MODE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "MAXPOP");
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "SIZE OF THE POPULATION FOR THE GA");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->MAXPOP);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OVERLAPING_GEN_FACTOR");
	(PARAM->P_LIST+PARAM->NP)->TYPE='F';
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "% OF INDIVIDUALS DYING AFTER EACH GENERATION (FLOAT)");
	(PARAM->P_LIST+PARAM->NP)->F=&((PARAM->BGA)->OVERLAPING_GEN_FACTOR);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "MAX_GEN");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "Maximum Number of Generations");
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->MAX_GEN);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "FACTOR_POP");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->FACTOR_POP);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "REPLACE_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->REPLACE_MODE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "RANDOM_SEED");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "Value for the Random Seed: -1 indicates a random value read from the clock \nAny positive value will make the random sequence deterministic");
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->RANDOM_SEED);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "ALLOW_DUPLICATE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->ALLOW_DUPLICATE);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "DYNAMIC_PARAM");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->DYNAMIC_PARAM);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "STABILISED");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->STABILISED);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "CREDIT_TO_REPORT");
	(PARAM->P_LIST+PARAM->NP)->TYPE='F';
	(PARAM->P_LIST+PARAM->NP)->F=&((PARAM->BGA)->CREDIT_TO_REPORT);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "FIT_TO_ADAPT");
	(PARAM->P_LIST+PARAM->NP)->TYPE='F';
	(PARAM->P_LIST+PARAM->NP)->F=&((PARAM->BGA)->FIT_TO_ADAPT);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "TOT_OP_FIT");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->TOT_OP_FIT);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "STEP_FIT_OPERATOR");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->STEP_FIT_OPERATOR);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "MAX_STAB");
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "Maximum Number of Generation That the GA can remain stabilised without being stopped");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->MAX_STAB);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "experience_name");
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "Name of your experience ( No space allowed)");
	(PARAM->P_LIST+PARAM->NP)->TYPE='S';
	(PARAM->P_LIST+PARAM->NP)->S=((PARAM->BGA)->experience_name);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "MAXIMISE");	
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->MAXIMISE);
	PARAM->NP++;
	

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "BIASE_POS");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->BIASE_POS);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "ANALYSE_OPERATOR");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->ANALYSE_OPERATOR);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_SCREEN");
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "0: No Output 1: Min Output 2: Verbose");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_SCREEN);
	PARAM->NP++;
	
	

	
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "PARALLEL_GA");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->PARALLEL_GA);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "PARALLEL_EXCHANGE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->PARALLEL_EXCHANGE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "PARALLEL_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->PARALLEL_MODE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "INTERACTIVE");
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "Set to 1 for a dose of interactivity");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->INTERACTIVE);
	PARAM->NP++;
	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_PARAM_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='S';
	(PARAM->P_LIST+PARAM->NP)->S=((PARAM->BGA)->OUTPUT_PARAM_FILE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_RESULT_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='S';
	(PARAM->P_LIST+PARAM->NP)->S=((PARAM->BGA)->OUTPUT_RESULT_FILE);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_ALN_FILE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='S';
	(PARAM->P_LIST+PARAM->NP)->S=((PARAM->BGA)->OUTPUT_ALN_FILE);
	PARAM->NP++;	
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_ALN");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_ALN);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_PARAM");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_PARAM);
	PARAM->NP++;
	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OUTPUT_RESULT");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OUTPUT_RESULT);
	PARAM->NP++;


	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "OF_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->BGA)->OF_MODE);
	PARAM->NP++;	

	
	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "SEEDING_MODE");
	(PARAM->P_LIST+PARAM->NP)->TYPE='D';
	(PARAM->P_LIST+PARAM->NP)->D=&((PARAM->SEED)->SEEDING_MODE);
	PARAM->NP++;

	sprintf ((PARAM->P_LIST+PARAM->NP)->name, "%s", "REF_DEC_CHROM");
	sprintf ((PARAM->P_LIST+PARAM->NP)->description, "File containing the sequences (PIR, FASTA, EMBL, SWISS_PROT)");
	(PARAM->P_LIST+PARAM->NP)->TYPE='S';
	(PARAM->P_LIST+PARAM->NP)->S=((PARAM->BGA)->REF_DEC_CHROM);
	PARAM->NP++;	  
	initialise_pb_parameter_list (PARAM);
	    
/*CHECK 1*/
	for ( a=0; a< PARAM->NP; a++)
		{
		if ( strcmp ((PARAM->P_LIST+a)->name, "")==0)
			crash ( "empty list for operator");
		
		if ((PARAM->P_LIST+a)->TYPE=='D' && (PARAM->P_LIST+a)->D==NULL)
			{
			sprintf (error, "Operator %s (%d) does not point to a %c",(PARAM->P_LIST+a)->name,a, (PARAM->P_LIST+a)->TYPE); 	
			crash ( error);
			}
		if ((PARAM->P_LIST+a)->TYPE=='F' && (PARAM->P_LIST+a)->F==NULL)
			{
			sprintf (error, "Operator %s(%d) does not point to a %c",(PARAM->P_LIST+a)->name,a, (PARAM->P_LIST+a)->TYPE); 	
			crash ( error);
			}
		if ((PARAM->P_LIST+a)->TYPE=='S' && (PARAM->P_LIST+a)->S==NULL)
			{
			sprintf (error, "Operator %s (%d) does not point to a %c",(PARAM->P_LIST+a)->name,a, (PARAM->P_LIST+a)->TYPE); 	
			crash ( error);
			}
		
		if ((PARAM->P_LIST+a)->TYPE=='O' && (PARAM->P_LIST+a)->O_1==NULL)
			{
			sprintf (error, "Operator %s does not point to a %c",(PARAM->P_LIST+a)->name, (PARAM->P_LIST+a)->TYPE); 	
			crash ( error);
			}
	
		}
/*CHECK 2*/
	for ( a=0; a< PARAM->NP; a++)
		{
		for ( b=a+1; b< PARAM->NP; b++)
			if ( strcmp ( (PARAM->P_LIST+a)->name,  (PARAM->P_LIST+b)->name)==0)
				{
				sprintf (error, "Parameter %d and %d have the same name: %s", a, b, (PARAM->P_LIST+a)->name);
				crash ( error);
				}
		}
	
	}	
void read_parameter_file ( char *fname, Parameter *PARAM)
	{
	FILE *fp;
	int a, flag;
	int c;
	char name[200];
	char buf[1000];
	int exit_flag=0;
	
	fp=vfopen (fname,"r");

	while ( (c=fgetc ( fp))!='$' && c!=EOF);
	while ( (c=fgetc (fp))!='$'  && c!=EOF)
		{
		if (c=='*')
			fgets ( buf, 999, fp);
		else if ( isalnum(c))
			{
			fscanf ( fp, "%s", name);
			
			for (flag=0, a=0; a< PARAM->NP; a++)
				{
				if ( strcmp (name, (PARAM->P_LIST+a)->name)==0)
					{
					
					flag=1;
					if ((PARAM->P_LIST+a)->TYPE!=c)
						{
						sprintf ( error, "%s is %c in FILE %s but declared as %c\n", name, c,fname,(PARAM->P_LIST+a)->TYPE);
						}
					else
						{
						if ( c=='D')
							fscanf (fp, "%d\n",(PARAM->P_LIST+a)->D);
						else if ( c=='F')
							fscanf (fp, "%f\n",(PARAM->P_LIST+a)->F);
						else if ( c=='S')
							fscanf (fp, "%s\n",(PARAM->P_LIST+a)->S);
						else if ( c=='O')
							fscanf (fp, "%d %d\n", (PARAM->P_LIST+a)->O_1,(PARAM->P_LIST+a)->O_2);
						(PARAM->P_LIST+a)->read=1;
						}
					}
				}
			if ( flag==0)
				{
			        fprintf (stderr, "\nCOULD NOT SET PARAMETER %s in FILE %s\n", name, fname);
				exit_flag=1;
				fgets ( buf, 999, fp);
				}
			}
		}

	if ( exit_flag)exit (0);
	fclose (fp);
	} 	 
	
void read_command_line ( int narg, char **list, Parameter *PARAM)
	{
	int a,b,c;
	int flag1, flag2;
	
	for ( b=1; b< narg; b++)
		{
		if ( !is_cl_flag ( list[b]))
			{
			for (flag1=0, flag2=0, a=0; a< PARAM->NP; a++)
				{
				if ( strcmp (list[b], (PARAM->P_LIST+a)->name)==0)
					{
					flag1=1;
					c=(PARAM->P_LIST+a)->TYPE;
					
					if ( c=='D')flag2=sscanf (list[++b], "%d",(PARAM->P_LIST+a)->D);
					else if ( c=='F')flag2=sscanf (list[++b], "%f",(PARAM->P_LIST+a)->F);
					else if ( c=='S')flag2=sscanf (list[++b], "%s",(PARAM->P_LIST+a)->S);
					else if ( c=='O')
						{
						sscanf (list[++b], "%d", (PARAM->P_LIST+a)->O_1);
						flag2=sscanf (list[++b], "%d",(PARAM->P_LIST+a)->O_2);
						(PARAM->P_LIST+a)->read=1;
						}
					}
				}
				
			if ( flag1==0)
				{
				fprintf ( stderr, "\n%s is not a parameter\n", list[b]);
				crash ( "");
				}
			if ( flag2==0 && c!='O')
				{
				fprintf ( stderr, "\ncould not read any value for %s\n", list[b-1]);
				crash ( "");
				}
			if ( flag2==0 && c=='O')
				{
				fprintf ( stderr, "\ncould not read the two values for %s\n", list[b-2]);
				crash ( "");
				}
			}
		else
			{
			b++;
			}
		}
	}	

void process_param ( Parameter *PARAM)
	{
	char *p; 
    	int a;
    	char pwd_fname[200];
    	FILE *fp;
    	char buf[1000];
    	char com[1000];
    	
    	
    	
    	
    	
    	if ( strcmp ( (PARAM->BGA)->experience_name, "DEFAULT")==0)
    		{
		sprintf ( (PARAM->BGA)->experience_name,"default");    	
    		}
    	
    	
    	if ( strcmp ((PARAM->BGA)->OUTPUT_RESULT_FILE, "DEFAULT")==0)sprintf ((PARAM->BGA)->OUTPUT_RESULT_FILE, "%s.ga_result",(PARAM->BGA)->experience_name);
    	if ( strcmp ((PARAM->BGA)->OUTPUT_PARAM_FILE, "DEFAULT")==0)sprintf ((PARAM->BGA)->OUTPUT_PARAM_FILE, "%s.saga_param",(PARAM->BGA)->experience_name);
    		 
    	input_parameter ( (PARAM->BGA)->experience_name,"THE SEARCH", (PARAM->BGA)->INTERACTIVE, PARAM);
        sprintf ( (PARAM->BGA)->NORMAL_END_FLAG_FILE,"end_%s.ga", (PARAM->BGA)->experience_name);
    	sprintf ( (PARAM->BGA)->FINISH_FLAG_FILE,"finish_%s.ga", (PARAM->BGA)->experience_name);
    	sprintf ( (PARAM->BGA)->POP_READ_FLAG_FILE, "population_read_%s.ga", (PARAM->BGA)->experience_name);
    	sprintf ( (PARAM->BGA)->POP_DUMP_FLAG_FILE, "population_dumped_%s.ga", (PARAM->BGA)->experience_name);
    	sprintf ( (PARAM->BGA)->ERROR_LOG_FILE, "error_%s.ga", (PARAM->BGA)->experience_name); 
    	(PARAM->BGA)->FP_ERROR=((PARAM->BGA)->OUTPUT_ERROR_FILE==0)?fopen( "/dev/null", "w"):fopen ( (PARAM->BGA)->ERROR_LOG_FILE,"w");
    
    
    	if ((PARAM->BGA)->PARALLEL_GA==1)
    		{ 
	   	set_flag_file ((PARAM->BGA)->NORMAL_END_FLAG_FILE , 0); 
    	   	set_flag_file ((PARAM->BGA)->FINISH_FLAG_FILE , 0); 
    	   	set_flag_file ((PARAM->BGA)->POP_READ_FLAG_FILE , 0);
    	   	set_flag_file ((PARAM->BGA)->POP_DUMP_FLAG_FILE , 0);
    		}
    	else 
    		{
    		if ((PARAM->BGA)->OUTPUT_FINISH_FLAG==0)sprintf ( (PARAM->BGA)->FINISH_FLAG_FILE, "/dev/null");
    		if ((PARAM->BGA)->OUTPUT_END_FLAG==0)sprintf ( (PARAM->BGA)->NORMAL_END_FLAG_FILE, "/dev/null");
    		if ((PARAM->BGA)->OUTPUT_POP_READ_FLAG==0)sprintf ((PARAM->BGA)->POP_READ_FLAG_FILE, "/dev/null");
    		if ((PARAM->BGA)->OUTPUT_POP_DUMP_FLAG==0)sprintf ((PARAM->BGA)->POP_DUMP_FLAG_FILE, "/dev/null");
    		}
    
    	(PARAM->BGA)->OVERLAPING_GEN_FACTOR=1-(PARAM->BGA)->OVERLAPING_GEN_FACTOR;
    	if ((PARAM->BGA)->OVERLAPING_GEN_FACTOR<=1)(PARAM->BGA)->OVERLAPING_GEN_FACTOR=(int)((PARAM->BGA)->OVERLAPING_GEN_FACTOR* (PARAM->BGA)->MAXPOP);
   	
    
    	a=((PARAM->BGA)->RANDOM_SEED<0)?1:(PARAM->BGA)->RANDOM_SEED;
    	sprintf ( pwd_fname, "pwd_%d",a);
    	sprintf ( com, "pwd > %s",pwd_fname);
    	system ( com); 
    	while ((fp=fopen (pwd_fname, "r"))==NULL)
    		{
    		system ( com);
    		}
    	fscanf ( fp, "%s",(PARAM->BGA)->PWD);
    	fclose ( fp);
    	sprintf ( com, "rm %s",pwd_fname);
    	system ( com); 	
    	
    	for ( (PARAM->DOS)->N_BI_PARENT_LIST=0,a=0; a< (PARAM->DOS)->TOT_OPERATOR; a++)
    		{
    		(PARAM->DOS)->N_BI_PARENT_LIST+=((PARAM->DOS)->BI_PARENT_LIST[a]==2)?1:0;
    		if ((PARAM->DOS)->OP_INIT_VAL[a][0]>0)
    			{
    			(PARAM->DOS)->OP_TO_USE[a][0]=1;
    			(PARAM->DOS)->N_OPERATOR++;
    			}
    		else
    			{
    			(PARAM->DOS)->OP_TO_USE[a][0]=0;
    			(PARAM->DOS)->OP_INIT_VAL[a][0]=0;
    			}
    		}
    	
    	if ( (PARAM->BGA)->OUTPUT_PARAM==1)output_parameters ( PARAM);
    	if ( (PARAM->BGA)->OUTPUT_RESULT==1)fprintf ( (PARAM->std1),      "\n            RESULT FILE: %s\n", (PARAM->BGA)-> OUTPUT_RESULT_FILE);
    		
        }
void output_parameters ( Parameter *PARAM)
	{
	FILE *fp;
	int a, b, c;
	
	fprintf ( (PARAM->std1), "\nYOUR PARAMETER FILE: %s", (PARAM->BGA)->OUTPUT_PARAM_FILE);
	
	fp=fopen ( (PARAM->BGA)->OUTPUT_PARAM_FILE, "w");
	fprintf (fp, "These are the parameters used for the %s run\n$\n", (PARAM->BGA)->experience_name);
	for ( a=0; a< (PARAM->NP); a++)
		{
		if ( (PARAM->P_LIST+a)->TYPE=='D')
			fprintf ( fp, "D	%s	%d\n", (PARAM->P_LIST+a)->name,*(PARAM->P_LIST+a)->D);
		else if (  (PARAM->P_LIST+a)->TYPE=='S')
			fprintf ( fp, "S	%s	%s\n", (PARAM->P_LIST+a)->name,(PARAM->P_LIST+a)->S);
		else if (  (PARAM->P_LIST+a)->TYPE=='F')
			fprintf ( fp, "F	%s	%f\n", (PARAM->P_LIST+a)->name,*(PARAM->P_LIST+a)->F);
		else if (  (PARAM->P_LIST+a)->TYPE=='O')
			fprintf ( fp, "O	%s	%d	%d\n", (PARAM->P_LIST+a)->name,*(PARAM->P_LIST+a)->O_1,*(PARAM->P_LIST+a)->O_2 );
		}
	fprintf (fp, "$\n");
	fclose (fp);
	}
			
		
int is_cl_flag (char *name)
	{
	if (name[0]=='-')return 1;
	else return 0;
	}

void option_help ( Parameter *PARAM, char *option)
	{
	int a;
	int l=strlen ( option);
	int flag=0;
	for ( a=0; a< PARAM->NP; a++)
		{
		if ( (strncmp (option, (PARAM->P_LIST+a)->name,l)==0) ||  (strcmp ( option, "ALL")==0) || (strstr((PARAM->P_LIST+a)->name, option)!=NULL))
			{
			flag=1;
			get_paramter_help ( (PARAM->P_LIST+a)->name, (PARAM->P_LIST+a));
			fprintf ( stdout, "\nNAME: %s\nTYPE: %c \nDESCRIPTION: \n%s\n",(PARAM->P_LIST+a)->name , (PARAM->P_LIST+a)->TYPE,(PARAM->P_LIST+a)->description);
			}
		}
	if ( flag==0)
		{
		fprintf ( stdout, "\nNo match to %s, try to shorten the name to get an approximate match", option);
		}
	}
void get_paramter_help ( char *name, P_list *P)
	{
	char buf[10000],c;
	FILE *fp;
	int x;
	fp=fopen ( "DOCUMENTATION", "r");
	if ((fp=find_token_in_file ("DOCUMENTATION" , fp,name))==NULL)
		sprintf ( buf, "\nDESCRIPTION NOT AVAILABLE");
	else
		{
		x=0;
		while ( (c=fgetc(fp))!='@')
			if ( c!='*')buf[x++]=c;
		buf[x]='\0';
		fclose (fp);
		}
	P->description=calloc ( strlen(buf)+1, sizeof (char));
	sprintf ( P->description, "%s", buf);
	}
	
void output_help ( Parameter *PARAM, char *exec)
	{
	fprintf ( stderr, "\n\nusage: %s <m: main_parameter_file> <s: secondary_parameter_file> <-PARAMETER value>", exec);
	fprintf (stderr, "\nenter: %s option: <parameter name> for some doc about a given parameter", exec);
	fprintf (stderr, "\nenter: %s option: ALL    (for some doc about ALL the parameters)", exec);
	return;
	}
	

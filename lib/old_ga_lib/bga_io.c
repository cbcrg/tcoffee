# define MAIN 1
# include "ga.h"

    
void do_command ( char*com, Parameter *PARAM)
    {
    FILE *fp;
    char string[1000];
    
    
	sprintf (string,"%s > /dev/null", com);
	system ( string);
    }
void get_ref_time (Parameter *PARAM)
	{
	struct tms time_buf[1];
	times ( time_buf);
	
	(PARAM->BGA)->ref_time=time_buf->tms_utime;
	} 
void save_individual (Chromosome *A,char *name, int generation, Parameter *PARAM)
       {
	 FILE *fp;
	 char best_file[100];
	 int a,c,b,d;
	 int line=60;
	 struct tms time_buf[1];
	 int tot;
	 int sec;
	 int min;
	 int hour;
	 char buf[100];
	 int tot_sec;
	 char fname[100];
	 char text [2000];
         char btext [2000];
	 
	 times (time_buf);
	 tot=(time_buf->tms_cutime+time_buf->tms_utime)+-(PARAM->BGA)->ref_time;

	 tot_sec=sec= tot/100;
	 min= sec/60;
	 sec= sec- min*60;
	 hour= min/60;
	 min= min - hour*60;
	
	 if ((PARAM->BGA)->OUTPUT_CPU_FILE==1 && name!=NULL)
	 	{
	 	sprintf ( fname, "%s_cpu.saga", name);
	 	fp=vfopen ( fname, "w");
	 	fprintf ( fp, "%d", tot_sec);
	 	fclose ( fp);
	 	}
	 if ((PARAM->BGA)->OUTPUT_SCORE_FILE==1 && name!=NULL)
	 	{	
	 	sprintf ( fname, "%s_score.saga", name);
	 	fp=vfopen ( fname, "w");
	 	fprintf( fp, "%d", (int)A->score);
	 	fclose ( fp);
		}
	 if ((PARAM->BGA)->OUTPUT_GENERATION_FILE==1 && name!=NULL)
	 	{ 
	 	sprintf ( fname, "%s_generation.saga", name);
	 	fp=vfopen ( fname, "w");
	 	fprintf( fp, "%d", generation);
	 	fclose ( fp);
   		}
   		
	
	 
	 sprintf ( text, "Total time : %d hour %d min %d sec\n", hour, min, sec);
	 sprintf (btext, "Random Seed= %d\n", (PARAM->BGA)->RANDOM_SEED);
	 strcat ( text, btext);
	 sprintf (btext, "Reference score=%d\n",((PARAM->BGA)->REFERENCE)?(int)((PARAM->BGA)->REFERENCE)->score:-1);
	 strcat ( text, btext);
	 sprintf (btext, "Score=%d\n", (int)A->score);
	 strcat ( text, btext);
	 sprintf (btext, "Generation=%d\n", generation);
	 	 
	 make_ga_output ( name, text,NULL,A, PARAM);
	 }

void print_result_ga ( Population *L_POP,int G, Statistic * L_STAT, int generation, Parameter *PARAM)    
     {
     static Chromosome *A;
     Decoded_chromosome *DA;
     static previous_score;
     static N_GEN;
     static int first;		
     FILE *fp;
     char fname[1000];
     static io_pld_score;
     static static_x;
     
     DA=get_mem_free_decoded_chrom(PARAM);
     if ( A==NULL)
	{A=(PARAM->MEM)->C1;
	 sprintf ( fname,"bilan_%s.saga",(PARAM->BGA)->experience_name);
	 if ( (PARAM->BGA)->OUTPUT_BILAN_FILE==1)
	 	{
	 	fp=vfopen ( fname, "w");
	 	fclose ( fp);
		}
	}
 
     sprintf (fname, "gbest_%s",(PARAM->BGA)->experience_name);
     
    
	
    first=1;
	   
    extract_chrom (L_POP,G,L_STAT->best_gen_chrom,A, PARAM);
    decode_chromosome(A, DA, PARAM);		    		    
    copy_chromosome ( A, L_STAT->BEST_GEN_CHROMOSOME, PARAM);
    copy_decoded_chromosome ( DA, L_STAT->BEST_GEN_DECODED_CHROMOSOME, PARAM);
	    				
    if ( (PARAM->BGA)->OUTPUT_GEN_BEST_FILE==1)save_individual (A,fname,generation, PARAM);
		
    if (a_better_than_b(L_STAT->score_gen_best,previous_score,(PARAM->BGA)->MAXIMISE )|| previous_score==0 || ((PARAM->BGA)->EVALUATE_ONLY==1))
		{
		save_individual (A,NULL,generation, PARAM);
		copy_decoded_chromosome ( DA, L_STAT->BEST_DECODED_CHROMOSOME, PARAM);
		copy_chromosome ( A, L_STAT->BEST_CHROMOSOME, PARAM);
		}	 	 
    previous_score= L_STAT->score_gen_best;
    N_GEN=0;

     sprintf ( fname,"bilan_%s.saga",(PARAM->BGA)->experience_name );
     if ( ((PARAM->BGA)->OUTPUT_BILAN_FILE)==1)
	{
	fp=vfopen (fname,"a");
	fprintf ( fp,"\nGEN=%d\tBEST=%d\tMEAN=%d\tSD=%d", generation,L_STAT->score_gen_best,L_STAT->mean, L_STAT->SD);
	fclose ( fp);
	}	
     
     fprintf ((PARAM->std2),"\nGEN=%d\tBEST=%d\tMEAN=%d\tSD=%d", generation,L_STAT->score_gen_best,L_STAT->mean, L_STAT->SD);    
     
     fprintf ((PARAM->std1),"%c", (io_pld_score==L_STAT->score_gen_best)?'=':'*');
     io_pld_score=L_STAT->score_gen_best;
     if ( static_x==30)
     	{
     	fprintf ((PARAM->std1),"[G%d %d]\n",generation,L_STAT->score_gen_best);
     	static_x=0;
     	}
     static_x++;
     
     free_mem_dc(DA, PARAM);
     }

void print_current_operator (Parameter *PARAM)
    {
    int a;
    for ( a=0; a< (PARAM->DOS)->TOT_OPERATOR; a++)
	if ( ((PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "\n#op %2d %10s : %d", a,(PARAM->DOS)->OP_NAME[a], (PARAM->DOS)->current_fitness[a][0]);

    }

void save_current_operator ( Parameter *PARAM, char *experience_name, int gen)
    {
    FILE *fp;
    static char *fname;
    int a;

    if ( (PARAM->BGA)->OUTPUT_CURRENT_OPERATOR_FILE==0)
	return;

    if ( fname==NULL)
	{
	fname=vmalloc ( sizeof(char)*100);
	sprintf ( fname, "current_operator_%s.saga", experience_name);
	fp=vfopen (fname, "w");
	fprintf ( fp, "#%d", (PARAM->DOS)->TOT_OPERATOR);
	fclose(fp);
	}

    fp=vfopen ( fname, "a");
    fprintf (fp, "\nG%d", gen);
    for ( a=0; a< (PARAM->DOS)->TOT_OPERATOR; a++)
	fprintf ( fp, "\n[%2d][%15s] : %d", a, (PARAM->DOS)->OP_NAME[a],(PARAM->DOS)->current_fitness[a][0]);
    fclose ( fp);
    }
    
void save_analyse_operator ( int gen, Parameter *PARAM)
    {
    FILE *fp;
    static int instant;
    char fname1[200];
    char fname2[200];
    int a;
    int tot_op;
    int tot_similar_p;
    int	tot_better_p;
    int	tot_integrated;
    int	P_tot_similar_p;
    int	P_tot_better_p;
    int	P_tot_integrated;    
    
    
    sprintf ( fname1, "IOA_%s.saga", (PARAM->BGA)->experience_name);
    sprintf ( fname2, "TOA_%s.saga", (PARAM->BGA)->experience_name);

    if (instant==0)
	{
	instant=1;
	fp=vfopen ( fname1, "w");
	fprintf ( fp, "\n1-number times op used\n2-number times child=parent\n3-number times child better than parent\n4-number times child integrated");
	fprintf ( fp, "\nGEN 0:");
	}
    else 
	{
	fp=vfopen ( fname1, "a");
	fprintf ( fp, "\nGEN %d:", gen);
	}

    for ( a=0; a< (PARAM->DOS)->TOT_OPERATOR; a++)
	{
	tot_op=(PARAM->DOS)->IOPUSE[a][0];
	tot_similar_p=(PARAM->DOS)->IOPUSE[a][1];
	tot_better_p=(PARAM->DOS)->IOPUSE[a][2];
	tot_integrated=(PARAM->DOS)->IOPUSE[a][3];
	P_tot_similar_p= (tot_op>0)?((tot_similar_p*100)/tot_op):100;
	P_tot_better_p=(tot_op>0)?((tot_better_p*100)/tot_op):100;
	P_tot_integrated= (tot_op>0)?((tot_integrated*100)/tot_op):100;

	fprintf ( fp, "\n%-25s%4d:[100]**%4d:[%-3d]**4%d:[%-3d]**%4d:[%-3d]",(PARAM->DOS)->OP_NAME[a], tot_op, tot_similar_p,P_tot_similar_p, tot_better_p, P_tot_better_p,tot_integrated, P_tot_integrated);      
	}
    fclose ( fp);

    fp=vfopen ( fname2, "w");
    fprintf ( fp, "\n1-number times op used\n2-number times child=parent\n3-number times child better than parent\n4-number times child integratedn\nGEN %d:", gen);
    for ( a=0; a< (PARAM->DOS)->TOT_OPERATOR; a++)
	{
	tot_op=(PARAM->DOS)->OPUSE[a][0];
	tot_similar_p=(PARAM->DOS)->OPUSE[a][1];
	tot_better_p=(PARAM->DOS)->OPUSE[a][2];
	tot_integrated=(PARAM->DOS)->OPUSE[a][3];
	P_tot_similar_p= (tot_op>0)?((tot_similar_p*100)/tot_op):100;
	P_tot_better_p=(tot_op>0)?((tot_better_p*100)/tot_op):100;
	P_tot_integrated=(tot_op>0)?((tot_integrated*100)/tot_op):100;

	fprintf ( fp, "\n%-25s%4d:[100]**%4d:[%-3d]**%4d:[%-3d]**%4d:[%-3d]",(PARAM->DOS)->OP_NAME[a], tot_op, tot_similar_p,P_tot_similar_p, tot_better_p, P_tot_better_p,tot_integrated, P_tot_integrated);      
	}
    fclose ( fp);

    }
void print_select_ga ( Population *L_POP, Parameter *L_PARAM)
    {
    int a;
    Decoded_chromosome *DC;
    
    DC=get_mem_free_decoded_chrom ( L_PARAM);
    if ( ((L_PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "************************************SELECTION***********************");
    for ( a=0; a<(L_PARAM->BGA)->MAXPOP; a++)
	{
	if ( ((L_PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "\n%d",a);
	if ( ((L_PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "\t%d", (int)(L_POP+0)->raw_fitness[a][0]);
	if ( ((L_PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "\t%d" , (int)(L_POP+0)->fitness[a][0]);
	if ( ((L_PARAM->BGA)->OUTPUT_SCREEN)==1)printf ( "\t%d" , (int)(L_POP+0)->expected_os[a][0]);
	}
    free_mem_dc ( DC, L_PARAM);	
    }

void make_ga_output ( char *name,char *text, Decoded_chromosome *ALN, Chromosome *C_ALN, Parameter *PARAM)
    {
    Decoded_chromosome *S, *A, *B;
    char score_name[200];
    char filter_name[200];
    char aln_name[200];
    int x1, x2, x3, x4;
    
    
    
    A=get_mem_free_decoded_chrom(PARAM);
    
    
    
    if ( C_ALN==NULL)A=copy_decoded_chromosome ( ALN, A, PARAM);
    else if ( ALN==NULL)A=decode_chromosome ( C_ALN, A, PARAM);
    
    
    if (name==NULL)
        {
        sprintf ( aln_name, "%s", (PARAM->BGA)->OUTPUT_ALN_FILE);
        
        }
    else
        {
        sprintf ( aln_name, "%s", name);
        sprintf ( score_name, "%s_rel", name);
        sprintf ( filter_name, "%s_filter", name);
        }
 
    A->score=evaluate_decoded_chromosome (A, (PARAM->BGA)->OF_MODE, PARAM);

    if ( (PARAM->BGA)->OUTPUT_ALN==1)output_decoded_chromosome (A, aln_name, text, PARAM);     
    free_mem_dc (A, PARAM);    
    }

void  print_chromosome (Chromosome *C, Parameter *PARAM)
      {
      static Decoded_chromosome *D;

      D=decode_chromosome ( C, D, PARAM);
      print_decoded_chromosome ( D,PARAM);
      }
      
 
void  print_decoded_chromosome ( Decoded_chromosome *B, Parameter *PARAM)
      {
      save_decoded_chromosome ( B, stdout, NULL, PARAM);
      }

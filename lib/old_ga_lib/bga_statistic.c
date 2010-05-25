#define MAIN 1
#include "ga.h"

void do_statistic ( Population *L_POP,int G, Statistic *L_STAT, int tot_pop, int generation, Parameter *PARAM)
    {
    float **dup_fitness;
    char fname[100];
    FILE *fp;
    int a,b;
    int B=0;
    int W=0;
    int max;
    
     max=(PARAM->BGA)->MAXIMISE;
    

    if (L_STAT->score_best==0)
	{
	 L_STAT->score_best=(int)(L_POP+G)->raw_fitness[0][0];
	}

 
    for ( a=0; a< tot_pop; a++)
	(L_POP+G)->raw_fitness[a][1]=(L_POP+G)->fitness[a][1]=a;

    dup_fitness=duplicate_float ( (L_POP+G)->raw_fitness, tot_pop, 2);
    sort_float ( dup_fitness, 2, 0, 0, tot_pop-1);
   
    if ( max==1)
	{B=(PARAM->BGA)->MAXPOP-1;
	 W=0;
	}   
    else 
	{B=0;
	 W=(PARAM->BGA)->MAXPOP-1;
	 }
    L_STAT->best_gen_chrom=(int)dup_fitness[B][1]; 
    L_STAT->score_gen_best=(int)dup_fitness[B][0];
    L_STAT->score_gen_worst=(int)dup_fitness[W][0];
  
    L_STAT->mean=(int) return_mean_float ( dup_fitness, tot_pop, 0);
    L_STAT->SD= (int) return_sd_float ( dup_fitness, tot_pop, 0, (float)L_STAT->mean);	    
    (PARAM->BGA)->T=return_mean_diff_float ( dup_fitness, tot_pop, 0, (float)L_STAT->mean);	
    
    for ( a=0, b=0; a< (PARAM->BGA)->MAXPOP; a++)
    	b+=((L_POP+G)->CHROMOSOME[a])->len;
    b=b/(PARAM->BGA)->MAXPOP;
    (PARAM->BGA)->T=(PARAM->BGA)->T/(float)b;
    
    
    if (L_STAT->score_best>L_STAT->score_gen_best)
	{
	L_STAT->score_best=L_STAT->score_gen_best;
	}
    free_float ( dup_fitness, tot_pop);
    }



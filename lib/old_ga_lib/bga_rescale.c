#define MAIN 1
#include "ga.h"

void rescale_fitness (Population *L_POP, int G, int tot_pop, int MODE, Parameter *L_PARAM)
    {
    int a,b;    

    
    if ( MODE==0)
	natural_fitness ( L_POP, G, tot_pop);

    else if ( MODE==1)
	ranked_fitness ( L_POP, G, tot_pop);

    else if ( MODE==2)
	affine_fitness ( L_POP, G, tot_pop, (L_PARAM->BGA)->RESCALING_C);
	
    else if ( MODE==3)
	sigma_fitness ( L_POP,G,tot_pop, (L_PARAM->BGA)->RESCALING_C);
    
    }

void natural_fitness (Population *L_POP, int G, int tot_pop)
    {
    int a;
    float sum;
    float b;

    sum=return_sum_float ((L_POP+G)->fitness, tot_pop, 0);
    for ( a=0; a< tot_pop; a++)
	{b=(L_POP+G)->scaled_fitness[a][0]=  ((L_POP+G)->fitness[a][0]*tot_pop)/sum;
	}
    }  
   

void ranked_fitness (Population *L_POP, int G, int tot_pop)
    {
    int a, b;
    float c;
    float sum;
    float bf1;
    float bf2;

    sort_float ( (L_POP+G)->ranked, 2, 0, 0, tot_pop-1);

    bf1=(L_POP+G)->ranked[0][0];
    (L_POP+G)->ranked[0][0]=1;
        
    for ( a=1; a<tot_pop; a++)
	{
	c=(L_POP+G)->ranked[a-1][0];
	bf2=(L_POP+G)->ranked[a][0];
	(L_POP+G)->ranked[a][0]=( bf1==(L_POP+G)->ranked[a][0])?c:c+1;
	bf1=bf2;
	}
    
    sum=return_sum_float ( (L_POP+G)->ranked, tot_pop,0);
    
    for ( a=0; a<tot_pop; a++)
	{
	(L_POP+G)->scaled_fitness[(int)(L_POP+G)->ranked[a][1]][0]=((L_POP+G)->ranked[a][0]*tot_pop)/sum;
	}
    }

void affine_fitness ( Population *L_POP, int G, int tot_pop, float C)
    {
    
    float min;
    float max;
    float mean;
    float a,b;
    int c;
    float sum; 

    min= return_min_float ( (L_POP+G)->fitness, tot_pop, 0);
    max= return_max_float ( (L_POP+G)->fitness, tot_pop, 0);
    mean=return_mean_float ( (L_POP+G)->fitness, tot_pop, 0);

    if ( min==mean && mean==max)
	{
	a=1;
	b=0;
	}
    else if ( min >((C*mean -max)/(C-1)))
	{
	a=(float)((C-1)*mean)/(float)(max-mean);
	b=(float)mean*(float)(max-c*mean)/(float)(max-mean);
	}
    else
	{
	a=mean/(mean-min);
	b= (-min*mean)/(mean-min);
	}

    if ( a==0)
	{   
	a=1;
	b=0;
	}

    
    for ( c=0; c< tot_pop; c++)
	(L_POP+G)->scaled_fitness[c][0]=a* (L_POP+G)->fitness[c][0]+b;

    sum=return_sum_float ( (L_POP+G)->scaled_fitness, tot_pop, 0);
    
    for ( c=0; c< tot_pop; c++)
	(L_POP+G)->scaled_fitness[c][0]=((L_POP+G)->scaled_fitness[c][0]*tot_pop)/sum;
    }

void sigma_fitness ( Population *L_POP, int G, int tot_pop, float C)
    {
    
    float min;
    float max;
    float mean;
    int a,b;
    int c;
    float sum; 
    float sd;

    mean=return_mean_float ( (L_POP+G)->fitness, tot_pop, 0);
    sd=return_sd_float ((L_POP+G)->fitness, tot_pop,0, mean);
     
    for ( c=0; c< tot_pop; c++)
	(L_POP+G)->scaled_fitness[c][0]= ((a=(L_POP+G)->fitness[c][0]-(mean -C*sd))>0)?a:0;

    sum=return_sum_float ( (L_POP+G)->scaled_fitness, tot_pop, 0);
    
    for ( c=0; c< tot_pop; c++)
	(L_POP+G)->scaled_fitness[c][0]=((L_POP+G)->scaled_fitness[c][0]*tot_pop)/sum;
    }


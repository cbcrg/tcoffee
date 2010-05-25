#define MAIN 1
#include "ga.h"

void preselect ( Population *L_POP,int G, int tot_pop, int MODE, Parameter *L_PARAM)
    {
    int a,b;
    
    
    if (MODE==0)
	for ( a=0; a< tot_pop; a++)
	    (L_POP+G)->expected_os[a][0]=(L_PARAM->BGA)->FACTOR_POP;	
    else if (MODE==1)
	no_preselect (L_POP, G, tot_pop, L_PARAM);
    else if ( MODE==2)
	remainder_preselect (L_POP, G, tot_pop, L_PARAM);
    
    /*if ( PARAM->DEBUG_3)
	{
	if ( PARAM->RESCALING_MODE!=1)
	    sort_float ( (L_POP+G)->ranked, 2, 0, 0, tot_pop-1);	    
	 (PARAM->BGA)->OUTPUT_SCREEN = val_para[23]
	 if ((PARAM->BGA)->OUTPUT_SCREEN==1) printf ("\nSUM_POP=%d MAXPOP=%d", return_sum_int ( (L_POP+G)->expected_os, tot_pop, 0),tot_pop*PARAM->FACTOR_POP);

	
	for ( a=0; a< tot_pop; a++)
	    {
	    b=(int)(L_POP+G)->ranked[a][1];
	    if ((PARAM->BGA)->OUTPUT_SCREEN==1)printf ("\n%d=%d\t%d=%d\t%f\tOS=%d",a,(int)(L_POP+G)->fitness[a][0],b, (int)(L_POP+G)->fitness[b][0],(L_POP+G)->scaled_fitness[b][0], (L_POP+G)->expected_os[b][0]);
	    }
	}
	*/
    }

void no_preselect (  Population *L_POP, int G, int tot_pop, Parameter *L_PARAM)
    {
    int a;

    for ( a=0; a<tot_pop; a++)
	(L_POP+G)->expected_os[a][0]=((L_PARAM->BGA)->FACTOR_POP* (L_POP+G)->scaled_fitness[a][0])+1;

    }
   
void remainder_preselect( Population *L_POP, int G, int tot_pop, Parameter *L_PARAM)
    {
    float **remainder;
    int a,b;
    int i;
    int total;

    remainder=declare_float ( tot_pop, 1);
    
    for ( a=0; a< tot_pop; a++)
	{
	(L_POP+G)->expected_os[a][0]=(int)(L_POP+G)->scaled_fitness[a][0];
	remainder[a][0]=(float)((float)(L_POP+G)->scaled_fitness[a][0]-(int)(L_POP+G)->expected_os[a][0]);
	}
    
    total=return_sum_int ( (L_POP+G)->expected_os, tot_pop, 0);
    
    i= -1;

   
    while( total<tot_pop)
	{
	i=( i==tot_pop-1)?0:i+1;
	
	if ( remainder[i][0]>0)
	    {
	    if ( (b=float_flip(remainder[i][0])))
		{
		total++;
		(L_POP+G)->expected_os[i][0]++;
		remainder[i][0]=0;
		}
	     }
	
	}

    for ( a=0; a< tot_pop; a++)
	{
	(L_POP+G)->expected_os[a][0]=(L_POP+G)->expected_os[a][0]*(L_PARAM->BGA)->FACTOR_POP;	
	}
    free_float (remainder, tot_pop); 	
	
    }		 
      
    		   

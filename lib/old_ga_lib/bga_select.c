#define MAIN 1
#include "ga.h"

void ga_select ( Population *L_POP,int G, int ** m_t,int tot_pop,int MODE, Parameter *PARAM)
    {
    int a,b;
    int sum;
    int **local_os;
    int max;
    int min;

    local_os=duplicate_int ( (L_POP+G)->expected_os, tot_pop, 1);
    
    sum=return_sum_int (local_os,  tot_pop,0);
 

    /*select the mates*/

    for ( a=0; a< tot_pop; a++)
	{
	    m_t[a][0]=wheel_select (local_os, tot_pop,sum,0,MODE);
	    sum=return_sum_int (local_os, tot_pop,0);
	    m_t[a][1]=wheel_select(local_os, tot_pop,sum,0,MODE);   
	    sum=return_sum_int (local_os, tot_pop,0);
	}

    
    		
    free_int ( local_os, tot_pop);

    /* select those who are going to be replaced*/

    local_os=duplicate_int ( (L_POP+G)->expected_os, tot_pop, 1);
    
    max=return_max_int (local_os, tot_pop,0);
    min=return_max_int (local_os, tot_pop,0);
    inverse_int ( local_os,tot_pop,0, max, min);

    MODE=(PARAM->BGA)->REPLACE_MODE;
    sum=return_sum_int (local_os,  tot_pop,0);
    
    b=0;
    for ( a=0; a< tot_pop; a++)
	{
	if ( MODE<2)
	   {
	   m_t[a][2]=wheel_select (local_os, tot_pop,sum,0,MODE);
	   local_os[m_t[a][2]][0]=0;
	   sum=return_sum_int (local_os, tot_pop,0);
	  
	   }  
	else if ( MODE==2)
	    {
	    m_t[a][2]=(L_POP+G)->ranked[a][1];
	    }
	else if ( MODE==3)
	    {
	    m_t[a][2]=a;
	    }
	}      	
	    
     free_int (local_os, tot_pop);
     }


/*select a random number with an exponential distribution between firdst and second, according
to factor 1 means uniform distribution*/
int exponential_selection ( int first_len, int second_len, float factor)
    {
    int **array;
    int T;
    float p;
    int choice;
    int sum;
    int a;
    int k;

    if ( second_len==0)
	return 0;
    else if ( first_len==second_len)
	return first_len;
    else if ( first_len>second_len)
	{
	printf ("\nWarning in chose_len maxlen<minlen");
	return 0;
	}
    else
	{
	
	T= int_fabs ( second_len-first_len) +1;
	if ( factor==1)
	    {
	    k=addrand ( (unsigned long) T);
	    choice= ( first_len < second_len)?(first_len+k):(first_len-k);
	    return choice;
	    }

	if ( factor< 0)
	    swap_int ( &second_len, &first_len, 1);
    
	array=declare_int (T, 1);
    
	for ( a=1;a<=T; a++)
	    {
	    p=1000* pow ( factor, a);
	    array[a-1][0]=(int)p;
	    }

	sum=return_sum_int ( array, T, 0);
	choice=wheel_select ( array,T, sum, 0, 0);

	free_int ( array, T);

	choice= ( first_len< second_len)?(first_len+choice):first_len-choice;
	return choice;
	}
    } 
        
        
int fatality_select ( int **fatality, Parameter *L_PARAM, int sum, int *x, int tot)
    {
    int c;
    if ((L_PARAM->BGA)->REPLACE_MODE==0 ||(L_PARAM->BGA)->REPLACE_MODE==1 )
	{
	c=wheel_select ( fatality, (L_PARAM->BGA)->MAXPOP, sum, 0, 0);
	x[0]=fatality[c][0];
	return c;
	}
    else if ( (L_PARAM->BGA)->REPLACE_MODE==2)
	{
	x[0]=0;
	return fatality [(L_PARAM->BGA)->MAXPOP-(tot+1)][1];
	}
     } 

int select_operator ( Parameter *L_PARAM)
    {
    return wheel_select ( (L_PARAM->DOS)->current_fitness, (L_PARAM->DOS)->TOT_OPERATOR,(L_PARAM->BGA)->TOT_OP_FIT, 0,0);
    }
 
int wheel_select ( int **array, int tot_pop, int sum, int field, int MODE)
    {
    int a;
    int b=0;
    int c=0;
	
    if (tot_pop==0)
    	return 0;
    else if (sum==0)
    	return addrand ((unsigned long)tot_pop);
    else
    	{	
    	a=addrand( (unsigned long) sum)+1;
    
    	while ( b<a)
		b+=array[c++][field];

    	return c-1;
    	}
    }	
    
/* generate an exponential distribution between first_len and second_len*/
int select_len ( int first_len, int second_len, float factor)
    {
    int **array;
    int T;
    float p;
    int choice;
    int sum;
    int a;
    int k;

    if ( second_len==0)
        return 0;
    else if ( first_len==second_len)
        return first_len;
    else if ( first_len>second_len)
        {
        printf ("\nWarning in chose_len maxlen<minlen");
        return 0;
        }
    else
        {
        T= (( (second_len-first_len)>0)?(second_len-first_len):(first_len-second_len)) +1;
        if ( factor==1 )
            {
            k=addrand ( (unsigned long) T);
 	    		
            choice= ( first_len < second_len)?(first_len+k):(first_len-k);
            return choice;
            }

        if ( factor< 0)
            swap_int ( &second_len, &first_len, 1);

        array=declare_int (T+1, 1);

        for ( a=1;a<=T; a++)
            {
            p=1000* pow ( factor, a);
            array[a-1][0]=(int)p;
            }

        sum=return_sum_int ( array, T, 0);
        choice=wheel_select ( array,T, sum, 0, 0);

        free_int ( array, T);

        choice= ( first_len< second_len)?(first_len+choice):first_len-choice;
        return choice;
        }
    }

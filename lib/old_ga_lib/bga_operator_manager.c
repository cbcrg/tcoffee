#define MAIN 1
#include "ga.h"

void report ( int a, int credit, Parameter *L_PARAM, int G)
    {
    int P_1=0;
    int P_2=1;
    int C=2;
    int OPER=3;
    int DELTA=4;
    float factor;
    int p1;
    int p2;
    int op_p1;
    int op_p2;
    float cp1;
    float cp2;

    
    factor=(float)(0.5*(float)(L_PARAM->BGA)->CREDIT_TO_REPORT);
    if ( G>0 && credit>500)
	{
	p1=(L_PARAM->DOS)->T[G][a][P_1];
	p2=(L_PARAM->DOS)->T[G][a][P_2];
	op_p1= (L_PARAM->DOS)->T[G-1][p1][OPER];
	op_p2= (L_PARAM->DOS)->T[G-1][p2][OPER];
	(L_PARAM->DOS)->T[G-1][p1][DELTA]+= (op_p1== -1)?0:(int)((float)credit*(float)factor); 
	(L_PARAM->DOS)->T[G-1][p2][DELTA]+= (op_p2== -1)?0:(int)((float)credit*(float)factor);
	
	if ( op_p1!= -1)
	    cp1=((float)credit*(float)factor);
	else
	    cp1=credit;

	if ( op_p2!= -1)
	    cp2=((float)credit*(float)factor);
	else
	    cp2=credit;

	report (p1,(int)cp1, L_PARAM, G-1);
	report (p2,(int)cp2, L_PARAM, G-1);
	}

    else 
	return;
    
    return;
    }

void init_op_fitness (Parameter *L_PARAM)
    {
    int f;
    int a;
    int sum=0;

    for ( a=0; a< (L_PARAM->DOS)->TOT_OPERATOR; a++)
	{
	sum+=(L_PARAM->DOS)->OP_INIT_VAL[a][0];
	}
    
	
    f= (int)((float)(L_PARAM->BGA)->TOT_OP_FIT/(float)sum);
    
    (L_PARAM->BGA)->TOT_OP_FIT=0;
    for ( a=0; a< (L_PARAM->DOS)->TOT_OPERATOR; a++)
	{
	if ( (L_PARAM->DOS)->OP_TO_USE[a][0]==1)
	    (L_PARAM->DOS)->current_fitness[a][0]=f*(L_PARAM->DOS)->OP_INIT_VAL[a][0];
	else
	    (L_PARAM->DOS)->current_fitness[a][0]=0;
	(L_PARAM->BGA)->TOT_OP_FIT+=(L_PARAM->DOS)->current_fitness[a][0];
	(L_PARAM->DOS)->OP_INIT_VAL[a][0]=(L_PARAM->DOS)->current_fitness[a][0];
	}

     save_current_operator (L_PARAM, (L_PARAM->BGA)->experience_name, 0);
     }

void evaluate_fitness_operator ( Parameter *L_PARAM)
    {
    int a,b,c;
    int sum_new_fit;
    float f1, f2;
    float fa, fb;

    for ( a=0; a< (L_PARAM->DOS)->TOT_OPERATOR; a++)
	for ( b=0; b< 3; b++)
	    (L_PARAM->DOS)->fitness_table[a][b]=0;

    for (a=0; a< (L_PARAM->BGA)->STEP_FIT_OPERATOR; a++)
	{
	for ( b=0; b<(L_PARAM->BGA)->MAXPOP; b++)
	    {
	     if ( (L_PARAM->DOS)->T[a][b][3]!=-1)
		{
		(L_PARAM->DOS)->fitness_table[(L_PARAM->DOS)->T[a][b][3]][1]+=(L_PARAM->DOS)->T[a][b][4];
		(L_PARAM->DOS)->fitness_table[(L_PARAM->DOS)->T[a][b][3]][0]++;
		}
	    }
	}
    
  

    sum_new_fit=0;
    for ( a=0; a<(L_PARAM->DOS)->TOT_OPERATOR; a++)
	{
	(L_PARAM->DOS)->fitness_table[a][2] = ((L_PARAM->DOS)->fitness_table[a][0]>0)?(L_PARAM->DOS)->fitness_table[a][1]/(L_PARAM->DOS)->fitness_table[a][0]:0;
	sum_new_fit+=(L_PARAM->DOS)->fitness_table[a][2];
	}

 
 
    if ( sum_new_fit>0)
	{
	f1=((L_PARAM->BGA)->TOT_OP_FIT*(L_PARAM->BGA)->FIT_TO_ADAPT)/sum_new_fit;
	f2= 1-(L_PARAM->BGA)->FIT_TO_ADAPT;

	(L_PARAM->BGA)->TOT_OP_FIT=0;
	for (a=0; a< (L_PARAM->DOS)->TOT_OPERATOR; a++)
	    {    
	    fa=(L_PARAM->DOS)->fitness_table[a][2]*f1;
	    fb=(L_PARAM->DOS)->OP_INIT_VAL[a][0]*f2;
	    (L_PARAM->DOS)->current_fitness[a][0]=fa+fb;
	    (L_PARAM->BGA)->TOT_OP_FIT+=(L_PARAM->DOS)->current_fitness[a][0];
	    }
	}
    else
	{
	(L_PARAM->BGA)->TOT_OP_FIT=0;
	for (a=0; a< (L_PARAM->DOS)->TOT_OPERATOR; a++)
	    {    
	    (L_PARAM->DOS)->current_fitness[a][0]=(L_PARAM->DOS)->OP_INIT_VAL[a][0];
	    (L_PARAM->BGA)->TOT_OP_FIT+=(L_PARAM->DOS)->current_fitness[a][0];
	    }
	}
  
    for ( a=0; a< (L_PARAM->BGA)->STEP_FIT_OPERATOR; a++)
	for ( b=0; b<(L_PARAM->BGA)->MAXPOP; b++)
	    for ( c=0; c< 5; c++)
		(L_PARAM->DOS)->T[a][b][c]=0;

    }

    			    	     

int is_bi_parent ( int x, Parameter *PARAM )
    {
    return ((PARAM->DOS)->BI_PARENT_LIST[x]==2)?1:0;
    }



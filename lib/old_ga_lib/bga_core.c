/*SAGA*/
#define MAIN 1
#include "ga.h"

static G;
static G_PARAL;
static *previous_score;
static  current_score;

Global_structure * do_ga (Population *L_POP, Parameter *L_PARAM ,Statistic *STAT)
    {
    int a;
    char end_file[100];
    FILE *fp;
    FILE *fp2;
    int gen=0, gen2=0;
    int r;
    int mem1, mem2;
    Global_structure *GS;
    
    
    previous_score= vmalloc ( sizeof (int)*(L_PARAM->BGA)->MAX_GEN);
    	

    fprintf (L_PARAM->std1, "\nSEEDING.....");
    seeding (L_POP,L_PARAM);
    if ( ((L_PARAM->BGA)->OUTPUT_SCREEN)==1) printf ( "DONE\n");
    
    get_fitness (L_POP,L_PARAM);
    init_op_fitness (L_PARAM);
   


    for ( gen=0; gen< (L_PARAM->BGA)->MAX_GEN; gen++)
	{	
	if ((r=do_generation(L_POP,L_PARAM,STAT, gen))== -1 || r== -2 )
		{
		gen2=gen;
		gen=(L_PARAM->BGA)->MAX_GEN+1;
		}
	}
  (L_PARAM->BGA)->EVALUATE_ONLY=1;
   print_result_ga (L_POP,0, STAT,gen2, L_PARAM);	
   
   set_flag_file ( (L_PARAM->BGA)->NORMAL_END_FLAG_FILE, 1);  
   
 
   G=0;
   free ( previous_score);
   current_score=0;
   
   GS= vcalloc ( 1, sizeof ( Global_structure));
   GS->PARAM= L_PARAM;
   GS->POP=L_POP;
   GS->STAT = STAT;
   return GS;
   }

int do_generation(Population *L_POP,Parameter *L_PARAM, Statistic *STAT, int gen)
    {
    
    static Chromosome *A1, *A2;
    Chromosome *R;
    int p1, p2,child, site_co2;
    int rep_pop;
    
    int P_1=0;
    int P_2=1;
    int Child=2;    
    int OPER=3;
    int DELTA=4;

    int **ofspring;
    int os_sum;
    int **fatality;
    int fat_sum;
    int os_p1, os_p2, os_c;
		
    int sum,tot;

    int max, min;
    int a,b,c,d,op;
    int s1, s2, s;

    static int **m_t;
    static int **temp;
    int current_best_score=0;
    int current_worst_score=0;
    int current_best_i;
    int current_worst_i;

    int dup;
    int new;

    int n;
    int x,z;
    int stab;

    float sp1, sp2, sp,ps1, ps2, ref;
    static int last_tot,last_n;	
    char command[1000];
    static pop_link;
    pop_link=-1;
     
    if ( A1==NULL)
	{
	A1=declare_chromosome (L_PARAM);
	A2=declare_chromosome (L_PARAM);
	m_t=declare_int ( (L_PARAM->BGA)->MAXPOP, 5);
	
	temp = vcalloc ( 52, sizeof ( int*));
	for ( a=0; a< 52; a++)
	    temp[a]= vcalloc ( 2, sizeof ( int));
	}
    
    if ( (L_PARAM->BGA)->ANALYSE_OPERATOR==1)
	for ( a=0; a< (L_PARAM->DOS)->TOT_OPERATOR; a++)
	    (L_PARAM->DOS)->IOPUSE[a][0]=(L_PARAM->DOS)->IOPUSE[a][1]=(L_PARAM->DOS)->IOPUSE[a][2]=(L_PARAM->DOS)->IOPUSE[a][3]=0;
    
    

    for ( a=0; a< (L_PARAM->BGA)->MAXPOP; a++)
	for ( b=0; b< 5; b++)
	    m_t[a][b]=0;


    rep_pop= (int)(L_PARAM->BGA)->OVERLAPING_GEN_FACTOR;    


    if ( G==(L_PARAM->BGA)->STEP_FIT_OPERATOR)
	{
	evaluate_fitness_operator ( L_PARAM);
	fprintf (L_PARAM->std2, "[R]");
	save_current_operator ( L_PARAM, (L_PARAM->BGA)->experience_name, gen);
	G=0;
	}
    
    
    process_fitness ( L_POP, 0, (L_PARAM->BGA)->MAXPOP, (L_PARAM->BGA)->MAXIMISE);   
    rescale_fitness ( L_POP,0,(L_PARAM->BGA)->MAXPOP, (L_PARAM->BGA)->RESCALING_MODE, L_PARAM);
    preselect ( L_POP,0,(L_PARAM->BGA)->MAXPOP, (L_PARAM->BGA)->PRESELECT_MODE, L_PARAM);
    



    ofspring=duplicate_int ( (L_POP+0)->expected_os,(L_PARAM->BGA)->MAXPOP, 1);
    os_sum=return_sum_int ( ofspring,(L_PARAM->BGA)->MAXPOP,0);

    fatality=generate_fatality ( L_POP, L_PARAM);
     
    fat_sum=return_sum_int ( fatality, (L_PARAM->BGA)->MAXPOP, 0);
    do_statistic(L_POP, 0, STAT, (L_PARAM->BGA)->MAXPOP, gen,L_PARAM);
    

    previous_score[gen]=STAT->score_gen_best;

    print_result_ga (L_POP,0, STAT,gen, L_PARAM);
    fprintf (L_PARAM->std2, " {%d->%d}", last_n,last_tot); 
    
    
/******TERMINATION*********/

    if ( gen== ( L_PARAM->BGA)->MAX_GEN)
	    {
	    free_int ( ofspring, (L_PARAM->BGA)->MAXPOP);
	    free_int ( fatality, (L_PARAM->BGA)->MAXPOP);
	    return 1;
	    }
     if ( gen>0)
     	{    
     	if ( (L_PARAM->BGA)->OF_MODE==8 && previous_score[gen-1]==100000)
     	   {
     	   return -1;
     	   }
     	}
     	
     	
    stab=1;
    for ( b=gen-1; b>=0; b--)
	{d= (previous_score[b]==STAT->score_gen_best)?1:0;
	 b= ( d==1)?b:0;
	 stab+=d;
	 } 

    if ( stab== (L_PARAM->BGA)->MAX_STAB+1)	
	{
	free_int ( ofspring, (L_PARAM->BGA)->MAXPOP);
	free_int ( fatality, (L_PARAM->BGA)->MAXPOP);
	return -1;
	}

    /*print_select_ga ( L_POP,L_PARAM);*/   
   
    
    tot=0;
    new=1;
    n=0;
    

    
    while ( tot<(rep_pop))
	{
	n++;
	if ( n== rep_pop*100)
	    {
	    free_int ( ofspring, (L_PARAM->BGA)->MAXPOP);
	    free_int ( fatality, (L_PARAM->BGA)->MAXPOP);
	    return -2;
	    }

	os_p1=os_p2=os_c=0;

	

	m_t[tot][OPER]=op=select_operator (L_PARAM);
	
	

	if ( rep_pop< (L_PARAM->BGA)->MAXPOP)
		{
		m_t[tot][Child]=child=fatality_select(fatality,L_PARAM, fat_sum, &os_c, tot);
		fat_sum-= os_c;
		fatality[child][0]-=os_c;
		}
	else
		m_t[tot][Child]=child=tot;
	


	m_t[tot][P_1]=p1= wheel_select (ofspring, (L_PARAM->BGA)->MAXPOP, os_sum,0,0);
	os_p1= ( ((L_PARAM->BGA)->FACTOR_POP/2)<ofspring[p1][0])? ((L_PARAM->BGA)->FACTOR_POP/2):ofspring[p1][0];
	

	
	ofspring[p1][0]-=((L_PARAM->BGA)->SELECT_MODE==1)?(os_p1):0;
	os_sum -= ((L_PARAM->BGA)->SELECT_MODE==1)?(os_p1):0; 

	if (is_bi_parent(op, L_PARAM) )
	    {
	    while (((m_t[tot][P_2]=p2= wheel_select (ofspring,(L_PARAM->BGA)->MAXPOP, os_sum,0,0))==p1));
	    
	    if ( (L_PARAM->BGA)->SELECT_MODE==1)
		{
		os_p2= ( ((L_PARAM->BGA)->FACTOR_POP/2)< ofspring[p2][0])? (L_PARAM->BGA)->FACTOR_POP/2: ofspring [p2][0];
		ofspring[p2][0]-=((L_PARAM->BGA)->SELECT_MODE==1)?os_p2:0;
		os_sum -= ((L_PARAM->BGA)->SELECT_MODE==1)?os_p2:0;
		}
	    
	    extract_chrom (L_POP,0,m_t[tot][0],A1,L_PARAM);
	    extract_chrom (L_POP,0,m_t[tot][1],A2,L_PARAM);

	    make_operation (A1, A2,op, L_PARAM);  
	    s1=(int)evaluate_chromosome (A1, (L_PARAM->BGA)->OF_MODE, L_PARAM);
	    s2=(int)evaluate_chromosome (A2, (L_PARAM->BGA)->OF_MODE, L_PARAM);
	    
	    R= (a_better_than_b (s1,s2,(L_PARAM->BGA)->MAXIMISE))?A1:A2;
	    s= (a_better_than_b (s1,s2,(L_PARAM->BGA)->MAXIMISE))?s1:s2;
	    

	    if ( (L_PARAM->BGA)->ANALYSE_OPERATOR == 1)
		{if ( is_same (A1,(L_POP+0)->CHROMOSOME[m_t[tot][0]],L_PARAM) ||  is_same (A2,(L_POP+0)->CHROMOSOME[m_t[tot][0]],L_PARAM)||is_same (A1,(L_POP+0)->CHROMOSOME[m_t[tot][1]],L_PARAM) ||is_same (A2,(L_POP+0)->CHROMOSOME[m_t[tot][1]],L_PARAM))
	 		{(L_PARAM->DOS)->OPUSE[op][1]++;
			 (L_PARAM->DOS)->IOPUSE[op][1]++;
			} 
		a=a_better_than_b ( s, ((L_POP+0)->CHROMOSOME[m_t[tot][0]])->score, (L_PARAM->BGA)->MAXIMISE);
		b=a_better_than_b ( s, ((L_POP+0)->CHROMOSOME[m_t[tot][1]])->score, (L_PARAM->BGA)->MAXIMISE);
		(L_PARAM->DOS)->OPUSE [op][2]+= ( a==1 && b==1)?1:0; 
		(L_PARAM->DOS)->IOPUSE [op][2]+= ( a==1 && b==1)?1:0; 
		   
		(L_PARAM->DOS)->OPUSE[op][0]++;
		(L_PARAM->DOS)->IOPUSE[op][0]++;                                            
		}
	    } 
	    
	else
	    {	
            m_t[tot][P_2]=p1;
	    if ( (L_PARAM->BGA)->SELECT_MODE==1)
		{
		m_t[tot][P_2]=p1;
		p2=p1;
		os_p2=( ((L_PARAM->BGA)->FACTOR_POP/2)< ofspring[p2][0])? (L_PARAM->BGA)->FACTOR_POP/2: ofspring [p2][0];
		os_sum-=os_p2;
		ofspring[p2][0]-=os_p2;
		}	    

	    extract_chrom  (L_POP,0, m_t[tot][P_1],A1,L_PARAM);
	    make_operation (A1,A1,op,L_PARAM);
	    s=(int)evaluate_chromosome (A1, (L_PARAM->BGA)->OF_MODE,L_PARAM); 
	    R=A1;    
	    
	    if ( (L_PARAM->BGA)->ANALYSE_OPERATOR == 1)
		{
		if ( is_same ( A1,(L_POP+0)->CHROMOSOME[m_t[tot][P_1]], L_PARAM) ) 
		    	{(L_PARAM->DOS)->OPUSE[op][1]++;
			 (L_PARAM->DOS)->IOPUSE[op][1]++;
			} 

		a=a_better_than_b ( s, ((L_POP+0)->CHROMOSOME[m_t[tot][0]])->score, (L_PARAM->BGA)->MAXIMISE);

		(L_PARAM->DOS)->OPUSE [op][2]+= (a==1)?1:0; 
		(L_PARAM->DOS)->IOPUSE [op][2]+=(a==1)?1:0; 
		   
		(L_PARAM->DOS)->OPUSE[op][0]++;
		(L_PARAM->DOS)->IOPUSE[op][0]++;                                            
		}
	    }
	R->score=s;
	

	dup=(rep_pop<(L_PARAM->BGA)->MAXPOP && (L_PARAM->BGA)->ALLOW_DUPLICATE==0)? is_dup(R, L_POP, m_t,tot,L_PARAM):0;	 
	
	if ( dup==0)
		{
		if ( (L_PARAM->BGA)->ANALYSE_OPERATOR==1)
		   {(L_PARAM->DOS)->OPUSE[op][3]++;
		    (L_PARAM->DOS)->IOPUSE[op][3]++;     
		   }
		insert_chrom (R,L_POP,1, m_t[tot][Child],L_PARAM);
		sp1= (L_POP+0)->raw_fitness[m_t[tot][P_1]][0];
		sp2= (L_POP+0)->raw_fitness[m_t[tot][P_2]][0];
		sp= (a_better_than_b (sp1,sp2,(L_PARAM->BGA)->MAXIMISE))?sp1:sp2; 
		ref=((L_PARAM->BGA)->DYNAMIC_PARAM==0)?(int)STAT->score_best:sp; 
			
		m_t[tot][DELTA]=((a_better_than_b( s,ref,(L_PARAM->BGA)->MAXIMISE)))?int_fabs ((int)(ref-s)):0;

		if ( current_best_score==current_worst_score==0)
		    current_best_score=current_worst_score=s;
		if ( a_better_than_b (s,current_best_score,(L_PARAM->BGA)->MAXIMISE))
		    {current_best_score=s;
		    current_best_i=m_t[tot][Child];
		    }
		if ( a_better_than_b (current_worst_score,s,(L_PARAM->BGA)->MAXIMISE))
		    {current_worst_score=s;
		    current_worst_i=m_t[tot][Child];
		    }
		tot++;
		new=1;
		}
	    else
		{
		os_sum+=os_p1+os_p2;
		fat_sum+=os_c;
		fatality[m_t[tot][Child]][0]+=os_c;
		ofspring[m_t[tot][P_1]][0]+= os_p1;
		ofspring[m_t[tot][P_2]][0]+= os_p2;
		}
	 }
	    	
   
    
   
    
    last_tot=tot;
    last_n =n;
    for ( a=0; a< tot; a++)
	copy_individual ( L_POP, 1, m_t[a][Child],0,m_t[a][Child],L_PARAM);
    
    
    check_end_signal ( (L_PARAM->BGA)->FINISH_FLAG_FILE, L_PARAM);
    	
    if ( (L_PARAM->BGA)->PARALLEL_GA==1)
    	{
    	G_PARAL++;
    	if ( G_PARAL==(L_PARAM->BGA)->PARALLEL_EXCHANGE)
    		{
    		G_PARAL=0;
    		save_population ( L_POP,0, L_PARAM);
    		wait_flag ((L_PARAM->BGA)->POP_READ_FLAG_FILE, L_PARAM);
    		read_population ( L_POP,0, L_PARAM);
    		}
    	}	 
    			    	
    for ( a=0; a< (L_PARAM->BGA)->MAXPOP; a++)
	{
	(L_PARAM->DOS)->T[G][a][0]=-1;
	}

    for ( a=0; a< tot; a++)
	for ( b=0; b< 5; b++)
	    (L_PARAM->DOS)->T[G][m_t[a][Child]][b]=m_t[a][b];

    for ( a=0; a< (L_PARAM->BGA)->MAXPOP; a++)
	if ( (L_PARAM->DOS)->T[G][a][0]==-1)
	    {
	    (L_PARAM->DOS)->T[G][a][P_1]=(L_PARAM->DOS)->T[G][a][P_2]=(L_PARAM->DOS)->T[G][a][Child]=a;
	    (L_PARAM->DOS)->T[G][a][OPER]=-1;
	    (L_PARAM->DOS)->T[G][a][DELTA]=0;
	    }	

    for ( a=0; a< (L_PARAM->BGA)->MAXPOP; a++)
	{
	if ( (L_PARAM->DOS)->T[G][a][OPER]>=0 && (L_PARAM->DOS)->T[G][a][DELTA]>0 )
	    report ( a, (L_PARAM->DOS)->T[G][a][DELTA],L_PARAM,G);
	
	}

    if ( (L_PARAM->BGA)->ANALYSE_OPERATOR==1)
	save_analyse_operator ( gen, L_PARAM);	
     		
    free_int ( ofspring, (L_PARAM->BGA)->MAXPOP);
    free_int ( fatality, (L_PARAM->BGA)->MAXPOP);
    G++;
    
    return 1;
    
    }

void process_fitness ( Population*L_POP, int G, int total_pop,int MODE)
    {
    /*MODE==0:invert, MODE=1:keep non inverted*/
    float min, max;
    int a;

    min=return_min_float ( (L_POP+G)->raw_fitness, total_pop, 0);
    max=return_max_float ( (L_POP+G)->raw_fitness, total_pop, 0); 
    
    for ( a=0; a< total_pop; a++)
	{
	if ( MODE==0)
	    (L_POP+G)->fitness[a][0]=max-(L_POP+G)->raw_fitness[a][0]+min;
	else
	    (L_POP+G)->fitness[a][0]=(L_POP+G)->raw_fitness[a][0];

	
	}
    
    min=return_min_float ( (L_POP+G)->fitness, total_pop, 0);
    for ( a=0; a< total_pop; a++)
	{
	if (min<0)
		{
		(L_POP+G)->fitness[a][0]-=min;
		}
	(L_POP+G)->ranked[a][0]=(L_POP+G)->fitness[a][0];
	(L_POP+G)->ranked[a][1]=a;
 	(L_POP+G)->age[a][0]++;
     	}
     
     }

void get_fitness(Population*L_POP,Parameter *L_PARAM)
    {
    int a,b;
    static Chromosome *A;
	
    for ( a=0; a< (L_PARAM->BGA)->MAXPOP; a++)
	{
	(L_POP)->raw_fitness [a][0]=((L_POP)->CHROMOSOME[a])->score=(int)evaluate_chromosome ((L_POP)->CHROMOSOME[a], (L_PARAM->BGA)->OF_MODE, L_PARAM);
	}
    }

int** generate_fatality ( Population *LPOP,Parameter *L_PARAM)
    {
    int a, b, min, max;
    int **fatality;

    fatality=declare_int ( (L_PARAM->BGA)->MAXPOP,2);
    
    for ( a=0; a< (L_PARAM->BGA)->MAXPOP; a++)
	    {
	    if ( (L_PARAM->BGA)->REPLACE_MODE==0)
	    	    fatality[a][0]=(LPOP+0)->expected_os[a][0];
	    else
		    fatality[a][0]=(int)(LPOP+0)->fitness[a][0];
	    fatality[a][1]=a;
	    }
	       

   
    max=return_max_int (fatality,(L_PARAM->BGA)->MAXPOP,0);
    min=return_min_int (fatality,(L_PARAM->BGA)->MAXPOP,0);
    inverse_int ( fatality,(L_PARAM->BGA)->MAXPOP,0, max, min);	  
	
     	       
    for  (a=0; a<(L_PARAM->BGA)->MAXPOP; a++)
	fatality[a][0]++;

    if ( (L_PARAM->BGA)->REPLACE_MODE ==2)
    	{
    	
	sort_int ( fatality, 2, 0, 0, (L_PARAM->BGA)->MAXPOP-1);
      
        }
    return fatality;    
    }		 



    
void check_end_signal ( char *fname, Parameter *PARAM)
	{
	FILE *fp;
	int a;
	
	
	if ( strcmp (fname, "/dev/null")==0)return;
	if ((fp=fopen( fname, "r"))==NULL)
		return;
	else
		{	
		fscanf ( fp, "%d", &a);
		fclose ( fp);
			if ( a==1)
				{
				set_flag_file ( (PARAM->BGA)->NORMAL_END_FLAG_FILE, 1);
				exit(0);
				}
			
		}		
	}

 
/*******************************************************************************************************/
/*                                                                                                     */
/*                                                                                                     */
/*                       OPERATORS                                                                     */
/*                                                                                                     */
/*                                                                                                     */
/*                                                                                                     */
/*******************************************************************************************************/
/* format of the file : 

D  OP_NAME number parents(1,2) proba).
*/

void make_operation ( Chromosome *A1, Chromosome *A2, int op, Parameter *L_PARAM)
    {
    int DEBUG=1;
    static int mem1;
    static int mem2;
    static int list[55];


#if MEMTRACK_OP
    mem1=(int)Mem_Used();
#endif

    (L_PARAM->DOS)->OP_FUNCTION[op](A1, A2, L_PARAM);

#if MEMTRACK_OP
    mem2=(int)Mem_Used();
    if ( mem1!=mem2)
    {
    printf ( "\n\tOPERATOR %s[%d] USED %d amount of memory",(L_PARAM->DOS)->OP_NAME[op], op, mem2-mem1);
    } 
#endif

   
    } 

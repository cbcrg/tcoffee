Alignment **  scan_for_repeats_2     (int **L, int ne,int scale,int gep,int maximise, Sequence *S,int sw_threshold,int sw_min_len, int sw_min_z_score )
           {
	   Alignment *A;
	   Alignment **aln_list=NULL;
	   int **L1;
	   int ne1;
	   int **start_list;
	  
	   L1=duplicate_int ( L, read_size_int ( L[-1],0),LIST_N_FIELDS);
	   ne1=ne;

	   while (1)
	         {
		 A=get_one_repeat_type ( L1, ne1, scale, gep,maximise, S,sw_threshold, sw_min_len, sw_min_z_score);
		 
		 if (A!=NULL && A->nseq>1 && A->len_aln > sw_min_len && A->score_aln> sw_min_z_score)
		    {
		    A->cpu=(int)get_time();
		    A=get_best_nol_local_aln(A, L1, ne1, scale, gep, sw_threshold, sw_min_len, sw_min_z_score, NON_GREEDY);		   
		    output_Alignment_with_res_number ( A,stderr);
		    L1=mask_list_with_aln (A,0, A->len_aln, L1, &ne1, UNDEFINED);
		    aln_list=realloc_aln_array (aln_list,1);
		    copy_aln (A,aln_list[(aln_list[-1])->nseq-1]);
		    free_Alignment (A);			 
		    /*get_one_repeat_type (NULL, 0,0,0,0, NULL,0,0,0);*/
		    }
		 else if ( A==NULL) break;
		 else
		     {
		     free_Alignment (A);			 		     
		     break;
		     }
		 }
	get_one_repeat_type (NULL, 0,0,0,0, NULL,0,0,0);
	L1=free_int (L1,-1);
        return aln_list;
	}



Alignment *  get_one_repeat_type_2  (int **L1, int n_elements_in,int scale,int penalty,int maximise, Sequence *S,int sw_t, int sw_l, int sw_z)
           {	
	   Alignment *A=NULL;
	   Alignment *RESULT=NULL;
	   static int series;
	   Alignment **aln_list;
	   Alignment **new_aln_list;
	   Alignment * ref_aln;
	   int *nseq_array;
	   int **L;
	   int n_elements;	   
	   int aln, n_aln;
	   int best_aln, best_score;	   
	   int a, b, c, d;
	   int cont=1;
	   int count=1;
	   int **start_list;
	   int **seq_coor;
	   int z;


	   
	   int iscale, ipenalty;
	   static int ** sub_matrix;
	   /*Returns:
	     A:The alignmnent containing the repeat 
	       NULL if no repeat was found
	   */

	   seq_coor=declare_int ( 10000, 3);
	   if ( L1==CLEAN_FUNCTION)
	      {
	      add_seq2aln (CLEAN_FUNCTION,0,0,0,NULL,0,NULL,0,0,0);	      
	      return NULL;
	      }
	   else if ( n_elements_in==0)
	      {
	      return NULL;
	      }
	   else
	      {
	      n_elements=n_elements_in;	   
	      nseq_array=vcalloc ( 100000, sizeof (int));	 	     
	      L=mirror_list (L1, &n_elements); 	
   	      start_list=extract_best_column ( L, n_elements, S, WE);
	      
	      for ( a=1; a<=start_list[0][0]; a++)
	          {
		  seq_coor[a-1][0]=0;
		  seq_coor[a-1][1]=start_list[a][0];
		  seq_coor[a-1][2]=strlen(S->seq[0])-start_list[a][0]+1;
		  A=build_progressive_nol_aln_with_seq_coor (L, n_elements,0,0,S, seq_coor,a );
		  
		  if ( A->nseq>1)
		     {
		     
		     z=evaluate_seq_z_score (A, 0, A->len_aln);
		     A->score_aln=z;
		     fprintf ( stderr, "\nSCORE=%d MIN=%d\n", z, sw_z);
		     /*  output_Alignment_with_res_number (A, stderr);*/

		     if (z< sw_z)break;
		     else
			 RESULT=copy_aln(A, RESULT);
		     }
		  else
		      RESULT=copy_aln(A, RESULT);
		  }
	      free_Alignment (A);
	      return RESULT;
	      }
	   }

Alignment ** add_seq2aln_2 (int **L,int n_elements,int scale,int penalty, Alignment *IN,int maximise,Sequence  *S,int sw_t, int sw_l, int sw_z)
           {	
	   int *n_groups;
	   int **group_list;
	   Alignment **aln_list;
	   static Alignment *A;
	   int n_aln=0;
	   int a,b;
	   int first=0;    
	   int len;
	   static int series=1;
	   
	   if ( A==NULL)
	       {
	       A=declare_Alignment (NULL);
	       A->S=S;
	       }
	   if ( L==CLEAN_FUNCTION)
	      {
	      get_best_local_aln ( CLEAN_FUNCTION, NULL, 0,0,0,0,0,0, GREEDY);
	      free_Alignment (A);
	      A=NULL;
	      return NULL;
	      }
	   
	   if ( IN==NULL)
	      {
	      
	      aln_list=declare_aln_array(1);    	      	    
	      
	      A=realloc_alignment2(A, 1, strlen (S->seq[0])+1);
	      A->S=S;
	      A->nseq=1;
	      sprintf ( A->seq_al[0], "%s", S->seq[0]);
	      sprintf (A->name[0], "%s_%d_1", S->name[0],series);
	     
	      A->len_aln=strlen ( A->seq_al[0]);
	      A->order[1][0]=0;
	      copy_aln ( A, aln_list[0]);
	      return aln_list;
	      }
	   else if (IN->finished==1)return NULL;
	   else
	      {
	      copy_aln ( IN, A);
	      A=realloc_alignment2 ( A, A->nseq+1,MAX(strlen ( S->seq[0])+1, A->len_aln+1));
	      n_groups=vcalloc ( 2, sizeof (int));
	      group_list=declare_int (2,A->nseq+1);
	      
	      n_groups[0]=A->nseq;
	      for ( a=0; a<A->nseq; a++)group_list[0][a]=a;
	      n_groups[1]=1;
	      group_list[1][0]=A->nseq;
	      sprintf (A->name[A->nseq], "%s_%d_%d",S->name[0],series,A->nseq+1);
	      sprintf (A->seq_al[A->nseq], "%s",S->seq[0]);
	      A->order[A->nseq][0]=0;
	      A->order[A->nseq][1]=0;
	      A->nseq++;
	      sw_pair_wise ( A, penalty, scale, n_groups, group_list, L, n_elements,maximise);
		      
	      /*		
	      if ( A->nseq==2)
		  {
		  for ( b=0; b< 100; b++)
	              {
		      sprintf (A->seq_al[0],  "%s",S->seq[0]);		     
		      sprintf (A->seq_al[1],  "%s",S->seq[0]);
		      A->order[0][0]=A->order[1][0]=A->order[0][1]=A->order[1][1]=0;
		      
		      sw_pair_wise ( A, penalty, scale*b, n_groups, group_list, L, n_elements,maximise);
		      fprintf( stderr, "\nSCORE=%.2f\n", (float)evaluate_seq_z_score(A, 0, A->len_aln));
		      output_Alignment_with_res_number(A,stderr);
		      }
		  }
	      */

	      free (n_groups);
	      free_int ( group_list,-1);	      
	      
	      aln_list=get_best_local_aln ( A, L, n_elements, scale, penalty,sw_t, sw_l, sw_z, NON_GREEDY);	     
	      
	      series+=(aln_list==NULL)?1:0;
	      return aln_list;
	      
	      }	   	   
	   }
Alignment * get_best_ol_local_aln ( Alignment *IN, int **L, int ne, int gop, int gep,int sw_t,int sw_l, int sw_z, int mode)
     {
     static Alignment *A;
     int a, b, c, d, e, f;
     int left, right;
     int best_a, score, best_score;

     for (best_a=0, best_score=0, score=0, a=0; a< IN->len_aln; a++)
         {
	 left =get_nol_aln_border(IN, a, GO_LEFT );
	 right=get_nol_aln_border(IN, a, GO_RIGHT);
	 
	 A=copy_aln(IN, A);
	 A=extract_aln (A, left, right);
	 score=evaluate_aln(A, 0, A->len_aln, L, ne, gop, gep, WE);
	 if ( score>best_score)
	    {
	    best_score=score;
	    best_a=a;
	    }
	 }
     left =get_nol_aln_border(IN, best_a, GO_LEFT );
     right=get_nol_aln_border(IN, best_a, GO_RIGHT);
	 
     A=copy_aln(IN, A);
     A=extract_aln (A, left, right);
     return A;
     }

Alignment * get_best_nol_local_aln ( Alignment *IN, int **L, int ne, int gop, int gep,int sw_t,int sw_l, int sw_z, int mode)
     {
     static Alignment *A;
     int a, score, best_score, best_a,ref;
     int end;

       
     return trunkate_local_aln(IN);
     if ( A->len_aln==IN->len_aln)return IN;
     else
         {
	 end =IN->len_aln;
	 for (best_score=0, best_a=0, a=0; a< IN->len_aln;)
             {

	     if ( a==IN->len_aln)break;
	     
	     A=copy_aln( IN, A);
	     A=extract_aln(A, a, A->len_aln);
	     A=trunkate_local_aln(A);
	     if (A->len_aln+a !=end)
	        {
		end=A->len_aln+a;	
		score=evaluate_aln(A, 0, A->len_aln, L, ne, gop, gep, WE);
		if ( score > best_score)
	           {
		   best_score=score;
		   best_a=a;
		   }
		}
	     a+=A->len_aln;
	     }
	 
	 fprintf ( stderr, "\nbest_a=%d",a);
	 IN=extract_aln (IN, best_a, IN->len_aln);
	 return trunkate_local_aln(IN);
	 }
     
     }
int old_pair_wise (Alignment *A, int gep, int gop,int*ns, int **l_s,int **L, int ne, int maximise )
	{
	/*********************************************************************
	*	makes DP between the the ns[0] sequences and the ns[1] sequences in A
	*
	*	for MODE, see the function  get_dp_cost
	**********************************************************************/
	
	int a, b, d, f, i, j;
	int *list;
	int *cc;
	int *dd, *ddg;
	int lenal[2], len;
	int nseq;
	int hD, gD, hI, gI,t, c,s, e,eg, k, ch,g,h;
	int sub;
	
	int fop;
	int score=0;
	
	static int **pos1;
	static int **pos0;
	
	static char **al, **aln;
	static int **trace;
	static max_len;

	int ala, alb,LEN;
	char *buffer;
	char *char_buf;
	
	WEIGHT_FIELD=WE;

	nseq=(A->S)->nseq;
	lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	
	len= ( lenal[0]>lenal[1])?lenal[0]:lenal[1];
	
	
	
/*DO MEMORY ALLOCATION FOR NW DP*/
nseq=A->nseq;	
           if ( len>=max_len)
		{
		if ( trace!=NULL)
			{
			free_int(trace,max_len);
			free_char (al, nseq);
			}
		trace=declare_int ( len+1, len+1);
		al=declare_char (nseq, len*2+1);
		max_len=len+1;
		
		}
	buffer=calloc ( 2*max_len+1, sizeof (char));
	char_buf= calloc ( 2*max_len+1, sizeof (char));	
	dd=calloc (max_len+1, sizeof (int));
	cc=calloc (max_len+1, sizeof (int));
	ddg=calloc (max_len+1, sizeof (int));
	
	pos0=aln2pos_simple ( A,-1, ns, l_s);
	
	
	
		/*
		0(s)   +(dd)
		  \      |
		   \     |
		    \    |
		     \   |
		      \  |
		       \ |
		        \|
		-(e)----O
		*/ 
		       
	for (i=0;i<=lenal[0];i++)
		{trace[i][0]= 1;
		}
	for (j=0;j<=lenal[1]; j++)
		{trace[0][j]=-1;
		}
	
	g=gop;
	h=gep;
	
	trace[0][0]=1;	
	cc[0]=0;
	
	t=g;
	for ( j=1; j<=lenal[1]; j++)
		{
		cc[j]=t=t+h;
		dd[j]=t+g;
		}
		
	
	t=g;			
	for (i=1; i<=lenal[0];i++)
			{
			s=cc[0];
			cc[0]=c=t=t+h;
			e=t+h;	
			for (eg=0,j=1; j<=lenal[1];j++)
				{
				sub=get_dp_cost (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,L, ne);
				if ( a_better_than_b ( e, c+g, maximise))
					eg++;
				else 
					eg=1;	
					
				e=best_of_a_b (e, c+g, maximise)+h;
				
				if ( a_better_than_b ( dd[j], cc[j]+g, maximise))
					ddg[j]++;
				else 
					ddg[j]=1;
					
				dd[j]=best_of_a_b( dd[j], cc[j]+g, maximise)+h;
				
				c=best_int(3,maximise,&fop, e, s+sub, dd[j], maximise);
			
				fop-=1;
				s=cc[j];
				cc[j]=c;
				
				if ( fop<0)
					{trace[i][j]=fop*eg;
					}
				else if ( fop>0)
					{trace[i][j]=fop*ddg[j];
					 }
				else 
					{trace[i][j]=0;	
					 }	
				fop= -2;
				}
			}


				
	score=c;	
	
	i=lenal[0];
	j=lenal[1];
	ala=alb=0;
		
	while (i>=0 && j>=0 && ((i+j)!=0))
			{
			if ( i==0)
				k=-1;
			else if ( j==0)
				k=1;
			else if ( j==0 && i==0)
				k=1;	
			else	
				k=trace[i][j];
				
				
				
			if (k==0)
				{
				al[0][ala++]=1;
				al[1][alb++]=1;
				
				i--;
				j--;
				}
			else if (k>0)
				{
				for ( a=0; a< k; a++)
					{
					al[0][ala++]=1;
					al[1][alb++]=0;
					i--;
					}
				}
			else if (k<0)
				{
				for ( a=0; a>k; a--)
					{
					al[0][ala++]=0;
					al[1][alb++]=1;
					j--;
					}
				}
			}
	
	
	LEN=ala;	
	c=LEN-1;  
	if ( A->declared_len<=LEN)realloc_alignment ( A, 2*LEN);  
	aln=A->seq_al;
	
	for ( c=0; c< 2; c++)
	    {for ( a=0; a< LEN; a++)
		buffer[a]=al[c][a];
	    b=0;
	    for ( a=LEN-1; a>=0; a--)
		al[c][b++]=buffer[a];
	    }
	
        for ( c=0; c< 2; c++)
	    {
	    for ( a=0; a< ns[c]; a++) 
		{
		ch=0;
		for ( b=0; b< LEN; b++)
		    {
		   
		    if (al[c][b]==1)
			char_buf[b]=aln[l_s[c][a]][ch++];
		    else
			char_buf[b]='-';
		   }
		char_buf[b]='\0';
		sprintf (aln[l_s[c][a]],"%s", char_buf);
	        }
	     }
	
	A->len_aln=LEN;
	A->nseq=ns[0]+ns[1];

	 
	free (cc);
	free (dd);		
	free (ddg);
	free (buffer);
	free (char_buf); 
	return score;
	}
Alignment *  get_one_repeat_type  (int **L1, int n_elements_in,int gop,int gep,int maximise, Sequence *S,int sw_t, int sw_l, int sw_z)
           {	
	   Alignment *A;
	   static int series;
	   Alignment **aln_list;
	   Alignment **new_aln_list;
	   Alignment * ref_aln;
	   int *nseq_array;
	   int **L;
	   int n_elements;	   
	   int aln, n_aln;
	   int best_aln, best_score;	   
	   int a, b, c, d;
	   int cont=1;
	   int count=1;
	   int igop, igep;
	   static int ** sub_matrix;
	   int ** residue_list;

	   if ( sub_matrix==NULL)sub_matrix=read_matrice ( "pam250mt");

	   if ( L1==CLEAN_FUNCTION)
	      {
	      add_seq2aln ( NULL,0,0,0,NULL,0,NULL,NULL,0,0,0);
	      
	      return NULL;
	      }
	   else if ( n_elements_in==0)
	      {
	      return NULL;
	      }
	   else
	      {
	      n_elements=n_elements_in;	   
	      nseq_array=vcalloc ( 100000, sizeof (int));	 	     
	      L=mirror_list (L1, &n_elements); 	
	      A=declare_Alignment (NULL);	   
	      A->S=S;
	      	      
	      igop=gop;
	      igep=gep;
	      
	      residue_list=extract_best_column(L, n_elements, S, WE);
	      
	      aln_list=add_seq2aln ( L,n_elements, gop, gep,NULL, maximise,S,residue_list,sw_t, sw_l, sw_z);	   
	      fprintf ( stderr, "\nREPEAT SET %d\n", (++series));
	      cont=1;
	      while (cont==1)
		   {	
		   cont=0;
		   ref_aln=((aln_list[0])->nseq>1)?aln_list[0]:NULL;
		   gop=igop;
		   gep=igep;		   
		   new_aln_list=add_seq2aln (L, n_elements,  gop, gep,aln_list[0], maximise,S,residue_list,sw_t, sw_l, sw_z);
		   if ( new_aln_list!=NULL)
		      {
		      cont=1;
		      copy_aln (new_aln_list[0], aln_list[0]);
		      fprintf ( local_stderr, "\t%s_%d [scale=%5d] [gep=%5d] [Score=%3d][start=%4d][End=%4d]\n", S->name[0],(aln_list[0])->nseq,gop, gep,(aln_list[0])->score_aln,(aln_list[0])->order[(aln_list[0])->nseq-1][1],(aln_list[0])->order[(aln_list[0])->nseq-1][1]+(aln_list[0])->len_aln);		      
		      }	
		   
		   else if ( new_aln_list==NULL)
		       (aln_list[0])->finished=1;			 
		   
		   count++;
		   }	     
	      copy_aln ( aln_list[0], A);
	      free ( nseq_array);
	      free_aln_array( aln_list);	      
	      return A;
	      }
	   }
Alignment ** add_seq2aln (int **L,int n_elements,int gop,int gep, Alignment *IN,int maximise,Sequence  *S,int **residue_list,int sw_t, int sw_l, int sw_z)
           {	
	   int *n_groups;
	   int **group_list;
	   Alignment **aln_list;
	   static Alignment *A;
	   int n_aln=0;
	   int a;
	   int first=0;    
	   int len;
	   static int series=1;
	   
	   if ( A==NULL)
	       {
	       A=declare_Alignment (NULL);
	       A->S=S;
	       }
	   if ( L==CLEAN_FUNCTION)
	      {
	      get_best_local_aln ( CLEAN_FUNCTION, NULL, 0,0,0,0,0,0, GREEDY);
	      free_Alignment (A);
	      A=NULL;
	      return NULL;
	      }
	   
	   if ( IN==NULL)
	      {
	      
	      aln_list=declare_aln_array(1);    	      	    
	      
	      A=realloc_alignment2(A, 1, strlen (S->seq[0])+1);
	      A->S=S;
	      A->nseq=1;
	      sprintf ( A->seq_al[0], "%s", S->seq[0]);
	      sprintf (A->name[0], "%s_%d_1", S->name[0],series);
	     
	      A->len_aln=strlen ( A->seq_al[0]);
	      A->order[1][0]=0;
	      copy_aln ( A, aln_list[0]);
	      return aln_list;
	      }
	   else if (IN->finished==1)return NULL;
	   else
	      {
	      copy_aln ( IN, A);
	      A=realloc_alignment2 ( A, A->nseq+1,MAX(strlen ( S->seq[0])+1, A->len_aln+1));
	      n_groups=vcalloc ( 2, sizeof (int));
	      group_list=declare_int (2,A->nseq+1);
	      
	      n_groups[0]=A->nseq;
	      for ( a=0; a<A->nseq; a++)group_list[0][a]=a;
	      n_groups[1]=1;
	      group_list[1][0]=A->nseq;
	      sprintf (A->name[A->nseq], "%s_%d_%d",S->name[0],series,A->nseq+1);
	      sprintf (A->seq_al[A->nseq], "%s",S->seq[0]);
	      A->order[A->nseq][0]=0;
	      A->order[A->nseq][1]=0;
	      A->nseq++;	      	     		     
	      
	      aln_list=trim_local_aln(A,L, n_elements, residue_list,S);
	      if ( aln_list==NULL)
		  {
		  fprintf ( stderr, "\nEND?");
		  
		  return NULL;
		  }
	      for ( a=0; a< 3; a++)compress_aln(aln_list[a]);
	      copy_aln(aln_list[0], A);
	      A=aln_cat (A, aln_list[1]);
	      A=aln_cat (A, aln_list[2]);

	      
	      
	     
	      output_Alignment_with_res_number(A, stderr);
	      
	      /*
	      pair_wise ( aln_list[0], gep, gop, n_groups, group_list, L, n_elements,maximise);
	      pair_wise ( aln_list[2], gep, gop, n_groups, group_list, L, n_elements,maximise);
	      */

	      aln_list[0]=aln_cat ( aln_list[0], aln_list[1]);
	      aln_list[0]=aln_cat ( aln_list[0], aln_list[2]);
	      A=copy_aln ( aln_list[0], A);
	      

	      
	      compress_aln(A);
	      A=remove_end(A);
	      pair_wise (A,0,0, n_groups, group_list, L, n_elements,maximise);
	      
	      
	      output_Alignment_with_res_number(A, stderr);
	      
	      verify_aln(A, S, "Final");
	      free_aln_array(aln_list);
	      
	      aln_list=declare_aln_array(1);
	      copy_aln (A, aln_list[0]);
	      free (n_groups);
	      free_int ( group_list,-1);	      
	      series+=(aln_list==NULL)?1:0;
	      return aln_list;
	      
	      }	   	   
	   }
Alignment ** get_best_local_aln_2 ( Alignment *IN,int **list,int ne,int gop, int gep, int sw_t, int sw_l, int sw_z, int greedy)
     {
     int a, b, c;
     int best_a,score;
     int len;
     static Alignment *A;
     static Alignment *B;     
     int prev_end, end;
     Alignment **a_list=NULL;
     int best_score;
     static int first;
     static int ref_score;
     static int tot_ref;
     static int n_ref;
     int first_p, last_p, best_p;
     int best_end;
     int * n_groups;
     int **group_list;
     int *start_list, start_score, start_best;
     int **pos,p1,p2;

     n_groups  =vcalloc ( 2, sizeof (int));
     group_list=declare_int ( 2, 1);

     
     if ( IN==NULL)
        {
	free_Alignment (A);
	A=NULL;
	free_Alignment (B);
	B=NULL;
	first=0;
	ref_score=tot_ref=n_ref=0;
	return NULL;
	}
     else
        {	
	A=copy_aln (IN,A );	
	if ( first==0)
	   {
	   start_list=do_analyse_list(list, ne);
	   pos=duplicate_int(aln2pos_simple (A, A->nseq), A->nseq, A->len_aln+1);
	   for (a=0, start_best=0, best_a=0; a< A->len_aln; a++)
	       {
	       p1=pos[0][a];
	       p2=pos[1][a];
	       
	       if ( p1>0 && p2>0)start_score=MIN(start_list[p1],start_list[p2]);
	       else start_score=0;
	       if ( start_score>start_best)
	          {
		  start_best=start_score;
		  best_a=a;
		  }
	       }
	   fprintf ( stderr, "\nBest_a=%d best_score=%d #1=%d #2=%d", best_a,start_best, pos[0][best_a], pos[1][best_a]);

	   A=extract_aln (A,best_a,A->len_aln);
	   A=trunkate_local_aln (A);
	   free_int (pos,-1);
	   free(start_list-1);
	   first=1;
	   }
        A=trunkate_local_aln(A);
	
	if ( a_list==NULL)a_list=realloc_aln_array ( a_list, 1);
	copy_aln ( A,a_list[(a_list[-1])->nseq-1]);

	score=evaluate_seq_z_score(a_list[0],0, (a_list[0])->len_aln);
	
	output_Alignment_with_res_number(A, stderr);
	fprintf ( stderr, "\nSCORE=%d", score);


	ref_score=(ref_score==0)?score:ref_score;
	if ((score*100)/ref_score>sw_t && (a_list[0])->len_aln > sw_l && score>sw_z )
	      {
	      (a_list[0])->score_aln=score;
	      copy_aln (a_list[0], B);	      
	      }
	else
	      {	
	     
	      free_aln_array(a_list);
	      a_list=NULL;
	      }
	ref_score=score;
	return a_list;
	}
     fprintf ( stderr, "\nUNEXPECTED CASE IN GET_BEST_LOCAL_ALN");
     exit(0);
     }     

Alignment ** get_best_local_aln ( Alignment *IN,int **list,int ne,int gop, int gep, int sw_t, int sw_l, int sw_z, int greedy)
     {
     int a, b, c;
     int best_a,score;
     int len;
     static Alignment *A;
     static Alignment *B;     
     int prev_end, end;
     Alignment **a_list=NULL;
     int best_score;
     static int first;
     static int ref_score;
     static int tot_ref;
     static int n_ref;
     int first_p, last_p, best_p;
     int best_end;
     int * n_groups;
     int **group_list;

     n_groups  =vcalloc ( 2, sizeof (int));
     group_list=declare_int ( 2, 1);

     if ( A==NULL)A=declare_Alignment (NULL);
     if ( B==NULL)B=declare_Alignment (NULL);

     if ( IN==NULL)
        {
	free_Alignment (A);
	A=NULL;
	free_Alignment (B);
	B=NULL;
	evaluate_aln (CLEAN_FUNCTION, 0,0,NULL, 0,0,0,0);	
	get_start_point (CLEAN_FUNCTION, NULL,0, 0,0,0,NULL, NULL);
	first=0;
	ref_score=tot_ref=n_ref=0;
	return NULL;
	}
     else
        {
	first=1;
	copy_aln (IN,A );	
	if ( first==0)
	   {
	   get_start_point(IN,list, ne, gop, gep,sw_t, &first_p, &last_p);
	   A=extract_aln (A,first_p,last_p);
	   first=1;
	   }
	trunkate_local_aln(A);
	if ( a_list==NULL)a_list=realloc_aln_array ( a_list, 1);		  
	copy_aln ( A,a_list[(a_list[-1])->nseq-1]);

        if( a_list==NULL)
	  {
	  score=0;
	  first=0;
	  n_ref=tot_ref=ref_score=0;
	  }
       else
          {
	   score=evaluate_seq_z_score(a_list[0],0, (a_list[0])->len_aln);
	   ref_score=(ref_score==0)?score:ref_score;
	   if ((score*100)/ref_score>sw_t && (a_list[0])->len_aln > sw_l && score>sw_z )
	      {
	      (a_list[0])->score_aln=score;
	      copy_aln (a_list[0], B);	      
	      }
	   else
	      {	
	     
	      free_aln_array(a_list);
	      a_list=NULL;
	      }
	   ref_score=score;
	   }
	return a_list;
	}
     fprintf ( stderr, "\nUNEXPECTED CASE IN GET_BEST_LOCAL_ALN");
     exit(0);
     }     
int  get_start_point ( Alignment *A,int **list, int ne,int gop, int gep, int T, int *first_p, int *last_p)
     {
     int a, b, c;
     static int *entry;
     static int **start_index;
     static int **end_index;
     int s1, s2, r1, r2;
     int **r;
     int **pos;
     int score=0;
     int tot=0;
     int ref;
     int **score_list;
     int **score_list2;
     

     
     int best_p, left_limit, right_limit;
     int n_it=10;
     int *list2;
     int right, left;
     int worst_p;
     
     
     if ( A==CLEAN_FUNCTION)
        {
	free_int ( start_index, -1);
	free_int ( end_index, -1);
	start_index=NULL;
	end_index=NULL;
	search_in_list_constraint (NULL,0,NULL,0, &start_index, &end_index);	       
	return -1;
	}
    
     score_list =declare_int ( A->len_aln+1,n_it+1);
     score_list2=declare_int ( A->len_aln+1,n_it+1);

     if ( entry==NULL) entry=calloc (LIST_N_FIELDS , sizeof (int));
     pos=aln2pos_simple(A, A->nseq);
     for ( c=0; c<A->len_aln; c++)
	 {
	 for (score=0,tot=0,a=0; a< A->nseq-1; a++)
	     {
	     b=A->nseq-1;	 
	     entry[SEQ1]=A->order[a][0];
	     entry[SEQ2]=A->order[b][0]; 
	     r1=entry[R1]=pos[a][c];
	     r2=entry[R2]=pos[b][c];
	     if ( r1<=0 || r2<=0);
 	     else if ((r=search_in_list_constraint ( entry, 3, list, ne, &start_index, &end_index))!=NULL)
 	          {
 		  tot++;
 		  if ( r[0][WE]!=UNDEFINED)score+=(r[0][CONS])*SCORE_K;
 		  }
 	     else tot++;
 		 
 	    }
 	 if ( tot!=0)score_list[c][0]=(int)((float)score/(float)tot);
 	 else score_list[c][0]=0;
 	 }  


     for ( a=1; a<=n_it; a++)
	 for ( b=1; b<A->len_aln-1; b++)
	     score_list[b][a]=(score_list[b-1][a-1]+score_list[b][a-1]+score_list[b+1][a-1])/3;

     return_2Dmax_coor_int(score_list,0,A->len_aln,n_it, n_it+1, &best_p, &b);

     left_limit =get_nol_aln_border (A, best_p, GO_LEFT );
     right_limit=get_nol_aln_border (A, best_p, GO_RIGHT);     
     first_p[0]=left_limit;

     last_p [0]=right_limit;
     free_int (score_list, -1);
     free_int (score_list2,-1);        
     return best_p;
     }

int **extract_best_column( int **in_L, int in_ne, Sequence *S, int WF)
        {
	int **L;
	int ne;
	int a, b, c,d;
	int max;
	int **seq;
	float **z_seq;
	int   **result; 
	int r1, r2;
	int n_it=1;
	double tot=0, sum=0, sum2=0, z;
	double x;
	double z_limit=0;
	double max_z, sum_z, tot_z;
	int res;
	int *entry;
	int **r;
	int score,last_score;
	int **end_index=NULL;
	int **start_index=NULL;

	in_L=undefine_list(in_L, in_ne);
	entry=vcalloc ( LIST_N_FIELDS, sizeof (int));
	L=duplicate_int (in_L, in_ne, LIST_N_FIELDS);
	ne=in_ne;

	
	max=MAX(return_max_int(L, ne, R1), return_max_int(L, ne, R2));
	seq=declare_int (n_it*2+1,max+1);
	result=declare_int(max+1, 2);

	z_seq=declare_float ( n_it*2+1,max+1);
	
	return_max_coor_int (L, ne,WF,&c);

	seq[0][L[c][R1]]+=L[c][WF];
	seq[0][L[c][R2]]+=L[c][WF];
	
	z_seq[0][L[c][R1]]+=100;
	z_seq[0][L[c][R2]]+=100;
	

	L[c][WE]=UNDEFINED;
	
	for ( a=1; a<= n_it; a++)
	    {
	    for (b=1; b<=max; b++)seq[a][b]=seq[a-1][b];
	    for (b=0; b< ne; b++)
	        {
		if ( L[b][WE]!=UNDEFINED)
		   {		   
		   r1=L[b][R1];
		   r2=L[b][R2];
		   if ( seq[a-1][r1]>0 || seq[a-1][r2]>0)
		      {
		      seq[a][r1]+=L[b][WF];
		      seq[a][r2]+=L[b][WF];
		      L[b][WE]=UNDEFINED;		   
		      }
		   }
		}
	    for (sum=0,sum2=0,tot=0, b=1; b<=max; b++)
	        {
		if ( seq[a][b]>0)
		   {
		   x=(float)seq[a][b]/1000;
		   sum+= x;
		   sum2+=x*x;
		   tot++;
		   }
		}
	    
	    for (c=0,sum_z=0, tot_z=0,max_z=0,b=1; b<=max; b++)
		if ( seq[a][b]>0)
		   {
		   z=z_seq[a][b]=(float)return_z_score ((double)seq[a][b]/1000, sum, sum2, tot);  
		   max_z=MAX(z, max_z);		   
		   sum_z+=MAX(0,z);
		   tot_z+=(z>0);
		   }

	    
	    fprintf ( stderr, "LIMIT=%.2f\n", z_limit);
	    }


	

	fprintf ( stderr, "AVG=%.2f=%.2f/%.2f",sum_z/tot_z, sum_z, tot_z);
	for ( c=1,a=1; a< max; a++)
	    {
	    z=(double)z_seq[n_it][a];
	    if (seq[n_it][a]>0 && z_seq[n_it][a]>=max_z)
		{
		result[c][0]=a;
		result[c][1]=z*100;
		c++;
		}
	    else
	        {
		seq[n_it][a]=z_seq[n_it][a]=0;
		}
		
	    }
	free_int (L, -1);
	L=duplicate_int (in_L, in_ne, LIST_N_FIELDS);
	ne=in_ne;
	
	sort_int_inv (result+1, 2, 1, 0, c-2);
	result[0][0]=c;
	
	
	for ( a=0; a< ne; a++)
	    if ( seq[n_it][L[a][R1]]>0 ||seq[n_it][L[a][R2]]>0) 
		{
		seq[n_it+1][L[a][R1]]+=L[a][WE];
		seq[n_it+1][L[a][R2]]+=L[a][WE];
		}


	for ( c=1,a=1; a<=max; a++)
	    if ( seq[n_it+1][a]>0)
	       {
	       result[c][0]=a;
	       result[c][1]=seq[n_it+1][a];
	       c++;
	       }


	sort_int_inv (result+1, 2, 1, 0, c-2);
	result[0][0]=c;
	last_score=0;
	for (d=1, a=1; a<=c; a++)
	    {
	    
	    for (score=0,b=1; b<d; b++)
	        {
		entry[SEQ1]=entry[SEQ2]=0;
		entry[R1]  =result[b][0];
		entry[R2]  =result[a][0];
		if ((r=search_in_list_constraint ( entry, 3, L, ne, &start_index, &end_index))!=NULL)
		    score+=r[0][CONS]*SCORE_K;
		}
	    score=(d==1)?10000:score/(d-1);
	    if (1)
	       {
	       result[d][0]=result[a][0];
	       result[d][1]=score;
	       d++;
	       }
	    else
		break;
	    }
	result[0][0]=d;
	/*sort_int_inv (result+1, 2, 1, 1, d-2);*/
	
	for ( a=1; a<=d; a++) fprintf ( stderr, "\n%4d %4d", result[a][0], result[a][1]);
	free_int (L, -1);
	
	return result;
	}
	

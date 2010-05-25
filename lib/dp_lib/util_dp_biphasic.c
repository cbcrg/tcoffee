int biphasic_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int ***t, ***m;
  int i,j, l1, l2, n, sub, trace,ntrace, a, b, c, score;
  int gop,rgop,tgop, gep, unmatch;
  int M1, M2, I1, D1, LEN;
  char **al, *char_buf, **aln;
  int **pos0;
  

  l1=strlen (A->seq_al[l_s[0][0]]);
  l2=strlen (A->seq_al[l_s[1][0]]);

  n=1;
  M1=n++;D1=n++;I1=n++;D2=n++;D2=n++;
  t=declare_arrayN(3, sizeof (int),n, l1+1, l2+1);
  m=declare_arrayN(3, sizeof (int),n, l1+1, l2+1);
  
  gop=CL->gop*SCORE_K;
  gep=CL->gep*SCORE_K;
  tgop=gop;
  unmatch=gep;
  pos0=aln2pos_simple ( A,-1, ns, l_s);
 
  
  for (j=1; j<=l2; j++)
    {
      m[D1][0][j]=gep1*j;
      m[D2][0][j]=gep2*j;
      m[M1][0][j]=2*gep1*j;
    }
  
  for (i=1; i<=l1; i++)
    {
      m[I1][i][0]=i*gep1;
      m[I2][i][0]=i*gep2;
      m[M1][i][0]=2*i*gep1;
      
      for ( j=1; j<=l2; j++)
	{
	  rgop1=(i==l1 || j==1)?0:gop1;
	  rgop2=(i==l1 || j==1)?0:gop2;
	  sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
	  
	  
	  
	  m[M1][i][j]=dp_max (&trace,4,M1,m[M1][i-1][j-1],I1, m[I1][i-1][j-1],D1,m[D1][i-1][j-1],M2,m[M2][i-1][j-1])+sub;
	  t[M1][i][j]=trace;
	  
	  m[D1][i][j]=dp_max (&trace,3, M1,m[M1][i][j-1]+rgop,D1, m[D1][i][j-1]+gep, M2, m[M2][i][j-1]);
	  t[D1][i][j]=trace;
	  
	  m[I1][i][j]=dp_max (&trace,3, M1,m[M1][i-1][j]+rgop, I1, m[I1][i-1][j]+gep, M2, m[M2][i-1][j]);
	  t[I1][i][j]=trace;

	  m[M2][i][j]=dp_max (&trace,4,M1,m[M1][i-1][j-1]+tgop,I1, m[I1][i-1][j-1]+tgop,D1,m[D1][i-1][j-1]+tgop,M2,m[M2][i-1][j-1])+unmatch;
	  t[M2][i][j]=trace;
	  
	}
	  
    }
  score=dp_max (&trace,4, M1,m[M1][l1][l2],D1,m[D1][l1][l2],I1, m[I1][l1][l2],M2,m[M2][l1][l2]);
  LEN=0;i=l1;j=l2;
  al=declare_char (2, l1+l2);
  

  trace=t[trace][i][j];
  while (!(i==0 &&j==0))
    {
  
      ntrace=t[trace][i][j];
      if (i==0)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( j==0)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      else if ( trace==M1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=1;
	  i--; j--;
	  LEN++;
	}
      else if ( trace==M2)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  LEN++;

	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  LEN++;

	  i--; j--;
	  
	}
      else if ( trace==D1)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( trace == I1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      trace=ntrace;     
      
    }
  
  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);	
  if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);	
  
  aln=A->seq_al;
  char_buf= vcalloc (LEN+1, sizeof (char));	
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++) 
	{		
	  int ch=0;
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
  free_arrayN((void *)m, 3);
  free_arrayN((void *)t, 3);
  vfree (char_buf);
  free_char (al, -1);
  return score;
}
    {
      m[I1][i][0]=i*gep;
      m[M2][i][0]=4*i*gep;
      m[M1][i][0]=2*i*gep;
                 
      for ( j=1; j<=l2; j++)
	{
	  rgop=(i==l1 || j==1)?0:gop;
	  rgop=gop;
	  sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
	  m[M1][i][j]=dp_max (&trace,4,M1,m[M1][i-1][j-1],I1, m[I1][i-1][j-1],D1,m[D1][i-1][j-1],M2,m[M2][i-1][j-1])+sub;
	  t[M1][i][j]=trace;
	  
	  m[D1][i][j]=dp_max (&trace,3, M1,m[M1][i][j-1]+rgop,D1, m[D1][i][j-1]+gep, M2, m[M2][i][j-1]);
	  t[D1][i][j]=trace;
	  
	  m[I1][i][j]=dp_max (&trace,3, M1,m[M1][i-1][j]+rgop, I1, m[I1][i-1][j]+gep, M2, m[M2][i-1][j]);
	  t[I1][i][j]=trace;

	  m[M2][i][j]=dp_max (&trace,4,M1,m[M1][i-1][j-1]+tgop,I1, m[I1][i-1][j-1]+tgop,D1,m[D1][i-1][j-1]+tgop,M2,m[M2][i-1][j-1])+unmatch;
	  t[M2][i][j]=trace;
	  
	}
	  
    }
  score=dp_max (&trace,4, M1,m[M1][l1][l2],D1,m[D1][l1][l2],I1, m[I1][l1][l2],M2,m[M2][l1][l2]);
  LEN=0;i=l1;j=l2;
  al=declare_char (2, l1+l2);
  

  trace=t[trace][i][j];
  while (!(i==0 &&j==0))
    {
  
      ntrace=t[trace][i][j];
      if (i==0)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( j==0)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      else if ( trace==M1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=1;
	  i--; j--;
	  LEN++;
	}
      else if ( trace==M2)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  LEN++;

	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  LEN++;

	  i--; j--;
	  
	}
      else if ( trace==D1)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( trace == I1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      trace=ntrace;     
      
    }
  
  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);	
  if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);	
  
  aln=A->seq_al;
  char_buf= vcalloc (LEN+1, sizeof (char));	
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++) 
	{		
	  int ch=0;
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
  free_arrayN((void *)m, 3);
  free_arrayN((void *)t, 3);
  vfree (char_buf);
  free_char (al, -1);
  return score;
}

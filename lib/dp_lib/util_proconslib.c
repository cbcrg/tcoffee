int procons_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
{
  char *seq1, *seq2;
  int l1, l2;
  int ***M;
  l1=strlen (  A->seq_al[l_s[0][0]]);
  seq1=vcalloc ( 1, sizeof (l1+1));
  sprintf ( seq1, "%s", A->seq_al[l_s[0][0]]);
  
  l2=strlen (  A->seq_al[l_s[1][0]]);
  seq2=vcalloc ( 1, sizeof (l2+1));
  sprintf ( seq1, "%s", A->seq_al[l_s[0][0]]);
  
  M=vcalloc (3, sizeof (int**));
  for ( a=0; a< 3; a++)M[a]=declare_int (l1+1, l2+1);
  
    


  forward_pair_wise (seq1, seq2, M);
  invert (seq1); invert(seq2);
  forward_pair_wise ( seq1, seq2, M);
  filter_matrix (M);
  for ( a=0; a<3; a++)free_int (M, -1);
  vfree(M);
}


int forward_pair_wise (char *seq1, char *seq2, char ***M)
{
  int l1, l2, i, j;
  int *matrix;
  
  for ( i=1; i<= l1; i++)
    for ( j=1; j<= l2; j++)
      {
	match=matrix[seq1[i-1]-'a'][seq2[i-1]-'a'];
	
	M[Match][i][j]+=match+M[Match][i-1][j-1]+M[D][i][j-1]+M[D][i-1][j]+M[LD][i-1][j]+M[LD][i][j-1];
	M[D][i][j]+=(M[i-1][j]+gep)+(M[i][j-1]+gep)+(M[i-1][j-1]+gop+gep);
	M[LD][i][j]+=(M[LD][i-1][j]+lgep)+(M[LD][i][j-1]+lgep) +(M[D][i-1][j]+lgep+lgop)+(M[D][i][j-1]+lgep+lgop);
      }
}
     



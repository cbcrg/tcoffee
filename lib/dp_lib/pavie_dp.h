Sequence *pavie_seq2random_seq (Sequence *S, char *subst); 
double **pavie_seq2pavie_aln(Sequence *S,char *mat, char *mode);
Alignment *pavie_seq2pavie_sort ( Sequence *S, char *mat, char *mode);
Alignment* pavie_seq2pavie_msa  ( Sequence *S, char *mat, char *mode);
NT_node pavie_seq2pavie_tree ( Sequence *S, char *mat, char *mode);
int **pavie_seq2trained_pavie_mat(Sequence *S, char *param);
int pavie_pair_wise (Alignment *A,int *ns, int **l_s,Constraint_list *CL );
Sequence *pavie_seq2noisy_seq (Sequence *S, int val,char *subst); 


  
      
  
 

int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int fasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
Dp_Model* initialize_dna_dp_model (Constraint_list *CL);
Dp_Result * make_fast_dp_pair_wise (Alignment *A,int*ns, int **l_s, Constraint_list *CL,Dp_Model *M);
int make_fasta_cdna_pair_wise (Alignment *B,Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);



int ** evaluate_diagonals_cdna ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup);
int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

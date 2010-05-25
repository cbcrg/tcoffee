Constraint_list * mask_list_with_aln (Alignment *A,int start, int len,Constraint_list *CL, int new_value);
Constraint_list* mask_list_with_aln_pair (Alignment *A,int start, int end,Constraint_list *CL,int new_value);
Constraint_list *mask_entry( Constraint_list *CL, int p, int new_value);
Constraint_list *prepare_list_and_seq4sw(Constraint_list *I, int n_seq, char **seq_name);
int ** get_undefined_list (Constraint_list *CL);
int      is_never_undefined (Constraint_list *CL,int r);
int* do_analyse_list ( Constraint_list *CL);



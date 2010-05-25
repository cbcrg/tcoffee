
Dp_Model * initialize_sseq_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL);
int ssec_pwaln_maln (Alignment *A, int *ns, int **ls, Constraint_list *CL);


int get_turn_gep_cost       (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_term_gep_cost  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_alpha_gep_cost      (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_alpha_start_gep_cost(Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_alpha_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_beta_gep_cost       (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_term_gep_cost  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_alpha_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_ssec_no_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

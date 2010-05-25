Alignment *clean_maln ( Alignment *A, Alignment *I, int T, int n_it);
Alignment *realign_segment (int seq, int start, int len,Alignment *A, Alignment *C);
Dp_Model * initialize_seg2prf_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL);

int get_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

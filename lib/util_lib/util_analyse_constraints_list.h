
void print_list(Constraint_list *CL);
void print_pair (Constraint_list *CL,int p);
int** bin_list (Constraint_list *CL,int field, int Threshold);
void   save_full_list (Constraint_list *CL, char*fname);
FILE * output_list ( Constraint_list *CL, FILE *fp);
FILE * output_pair (Constraint_list *CL,int p, FILE *fp);

void make_ga_output ( char *name,char *text, Decoded_chromosome *ALN, Chromosome *C_ALN, Parameter *PARAM);
void do_command ( char*com, Parameter *PARAM);
void save_individual (Chromosome *A,char *name, int generation, Parameter *PARAM);
void print_result_ga ( Population *L_POP,int G, Statistic * L_STAT, int generation, Parameter *PARAM);
void print_current_operator (Parameter *PARAM);
void save_current_operator ( Parameter *PARAM, char *experience_name, int gen);
void save_analyse_operator ( int gen, Parameter *PARAM);
void print_select_ga ( Population *L_POP, Parameter *L_PARAM);
void read_secondary_param ( char *fname, Parameter *PARAM);
void read_tertiary_param ( char *fname, Parameter *PARAM);
void get_ref_time (Parameter *PARAM);

void  print_chromosome (Chromosome *C, Parameter *PARAM);
void  print_decoded_chromosome ( Decoded_chromosome *B, Parameter *PARAM);


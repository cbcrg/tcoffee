void ga_select ( Population *L_POP,int G, int ** m_t,int tot_pop,int MODE, Parameter *PARAM);
int exponential_selection ( int first_len, int second_len, float factor);
int fatality_select ( int **fatality, Parameter *L_PARAM, int sum, int *x, int tot);
int select_operator ( Parameter *L_PARAM);
int wheel_select ( int **array, int tot_pop, int sum, int field, int MODE);
int select_len ( int first_len, int second_len, float factor);

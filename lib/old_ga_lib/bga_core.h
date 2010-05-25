Global_structure* do_ga (Population *L_POP, Parameter *L_PARAM, Statistic *L_STAT);
int do_generation(Population *L_POP,Parameter *L_PARAM, Statistic *STAT, int gen);
void process_fitness ( Population*L_POP, int G, int total_pop,int MODE);
void get_fitness(Population*L_POP,Parameter *L_PARAM);
int** generate_fatality ( Population *LPOP,Parameter *L_PARAM);
void check_end_signal ( char *fname, Parameter *PARAM);

void rescale_fitness (Population *POP, int G, int tot_pop, int MODE, Parameter *PARAM);
void natural_fitness (Population *POP, int G, int tot_pop);
void ranked_fitness (Population *POP, int G, int tot_pop);
void affine_fitness ( Population *POP, int G, int tot_pop, float C);
void sigma_fitness ( Population *POP, int G, int tot_pop, float C);

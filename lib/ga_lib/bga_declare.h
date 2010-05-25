/*insert here structure*/


Parameter * declare_parameter ();
void declare_bga (Parameter *PARAM);
void proto_mut ( Chromosome *, Chromosome *, Parameter *);
void declare_dos_1 ( Parameter *L_PARAM);
void declare_dos_2 ( Parameter *L_PARAM);
void reset_dos ( Parameter *PARAM);
Statistic * declare_statistic (Parameter *PARAM);
void reset_statistic (Statistic *STAT);
Population* declare_population (int tot_pop, Parameter *PARAM);

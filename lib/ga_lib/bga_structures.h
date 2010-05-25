
typedef struct
    {
    char *name;
    char *description;
    int used;
    char TYPE;
    char *S;
    int *D;
    float *F;
    int *O_1;
    int *O_2;
    int read;
    
    }P_list;

typedef struct
    {
    int RESCALING_MODE ;
    float RESCALING_C ;
    int PRESELECT_MODE;
    int SELECT_MODE;
    int MAXPOP;
    float OVERLAPING_GEN_FACTOR ;
    int MAX_GEN;
    int FACTOR_POP;
    int REPLACE_MODE;
    int RANDOM_SEED;
    int ALLOW_DUPLICATE;
    int DYNAMIC_PARAM;
    int STABILISED;  
    float CREDIT_TO_REPORT;
    float FIT_TO_ADAPT;
    int TOT_OP_FIT;
    int STEP_FIT_OPERATOR;
    int MAX_STAB;
    char experience_name[100];
    int GA_MAXIMISE;
    int BIASE_POS;
    int ANALYSE_OPERATOR;
    int OUTPUT_SCREEN;
    int OUTPUT_FILE;
    int INTERACTIVE;
    int ANALYSE_CONSISTENCY;
    int PARALLEL_GA;
    int PARALLEL_EXCHANGE;
    int PARALLEL_MODE;
    char TYPE[20];
    char NORMAL_END_FLAG_FILE [500];
    char POP_READ_FLAG_FILE[500];
    char POP_DUMP_FLAG_FILE[500];
    char FINISH_FLAG_FILE[500];
    char ERROR_LOG_FILE[500];
    FILE *FP_ERROR;     
    char PWD[500];
    char OUTPUT_PARAM_FILE[200];
    int ref_time;
    int EVALUATE_ONLY;
    int OF_MODE;
    int  OUTPUT_PARAM;
    int OF_CHECK;
    int OUTPUT_ERROR_FILE;
    int OUTPUT_FINISH_FLAG;
    int OUTPUT_END_FLAG;
    int OUTPUT_POP_READ_FLAG;
    int OUTPUT_POP_DUMP_FLAG;
    int OUTPUT_GENERATION_FILE;
    int OUTPUT_SCORE_FILE;
    int  OUTPUT_CPU_FILE;
    int  OUTPUT_GEN_BEST_FILE;
    int  OUTPUT_BILAN_FILE;
    int  OUTPUT_CURRENT_OPERATOR_FILE;
    char OUTPUT_RESULT_FILE[200];
    char OUTPUT_ALN_FILE[200];
    int OUTPUT_ALN;
    int OUTPUT_ERROR_LOG_FILE;
    int OUTPUT_RESULT;
    char REF_DEC_CHROM[200];
    Decoded_chromosome * REFERENCE;
    float REF_SCORE;
    float T;
    }Bga;

typedef struct
    {   
    int N_BI_PARENT_LIST;
    int *BI_PARENT_LIST;
    int TOT_OPERATOR;
    int ***T;
    int **fitness_table;
    int **current_fitness;
    int **OP_TO_USE;
    int **OP_INIT_VAL;
    int N_OPERATOR;
    int **OPUSE;
    int **IOPUSE;
    char **OP_NAME;
    Chromosome* (*OP_FUNCTION[200])(Chromosome *, Chromosome *, void*);
    }Dos;

typedef struct 
    {
    float DE;
    Chromosome ** CHROMOSOME;
    float **raw_fitness;
    float **threshold;
    float **fitness;
    float **scaled_fitness;
    int **expected_os;
    int **age;
    float **ranked;
    }Population;

typedef struct
    {
    int mean;
    int SD;
    int best_gen_chrom;
    int score_best;
    int score_gen_worst;
    int score_gen_best;
    Chromosome *BEST_CHROMOSOME;
    Chromosome *BEST_GEN_CHROMOSOME;
    Decoded_chromosome *BEST_DECODED_CHROMOSOME;
    Decoded_chromosome *BEST_GEN_DECODED_CHROMOSOME;
    }Statistic;


typedef struct
    {
     Decoded_chromosome *SD1;
     Decoded_chromosome *SD2;
     Chromosome *C1;
     Chromosome *C2;
     Decoded_chromosome **C_BUF;
     int NC;
     }Mem;

typedef struct
    {
    Bga *BGA;
    Dos *DOS;
    Opp *OPP;
    Data *DATA;
    Ofp *OFP;
    Seed *SEED;
    P_list *P_LIST;
    Population *POPULATION;

    Mem *MEM;
    char OS[100];
    int START;
    int NP;
    int LAST_OPERATOR;
    int FIRST_OPERATOR;
    FILE *std1;
    FILE *std2;
    FILE *std0;
    }Parameter;


typedef struct {
    Parameter *PARAM;
    Population *POP;
    Statistic *STAT;
    }Global_structure;

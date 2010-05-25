struct Dps_result
    {
      int njobs;
      struct Dps_job **dps_job;
};
typedef struct Dps_result Dps_result;

struct Dps_job
    {
      int JobId;
      struct Constraint_list *CL;
      char *input_file;
      char *output_file;
};
typedef struct Dps_job Dps_job;

struct Dps_result *seq2list_DPS (struct Constraint_list *CL,char *method, char *aln_command, char *seq_command, char *weight, Dps_result *dps_result);
struct Constraint_list * gather_results_DPS ( Dps_result *DPS, struct Constraint_list *CL);
Dps_result *declare_dps_result ( int naln, Dps_result *dps);

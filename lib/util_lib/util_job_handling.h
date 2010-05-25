struct Job_TC
    {
      int jobid;
      int status;
      
      struct  Job_TC *c;
      struct  Job_TC *p;
      struct  Job_io_TC *io;
      struct  Job_control_TC *control;
      
      struct  Job_param_TC *param;
        
      /*memory mangement*/
      char **pl;
      int np;
};
typedef struct Job_TC Job_TC;

struct Job_control_TC
    {
      
      struct Job_TC* (*submitF) (struct Job_TC*);
      struct Job_TC* (*retrieveF)(struct Job_TC*);
      char *mode;
};
typedef struct Job_control_TC Job_control_TC;

struct Job_io_TC
    {
      char *in;
      char *out;
      struct Constraint_list *CL;
      struct Alignment *A;
};
typedef struct Job_io_TC Job_io_TC;

struct Job_param_TC
{
  char *method;
  struct TC_method *TCM; 
  char *temp_c;
  char *aln_c;
  char *seq_c;
  char *aln_mode;
};
typedef struct Job_param_TC Job_param_TC;

Job_TC* print_lib_job ( Job_TC *job,char *string, ...);
Job_TC *print_lib_job2 ( Job_TC* job, int n, char **name, char **value);


/*Stack Manipulation*/
Job_TC *free_queue  (Job_TC *job);
Job_TC *free_job  (Job_TC *job);
Job_TC * queue2heap (Job_TC*job);
Job_TC * queue2last (Job_TC*job);
int queue2n (Job_TC*job);
Job_TC * descend_queue (Job_TC*job);
Job_TC *queue_cat  (Job_TC *P, Job_TC *C);
Job_TC *delete_job (Job_TC *job);
/*Job Control*/
struct Job_TC* submit_job ( Job_TC *job);
struct Job_TC* retrieve_job ( Job_TC *job);
Job_TC*** split_job_list (Job_TC *job, int ns);

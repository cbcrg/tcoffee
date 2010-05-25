


void process_command_line ( int nargc,char **list, Parameter *PARAM);
void input_main_parameter  (char *fname, Parameter *PARAM);
void initialise_parameter_list (Parameter *PARAM);
void read_parameter_file ( char *fname, Parameter *PARAM);
void read_command_line ( int narg, char **list, Parameter *PARAM);
void process_param ( Parameter *PARAM);
void option_help ( Parameter *PARAM, char *option);
void output_help ( Parameter *PARAM, char *exec);
int is_cl_flag (char *name);
void get_paramter_help ( char *name, P_list *P);
void output_parameters ( Parameter *PARAM);

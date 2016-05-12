char * pdb2lib (Sequence *S, char *mode,float max, char *name);

int apdb (int argc, char *argv[]);

Constraint_list * set_constraint_list4align_pdb (Constraint_list *inCL,int seq, char *dp_mode, char *hasch_mode, char *param_file);
int evaluate_ca_trace_sap2_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_ca_trace_nb (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_3D_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_transversal (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble_2 (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble_3 (Constraint_list *CL, int s1, int s2, int r1, int r2);


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:START                                   */
/*                                                                                           */
/*********************************************************************************************/
float matrix_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float transversal_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2,Struct_nb *nbs1, Struct_nb *nbs2);
float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float sap2_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:END                                     */
/*                                                                                           */
/*********************************************************************************************/
Alignment * analyse_pdb ( Alignment *A, Alignment *ST, char *filename);
Alignment * msa2struc_dist ( Alignment *A, Alignment *ST, char *prefix, char *output,int gapped, int min_ncol);
float **** analyse_pdb_residues ( Alignment *A, Constraint_list *CL, Pdb_param *pdb_param);

float square_atom ( Atom *X);
Atom* reframe_atom ( Atom *X, Atom*Y, Atom *Z, Atom *IN, Atom *R);
Atom* add_atom ( Atom *A, Atom*B, Atom *R);
Atom* diff_atom ( Atom *A, Atom*B, Atom *R);
Atom * copy_atom ( Atom *A, Atom*R);
/************************************************************************/
/*                                                                      */
/*                            NUSSINOV                                  */
/*                                                                      */
/************************************************************************/
char *nussinov (char *S, int min_dist);
char * rna_struc2rna_lib ( char *seq_name, char *seq, char *name);
int display_rna_ss ( int pos, char *seq, char *struc, FILE *fp);

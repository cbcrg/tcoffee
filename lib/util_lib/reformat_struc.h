char * process_repeat (char *aln, char *seq, char *pdb);

char     * normalize_pdb_file  (char *name, char *seq,char *out_file);
Ca_trace * trim_ca_trace (Ca_trace *st, char *seq );
void bfactor2x_in_pdb (char *infile, char *outfile, float *v);
Ca_trace * read_ca_trace (char *file, char *seq_field );
Ca_trace * simple_read_ca_trace (char *file );
Ca_trace * hasch_ca_trace             ( Ca_trace *T);
Ca_trace * hasch_ca_trace_nb          ( Ca_trace *T);
Ca_trace * hasch_ca_trace_bubble      ( Ca_trace *T);
Ca_trace * hasch_ca_trace_transversal ( Ca_trace *TRACE);

float get_closest_vdw_contact ( Atom *A, Atom*B);
float get_atomic_distance ( Atom *A, Atom*B);
int   residues2contacts (Atom *A, Atom*B, float probe);
float ** trace2ca_distances(Ca_trace *T, float max);
float ** trace2contacts(Ca_trace *T, float probe);
float ** trace2best_contacts(Ca_trace *T, float probe );
float ** trace2closest_contacts (Ca_trace *T, float probe);
float ** trace2count_contacts(Ca_trace *T, float probe);


float ** traces2ca_distances    (Ca_trace *T1,Ca_trace *T2, float max);
float ** traces2contacts        (Ca_trace *T1,Ca_trace *T2, float max);
float ** traces2best_contacts   (Ca_trace *T1,Ca_trace *T2, float max);
float ** traces2closest_contacts(Ca_trace *T1,Ca_trace *T2, float max);
float ** traces2count_contacts  (Ca_trace *T1,Ca_trace *T2, float max);

float** print_contacts ( char  *file1, char *file2, float T);
char *  map_contacts ( char  *file1, char *file2, float T);
int * identify_contacts (Ca_trace *ST1,Ca_trace *ST2, float T);
Sequence *seq2contacts4lignads ( Sequence *S, float T);
char *string2contacts (char *seq,char *name,char *comment, float T);
char **struc2nb (char *name,char *seq, char *comment, float Threshold, char *atom_list, char *output);
char **struclist2nb (char *name,char *seq, char *comment, float Threshold, char *atom_list, char *output);

float atom2radius (char *t);

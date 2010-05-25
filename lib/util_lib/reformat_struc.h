char * process_repeat (char *aln, char *seq, char *pdb);
char     * normalize_pdb_file  (char *name, char *seq,char *out_file);
Ca_trace * trim_ca_trace (Ca_trace *st, char *seq );

Ca_trace * read_ca_trace (char *file, char *seq_field );
Ca_trace * simple_read_ca_trace (char *file );
Ca_trace * hasch_ca_trace             ( Ca_trace *T);
Ca_trace * hasch_ca_trace_nb          ( Ca_trace *T);
Ca_trace * hasch_ca_trace_bubble      ( Ca_trace *T);
Ca_trace * hasch_ca_trace_transversal ( Ca_trace *TRACE);

float get_atomic_distance ( Atom *A, Atom*B);
float ** measure_ca_distances(Ca_trace *T);

float** print_contacts ( char  *file1, char *file2, float T);
char *  map_contacts ( char  *file1, char *file2, float T);
int * identify_contacts (Ca_trace *ST1,Ca_trace *ST2, float T);
Sequence *seq2contacts ( Sequence *S, float T);
char *string2contacts (char *seq,char *name,char *comment, float T);
char **struc2nb (char *name,char *seq, char *comment, float Threshold, char *atom_list, char *output);
char **struclist2nb (char *name,char *seq, char *comment, float Threshold, char *atom_list, char *output);


int **pavie_seq2trained_pavie_mat(Sequence *S);
double delta_matrix ( int **mat1, **mat2, char *alp);
double **pavie_seq2pavie_mat (Sequence *S, char *mat  );
int **pavie_fmat2pavie_log_mat (double **fmat, char *alphabet);
double **pavie_aln2fmat(Alignment *A, double **fmat);
Alignment* pavie_seq2pavie_msa ( Sequence *S, char *mat, char *mode);
float **pavie_seq2pavie_aln(Sequence *S,char *mat, char *mode);

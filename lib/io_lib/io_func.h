#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

typedef struct
    {
    int in_seq;
    FILE *fp;
    int font;
    int x0;
    int y0;
    int x;
    int y;
    int n_pages;
    int max_line_ppage;
    int n_line;
    int line;
    int eop;
    int in_html_span;
    char previous_html_color[100];

    }
FILE_format;

typedef struct
    {
    float r;
    float g;
    float b;
    char html_color[30];
    char html_color_class[30];
    int ascii_value;
    }
Color;


Sequence * fill_sequence_struc ( int nseq, char **sequences, char **seq_name, Genomic_info *genome_co);
Sequence * cw_read_sequences ( char *seq_name);
Sequence * get_sequence_type (Sequence *S);
char     * get_array_type (int n, char **s);
Alignment* get_aln_type (Alignment *A);

char     * get_string_type   (char *string);

char *store_mode (char *val);
char *retrieve_mode ();
char *unset_mode ();
char *set_mode (int mode, char *val);

char *store_seq_type (char *val);
char *retrieve_seq_type ();
char *unset_seq_type ();
char *set_seq_type (int mode, char *val);

void get_sequence (char *seq_file,int *NSEQ, char ***SEQ, char ***SN, int **sl, int *min, int *max);

int ** get_matrix   ( char *name, char *format);
int ** read_matrice (char *mat_name);
int **neg_matrix2pos_matrix ( int **matrix);


void   print_aln ( Alignment *B);

int       output_reliability_ps     ( Alignment *B,Alignment *S, char *name);
int       output_reliability_pdf    ( Alignment *B,Alignment *S, char *name);
int       output_reliability_html   ( Alignment *B,Alignment *S, char *name);
int       output_color_ps     ( Alignment *B,Alignment *S, char *name);
int       output_color_pdf    ( Alignment *B,Alignment *S, char *name);
int       output_color_html   ( Alignment *B,Alignment *S, char *name);
int       output_hit_color_html   (Alignment *B, float **ffPScoreTable, int nl, char *name);	//JM_ADD
void 	  output_hit_matrix(char *fileName, float **ffpHitScoreMatrix, int nl);		//JM_ADD
void      get_rgb_values(int val, Color *C);
int       output_reliability_format ( Alignment *B,Alignment *S, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));
int       output_score_format ( Alignment *B,Alignment *S, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));


FILE_format * print_ps_string      ( char *s , Color *box, Color *ink, FILE_format *f);
FILE_format * print_ps_char        ( int   c,    Color *box, Color *ink, FILE_format *f);



void get_rgb_values_ps ( int val, Color *C);
FILE_format* vfopen_ps ( char *name);
FILE_format* vfclose_ps ( FILE_format *fps);

FILE_format *print_html_string( char *s, Color *box, Color *ink, FILE_format *fhtml);
FILE_format * print_html_char ( int c, Color *box, Color *ink, FILE_format *f);
void get_rgb_values_html ( int val, Color *C);
FILE_format* vfopen_html ( char *name);
FILE_format* vfclose_html ( FILE_format *fhtml);

int       output_reliability_ascii     ( Alignment *B,Alignment *S, char *name);
int       output_reliability_fasta     ( Alignment *B,Alignment *S, char *name);
int       output_color_ascii           ( Alignment *B,Alignment *S, char *name);

FILE_format *print_ascii_string( char *s, Color *box, Color *ink, FILE_format *fascii);
FILE_format * print_ascii_char ( int c, Color *box, Color *ink, FILE_format *f);
void get_rgb_values_ascii ( int val, Color *C);

FILE_format* vfopen_ascii ( char *name);
FILE_format* vfclose_ascii ( FILE_format *fascii);
int       output_seq_reliability_ascii     ( Alignment *B,Alignment *S, char *name);

#include <stdio.h>
struct bin_node
      {
      struct bin_node *parent;
      struct bin_node *left_child;
      struct bin_node *right_child;
      int n;
      int score;
      int used;
      int *list;
      int index;
      int max_list;
      };
typedef struct bin_node Bin_node;
struct mot
      {
      struct mot  * parent;
      struct mot **child;
      struct mot **leaf_list;
      struct bin_node *bin_tree_val;
      char       *motif;
      int val;
      int *bin_val;
      int weight;
      int motif_offset;
      int n_motifs;
      short *cache;
      int n_children;
      };
typedef struct mot Motif;
struct oligo
      {
      int NSEQ;
      int LEN;
      int WSIZE;
      char *ALPHABET;
      char *EALPHABET;
      char *AMBIGUITIES;
      char **seq;
      Motif  *motif;
      Motif *emotif;
      };
typedef struct oligo Oligo;

Motif*	build_motif_tree ( Motif *motif_tree, int len, int n, char *al,int ns, int n_motifs);
Motif* find_motif_in_tree ( char *motif,int p, int len, Motif *mtree);
int bin2int ( int *in, int len);
Motif * declare_motif (int len, char *al, int ns, int n_motif);
int same_motif ( char *m1, char *m2);
Oligo * read_oligo_list ( char *fname);
int motif_is_in_seq ( char *seq, char *motif);
int weight_motif ( char *motif);
FILE *print_motif_and_bin ( FILE *fp,  char *p, Motif *M, int size);
int get_bin_edit_distance ( int *b1, int *b2, int l);
int get_bin_edit_distance_fp ( int *b1, int *b2, int l);
int get_bin_edit_distance_fn ( int *b1, int *b2, int l);



int bin_is_included_in_bin ( int *b1, int *b2, int l);
int count_bits (int *b, int l);
Bin_node * find_val_in_bin_tree ( Bin_node ** Node_list,Bin_node *N,int index, int *bin_val, int n, int len, int *used);

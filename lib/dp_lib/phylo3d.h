#include "define_header.h"
struct p3D
    {
      int N;
      Alignment *A;
      Constraint_list *CL;
      Sequence *S;
      
      float max_gap;
      char *tree_mode;
      
      int replicates;
      int enb;
      double maxd;
      double extremed;
      
      int **pos;
      int ***dm3d;
      double **dm;
      int **col;
      int **colrep;
      int *used_site;
      int nsites;
      int **used_site_pair;
      int nsitepairs;
      
      

    };
typedef struct p3D p3D;
Alignment *phylo3d (Alignment *inA, Constraint_list *CL);
Alignment *phylo3d_gt (Alignment *inA, Constraint_list *CL);
double scan_maxd (p3D *D);

int free_3dD (p3D*D);
p3D* fill_p3D (Alignment *A, Constraint_list *CL);
int col2n (int **col);
p3D * makerep (p3D *D, int mode);
int **col2colrep (int **colin,int **colout, int ni, int mode);

int **col2rep   (int **colin,int **colout, int ni);
int **col2bsrep (int **colin, int **colout, int ni);

Alignment * addtree (p3D *D,Alignment *A);
int filter_columns_with_dist(Alignment *B, int **pos,int **col,int ***dm, double maxd);
int filter_columns_with_dist_strict(Alignment *B, int **pos,int **col,int ***dm, double maxd);
int filter_columns_with_dist_relaxed(Alignment *B, int **pos,int **col,int ***dm, double maxd);

int filter_columns_with_gap (int **col, Alignment *B, float max_gap);

int**  msa2column_list (Alignment *B   , int **col);
int** file2column_list (char      *file, int **col);

int*** aln2dm3d (Alignment *A, Constraint_list*CL, double *extremed);
Alignment *aln2trim3d (Alignment *A, Constraint_list *CL);
double pair2dist(p3D *D, int s1, int s2);
int aln2dm (p3D *D, Alignment *A);
double phylo3d2score (double w1, double w2, double *rscore, double *rmax);
 

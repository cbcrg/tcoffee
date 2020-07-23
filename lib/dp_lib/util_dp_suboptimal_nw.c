#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

float * backward_proba_pair_wise_test ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *matchProb, float ***insProb, float **transProb);
float * forward_proba_pair_wise_test ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *matchProb, float ***insProb, float **transProb);
Constraint_list *ProbaMatrix2CL_test (Alignment *A, int *ns, int **ls, int NumMatrixTypes, int NumInsertStates, float *forward, float *backward, float thr, Constraint_list *CL);
float ComputeTotalProbability_test (int seq1Length, int seq2Length,int NumMatrixTypes, int NumInsertStates,float *forward, float *backward);
float * backward_proba_pair_wise_test_old ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *matchProb, float ***insProb, float **transProb);

//Values as provided in Probcons V1.1
static float EXP_UNDERFLOW_THRESHOLD = -4.60f;
static float LOG_UNDERFLOW_THRESHOLD = 7.50f;
//static float LOG_ZERO = -FLT_MAX;
static float LOG_ZERO=-200000004008175468544.000000;
static float LOG_ONE = 0.0f;
//DNA Alignment Models
static float DNAinitDistrib2Default[] ={ 0.9588437676f, 0.0205782652f, 0.0205782652f };
static float DNAgapOpen2Default[] = { 0.0190259293f, 0.0190259293f };
static float DNAgapExtend2Default[] = { 0.3269913495f, 0.3269913495f };

static char  DNAalphabetDefault[] = "ACGUTN";
static float DNAemitSingleDefault[6] = {0.2270790040f, 0.2422080040f, 0.2839320004f, 0.2464679927f, 0.2464679927f, 0.0003124650f};

static float DNAemitPairsDefault[6][6] = {
  { 0.1487240046f, 0.0184142999f, 0.0361397006f, 0.0238473993f, 0.0238473993f, 0.0000375308f },
  { 0.0184142999f, 0.1583919972f, 0.0275536999f, 0.0389291011f, 0.0389291011f, 0.0000815823f },
  { 0.0361397006f, 0.0275536999f, 0.1979320049f, 0.0244289003f, 0.0244289003f, 0.0000824765f },
  { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f, 0.1557479948f, 0.0000743985f },
  { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f, 0.1557479948f, 0.0000743985f },
  { 0.0000375308f, 0.0000815823f, 0.0000824765f, 0.0000743985f, 0.0000743985f, 0.0000263252f }
};
//RNA Alignment Models

static float RNAinitDistrib2Default[] = { 0.9615409374f, 0.0000004538f, 0.0000004538f, 0.0192291681f, 0.0192291681f };
static float RNAgapOpen2Default[] = { 0.0082473317f, 0.0082473317f, 0.0107844425f, 0.0107844425f };
static float RNAgapExtend2Default[] = { 0.3210460842f, 0.3210460842f, 0.3298229277f, 0.3298229277f };

static char RNAalphabetDefault[] = "ACGUTN";
static float RNAemitSingleDefault[6] = {0.2270790040f, 0.2422080040f, 0.2839320004f, 0.2464679927f, 0.2464679927f, 0.0003124650f};

static float RNAemitPairsDefault[6][6] = {
  { 0.1487240046f, 0.0184142999f, 0.0361397006f, 0.0238473993f, 0.0238473993f, 0.0000375308f },
  { 0.0184142999f, 0.1583919972f, 0.0275536999f, 0.0389291011f, 0.0389291011f, 0.0000815823f },
  { 0.0361397006f, 0.0275536999f, 0.1979320049f, 0.0244289003f, 0.0244289003f, 0.0000824765f },
  { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f, 0.1557479948f, 0.0000743985f },
  { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f, 0.1557479948f, 0.0000743985f },
  { 0.0000375308f, 0.0000815823f, 0.0000824765f, 0.0000743985f, 0.0000743985f, 0.0000263252f }
};

//Protein Alignment Gap Model: Monophasic
float initDistrib1Default[] = { 0.6080327034f, 0.1959836632f, 0.1959836632f };
float gapOpen1Default[] = { 0.01993141696f, 0.01993141696f };
float gapExtend1Default[] = { 0.7943345308f, 0.7943345308f };

//Protein Alignment Models: bi-phasic
static float initDistrib2Default[] = { 0.6814756989f, 8.615339902e-05f, 8.615339902e-05f, 0.1591759622f, 0.1591759622f };
static float gapOpen2Default[] = { 0.0119511066f, 0.0119511066f, 0.008008334786f, 0.008008334786f };
static float gapExtend2Default[] = { 0.3965826333f, 0.3965826333f, 0.8988758326f, 0.8988758326f };

static char alphabetDefault[] = "ARNDCQEGHILKMFPSTWYV";
static float emitSingleDefault[20] = {
  0.07831005f, 0.05246024f, 0.04433257f, 0.05130349f, 0.02189704f,
  0.03585766f, 0.05615771f, 0.07783433f, 0.02601093f, 0.06511648f,
  0.09716489f, 0.05877077f, 0.02438117f, 0.04463228f, 0.03940142f,
  0.05849916f, 0.05115306f, 0.01203523f, 0.03124726f, 0.07343426f
};

static float emitPairsDefault[20][20] = {
  {0.02373072f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00244502f, 0.01775118f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00210228f, 0.00207782f, 0.01281864f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00223549f, 0.00161657f, 0.00353540f, 0.01911178f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00145515f, 0.00044701f, 0.00042479f, 0.00036798f, 0.01013470f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00219102f, 0.00253532f, 0.00158223f, 0.00176784f, 0.00032102f, 0.00756604f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00332218f, 0.00268865f, 0.00224738f, 0.00496800f, 0.00037956f, 0.00345128f, 0.01676565f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00597898f, 0.00194865f, 0.00288882f, 0.00235249f, 0.00071206f, 0.00142432f, 0.00214860f, 0.04062876f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00114353f, 0.00132105f, 0.00141205f, 0.00097077f, 0.00026421f, 0.00113901f, 0.00131767f, 0.00103704f, 0.00867996f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00318853f, 0.00138145f, 0.00104273f, 0.00105355f, 0.00094040f, 0.00100883f, 0.00124207f, 0.00142520f, 0.00059716f, 0.01778263f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00449576f, 0.00246811f, 0.00160275f, 0.00161966f, 0.00138494f, 0.00180553f, 0.00222063f, 0.00212853f, 0.00111754f, 0.01071834f, 0.03583921f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00331693f, 0.00595650f, 0.00257310f, 0.00252518f, 0.00046951f, 0.00312308f, 0.00428420f, 0.00259311f, 0.00121376f, 0.00157852f, 0.00259626f, 0.01612228f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00148878f, 0.00076734f, 0.00063401f, 0.00047808f, 0.00037421f, 0.00075546f, 0.00076105f, 0.00066504f, 0.00042237f, 0.00224097f, 0.00461939f, 0.00096120f, 0.00409522f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00165004f, 0.00090768f, 0.00084658f, 0.00069041f, 0.00052274f, 0.00059248f, 0.00078814f, 0.00115204f, 0.00072545f, 0.00279948f, 0.00533369f, 0.00087222f, 0.00116111f, 0.01661038f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00230618f, 0.00106268f, 0.00100282f, 0.00125381f, 0.00034766f, 0.00090111f, 0.00151550f, 0.00155601f, 0.00049078f, 0.00103767f, 0.00157310f, 0.00154836f, 0.00046718f, 0.00060701f, 0.01846071f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00631752f, 0.00224540f, 0.00301397f, 0.00285226f, 0.00094867f, 0.00191155f, 0.00293898f, 0.00381962f, 0.00116422f, 0.00173565f, 0.00250962f, 0.00312633f, 0.00087787f, 0.00119036f, 0.00180037f, 0.01346609f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00389995f, 0.00186053f, 0.00220144f, 0.00180488f, 0.00073798f, 0.00154526f, 0.00216760f, 0.00214841f, 0.00077747f, 0.00248968f, 0.00302273f, 0.00250862f, 0.00093371f, 0.00107595f, 0.00147982f, 0.00487295f, 0.01299436f, 0.0f, 0.0f, 0.0f},
  {0.00039119f, 0.00029139f, 0.00021006f, 0.00016015f, 0.00010666f, 0.00020592f, 0.00023815f, 0.00038786f, 0.00019097f, 0.00039549f, 0.00076736f, 0.00028448f, 0.00016253f, 0.00085751f, 0.00015674f, 0.00026525f, 0.00024961f, 0.00563625f, 0.0f, 0.0f},
  {0.00131840f, 0.00099430f, 0.00074960f, 0.00066005f, 0.00036626f, 0.00070192f, 0.00092548f, 0.00089301f, 0.00131038f, 0.00127857f, 0.00219713f, 0.00100817f, 0.00054105f, 0.00368739f, 0.00047608f, 0.00102648f, 0.00094759f, 0.00069226f, 0.00999315f, 0.0f},
  {0.00533241f, 0.00169359f, 0.00136609f, 0.00127915f, 0.00119152f, 0.00132844f, 0.00178697f, 0.00194579f, 0.00071553f, 0.01117956f, 0.00914460f, 0.00210897f, 0.00197461f, 0.00256159f, 0.00135781f, 0.00241601f, 0.00343452f, 0.00038538f, 0.00148001f, 0.02075171f}
};


static int suboptimal_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int mode);
static int *** forward_so_dp ( Alignment *A, int *ns, int **ls, int **pos,int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL);
static int *** backward_so_dp ( Alignment *A, int *ns, int **ls,int **pos,int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL);
static int *** forward_so_dp_biphasic ( Alignment *A, int *ns, int **ls, int **pos,int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL);
static int *** backward_so_dp_biphasic ( Alignment *A, int *ns, int **ls,int **pos,int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL);
static int *** forward_so_dp_glocal ( Alignment *A, int *ns, int **ls, int **pos,int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL);
static int *** backward_so_dp_glocal ( Alignment *A, int *ns, int **ls,int **pos,int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL);

static int match=0;
static int ins=1;
static int del=2;
static int umatch=3;
static int ins2=3;
static int del2=4;
float **    get_emitPairs (char *mat, char *alp, float **p, float *s);
int subop1_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  return suboptimal_pair_wise ( A, ns, ls, CL, 1);
}

int subop2_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  return suboptimal_pair_wise ( A, ns, ls, CL, 3);
}



int suboptimal_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int mode)
{
  int ***F=NULL;
  int ***B=NULL;
  int **pos0;
  int gop, gep,gop2, gep2;
  int i, I, j, J, n, s1, s2;
  char *seqI, *seqJ;
  int id;
  int *entry;
  float opt, min, score, nscore, thres;
  int l1, l2, set;


  gop=CL->gop*SCORE_K;
  gep=CL->gep*SCORE_K;

  /*gop2=CL->gop*10*SCORE_K;*/
  gop2=CL->gop*2*SCORE_K;
  gep2=0;

  //Values Adapted from Probcons 1.1
  gop=-132;
  gep=-27;

  gop2=-144;
  gep2=-3;

  ungap(A->seq_al[ls[0][0]]);
  ungap(A->seq_al[ls[1][0]]);

  seqI=A->seq_al[ls[0][0]];
  seqJ=A->seq_al[ls[1][0]];

  I=strlen (seqI); J=strlen (seqJ);
  pos0=aln2pos_simple ( A,-1, ns, ls);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);

  if ( mode==1)
    {
      F=forward_so_dp (A, ns, ls, pos0,I, J,gop, gep,gop2, gep2,CL);
      B=backward_so_dp (A, ns, ls, pos0,I, J,gop, gep,gop2, gep2, CL);
    }
  else if ( mode ==2)
    {
      F=forward_so_dp_glocal (A, ns, ls, pos0,I, J,gop, gep,gop2, gep2,CL);
      B=backward_so_dp_glocal (A, ns, ls, pos0,I, J,gop, gep,gop2, gep2, CL);
    }
  else if ( mode ==3)
    {
      F=forward_so_dp_biphasic (A, ns, ls, pos0,I, J,gop, gep,gop2, gep2,CL);
      B=backward_so_dp_biphasic (A, ns, ls, pos0,I, J,gop, gep,gop2, gep2, CL);
    }
  if ( MAX5(F[match][l1][l2], F[ins][l1][l2], F[del][l1][l2],F[ins2][l1][l2], F[del2][l1][l2] )!=MAX5( B[match][1][1], B[ins][1][1], B[del][1][1], B[ins2][1][1], B[del2][1][1]))
    {
      HERE ("ERROR in subop_pair");
      fprintf ( stdout, "\nForward:  %d", MAX3(F[match][l1][l2], F[ins][l1][l2], F[del][l1][l2]));
      fprintf ( stdout, "\nBackWard: %d \n\n",MAX3( B[match][1][1], B[ins][1][1], B[del][1][1]));
    }


  for (opt=0,min=0, set=0, i=1; i<=I; i++)
    for (j=1; j<=J; j++)
      {
	if ( F[match][i][j]==UNDEFINED)continue;
	F[match][i][j]+=B[match][i][j]-(CL->get_dp_cost) (A, pos0, ns[0], ls[0], i-1, pos0, ns[1], ls[1],j-1,CL);
	if (set==0)
	  {set=1; opt=F[match][i][j];min=F[match][i][j];}
	opt=MAX(F[match][i][j],opt);
	min=MIN(F[match][i][j],min);
      }


  s1=name_is_in_list (A->name[ls[0][0]], (CL->S)->name, (CL->S)->nseq, 100);
  s2=name_is_in_list (A->name[ls[1][0]], (CL->S)->name, (CL->S)->nseq, 100);

  id=idscore_pairseq(seqI,seqJ,-12, -1, CL->M, "idmat");
  
  entry=(int*)vcalloc ( CL->entry_len+1, CL->el_size);
  entry[SEQ1]=s1;entry[SEQ2]=s2;

  thres=opt;
  for ( n=0,i=1; i<=I; i++)
    {
      for (j=1; j<=J; j++)
	{
	  score=F[0][i][j];
	  nscore=((score-min))/(opt-min);

	  if (score==opt)
	    {
	      n++;
	      entry[R1]=i;entry[R2]=j;
	      entry[WE]=id;
	      entry[CONS]=1;

	      add_entry2list (entry,A->CL);
	    }
	}
    }

  vfree (entry);
  free_int (pos0, -1);
  free_arrayN (F, 3);
  free_arrayN (B, 3);

  return A->score_aln;
}
/************************************************************************************************************************/
/*                                                                                                                      */
/*                                                                                                                      */
/*                                                     GLOCAL                                                           */
/*                                                                                                                      */
/*                                                                                                                      */
/************************************************************************************************************************/
int *** forward_so_dp_glocal ( Alignment *A, int *ns, int **ls, int **pos0,int I, int J,int gop, int gep,int gop2, int gep2,Constraint_list *CL)
{
  int i,j;
  int c;
  int sub;
  int ***M;
  int match=0, del=1, ins=2;

  M=(int***)declare_arrayN (3, sizeof (int), 5, I+1, J+1);

  for ( i=0; i<=I; i++)for (j=0; j<=J; j++)for (c=0; c<5; c++)M[c][i][j]=-999999;

  M[match][0][0]=0;

  for (i=1; i<=I; i++){M[del]  [i][0]=i*gep;M[umatch][i][0]=i*gep2+gop2;}
  for (j=1; j<=J; j++){M[ins]  [0][j]=j*gep;M[umatch][0][j]=j*gep2+gop2;}


  for (i=1; i<=I; i++)
    {
      for ( j=1; j<=J; j++)
	{
	sub=(CL->get_dp_cost) (A, pos0, ns[0], ls[0], i-1, pos0, ns[1], ls[1],j-1,CL);

	M[match][i][j] =MAX4  (M[match][i-1][j-1],M[del][i-1][j-1], M[ins][i-1][j-1],M[umatch][i-1][j-1])+sub;
	M[del][i][j]   =MAX2       ((M[match][i-1][j]+gop), M[del][i-1][j])+gep;
	M[ins][i][j]   =MAX2       ((M[match][i][j-1]+gop), M[ins][i][j-1])+gep;
	M[umatch][i][j]=MAX6 (M[match][i-1][j-1]+gop2, M[match][i][j-1]+gop2, M[match][i-1][j]+gop2,M[umatch][i-1][j-1], M[umatch][i-1][j], M[umatch][i][j-1])+gep2;
	}
    }
  return M;
}
int *** backward_so_dp_glocal ( Alignment *A, int *ns, int **ls, int **pos0, int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL)
{
  int i,j;
  int c;
  int sub;
  int ***M;



  M=(int***)declare_arrayN (3, sizeof (int), 5, I+2, J+2);
  for ( i=I+1; i>=0; i--)for (j=J+1; j>=0; j--)for (c=0; c<5; c++)M[c][i][j]=-999999;
  M[match][I+1][J+1]=0;

  for (i=I; i>0; i--){M[ins]  [i][J+1]=i*gep;M[umatch]  [i][J+1]=i*gep2+gop2;}
  for (j=J; j>0; j--){M[del]  [I+1][j]=j*gep;M[umatch]  [I+1][j]=j*gep2+gop2;}

  for (i=I; i>0; i--)
    {
      for ( j=J; j>0; j--)
	{
	sub=(CL->get_dp_cost) (A, pos0, ns[0], ls[0], i-1, pos0, ns[1], ls[1],j-1,CL);

	M[match ][i][j]  =MAX4 ((M[del][i+1][j+1]+gop), (M[ins][i+1][j+1]+gop), M[match][i+1][j+1], M[umatch][i+1][j+1]+gop2)+sub;
	M[del   ][i][j]  =MAX2 (M[match][i+1][j], M[del][i+1][j])+gep;
	M[ins   ][i][j]  =MAX2 (M[match][i][j+1], M[ins][i][j+1])+gep;
	M[umatch][i][j]  =MAX6 (M[match][i+1][j+1], M[match][i+1][j],M[match][i][j+1], M[umatch][i+1][j+1], M[umatch][i+1][j], M[umatch][i][j+1])+gep2;

	}
    }
  return M;
}




/************************************************************************************************************************/
/*                                                                                                                      */
/*                                                                                                                      */
/*                                                     SIMPLE                                                           */
/*                                                                                                                      */
/*                                                                                                                      */
/************************************************************************************************************************/

int *** forward_so_dp ( Alignment *A, int *ns, int **ls, int **pos0,int I, int J,int gop, int gep,int gop2, int gep2,Constraint_list *CL)
{
  int i,j;
  int c;
  int sub;
  int ***M;
  int lgop;



  M=(int***)declare_arrayN (3, sizeof (int), 5, I+1, J+1);
  for ( i=0; i<=I; i++)for (j=0; j<=J; j++)for (c=0; c<3; c++)M[c][i][j]=-999999;

  M[match][0][0]=0;
  for (i=1; i<=I; i++){M[del]  [i][0]=i*gep;}
  for (j=1; j<=J; j++){M[ins]  [0][j]=j*gep;}



  for (i=1; i<=I; i++)
    {
      for ( j=1; j<=J; j++)
	{
	  lgop=(i==I || j==J)?0:gop;
	  sub=(CL->get_dp_cost) (A, pos0, ns[0], ls[0], i-1, pos0, ns[1], ls[1],j-1,CL);

	  M[match][i][j]=MAX3  (M[del][i-1][j-1], M[ins][i-1][j-1], M[match][i-1][j-1])+sub;
	  M[del][i][j]  =MAX ((M[match][i-1][j]+lgop),M[del][i-1][j])+gep;
	  M[ins][i][j]  =MAX ((M[match][i][j-1]+lgop),  M[ins][i][j-1])+gep;
	}

    }

  return M;
  }
int *** backward_so_dp ( Alignment *A, int *ns, int **ls, int **pos0, int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL)
{
  int i,j, a, b;


  int ***M, ***T;


  for (a=0; a<2; a++)
    for (b=0; b<ns[a]; b++)
      {
	invert_string2(A->seq_al[ls[a][b]]);
	invert_string2((CL->S)->seq[A->order[ls[a][b]][0]]);
      }
  T=forward_so_dp(A,ns,ls,pos0, I, J, gop, gep, gop2, gep2, CL);
  for (a=0; a<2; a++)
    for (b=0; b<ns[a]; b++)
      {
	invert_string2(A->seq_al[ls[a][b]]);
	invert_string2((CL->S)->seq[A->order[ls[a][b]][0]]);
      }

  M=(int***)declare_arrayN (3, sizeof (int), 5, I+2, J+2);


  for (i=0; i<=I; i++)
    for (j=0; j<=J; j++)
      {
	M[match][i+1][j+1]=T[match][I-i][J-j];
	M[ins][i+1][j+1]=T[ins][I-i][J-j];
	M[del][i+1][j+1]=T[del][I-i][J-j];
      }
  return M;
}

/************************************************************************************************************************/
/*                                                                                                                      */
/*                                                                                                                      */
/*                                                     BI-PHASIC                                                        */
/*                                                                                                                      */
/*                                                                                                                      */
/************************************************************************************************************************/
int biphasic_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  int i,j,a,b;
  int c;
  int sub;
  int ***m, ***t;
  int M1, D1, D2, I1, I2, LEN;
  int I, J;
  int n=1;
  char **al, **aln, *char_buf;
  int gop1, gop2, gep1, gep2;
  int **pos0;
  int score, trace, ntrace;
  M1=n++; D1=n++; D2=n++; I1=n++, I2=n++;

  I=strlen (A->seq_al[ls[0][0]]);
  J=strlen (A->seq_al[ls[1][0]]);
  m=(int***)declare_arrayN (3, sizeof (int),n, I+1, J+1);
  t=(int***)declare_arrayN (3, sizeof (int),n, I+1, J+1);
  pos0=aln2pos_simple ( A,-1, ns, ls);
  al=declare_char (2, I+J+1);
  for ( i=0; i<=I; i++)for (j=0; j<=J; j++)for (c=0; c<n; c++)m[c][i][j]=-999999;

  gop1=CL->gop*SCORE_K*2;
  gep1=CL->gep*SCORE_K/2;

  gop2=CL->gop*SCORE_K/2;
  gep2=CL->gep*SCORE_K*2;

  m[M1][0][0]=0;
  for (i=1; i<=I; i++){m[I1][i][0]=gep1*i;}
  for (j=1; j<=J; j++){m[D1][0][j]=gep1*j;}

  for (i=1; i<=I; i++){m[I2]  [i][0]=gep2*i;}
  for (j=1; j<=J; j++){m[D2]  [0][j]=gep2*j;}

  for (i=1; i<=I; i++)
    {
      for ( j=1; j<=J; j++)
	{
	  sub=(CL->get_dp_cost) (A, pos0, ns[0], ls[0], i-1, pos0, ns[1], ls[1],j-1,CL);
	  m[M1][i][j]=max_int  (&t[M1][i][j],D1,m[D1][i-1][j-1],I1,m[I1][i-1][j-1], M1, m[M1][i-1][j-1],D2,m[D2][i-1][j-1],I2,m[I2][i-1][j-1], -1)+sub;

	  m[D1][i][j]=max_int  (&t[D1][i][j],M1,(m[M1][i][j-1]+gop1),D1,m[D1][i][j-1], -1)+gep1;
	  m[I1][i][j]=max_int  (&t[I1][i][j],M1,(m[M1][i-1][j]+gop1),I1,m[I1][i-1][j], -1)+gep1;

	  m[D2][i][j]=max_int  (&t[D2][i][j],M1,(m[M1][i][j-1]+gop2),D2,m[D2][i][j-1], -1)+gep2;
	  m[I2][i][j]=max_int  (&t[I2][i][j],M1,(m[M1][i-1][j]+gop2),I2,m[I2][i-1][j], -1)+gep2;
	}
    }

  score=max_int (&trace,M1,m[M1][I][J],D1,m[D1][I][J],I1, m[I1][I][J],D2,m[D2][I][J],I2,m[I2][I][J], -1);
  LEN=0;i=I;j=J;


  trace=t[trace][i][j];
  while (!(i==0 &&j==0))
    {

      ntrace=t[trace][i][j];
      if (i==0)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( j==0)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      else if ( trace==M1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=1;
	  i--; j--;
	  LEN++;
	}

      else if ( trace==D1 || trace==D2)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( trace == I1 || trace==I2)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      trace=ntrace;

    }

  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);
  if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);

  aln=A->seq_al;
  char_buf=(char*) vcalloc (LEN+1, sizeof (char));
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++)
	{
	  int ch=0;
	  for ( b=0; b< LEN; b++)
	    {
	      if (al[c][b]==1)
		char_buf[b]=aln[ls[c][a]][ch++];
	      else
		char_buf[b]='-';
	    }
	  char_buf[b]='\0';
	  sprintf (aln[ls[c][a]],"%s", char_buf);
	}
    }


  A->len_aln=LEN;
  A->nseq=ns[0]+ns[1];
  free_arrayN((void *)m, 3);
  free_arrayN((void *)t, 3);
  vfree (char_buf);
  free_char (al, -1);
  return score;
  }
int *** forward_so_dp_biphasic ( Alignment *A, int *ns, int **ls, int **pos0,int I, int J,int gop1, int gep1,int gop2, int gep2,Constraint_list *CL)
{
  int i,j;
  int c;
  int sub;
  int ***M;
  int match=0, del=1, ins=2;
  int lgop1, lgop2, lgep1, lgep2;

  M=(int***)declare_arrayN (3, sizeof (int), 5, I+1, J+1);

  for ( i=0; i<=I; i++)for (j=0; j<=J; j++)for (c=0; c<5; c++)M[c][i][j]=-999999;

  M[match][0][0]=0;

  for (i=1; i<=I; i++){M[del]  [i][0]=gep1*i+gop1;}
  for (j=1; j<=J; j++){M[ins]  [0][j]=gep1*j+gop1;}

  for (i=1; i<=I; i++){M[del2]  [i][0]=gep2*i+gop2;}
  for (j=1; j<=J; j++){M[ins2]  [0][j]=gep2*j+gop2;}

  for (i=1; i<=I; i++)
    {
      for ( j=1; j<=J; j++)
	{
	  lgop1=(i==I || j==J)?gop1:gop1;
	  lgop2=(i==I || j==J)?gop2:gop2;
	  lgep1=gep1;
	  lgep2=gep2;

	  sub=(CL->get_dp_cost) (A, pos0, ns[0], ls[0], i-1, pos0, ns[1], ls[1],j-1,CL);
	  M[match][i][j]=MAX5  (M[del][i-1][j-1], M[ins][i-1][j-1], M[match][i-1][j-1], M[ins2][i-1][j-1], M[del2][i-1][j-1])+sub;

	  M[del ][i][j] =MAX2 ((M[match][i-1][j]+lgop1), M[del ][i-1][j])+lgep1;
	  M[del2][i][j] =MAX2 ((M[match][i-1][j]+lgop2), M[del2][i-1][j])+lgep2;

	  M[ins ][i][j] =MAX2 ((M[match][i][j-1]+lgop1), M[ins ][i][j-1] )+lgep1;
	  M[ins2][i][j] =MAX2 ((M[match][i][j-1]+lgop2), M[ins2][i][j-1] )+lgep2;
	}
    }
  return M;
  }
int *** backward_so_dp_biphasic ( Alignment *A, int *ns, int **ls, int **pos0, int I, int J, int gop, int gep,int gop2, int gep2,Constraint_list *CL)
{
  int i,j, a, b;


  int ***M, ***T;


  for (a=0; a<2; a++)
    for (b=0; b<ns[a]; b++)
      {
	invert_string2(A->seq_al[ls[a][b]]);
	invert_string2((CL->S)->seq[A->order[ls[a][b]][0]]);
      }
  T=forward_so_dp_biphasic(A,ns,ls,pos0, I, J, gop, gep, gop2, gep2, CL);
  for (a=0; a<2; a++)
    for (b=0; b<ns[a]; b++)
      {
	invert_string2(A->seq_al[ls[a][b]]);
	invert_string2((CL->S)->seq[A->order[ls[a][b]][0]]);
      }

  M=(int***)declare_arrayN (3, sizeof (int), 5, I+2, J+2);


  for (i=0; i<=I; i++)
    for (j=0; j<=J; j++)
      {
	M[match][i+1][j+1]=T[match][I-i][J-j];
	M[ins][i+1][j+1]=T[ins][I-i][J-j];
	M[del][i+1][j+1]=T[del][I-i][J-j];
	M[ins2][i+1][j+1]=T[ins2][I-i][J-j];
	M[del2][i+1][j+1]=T[del2][I-i][J-j];
      }
  free_arrayN(T,3);
  return M;
}


int get_tot_prob (Alignment *A1,Alignment *A2, int *ns, int **ls, int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, Constraint_list *CL,int mode);

int get_tot_prob2 (Alignment *A1,Alignment *A2, int *ns, int **ls, int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, Constraint_list *CL,int mode);

int get_tot_prob3 (Alignment *A1,Alignment *A2, int *ns, int **ls, int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, Constraint_list *CL,int mode);


float * forward_proba_pair_wise  ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *TmatchProb, float ***TinsProb, float **transProb);
float * backward_proba_pair_wise ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *TmatchProb, float ***TinsProb,float **transProb);
float ComputeTotalProbability (int seq1Length, int seq2Length,int NumMatrixTypes, int NumInsertStates,float *forward, float *backward) ;
int ProbabilisticModel (int NumMatrixTypes, int NumInsertStates,float *initDistribMat,float *emitSingle,  float** emitPairs, float *gapOpen, float *gapExtend, float **transMat, float *initialDistribution, float **matchProb, float **insProb, float **transProb);

Constraint_list *ProbaMatrix2CL (Alignment *A, int *ns, int **ls, int NumMatrixTypes, int NumInsertStates, float *forward, float *backward, float thr, Constraint_list *CL);


int proba_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
   static int NumMatrixTypes;
   static int NumInsertStates;
   static float **transMat, **insProb, **matchProb, *initialDistribution, **transProb, **emitPairs, *emitSingle, ***TinsProb, *TmatchProb;
   static int TinsProb_ml, TmatchProb_ml;
   int i, j,I, J;
   float *F, *B;

   int l;
   float thr=0.01;//ProbCons Default
   char *alphabet;


   //Free all the memory
   if (A==NULL)
     {
       free_float (transMat, -1);transMat=NULL;
       free_float (insProb, -1);insProb=NULL;
       free_float (matchProb, -1);matchProb=NULL;
       vfree (initialDistribution); initialDistribution=NULL;
       free_float (transProb, -1);transProb=NULL;
       free_float (emitPairs, -1);emitPairs=NULL;
       vfree (emitSingle);emitSingle=NULL;


       free_arrayN((void***)TinsProb, 3);TinsProb=NULL;
       vfree (TmatchProb);TmatchProb=NULL;
       TinsProb_ml=0; TmatchProb_ml=0;

       forward_proba_pair_wise (NULL, NULL, 0,0,NULL,NULL,NULL,NULL,NULL);
       backward_proba_pair_wise (NULL, NULL, 0,0,NULL,NULL,NULL,NULL,NULL);
       ProbaMatrix2CL(NULL, NULL, NULL, 0, 0, NULL, NULL, 0, NULL);
       return 0;
     }

   if (!transMat && (strm (retrieve_seq_type(), "DNA")))
     {
     static float **p;
     static float *s;
     NumInsertStates=1;
     NumMatrixTypes=3;
     if (!p)
       {
	 int l,a,b;
	 l=strlen (DNAalphabetDefault);
	 p=declare_float (l,l);
	 s=(float*)vcalloc (l, sizeof (float));
	 for (a=0; a<l; a++)
	   {
	     s[a]=DNAemitSingleDefault[a];
	     for (b=0; b<l; b++)
	       p[a][b]=RNAemitPairsDefault[a][b];
	   }
       }
     p=get_emitPairs (CL->method_matrix, DNAalphabetDefault,p,s);
     alphabet=RNAalphabetDefault;
     emitPairs=declare_float (256, 256);
     emitSingle=(float*)vcalloc (256, sizeof (float));
     for (i=0; i<256; i++)
       {
	 emitSingle[i]=1e-5;
	 for (j=0; j<256; j++)
	   emitPairs[i][j]=1e-10;
       }
     l=strlen (alphabet);

     for (i=0; i<l; i++)
       {
	 int C1,c1, C2,c2;
	 c1=tolower(alphabet[i]);
	 C1=toupper(alphabet[i]);
	 emitSingle[c1]=s[i];
	 emitSingle[C1]=s[i];
	 for (j=0; j<=i; j++)
	   {
	     c2=tolower(alphabet[j]);
	     C2=toupper(alphabet[j]);

	     emitPairs[c1][c2]=p[i][j];
	     emitPairs[C1][c2]=p[i][j];
	     emitPairs[C1][C2]=p[i][j];
	     emitPairs[c1][C2]=p[i][j];
	     emitPairs[c2][c1]=p[i][j];
	     emitPairs[C2][c1]=p[i][j];
	     emitPairs[C2][C1]=p[i][j];
	     emitPairs[c2][C1]=p[i][j];
	   }
       }


     transMat=declare_float (2*NumInsertStates+1, 2*NumInsertStates+1);
     transProb=declare_float (2*NumInsertStates+1,2* NumInsertStates+1);
     insProb=declare_float (256,NumMatrixTypes);
     matchProb=declare_float (256, 256);
     initialDistribution=(float*)vcalloc (2*NumMatrixTypes+1, sizeof (float));

     ProbabilisticModel (NumMatrixTypes,NumInsertStates,initDistrib2Default, emitSingle,emitPairs,DNAgapOpen2Default,DNAgapExtend2Default, transMat,initialDistribution,matchProb, insProb,transProb);

     }
   else if (!transMat && (strm (retrieve_seq_type(), "RNA")))
     {
       static float **p;
       static float *s;
       NumInsertStates=2;
       NumMatrixTypes=5;

       if (!p)
	 {
	   int l,a,b;
	   l=strlen (RNAalphabetDefault);
	   p=declare_float (l,l);
	   s=(float*)vcalloc (l, sizeof (float));
	   for (a=0; a<l; a++)
	     {
	       s[a]=RNAemitSingleDefault[a];
	       for (b=0; b<l; b++)
		 p[a][b]=RNAemitPairsDefault[a][b];
	     }
	 }
       p=get_emitPairs (CL->method_matrix, RNAalphabetDefault,p,s);
       alphabet=RNAalphabetDefault;
       emitPairs=declare_float (256, 256);
       emitSingle=(float*)vcalloc (256, sizeof (float));
       for (i=0; i<256; i++)
	 {
	   emitSingle[i]=1e-5;
	   for (j=0; j<256; j++)
	     emitPairs[i][j]=1e-10;
	 }
       l=strlen (alphabet);

       for (i=0; i<l; i++)
	 {
	   int C1,c1, C2,c2;
	   c1=tolower(alphabet[i]);
	   C1=toupper(alphabet[i]);
	   emitSingle[c1]=s[i];
	   emitSingle[C1]=s[i];
	   for (j=0; j<=i; j++)
	     {
	       c2=tolower(alphabet[j]);
	       C2=toupper(alphabet[j]);

	       emitPairs[c1][c2]=p[i][j];
	       emitPairs[C1][c2]=p[i][j];
	       emitPairs[C1][C2]=p[i][j];
	       emitPairs[c1][C2]=p[i][j];
	       emitPairs[c2][c1]=p[i][j];
	       emitPairs[C2][c1]=p[i][j];
	       emitPairs[C2][C1]=p[i][j];
	       emitPairs[c2][C1]=p[i][j];
	     }
	 }


       transMat=declare_float (2*NumInsertStates+1, 2*NumInsertStates+1);
       transProb=declare_float (2*NumInsertStates+1,2* NumInsertStates+1);
       insProb=declare_float (256,NumMatrixTypes);
       matchProb=declare_float (256, 256);
       initialDistribution=(float*)vcalloc (2*NumMatrixTypes+1, sizeof (float));

       ProbabilisticModel (NumMatrixTypes,NumInsertStates,initDistrib2Default, emitSingle,emitPairs,RNAgapOpen2Default,RNAgapExtend2Default, transMat,initialDistribution,matchProb, insProb,transProb);
     }
   else if ( !transMat && strm (retrieve_seq_type(), "PROTEIN"))
     {
       static float **p;
       static float *s;
       NumInsertStates=2;
       NumMatrixTypes=5;
       if (atoigetenv ("NOBIPHASIC"))
	 {
	   NumInsertStates=1;
	   NumMatrixTypes=3;
	 }
       if (!p)
	 {
	   int l,a,b;
	   l=strlen (alphabetDefault);
	   p=declare_float (l,l);
	   s=(float*)vcalloc (l, sizeof (float));
	   for (a=0; a<l; a++)
	     {
	       s[a]=emitSingleDefault[a];
	       for (b=0; b<l; b++)
		 p[a][b]=emitPairsDefault[a][b];
	     }
	 }
       p=get_emitPairs (CL->method_matrix, alphabetDefault,p,s);
       alphabet=alphabetDefault;
       emitPairs=declare_float (256, 256);
       emitSingle=(float*)vcalloc (256, sizeof (float));
       for (i=0; i<256; i++)
	 {
	   //emitSingle[i]=1e-5;
	   emitSingle[i]=1;
	   for (j=0; j<256; j++)
	     //emitPairs[i][j]=1e-10;
	     emitPairs[i][j]=1;

	 }
       l=strlen (alphabet);

       for (i=0; i<l; i++)
	 {
	   int C1,c1, C2,c2;
	   c1=tolower(alphabet[i]);
	   C1=toupper(alphabet[i]);
	   emitSingle[c1]=s[i];
	   emitSingle[C1]=s[i];
	   for (j=0; j<=i; j++)
	     {
	       c2=tolower(alphabet[j]);
	       C2=toupper(alphabet[j]);

	       emitPairs[c1][c2]=p[i][j];
	       emitPairs[C1][c2]=p[i][j];
	       emitPairs[C1][C2]=p[i][j];
	       emitPairs[c1][C2]=p[i][j];
	       emitPairs[c2][c1]=p[i][j];
	       emitPairs[C2][c1]=p[i][j];
	       emitPairs[C2][C1]=p[i][j];
	       emitPairs[c2][C1]=p[i][j];

	     }
	 }


       transMat=declare_float (2*NumInsertStates+1, 2*NumInsertStates+1);
       transProb=declare_float (2*NumInsertStates+1,2* NumInsertStates+1);
       insProb=declare_float (256,NumMatrixTypes);
       matchProb=declare_float (256, 256);
       initialDistribution=(float*)vcalloc (2*NumMatrixTypes+1, sizeof (float));
       if (atoigetenv ("NOBIPHASIC"))
	 ProbabilisticModel (NumMatrixTypes,NumInsertStates,initDistrib1Default, emitSingle,emitPairs,gapOpen1Default,gapExtend1Default, transMat,initialDistribution,matchProb, insProb,transProb);
       else
	 ProbabilisticModel (NumMatrixTypes,NumInsertStates,initDistrib2Default, emitSingle,emitPairs,gapOpen2Default,gapExtend2Default, transMat,initialDistribution,matchProb, insProb,transProb);
     }

   I=strlen (A->seq_al[ls[0][0]]);
   J=strlen (A->seq_al[ls[1][0]]);
   //TmatchProb=vcalloc ((I+1)*(J+1), sizeof (float));
   //TinsProb=declare_arrayN (3, sizeof (float),2,NumMatrixTypes,MAX(I,J)+1);

   l=(I+1)*(J+1);
   if (l>TmatchProb_ml)
     {
       TmatchProb_ml=l;
       if (TmatchProb)TmatchProb=(float*)vrealloc(TmatchProb,TmatchProb_ml*sizeof (float));
       else TmatchProb=(float*)vcalloc ( l, sizeof (float));
     }
   l=MAX(I,J)+1;
   if ( l>TinsProb_ml)
     {
       TinsProb_ml=l;
       if (TinsProb)free_arrayN (TinsProb, 3);
       TinsProb=(float***)declare_arrayN (3, sizeof (float),2,NumMatrixTypes,TinsProb_ml);
     }
   if (strm (retrieve_seq_type(), "RNA"))
     get_tot_prob (A,A, ns,ls,NumMatrixTypes, matchProb, insProb,TmatchProb,TinsProb, CL, SEQUENCE);
   else
     get_tot_prob2 (A,A, ns,ls,NumMatrixTypes, matchProb, insProb,TmatchProb,TinsProb, CL, SEQUENCE);
   F=forward_proba_pair_wise (A->seq_al[ls[0][0]], A->seq_al[ls[1][0]], NumMatrixTypes,NumInsertStates,transMat, initialDistribution,TmatchProb,TinsProb, transProb);
   B=backward_proba_pair_wise (A->seq_al[ls[0][0]], A->seq_al[ls[1][0]], NumMatrixTypes,NumInsertStates,transMat, initialDistribution,TmatchProb,TinsProb, transProb);
   A->CL=ProbaMatrix2CL(A,ns, ls,NumMatrixTypes,NumInsertStates, F, B, thr,CL);

   //free_proba_pair_wise();
   return 1;
   }

void free_proba_pair_wise ()
{
  proba_pair_wise (NULL, NULL, NULL, NULL);
}

int get_tot_prob (Alignment *A1,Alignment *A2, int *ns, int **ls, int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, Constraint_list *CL, int mode)
{
  int i, j, a, b, c,d, k, n,n1,n2, ij;
  int  c1, c2;
  int I, J;
  int ***VA1,***VA2, *observed, index;
  char *ss1=NULL;
  char *ss2=NULL;
  int uss=0;
  static int gtp=0;
  //Pre-computation of the pairwise scores in order to use potential profiles
  //The profiles are vectorized AND Compressed so that the actual alphabet size (proteins/DNA) does not need to be considered
  int use_cons=atoigetenv ("KM_COFFEE_PRF_CONS");
  use_cons=0;
  if (mode==SEQUENCE)
    {
      int s1, s2;
      int *nns, **nls;
      Alignment *NA1, *NA2;
      char *sst1;
      char *sst2;


      nns=(int*)vcalloc ( 2, sizeof (int));
      nls=(int**)vcalloc (2, sizeof (int*));

      s1=A1->order[ls[0][0]][0];
      s2=A2->order[ls[1][0]][0];
      NA1=seq2R_template_profile (CL->S,s1);
      NA2=seq2R_template_profile (CL->S,s2);

      sst1=seq2T_template_string((CL->S),s1);
      sst2=seq2T_template_string((CL->S),s2);


      if (NA1 || NA2)
	{
	  if (NA1)
	    {
	      nns[0]=NA1->nseq;
	      nls[0]=(int*)vcalloc (NA1->nseq, sizeof (int));
	      for (a=0; a<NA1->nseq; a++)
		nls[0][a]=a;
	      NA1->seq_al[NA1->nseq]=sst1;
	      sprintf (NA1->name[NA1->nseq], "sst1");
	    }
	  else
	    {
	    NA1=A1;
	    nns[0]=ns[0];
	    nls[0]=(int*)vcalloc (ns[0], sizeof (int));
	    for (a=0; a<ns[0]; a++)
	      nls[0][a]=ls[0][a];
	    }

	  if (NA2)
	    {
	      nns[1]=NA2->nseq;
	      nls[1]=(int*)vcalloc (NA2->nseq, sizeof (int));
	      for (a=0; a<NA2->nseq; a++)
		nls[1][a]=a;
	      NA2->seq_al[NA2->nseq]=sst2;
	      sprintf (NA2->name[NA2->nseq], "sst2");
	    }
	  else
	    {
	      NA2=A2;
	      nns[1]=ns[1];
	      nls[1]=(int*)vcalloc (ns[1], sizeof (int));
	      for (a=0; a<ns[1]; a++)
		nls[1][a]=ls[1][a];
	    }

	  get_tot_prob (NA1, NA2, nns, nls, nstates, matchProb, insProb, TmatchProb, TinsProb, CL,PROFILE);
	  vfree (nns); free_int (nls,-1);
	  return 1;
	}
    }
  if ( A1!=A2)
    {
      if (strm (A1->name[A1->nseq], "sst1"))ss1=A1->seq_al[A1->nseq];
      if (strm (A2->name[A2->nseq], "sst2"))ss2=A2->seq_al[A2->nseq];
      uss=(ss1&&ss2)?1:0;
    }
  else
    uss=0;

  I=strlen (A1->seq_al[ls[0][0]]);
  J=strlen (A2->seq_al[ls[1][0]]);



  //get Ins for I
  for (i=1; i<=I; i++)
    {
      for (k=0; k<nstates; k++)
	{
	  TinsProb[0][k][i]=0;
	  if (use_cons)
	    {
	      c1=A1->seq_al[A1->nseq][i-1];
	      if (c1!='-') TinsProb[0][k][i]+=insProb[c1][k];
	    }
	  else
	    {
	      for (n=0,b=0; b<ns[0]; b++)
		{
		  c1=A1->seq_al[ls[0][b]][i-1];
		  if (c1!='-')
		    {
		      TinsProb[0][k][i]+=insProb[c1][k];
		      n++;
		    }
		}
	      if (n)TinsProb[0][k][i]/=n;
	    }
	}
    }
  //Get Ins for J
  for (j=1; j<=J; j++)
    {
      for (k=0; k<nstates; k++)
	{
	  TinsProb[1][k][j]=0;
	  if (use_cons)
	    {
	      c2=A2->seq_al[A2->nseq][j-1];
	      if (c2!='-') TinsProb[1][k][j]+=insProb[c2][k];
	    }
	  else
	    {
	      for (n=0,b=0; b<ns[1]; b++)
		{
		  c2=A2->seq_al[ls[1][b]][j-1];
		  if (c2!='-')
		    {
		      TinsProb[1][k][j]+=insProb[c2][k];
		      n++;
		    }
		}
	      if (n)TinsProb[1][k][j]/=n;
	    }
	}
    }
  
  observed=(int*)vcalloc ( 26, sizeof (int));
  VA1=(int***)declare_arrayN (3, sizeof (int),2,26,I);
  for (i=0; i<I; i++)
    {
      if (use_cons)
	{
	  c1=tolower(A1->seq_al[A1->nseq][i]);
	  if ( c1=='-' || c1=='.' || c1=='~') VA1[0][0][i]=-1;
	  else
	    {
	      c1-='a';
	      VA1[0][0][i]=c1;
	      VA1[1][0][i]++;
	      VA1[0][1][i]=-1;
	    }
	}
      else
	{
	  for (index=0, b=0; b<ns[0]; b++)
	    {
	      int in;
	      c1=tolower(A1->seq_al[ls[0][b]][i]);
	      if ( c1=='-' || c1=='.' || c1=='~')continue;
	      c1-='a';
	      
	      if (!(in=observed[c1])){in=observed[c1]=++index;}
	      
	      VA1[0][in-1][i]=c1;
	      VA1[1][in-1][i]++;
	    }
	  
      
	  VA1[0][index][i]=-1;
	  for (b=0; b<26; b++)observed[b]=0;
	}
    }

  VA2=(int***)declare_arrayN (3, sizeof (int),2,26,J);
  for (i=0; i<J; i++)
    {
      if (use_cons)
	{
	  c1=tolower(A2->seq_al[A2->nseq][i]);
	  if ( c1=='-' || c1=='.' || c1=='~')VA2[0][0][i]=-1;
	  else
	    {
	      c1-='a';
	      VA2[0][0][i]=c1;
	      VA2[1][0][i]++;
	      VA2[0][1][i]=-1;
	    }
	}
      else
	{
	  for (index=0, b=0; b<ns[1]; b++)
	    {
	      int in;
	      
	      c1=tolower(A2->seq_al[ls[1][b]][i]);
	      if ( c1=='-')continue;
	      c1-='a';
	      
	      if (!(in=observed[c1])){in=observed[c1]=++index;}
	      
	      VA2[0][in-1][i]=c1;
	      VA2[1][in-1][i]++;
	    }
	
	  VA2[0][index][i]=-1;
	  for (b=0; b<26; b++)observed[b]=0;
	}
    }
  vfree (observed);

  for ( ij=0,i=0; i<=I; i++)
    {
      for ( j=0; j<=J ; j++, ij++)
	{
	  n=0;
	  TmatchProb[ij]=0;
	  if (i==0 || j==0);
	  else
	    {
	      float sfac;
	      if (!uss)sfac=1;
	      else if (ss1[i-1]!=ss2[j-1])sfac=1;
	      else if (ss1[i-1]==ss2[j-1])sfac=1;
	      else sfac=1;


	      c=0;
	      while (VA1[0][c][i-1]!=-1)
		{
		  c1=VA1[0][c][i-1]+'a';
		  n1=VA1[1][c][i-1];
		  d=0;
		  while (VA2[0][d][j-1]!=-1)
		    {
		      c2=VA2[0][d][j-1]+'a';
		      n2=VA2[1][d][j-1];
		      TmatchProb[ij]+=matchProb[c1][c2]*(double)n1*(double)n2*sfac;
		      n+=n1*n2;
		      d++;
		    }
		  c++;
		}
	    }
	  if (n)TmatchProb[ij]/=n;
	}
    }

  free_arrayN ((void **)VA1, 3);
  free_arrayN ((void **)VA2, 3);
  return 1;
}
int get_tot_prob2 (Alignment *A1,Alignment *A2, int *ns, int **ls, int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, Constraint_list *CL, int mode)
{
  static double **prf1;
  static double **prf2;
  int i, I, j, J, k, ij,r,r1,r2; 
  int *lu;
 
  if (mode==SEQUENCE)
    {
      int s1, s2, a;
      int *nns, **nls;
      Alignment *NA1, *NA2;
      char *sst1;
      char *sst2;

      
      nns=(int*)vcalloc ( 2, sizeof (int));
      nls=(int**)vcalloc (2, sizeof (int*));

      s1=A1->order[ls[0][0]][0];
      s2=A2->order[ls[1][0]][0];
      NA1=seq2R_template_profile (CL->S,s1);
      NA2=seq2R_template_profile (CL->S,s2);

      sst1=seq2T_template_string((CL->S),s1);
      sst2=seq2T_template_string((CL->S),s2);


      if (NA1 || NA2)
	{
	  if (NA1)
	    {
	      nns[0]=NA1->nseq;
	      nls[0]=(int*)vcalloc (NA1->nseq, sizeof (int));
	      for (a=0; a<NA1->nseq; a++)
		nls[0][a]=a;
	      NA1->seq_al[NA1->nseq]=sst1;
	      sprintf (NA1->name[NA1->nseq], "sst1");
	    }
	  else
	    {
	    NA1=A1;
	    nns[0]=ns[0];
	    nls[0]=(int*)vcalloc (ns[0], sizeof (int));
	    for (a=0; a<ns[0]; a++)
	      nls[0][a]=ls[0][a];
	    }

	  if (NA2)
	    {
	      nns[1]=NA2->nseq;
	      nls[1]=(int*)vcalloc (NA2->nseq, sizeof (int));
	      for (a=0; a<NA2->nseq; a++)
		nls[1][a]=a;
	      NA2->seq_al[NA2->nseq]=sst2;
	      sprintf (NA2->name[NA2->nseq], "sst2");
	    }
	  else
	    {
	      NA2=A2;
	      nns[1]=ns[1];
	      nls[1]=(int*)vcalloc (ns[1], sizeof (int));
	      for (a=0; a<ns[1]; a++)
		nls[1][a]=ls[1][a];
	    }
	 
	  get_tot_prob2 (NA1, NA2, nns, nls, nstates, matchProb, insProb, TmatchProb, TinsProb, CL,PROFILE);
	  vfree (nns); free_int (nls,-1);
	  return 1;
	}
    }
  else if (mode == PROFILE)
    {
      static int display_mode;
      if (!display_mode)
	{
	  fprintf ( stderr, "\n! profile/profile alignment --- use proba_pair\n");
	  display_mode=1;
	}
    }
        
  I=strlen (A1->seq_al[ls[0][0]]);
  J=strlen (A2->seq_al[ls[1][0]]);
     
  lu=dirichlet_code2aa_lu();
  
  prf1=aln2prf (A1, ns[0], ls[0], I, prf1);
  prf2=aln2prf (A2, ns[1], ls[1], J, prf2);

    
  //get Ins for I
  for (i=1; i<=I; i++)
    {
      for (k=0; k<nstates; k++)
	{
	  TinsProb[0][k][i]=0;
	  for (r=0; r<20; r++)
	    {
	      TinsProb[0][k][i]+=(float)prf1[i-1][r]*insProb[lu[r]][k];
	    }
	}
    }
  

  
  
  //Get Ins for J
  for (j=1; j<=J; j++)
    {
      for (k=0; k<nstates; k++)
	{
	  TinsProb[1][k][j]=0;
	  for (r=0; r<20; r++)
	    {
	      TinsProb[1][k][j]+=(float)prf2[j-1][r]*insProb[lu[r]][k]; 
	    }
	}
    }
  for (ij=0,i=0; i<=I; i++)
    for (j=0; j<=J; j++, ij++)
      {
	float tot=0,f;
	if (i==0 || j==0)continue;
	TmatchProb[ij]=0;
	for (tot=0,r1=0; r1<20; r1++)
	  {
	    for (r2=0; r2<20; r2++)
	      {
		f=(float)prf1[i-1][r1]*(float)prf2[j-1][r2];
		TmatchProb[ij]+=matchProb[lu[r1]][lu[r2]]*f;
	      }
	  }
      }
  
  return 1;
}
int get_tot_prob3 (Alignment *A1,Alignment *A2, int *ns, int **ls, int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, Constraint_list *CL, int mode)
{
  static double **prf1;
  static double **prf2;
  static double **dmx1;
  static double **dmx2;
  
  int i, I, j, J, k, ij,r,r1,r2; 
  int *lu;
   
  if (mode==SEQUENCE)
    {
      int s1, s2, a;
      int *nns, **nls;
      Alignment *NA1, *NA2;
      char *sst1;
      char *sst2;


      nns=(int*)vcalloc ( 2, sizeof (int));
      nls=(int**)vcalloc (2, sizeof (int*));

      s1=A1->order[ls[0][0]][0];
      s2=A2->order[ls[1][0]][0];
      NA1=seq2R_template_profile (CL->S,s1);
      NA2=seq2R_template_profile (CL->S,s2);

      sst1=seq2T_template_string((CL->S),s1);
      sst2=seq2T_template_string((CL->S),s2);


      if (NA1 || NA2)
	{
	  if (NA1)
	    {
	      nns[0]=NA1->nseq;
	      nls[0]=(int*)vcalloc (NA1->nseq, sizeof (int));
	      for (a=0; a<NA1->nseq; a++)
		nls[0][a]=a;
	      NA1->seq_al[NA1->nseq]=sst1;
	      sprintf (NA1->name[NA1->nseq], "sst1");
	    }
	  else
	    {
	    NA1=A1;
	    nns[0]=ns[0];
	    nls[0]=(int*)vcalloc (ns[0], sizeof (int));
	    for (a=0; a<ns[0]; a++)
	      nls[0][a]=ls[0][a];
	    }

	  if (NA2)
	    {
	      nns[1]=NA2->nseq;
	      nls[1]=(int*)vcalloc (NA2->nseq, sizeof (int));
	      for (a=0; a<NA2->nseq; a++)
		nls[1][a]=a;
	      NA2->seq_al[NA2->nseq]=sst2;
	      sprintf (NA2->name[NA2->nseq], "sst2");
	    }
	  else
	    {
	      NA2=A2;
	      nns[1]=ns[1];
	      nls[1]=(int*)vcalloc (ns[1], sizeof (int));
	      for (a=0; a<ns[1]; a++)
		nls[1][a]=ls[1][a];
	    }

	  get_tot_prob3 (NA1, NA2, nns, nls, nstates, matchProb, insProb, TmatchProb, TinsProb, CL,PROFILE);
	  vfree (nns); free_int (nls,-1);
	  return 1;
	}
    }
  
  I=strlen (A1->seq_al[ls[0][0]]);
  J=strlen (A2->seq_al[ls[1][0]]);
     
  lu=dirichlet_code2aa_lu();
  
  prf1=aln2prf (A1, ns[0], ls[0], I, prf1);
  dmx1=prf2dmx (prf1, dmx1, I);
  prf2=aln2prf (A2, ns[1], ls[1], J, prf2);
  dmx2=prf2dmx (prf2, dmx2, J);
  
    
  //get Ins for I
  for (i=1; i<=I; i++)
    {
      for (k=0; k<nstates; k++)
	{
	  TinsProb[0][k][i]=0;
	  for (r=0; r<20; r++)
	    {
	      TinsProb[0][k][i]+=(float)prf1[i-1][r]*insProb[lu[r]][k];
	    }
	}
    }
  
  //Get Ins for J
  for (j=1; j<=J; j++)
    {
      for (k=0; k<nstates; k++)
	{
	  TinsProb[1][k][j]=0;
	  for (r=0; r<20; r++)
	    {
	      TinsProb[1][k][j]+=(float)prf2[j-1][r]*insProb[lu[r]][k];
	    }
	}
    }
  
  for (ij=0,i=0; i<=I; i++)
    for (j=0; j<=J; j++, ij++)
      {
	if (i==0 || j==0)continue;
	TmatchProb[ij]=0;
	for (r1=0; r1<20; r1++)
	  for (r2=0; r2<20; r2++)
	    {
	      TmatchProb[ij]+=(prf1[i-1][r1]*prf2[j-1][r2])*(dmx2[j-1][r1]+dmx1[i-1][r2]);
	    }
	
      }
  
  return 1;
}



Constraint_list *ProbaMatrix2CL (Alignment *A, int *ns, int **ls, int NumMatrixTypes, int NumInsertStates, float *forward, float *backward, float thr, Constraint_list *CL)
{
  float totalProb;
  int ij, i, j,k, I, J, s1, s2;
  static int *entry;
  static int **list;
  static int list_max;
  int sim;
  int list_size;
  int list_n;
  int old_n=0;
  double v;
  int a;
  static float F=4; //potential number of full suboptimal alignmnents incorporated in the library
  static int tot_old, tot_new;
  
  if (!A)
    {
      free_int (list, -1);list=NULL;
      list_max=0;

      vfree(entry); entry=NULL;
      return NULL;
    }
  
  I=strlen (A->seq_al[ls[0][0]]);
  J=strlen (A->seq_al[ls[1][0]]);

  

  s1=name_is_in_list (A->name[ls[0][0]], (CL->S)->name, (CL->S)->nseq, 100);
  s2=name_is_in_list (A->name[ls[1][0]], (CL->S)->name, (CL->S)->nseq, 100);
  
  list_size=I*J;

  if ( list_max<list_size)
    {
      free_int (list, -1);
      list_max=list_size;
      list=declare_int (list_max, 3);
    }


  totalProb = ComputeTotalProbability (I,J,NumMatrixTypes, NumInsertStates,forward, backward);

  ij = 0;
  for (list_n=0,ij=0,i =0; i <= I; i++)
    {
      for (j =0; j <= J; j++, ij+=NumMatrixTypes)
	{
	  v= EXP (MIN(LOG_ONE,(forward[ij] + backward[ij] - totalProb)));
	  if (v>thr)//Conservative reduction of the list size to speed up the sorting
	    {
	      list[list_n][0]=i;
	      list[list_n][1]=j;
	      list[list_n][2]=(int)((float)v*(float)NORM_F);
	      list_n++;
	    }
	  if (v>0.01)old_n++;
	}
    }

  sort_int_inv (list, 3, 2, 0, list_n-1);
  if (!entry)entry=(int*)vcalloc ( CL->entry_len+1, CL->el_size);

  list_n=MIN(list_n,(F*MIN(I,J)));
  for (i=0; i<list_n; i++)
    {
       entry[SEQ1]=s1;
       entry[SEQ2]=s2;
       entry[R1]  =list[i][0];
       entry[R2]  =list[i][1];
       entry[WE]  =list[i][2];
       entry[CONS]=1;
       add_entry2list (entry,A->CL);
    }
  tot_new+=list_n;
  tot_old+=old_n;
  // HERE ("LIB_SIZE NEW: %d (new) %d (old) [%.2f]", list_n, old_n, (float)tot_new/(float)tot_old);
  return A->CL;
}



float ComputeTotalProbability (int seq1Length, int seq2Length,int NumMatrixTypes, int NumInsertStates,float *forward, float *backward)
{

    float totalForwardProb = LOG_ZERO;
    float totalBackwardProb = LOG_ZERO;
    int k;

    for (k = 0; k < NumMatrixTypes; k++)
      {
      LOG_PLUS_EQUALS (&totalForwardProb,forward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + backward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
      }

    totalBackwardProb =forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)];

    for (k = 0; k < NumInsertStates; k++)
      {
      LOG_PLUS_EQUALS (&totalBackwardProb,forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)]);
      LOG_PLUS_EQUALS (&totalBackwardProb,forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)]);
      }
    
    return (totalForwardProb + totalBackwardProb) / 2;
  }


float * backward_proba_pair_wise ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *matchProb, float ***insProb, float **transProb)
{
  static float *backward;
  static int max_l;


  int k, i, j,ij, i1j1, i1j, ij1,a, l, seq1Length, seq2Length, m;
  char c1, c2;
  char *iter1, *iter2;

  if (!seq1)
    {
      vfree (backward);
      backward=NULL; max_l=0;
      return NULL;
    }

  iter1=seq1-1;
  iter2=seq2-1;
  seq1Length=strlen (seq1);
  seq2Length=strlen (seq2);
  l=(seq1Length+1)*(seq2Length+1)*NumMatrixTypes;

  if (!backward)
    {
      backward=(float*)vcalloc (l, sizeof (float));
      max_l=l;
    }
  else if (max_l<l)
    {
      backward=(float*)vrealloc (backward, l*sizeof(float));
      max_l=l;
    }

  for (a=0; a<l; a++)backward[a]=LOG_ZERO;

  for (k = 0; k < NumMatrixTypes; k++)
    backward[NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1) + k] = initialDistribution[k];
  
  //Difference with Probcons: this emission is not added to the bward
  backward[NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1) + 0]+=matchProb[(seq1Length+1) * (seq2Length+1) - 1];
  // remember offset for each index combination
  ij = (seq1Length+1) * (seq2Length+1) - 1;

  i1j = ij + seq2Length + 1;
  ij1 = ij + 1;
  i1j1 = ij + seq2Length + 2;
  ij *= NumMatrixTypes;
  i1j *= NumMatrixTypes;
  ij1 *= NumMatrixTypes;
  i1j1 *= NumMatrixTypes;

  // compute backward scores
  for (i = seq1Length; i >= 0; i--)
    {
      c1 = (i == seq1Length) ? '~' : (unsigned char) iter1[i+1];
      for (j = seq2Length; j >= 0; j--)
	{
	  c2 = (j == seq2Length) ? '~' : (unsigned char) iter2[j+1];

	  if (i < seq1Length && j < seq2Length)
	    {
	      m=((i+1)*(seq2Length+1))+j+1;//The backward and the forward are offset by 1
	      float ProbXY = backward[0 + i1j1] + matchProb[m];


	      for (k = 0; k < NumMatrixTypes; k++)
		{
		  LOG_PLUS_EQUALS (&backward[k + ij], ProbXY + transProb[k][0]);
		}
	    }
	  if (i < seq1Length)
	    {
	      for (k = 0; k < NumInsertStates; k++)
		{
		LOG_PLUS_EQUALS (&backward[0 + ij], backward[2*k+1 + i1j] + insProb[0][k][i+1] + transProb[0][2*k+1]);
		LOG_PLUS_EQUALS (&backward[2*k+1 + ij], backward[2*k+1 + i1j] + insProb[0][k][i+1] + transProb[2*k+1][2*k+1]);
		}
	    }
        if (j < seq2Length)
	  {
	    for (k = 0; k < NumInsertStates; k++)
	      {
		//+1 because the backward and the forward are offset by 1
		LOG_PLUS_EQUALS (&backward[0 + ij], backward[2*k+2 + ij1] + insProb[1][k][j+1] + transProb[0][2*k+2]);
		LOG_PLUS_EQUALS (&backward[2*k+2 + ij], backward[2*k+2 + ij1] + insProb[1][k][j+1] + transProb[2*k+2][2*k+2]);
	      }
	  }

        ij -= NumMatrixTypes;
        i1j -= NumMatrixTypes;
        ij1 -= NumMatrixTypes;
        i1j1 -= NumMatrixTypes;
	}
    }
 
  return backward;
}

float * forward_proba_pair_wise ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *matchProb, float ***insProb, float **transProb)
{
  static float *forward;
  static int max_l;
  int k, i, j,ij, i1j1, i1j, ij1, seq1Length, seq2Length, m;
  char *iter1, *iter2;
  int l,a;

  if (!seq1)
    {
      vfree (forward);
      forward=NULL; max_l=0;
      return NULL;
    }
  iter1=seq1-1;
  iter2=seq2-1;
  seq1Length=strlen (seq1);
  seq2Length=strlen (seq2);
  l=(seq1Length+1)*(seq2Length+1)*NumMatrixTypes;

  if (!forward)
    {
      forward=(float*)vcalloc (l, sizeof (float));
      max_l=l;
    }
  else if (max_l<l)
    {
      forward=(float*)vrealloc (forward, l*sizeof(float));
      max_l=l;
    }
  for (a=0; a<l; a++)forward[a]=LOG_ZERO;


  forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] = initialDistribution[0] + matchProb[seq2Length+2];

  for (k = 0; k < NumInsertStates; k++)
    {
      forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] = initialDistribution[2*k+1] + insProb[0][k][1];
      forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] = initialDistribution[2*k+2] + insProb[1][k][1];
    }

  // remember offset for each index combination
    ij = 0;
    i1j = -seq2Length - 1;
    ij1 = -1;
    i1j1 = -seq2Length - 2;

    ij *= NumMatrixTypes;
    i1j *= NumMatrixTypes;
    ij1 *= NumMatrixTypes;
    i1j1 *= NumMatrixTypes;


    // compute forward scores
    for (m=0,i = 0; i <= seq1Length; i++)
      {
	for (j = 0; j <= seq2Length; j++, m++)
	  {
	  if (i > 1 || j > 1)
	    {
	    if (i > 0 && j > 0)
	      {
		//Sum over all possible alignments
		forward[0 + ij] = forward[0 + i1j1] + transProb[0][0];
		for (k = 1; k < NumMatrixTypes; k++)
		  {
		    LOG_PLUS_EQUALS (&forward[0 + ij], forward[k + i1j1] + transProb[k][0]);
		  }
		forward[0 + ij] += matchProb[m];
	      }
	    if ( i > 0)
	      {
	      for (k = 0; k < NumInsertStates; k++)
		{
		  forward[2*k+1 + ij] = insProb[0][k][i] + LOG_ADD (forward[0 + i1j] + transProb[0][2*k+1],forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1]);
		}
	      }
	  if (j > 0)
	    {
	    for (k = 0; k < NumInsertStates; k++)
	      {
		forward[2*k+2 + ij] = insProb[1][k][j] +LOG_ADD (forward[0 + ij1] + transProb[0][2*k+2],forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2]);
	      }
	    }
	  }

        ij += NumMatrixTypes;
        i1j += NumMatrixTypes;
        ij1 += NumMatrixTypes;
        i1j1 += NumMatrixTypes;
      }

    }
    return forward;
  }
int ProbabilisticModel (int NumMatrixTypes, int NumInsertStates,float *initDistribMat,float *emitSingle,  float **emitPairs, float *gapOpen, float *gapExtend, float **transMat, float *initialDistribution, float **matchProb, float **insProb, float **transProb)
{


    // build transition matrix
  int i, j;

  //Maybe an Issue with this topology

  transMat[0][0] = 1;
  for (i = 0; i < NumInsertStates; i++)
    {
    transMat[0][2*i+1] = gapOpen[2*i];
    transMat[0][2*i+2] = gapOpen[2*i+1];
    transMat[0][0] -= (gapOpen[2*i] + gapOpen[2*i+1]);

    transMat[2*i+1][2*i+1] = gapExtend[2*i];
    transMat[2*i+2][2*i+2] = gapExtend[2*i+1];
    transMat[2*i+1][2*i+2] = 0;
    transMat[2*i+2][2*i+1] = 0;
    transMat[2*i+1][0] = 1 - gapExtend[2*i];
    transMat[2*i+2][0] = 1 - gapExtend[2*i+1];
    }



  // create initial and transition probability matrices
  for (i = 0; i < NumMatrixTypes; i++){
    initialDistribution[i] = (float)log ((float)initDistribMat[i]);
    for (j = 0; j < NumMatrixTypes; j++)
      transProb[i][j] = (float)log ((float)transMat[i][j]);
  }

  // create insertion and match probability matrices
  for (i = 0; i < 256; i++)
    {
      for (j = 0; j < NumMatrixTypes; j++)
	{
	  insProb[i][j] = (float)log((float)emitSingle[i]);
	}
      for (j = 0; j < 256; j++)
	{
	  matchProb[i][j] = (float)log((float)emitPairs[i][j]);
	}
    }
  return 1;
}


int viterbi_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  int C1,c1, C2,c2;
  char *alphabet, *char_buf;
  char **al, **aln;
  int seq1Length, seq2Length, I, J;
  int i, j,ij, i1j1, i1j, ij1, k, a, b,l, LEN, r, c, m, state;
  int NumMatrixTypes=5;
  int NumInsertStates=2;
  int *traceback;
  float bestProb;
  static float **transMat, **insProb, **matchProb, *initialDistribution, **transProb, **emitPairs, *emitSingle, *TmatchProb, ***TinsProb;
  float *viterbi;

  ungap_sub_aln (A, ns[0],ls[0]);
  ungap_sub_aln (A, ns[1],ls[1]);

  seq1Length=I=strlen (A->seq_al[ls[0][0]]);
  seq2Length=J=strlen (A->seq_al[ls[1][0]]);


  if (!transMat)
    {
       alphabet=alphabetDefault;
       emitPairs=declare_float (256, 256);
       emitSingle=(float*)vcalloc (256, sizeof (float));
       for (i=0; i<256; i++)
	 {
	   emitSingle[i]=1e-5;
	   for (j=0; j<256; j++)
	     emitPairs[i][j]=1e-10;
	 }
       l=strlen (alphabet);

       for (i=0; i<l; i++)
	 {

	   c1=tolower(alphabet[i]);
	   C1=toupper(alphabet[i]);
	   emitSingle[c1]=emitSingleDefault[i];
	   emitSingle[C1]=emitSingleDefault[i];
	   for (j=0; j<=i; j++)
	     {
	       c2=tolower(alphabet[j]);
	       C2=toupper(alphabet[j]);

	       emitPairs[c1][c2]=emitPairsDefault[i][j];
	       emitPairs[C1][c2]=emitPairsDefault[i][j];
	       emitPairs[C1][C2]=emitPairsDefault[i][j];
	       emitPairs[c1][C2]=emitPairsDefault[i][j];
	       emitPairs[c2][c1]=emitPairsDefault[i][j];
	       emitPairs[C2][c1]=emitPairsDefault[i][j];
	       emitPairs[C2][C1]=emitPairsDefault[i][j];
	       emitPairs[c2][C1]=emitPairsDefault[i][j];
	     }
	 }


       transMat=declare_float (2*NumInsertStates+1, 2*NumInsertStates+1);
       transProb=declare_float (2*NumInsertStates+1,2* NumInsertStates+1);
       insProb=declare_float (256,NumMatrixTypes);
       matchProb=declare_float (256, 256);
       initialDistribution=(float*)vcalloc (2*NumMatrixTypes+1, sizeof (float));

       ProbabilisticModel (NumMatrixTypes,NumInsertStates,initDistrib2Default, emitSingle,emitPairs,gapOpen2Default,gapExtend2Default, transMat,initialDistribution,matchProb, insProb,transProb);
     }


   TmatchProb=(float*)vcalloc ((I+1)*(J+1), sizeof (float));
   TinsProb=(float***)declare_arrayN (3, sizeof (float),2,NumMatrixTypes,MAX(I,J)+1);
   get_tot_prob2 (A,A, ns,ls,NumMatrixTypes, matchProb, insProb,TmatchProb,TinsProb, CL,SEQUENCE);

   // create viterbi matrix
   l=NumMatrixTypes * (seq1Length+1) * (seq2Length+1);
   viterbi =(float*)vcalloc (l, sizeof (float));
   for (a=0; a<l; a++)viterbi[a]=LOG_ZERO;
   traceback=(int*)vcalloc (l, sizeof (int));
   for (a=0; a<l; a++)traceback[a]=-1;

   // initialization condition
   for (k = 0; k < NumMatrixTypes; k++)
     viterbi[k] = initialDistribution[k];

   // remember offset for each index combination
   ij = 0;
   i1j = -seq2Length - 1;
   ij1 = -1;
   i1j1 = -seq2Length - 2;

   ij *= NumMatrixTypes;
   i1j *= NumMatrixTypes;
   ij1 *= NumMatrixTypes;
   i1j1 *= NumMatrixTypes;

   // compute viterbi scores
   for (m=0,i = 0; i <= seq1Length; i++)
     {
       for ( j = 0; j <= seq2Length; j++, m++)
	 {
	   if (i > 0 && j > 0)
	     {
	       for (k = 0; k < NumMatrixTypes; k++)
		 {
		   float newVal = viterbi[k + i1j1] + transProb[k][0] + TmatchProb[m];
		   if (viterbi[0 + ij] < newVal)
		     {
		       viterbi[0 + ij] = newVal;
		       traceback[0 + ij] = k;
		     }
		 }
	     }
	   if (i > 0)
	     {
	       for (k = 0; k < NumInsertStates; k++)
		 {
		   float valFromMatch = TinsProb[0][k][i] + viterbi[0 + i1j] + transProb[0][2*k+1];
		   float valFromIns = TinsProb[0][k][i] + viterbi[2*k+1 + i1j] + transProb[2*k+1][2*k+1];
		   if (valFromMatch >= valFromIns){
		     viterbi[2*k+1 + ij] = valFromMatch;
		     traceback[2*k+1 + ij] = 0;
		   }
		   else {
		     viterbi[2*k+1 + ij] = valFromIns;
		     traceback[2*k+1 + ij] = 2*k+1;
		   }
		 }
	     }
	   if (j > 0)
	     {
	       for (k = 0; k < NumInsertStates; k++){
		 float valFromMatch = TinsProb[1][k][j] + viterbi[0 + ij1] + transProb[0][2*k+2];
		 float valFromIns = TinsProb[1][k][j] + viterbi[2*k+2 + ij1] + transProb[2*k+2][2*k+2];
		 if (valFromMatch >= valFromIns){
		   viterbi[2*k+2 + ij] = valFromMatch;
		   traceback[2*k+2 + ij] = 0;
		 }
		 else
		   {
		     viterbi[2*k+2 + ij] = valFromIns;
		     traceback[2*k+2 + ij] = 2*k+2;
		   }
	       }
	     }

	   ij += NumMatrixTypes;
	   i1j += NumMatrixTypes;
	   ij1 += NumMatrixTypes;
	   i1j1 += NumMatrixTypes;
	 }
     }

   // figure out best terminating cell
   bestProb = LOG_ZERO;
   state = -1;
   for (k = 0; k < NumMatrixTypes; k++)
     {
       float thisProb = viterbi[k + NumMatrixTypes * ((seq1Length+1)*(seq2Length+1) - 1)] + initialDistribution[k];
       if (bestProb < thisProb)
	 {
	   bestProb = thisProb;
	   state = k;
	 }
     }



   // compute traceback
   al=declare_char(2,seq1Length+seq2Length);
   LEN=0;
   r = seq1Length, c = seq2Length;
   while (r != 0 || c != 0)
     {
       int newState = traceback[state + NumMatrixTypes * (r * (seq2Length+1) + c)];

       if (state == 0){ c--; r--; al[0][LEN]=1;al[1][LEN]=1;}
       else if (state % 2 == 1) {r--; al[0][LEN]=1;al[1][LEN]=0;}
       else { c--; al[0][LEN]=0;al[1][LEN]=1;}
       LEN++;
       state = newState;
     }


   invert_list_char ( al[0], LEN);
   invert_list_char ( al[1], LEN);
   if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);
   aln=A->seq_al;
   char_buf=(char*) vcalloc (LEN+1, sizeof (char));
   for ( c=0; c< 2; c++)
     {
       for ( a=0; a< ns[c]; a++)
	 {
	   int ch=0;
	   for ( b=0; b< LEN; b++)
	     {
	       if (al[c][b]==1)
		 char_buf[b]=aln[ls[c][a]][ch++];
	       else
		 char_buf[b]='-';
	     }
	   char_buf[b]='\0';
	   sprintf (aln[ls[c][a]],"%s", char_buf);
	 }
     }


   A->len_aln=LEN;
   A->nseq=ns[0]+ns[1];
   vfree (char_buf);
   free_char (al, -1);





   return (int)(bestProb*(float)1000);
}

float ** get_emitPairs (char *mat, char *alp, float **p, float *s)
  {
    static char *rmat;
    float k=0, t=0;
    int a, b, c, l;
    int **M;

    if (!rmat)rmat=(char*)vcalloc (100, sizeof (char));

    if (!mat || !mat[0] || strm (mat, "default"))return p;
    else if (strm (rmat, mat))return p;

    sprintf (rmat,"%s", mat);

    M=read_matrice (mat);
    l=strlen (alp);

    k=log (2)/2;
    for (a=0; a<l; a++)
      for (b=0; b<l; b++)
	{
	  int sc;
	  float e;
	  e=s[a]*s[b];
	  sc=M[alp[a]-'A'][alp[b]-'A'];
	  p[a][b]=e*exp ((double)sc*k);
	}

    for (a=0; a<l; a++)
      for (b=0; b<l; b++)
	t+=p[a][b];

    for (a=0; a<l; a++)
      for (b=0; b<l; b++)
	p[a][b]=p[a][b]/t;

    t=0;

    for (a=0; a<l; a++)
      for (b=0; b<l; b++)
	t+=p[a][b];

    return p;
  }

float * forward_proba_pair_wise_test ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *matchProb, float ***insProb, float **transProb)
{
  static float *forward;
  static int max_l;
  int k, i, j,ij, i1j1, i1j, ij1, seq1Length, seq2Length, m;
  char *iter1, *iter2;
  int l,a;

 

  if (!seq1)
    {
      vfree (forward);
      forward=NULL; max_l=0;
      return NULL;
    }
  iter1=seq1-1;
  iter2=seq2-1;
  seq1Length=strlen (seq1);
  seq2Length=strlen (seq2);
  l=(seq1Length+1)*(seq2Length+1)*NumMatrixTypes;

  if (!forward)
    {
      forward=(float*)vcalloc (l, sizeof (float));
      max_l=l;
    }
  else if (max_l<l)
    {
      forward=(float*)vrealloc (forward, l*sizeof(float));
      max_l=l;
    }
  for (a=0; a<l; a++)forward[a]=LOG_ZERO;

  //f[1][1][0]
  forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] = initialDistribution[0] + matchProb[seq2Length+2];
  for (k = 0; k < NumInsertStates; k++)
    {
      //forward[1][0][k]
      forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] = initialDistribution[2*k+1] + insProb[0][k][1];
      //forward[1][0][k]
      forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] = initialDistribution[2*k+2] + insProb[1][k][1];
    }

  // remember offset for each index combination
  ij = 0;                //f[i][j]
  i1j = -seq2Length - 1; //f[i-1][j]
  ij1 = -1;              //f[i][j-1]
  i1j1 = -seq2Length - 2;//f[i-1][j-1]
  
  ij *= NumMatrixTypes;
  i1j *= NumMatrixTypes;
  ij1 *= NumMatrixTypes;
  i1j1 *= NumMatrixTypes;
  

  // compute forward scores
  for (m=0,i = 0; i <= seq1Length; i++)
    {
      for (j = 0; j <= seq2Length; j++, m++)
	{
	  if (i > 1 || j > 1)
	    {
	      if (i > 0 && j > 0)
		{
		  //Sum over all possible alignments
		  forward[0 + ij] = forward[0 + i1j1] + transProb[0][0];
		  for (k = 1; k < NumMatrixTypes; k++)
		    {
		      LOG_PLUS_EQUALS (&forward[0 + ij], forward[k + i1j1] + transProb[k][0]);
		    }
		  forward[0 + ij] += matchProb[m];
		}
	      if ( i > 0)
		{
		  for (k = 0; k < NumInsertStates; k++)
		    {
		      forward[2*k+1 + ij] = insProb[0][k][i] + LOG_ADD (forward[0 + i1j] + transProb[0][2*k+1],forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1]);
		    }
		}
	      if (j > 0)
		{
		  for (k = 0; k < NumInsertStates; k++)
		    {
		      forward[2*k+2 + ij] = insProb[1][k][j] +LOG_ADD (forward[0 + ij1] + transProb[0][2*k+2],forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2]);
		    }
		}
	    }
	  
	  ij += NumMatrixTypes;
	  i1j += NumMatrixTypes;
	  ij1 += NumMatrixTypes;
	  i1j1 += NumMatrixTypes;
	}
      
    }
  return forward;
}
//old




float * backward_proba_pair_wise_test ( char *seq1, char *seq2, int NumMatrixTypes, int NumInsertStates, float **transMat, float *initialDistribution,float *imatchProb, float ***insProb, float **transProb)
{
  static float *backward;
  static int max_l;
  float ***bw;
  float **matchProb;
  int l1, l2, ns;
  
  int k, i, j,ij, i1j1, i1j, ij1,a, l, seq1Length, seq2Length, m;
  char c1, c2;
  char *iter1, *iter2;

  if (!seq1)
    {
      vfree (backward);
      backward=NULL; max_l=0;
      return NULL;
    }

  iter1=seq1-1;
  iter2=seq2-1;
  l1=seq1Length=strlen (seq1);
  l2=seq2Length=strlen (seq2);
  l=(seq1Length+1)*(seq2Length+1)*NumMatrixTypes;
  ns=NumMatrixTypes;
  
 
  if (!backward)
    {
      backward=(float*)vcalloc (l, sizeof (float));
      max_l=l;
    }
  else if (max_l<l)
    {
      backward=(float*)vrealloc (backward, l*sizeof(float));
      max_l=l;
    }
  for (a=0; a<l; a++)backward[a]=LOG_ZERO;
  
  bw=(float***)vcalloc (ns, sizeof (float**));
  for (k=0; k<ns; k++)
    {
      bw[k]=declare_float (l1+2, l2+2);
      for (i=0; i<(l1+2); i++)
	for (j=0; j<(l2+2); j++)
	  bw[k][i][j]=LOG_ZERO;
    }
  matchProb=declare_float (l1+1, l2+1);
  for (m=0,i=0; i<=l1; i++)
    for (j=0; j<=l2; j++)
      matchProb[i][j]=imatchProb[m++];
  
  
  bw[0][l1][l2]=initialDistribution[0]+matchProb[l1][l2];
  for (k = 0; k < NumInsertStates; k++)
    {
      bw[2*k+1][l1][l2+1]= initialDistribution[2*k+1] + insProb[0][k][1];
      bw[2*k+2][l1+1][l2]= initialDistribution[2*k+2] + insProb[1][k][1];
    }
  for (i=l1+1; i > 0; i--)
    {
      for (j =l2+1; j > 0; j--)
	{
	  if (i<l1 || j<l2 )
	    {
	      if (i <(l1+1) && j < (l2+1))
		{
		  bw[0][i][j]=bw[0][i+1][j+1]+transProb[0][0];
		  for (k = 1; k < NumMatrixTypes; k++)
		    {
		      LOG_PLUS_EQUALS (&bw[0][i][j], bw[k][i+1][j+1] + transProb[k][0]);
		    }
		  bw[0][i][j] += matchProb[i][j];
		}
	      if (i <l1+1)
		{
		  for (k = 0; k < NumInsertStates; k++)
		    {
		      bw[2*k+1][i][j] = insProb[0][k][i] + LOG_ADD (bw[0][i+1][j] + transProb[0][2*k+1],bw[2*k+1][i+1][j] + transProb[2*k+1][2*k+1]);
		    }
		}
	      
	      if (j <l2+1)
		{
		  for (k = 0; k < NumInsertStates; k++)
		    {
		      bw[2*k+2][i][j] = insProb[1][k][j] +LOG_ADD (bw[0][i][j+1] + transProb[0][2*k+2],bw[2*k+2][i][j+1] + transProb[2*k+2][2*k+2]);
		    }
		}
	    }
	}
    }
 
  //resize backward to make it compatible with forward
  for (m=0,i=0; i<=l1; i++)
    for (j=0; j<=l2; j++)
      {
	for (k=0; k<ns; k++)
	  {
	    backward[m++]=bw[k][i][j];
	  }
      }
  
  return backward;
}


Constraint_list *ProbaMatrix2CL_test (Alignment *A, int *ns, int **ls, int NumMatrixTypes, int NumInsertStates, float *forward, float *backward, float thr, Constraint_list *CL)
{
  float totalProb;
  int ij, i, j,k, I, J, s1, s2;
  static int *entry;
  static int **list;
  static int list_max;
  int sim;
  int list_size;
  int list_n;
  int old_n=0;
  double v;
  static float F=4; //potential number of full suboptimal alignmnents incorporated in the library
  static int tot_old, tot_new;
  float ***fw, tfw;
  float ***bw, tbw;

  
  
  
  if (!A)
    {
      free_int (list, -1);list=NULL;
      list_max=0;

      vfree(entry); entry=NULL;
      return NULL;
    }

  
  
  I=strlen (A->seq_al[ls[0][0]]);
  J=strlen (A->seq_al[ls[1][0]]);
  s1=name_is_in_list (A->name[ls[0][0]], (CL->S)->name, (CL->S)->nseq, 100);
  s2=name_is_in_list (A->name[ls[1][0]], (CL->S)->name, (CL->S)->nseq, 100);

  list_size=I*J;
  
  fw=(float***)vcalloc (NumMatrixTypes, sizeof (float**));
  bw=(float***)vcalloc (NumMatrixTypes, sizeof (float**));
  
  for (k=0; k<NumMatrixTypes; k++)
    {
      fw[k]=declare_float (I+1, J+1);
      bw[k]=declare_float (I+1, J+1);
    }
  
  
  if ( list_max<list_size)
     {
      free_int (list, -1);
      list_max=list_size;
      list=declare_int (list_max, 3);
    }


  totalProb = ComputeTotalProbability_test (I,J,NumMatrixTypes, NumInsertStates,forward, backward);

  ij = 0;
  for (list_n=0,ij=0,i =0; i <= I; i++)
    {
      for (j =0; j <= J; j++, ij+=NumMatrixTypes)
	{
	  for (k=0; k<NumMatrixTypes; k++)
	    {
	      fw[k][i][j]=forward [k+ij];
	      bw[k][i][j]=backward[k+ij];
	    }
	 
	  v= EXP (MIN(LOG_ONE,(forward[ij] + backward[ij] - totalProb)));
	  if (v>thr)//Conservative reduction of the list size to speed up the sorting
	    {
	      list[list_n][0]=i;
	      list[list_n][1]=j;
	      list[list_n][2]=(int)((float)v*(float)NORM_F);
	      list_n++;
	    }
	  if (v>0.01)old_n++;
	}
    }

  if (2==1)
    {
      for (i=0; i<=I; i++)
	for (j=0; j<=J; j++)
	  {
	    int ns=NumMatrixTypes;
	    for (k=0; k<ns; k++)
	      {
		fprintf ( stderr, "\ni=%2d j=%2d K=%d F=%f B=%f", i, j, k, fw[k][i][j], bw[k][i][j]);
	      }
	  }
    }
  tbw=tfw=LOG_ZERO;
  for (k=0; k<NumMatrixTypes; k++)
    {
      //LOG_PLUS_EQUALS(&tbw, bw[k][1][1]);
      //LOG_PLUS_EQUALS(&tbw, fw[k][1][1]);
      
      LOG_PLUS_EQUALS(&tfw, fw[k][I][J]);
      LOG_PLUS_EQUALS(&tfw, bw[k][I][J]);
    }
      
    
  sort_int_inv (list, 3, 2, 0, list_n-1);
  if (!entry)entry=(int*)vcalloc ( CL->entry_len+1, CL->el_size);

  list_n=MIN(list_n,(F*MIN(I,J)));
  for (i=0; i<list_n; i++)
    {
       entry[SEQ1]=s1;
       entry[SEQ2]=s2;
       entry[R1]  =list[i][0];
       entry[R2]  =list[i][1];
       entry[WE]  =list[i][2];
       entry[CONS]=1;
       add_entry2list (entry,A->CL);
    }
  tot_new+=list_n;
  tot_old+=old_n;
  // HERE ("LIB_SIZE NEW: %d (new) %d (old) [%.2f]", list_n, old_n, (float)tot_new/(float)tot_old);
  return A->CL;
}

float ComputeTotalProbability_test (int seq1Length, int seq2Length,int NumMatrixTypes, int NumInsertStates,float *forward, float *backward)
{

    float totalForwardProb = LOG_ZERO;
    float totalBackwardProb= LOG_ZERO;
    int k;
   
    for (k = 0; k < NumMatrixTypes; k++)
      {
      LOG_PLUS_EQUALS (&totalForwardProb,forward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + backward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
      }

    totalBackwardProb =forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)];

    for (k = 0; k < NumInsertStates; k++)
      {
      LOG_PLUS_EQUALS (&totalBackwardProb,forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)]);
      LOG_PLUS_EQUALS (&totalBackwardProb,forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)]);
      }
    return (totalForwardProb + totalBackwardProb) / 2;
  }


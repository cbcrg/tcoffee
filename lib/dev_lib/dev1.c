#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
#include "dev1.h"

//Insert functions here
void aln2hitMat_help()
{
	fprintf ( stdout, "\n+aln2hitMat| _MODE_    : how to compare the two positions of the alignment (default:id)");
	fprintf ( stdout, "\n.................id    : the sequence identity of those two positions");
	fprintf ( stdout, "\n.............pairscore : the pairwise score of the residues of those two positions");
	fprintf ( stdout, "\n+aln2hitMat| _MATRIX_  : matrix used for the comparison (idmat, blosum62mt, pam250mt.. default:blosum62mt)\n");
	exit (EXIT_SUCCESS);
}

void aln2hitMat (Alignment *A, char *phitmat)
{
	float **ffpHitScoreMatrix;
	int i, j, k, l, s;
	int nl = A->len_aln;
	int inseq = A->nseq;
	int itmpScore;
  	char matrix[100];
	char mode[100];
	int isim_count, itotal_count, r1, r2;

//Initialization for files	
	char *pcFileName = A->file[0];
	char prefix[200] ={0};
	char *hit_matrix_file = vcalloc(200, sizeof (char));
	char *hit_html_file = vcalloc(200, sizeof (char));
	int len = (strrchr(pcFileName,'.')?strrchr(pcFileName,'.')-pcFileName:strlen(pcFileName));

	strncpy(prefix, pcFileName, len);	
	sprintf(hit_matrix_file, "%s%s", prefix, "_aln.hit_matrix");
  	sprintf(hit_html_file, "%s%s", prefix, ".alnhit_html");

	if ( phitmat && strstr ( phitmat, "help"))
		aln2hitMat_help();

	if(phitmat == NULL) phitmat = vcalloc(1, sizeof(char)); //such that program could get default value

  	strget_param (phitmat, "_MODE_", "id", "%s", mode);
	strget_param (phitmat, "_MATRIX_", "blosum62mt", "%s", matrix);

	fprintf ( stdout, "[START] aln to hit matrix\n");
	fprintf ( stdout, "	   Mode:%s\n", mode);
	fprintf ( stdout, "	   Matrix:%s\n", matrix);

	int **mat = read_matrice(matrix);

	ffpHitScoreMatrix=vcalloc (nl, sizeof (float*));
	for(i = 0; i < nl; i++)
      		ffpHitScoreMatrix[i]=vcalloc (nl-i, sizeof (float));

	fprintf (stdout, "Process positions\n", i);
	for(i = 0; i < nl; i++)
	{
		fprintf (stdout, "%d, ", i);
		for(j = i; j < nl; j++)
		{
			if(strm (mode, "id"))
				ffpHitScoreMatrix[i][j-i]=generic_get_seq_sim (aln_column2string(A, i), aln_column2string(A, j), (A->cdna_cache)?A->cdna_cache[0]:NULL, matrix);
			else if(strm (mode, "pairscore"))
			{
				isim_count = itotal_count = 0;
				for (k=0; k< inseq; k++)
				{
					r1=tolower(A->seq_al[k][i]);
					if (is_gap(r1))continue;
					for (l=0; l< inseq; l++)
					{
						r2=tolower(A->seq_al[l][j]);
						if (is_gap (r2))continue;						
						s=mat[r2-'A'][r1-'A'];
						s=(s<=0)?0:1;
						isim_count += s;
						itotal_count++;
					}
				}
				r1=(isim_count*100)/itotal_count;
				ffpHitScoreMatrix[i][j-i] = r1;
			}
			else
				aln2hitMat_help();
		}
	}
	fprintf (stdout, "\n");
	output_hit_matrix(hit_matrix_file, ffpHitScoreMatrix, nl);

//Output Hit Score into color html
	output_hit_color_html  (A, ffpHitScoreMatrix, nl, hit_html_file);
 	vfree(ffpHitScoreMatrix);
	vfree(hit_matrix_file);
	vfree(hit_html_file);
	fprintf ( stdout, "[END] aln to hit matrix\n");
}

/*
 * 5_state_hmm.cpp
 *
 *  Created on: Apr 12, 2012
 *      Author: Carsten Kemena
 *
 * This file is part of BioTools++.
 *
 * BioTools++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "HMM.h"

namespace BioTools {
namespace HMM {


using namespace std;
using namespace BioTools::Utils;
using namespace BioTools::Seq;

//Protein Alignment Models
static float prot_initDistrib2Default[] = { 0.6814756989f, 8.615339902e-05f, 8.615339902e-05f, 0.1591759622f, 0.1591759622 };
static float prot_gapOpen2Default[] = { 0.0119511066f, 0.0119511066f, 0.008008334786f, 0.008008334786 };
static float prot_gapExtend2Default[] = { 0.3965826333f, 0.3965826333f, 0.8988758326f, 0.8988758326 };

static char prot_alphabetDefault[] = "ARNDCQEGHILKMFPSTWYV";
static float prot_emitSingleDefault[20] = {
  0.07831005f, 0.05246024f, 0.04433257f, 0.05130349f, 0.02189704f,
  0.03585766f, 0.05615771f, 0.07783433f, 0.02601093f, 0.06511648f,
  0.09716489f, 0.05877077f, 0.02438117f, 0.04463228f, 0.03940142f,
  0.05849916f, 0.05115306f, 0.01203523f, 0.03124726f, 0.07343426f
};

static float prot_emitPairsDefault[20][20] = {
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




HMM::HMM(char type)
{
	float *initDistribDef, *gapOpen, *gapExtend, *emitSingle, (*emitPair)[20];
	char *alphabet;
	short alpha_len;
	_num_states = 5;
	_num_ins_states = 2;
	if (type=='P')
	{
		initDistribDef= &prot_initDistrib2Default[0];
		gapOpen = &prot_gapOpen2Default[0];
		gapExtend= &prot_gapExtend2Default[0];
		emitSingle =&prot_emitSingleDefault[0];
		emitPair = &prot_emitPairsDefault[0];
		alphabet = &prot_alphabetDefault[0];
		alpha_len =20;//sizeof(prot_alphabetDefault);
	}

	_initDistr= new float[_num_states];
	_transProb.resize(2*_num_ins_states+1, 2*_num_ins_states+1);
	_transProb[0][0]=1;
	for (short i = 0; i < _num_ins_states; i++)
	{
		_transProb[0][2*i+1] = log(gapOpen[2*i]);
		_transProb[0][2*i+2] = log(gapOpen[2*i+1]);
		_transProb[0][0] -= (gapOpen[2*i] + gapOpen[2*i+1]);

		_transProb[2*i+1][2*i+1] = log(gapExtend[2*i]);
		_transProb[2*i+2][2*i+2] = log(gapExtend[2*i+1]);
		_transProb[2*i+1][2*i+2] = FLT_MIN;
		_transProb[2*i+2][2*i+1] = FLT_MIN;
		_transProb[2*i+1][0] = log(1 - gapExtend[2*i]);
		_transProb[2*i+2][0] = log(1 - gapExtend[2*i+1]);
	}
	_transProb[0][0]=log(_transProb[0][0]);

	short i,j;
	for (i = 0; i < _num_states; i++)
		_initDistr[i] = log(initDistribDef[i]);

	_matchProb.resize(256,256);
	_insProb.resize(256,_num_states);
	int c1,C1,c2,C2;

	for (i=0; i<256; ++i)
	{
		for (j=0; j<_num_states; ++j)
			_insProb[i][j] = log(1e-5);
		for (j=0; j<256; j++)
			_matchProb[i][j]=log(1e-10);
	}

	for (i=0; i<alpha_len; ++i)
	{

		c1=static_cast<int>(tolower(alphabet[i]));
		C1=static_cast<int>(toupper(alphabet[i]));
		for (j=0; j<_num_states; ++j)
		{
			_insProb[c1][j]=log(emitSingle[i]);
			_insProb[C1][j]=log(emitSingle[i]);
		}
		for (j=0; j<=i; j++)
		{
			c2=static_cast<int>(tolower(alphabet[j]));
			C2=static_cast<int>(toupper(alphabet[j]));
			_matchProb[c1][c2] = _matchProb[C1][c2] = _matchProb[C1][C2] = _matchProb[c1][C2] = _matchProb[c2][c1] = _matchProb[C2][c1] = _matchProb[C2][C1] = _matchProb[c2][C1] = log(emitPair[i][j]);
		}
	}
}


HMM::~HMM()
{

}


void
HMM::aln_probs(const BioTools::Seq::Alignment &aln1, const BioTools::Seq::Alignment &aln2, BioTools::Utils::Matrix<double> &ins_probs1, BioTools::Utils::Matrix<double> &ins_probs2, BioTools::Utils::Matrix<double> &match_probs)
{

	size_t aln_len1 = aln1.size();
	size_t aln_len2 = aln2.size();
	size_t i, j, k, l;

	//calculate new insertion scores
	size_t n_ins_states=num_ins_states();
	ins_probs1.resize(aln_len1, n_ins_states);
	ins_probs2.resize(aln_len2, n_ins_states);

	size_t n_seq1 = aln1.n_seqs();
	char c;
	size_t n_no_gaps;
	for (i=0; i<aln_len1; ++i)
	{
		n_no_gaps=0;
		for (j=0; j<n_seq1; ++j)
		{
			c=aln1[j][i];
			if ((c!='-') || (c!='.'))
			{	++n_no_gaps;
				for (k=0; k<n_ins_states; ++k)
					ins_probs1[i][k] += _insProb[c][k];
			}
		}
		for (k=0; k<n_ins_states; ++k)
			ins_probs1[i][k] /= n_no_gaps;
	}

	size_t n_seq2 = aln2.n_seqs();
	for (i=0; i<aln_len2; ++i)
	{
		n_no_gaps=0;
		for (j=0; j<n_seq2; ++j)
		{
			c=aln2[j][i];
			if ((c!='-') || (c!='.'))
			{	++n_no_gaps;
				for (k=0; k<n_ins_states; ++k)
					ins_probs2[i][k] += _insProb[c][k];
			}
		}
		for (k=0; k<n_ins_states; ++k)
			ins_probs2[i][k] /= n_no_gaps;
	}

	// calculate match_scores
	Matrix<double> prof1(aln_len1, 26, 0);
	Matrix<double> prof2(aln_len2, 26, 0);

	for (i=0; i<n_seq1; ++i)
	{
		const Sequence &seq = aln1[i];
		for (j=0; j<aln_len1; ++j)
		{
			c=seq[j];
			if ((c!='-') && (c!='.'))
				++prof1[j][tolower(c)-'a'];
		}
	}

	for (i=0; i<n_seq2; ++i)
	{
		const Sequence &seq = aln2[i];
		for (j=0; j<aln_len2; ++j)
		{
			c=seq[j];
			if ((c!='-') && (c!='.'))
				++prof2[j][tolower(c)-'a'];
		}
	}

	match_probs.resize(aln_len1, aln_len2);
	double tmp;
	for (i=0; i<aln_len1; ++i)
	{
		for (j=0; j<aln_len2; ++j)
		{
			tmp=0;
			for (k=0; k<26; ++k)
			{
				for (l=0; l<26; ++l)
				{
			//		cout << i << " "<< j << " "<< k << " "<< l << endl;
					match_probs[i][j] += _matchProb[k+65][l+65] * prof1[i][k] * prof2[j][l];
					tmp += prof1[i][k] * prof2[j][l];
				}
			}
			match_probs[i][j] /= tmp;
		}
	}
}



}
}



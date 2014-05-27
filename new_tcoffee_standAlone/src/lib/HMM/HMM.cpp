/*
 * 5_state_hmm.cpp
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

static char prot_alphabetDefault[] = "ARNDCQEGHILKMFPSTWYVBZX";
static float prot_emitSingleDefault[23] = {
  0.07831005f, 0.05246024f, 0.04433257f, 0.05130349f, 0.02189704f,
  0.03585766f, 0.05615771f, 0.07783433f, 0.02601093f, 0.06511648f,
  0.09716489f, 0.05877077f, 0.02438117f, 0.04463228f, 0.03940142f,
  0.05849916f, 0.05115306f, 0.01203523f, 0.03124726f, 0.07343426f, 0.04781803f, 0.046007685f, 0.0001
};

static float prot_emitPairsDefault[23][23] = {
  {0.02373072f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00244502f, 0.01775118f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00210228f, 0.00207782f, 0.01281864f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00223549f, 0.00161657f, 0.00353540f, 0.01911178f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00145515f, 0.00044701f, 0.00042479f, 0.00036798f, 0.01013470f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00219102f, 0.00253532f, 0.00158223f, 0.00176784f, 0.00032102f, 0.00756604f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00332218f, 0.00268865f, 0.00224738f, 0.00496800f, 0.00037956f, 0.00345128f, 0.01676565f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00597898f, 0.00194865f, 0.00288882f, 0.00235249f, 0.00071206f, 0.00142432f, 0.00214860f, 0.04062876f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00114353f, 0.00132105f, 0.00141205f, 0.00097077f, 0.00026421f, 0.00113901f, 0.00131767f, 0.00103704f, 0.00867996f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00318853f, 0.00138145f, 0.00104273f, 0.00105355f, 0.00094040f, 0.00100883f, 0.00124207f, 0.00142520f, 0.00059716f, 0.01778263f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00449576f, 0.00246811f, 0.00160275f, 0.00161966f, 0.00138494f, 0.00180553f, 0.00222063f, 0.00212853f, 0.00111754f, 0.01071834f, 0.03583921f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00331693f, 0.00595650f, 0.00257310f, 0.00252518f, 0.00046951f, 0.00312308f, 0.00428420f, 0.00259311f, 0.00121376f, 0.00157852f, 0.00259626f, 0.01612228f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00148878f, 0.00076734f, 0.00063401f, 0.00047808f, 0.00037421f, 0.00075546f, 0.00076105f, 0.00066504f, 0.00042237f, 0.00224097f, 0.00461939f, 0.00096120f, 0.00409522f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00165004f, 0.00090768f, 0.00084658f, 0.00069041f, 0.00052274f, 0.00059248f, 0.00078814f, 0.00115204f, 0.00072545f, 0.00279948f, 0.00533369f, 0.00087222f, 0.00116111f, 0.01661038f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00230618f, 0.00106268f, 0.00100282f, 0.00125381f, 0.00034766f, 0.00090111f, 0.00151550f, 0.00155601f, 0.00049078f, 0.00103767f, 0.00157310f, 0.00154836f, 0.00046718f, 0.00060701f, 0.01846071f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00631752f, 0.00224540f, 0.00301397f, 0.00285226f, 0.00094867f, 0.00191155f, 0.00293898f, 0.00381962f, 0.00116422f, 0.00173565f, 0.00250962f, 0.00312633f, 0.00087787f, 0.00119036f, 0.00180037f, 0.01346609f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00389995f, 0.00186053f, 0.00220144f, 0.00180488f, 0.00073798f, 0.00154526f, 0.00216760f, 0.00214841f, 0.00077747f, 0.00248968f, 0.00302273f, 0.00250862f, 0.00093371f, 0.00107595f, 0.00147982f, 0.00487295f, 0.01299436f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00039119f, 0.00029139f, 0.00021006f, 0.00016015f, 0.00010666f, 0.00020592f, 0.00023815f, 0.00038786f, 0.00019097f, 0.00039549f, 0.00076736f, 0.00028448f, 0.00016253f, 0.00085751f, 0.00015674f, 0.00026525f, 0.00024961f, 0.00563625f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00131840f, 0.00099430f, 0.00074960f, 0.00066005f, 0.00036626f, 0.00070192f, 0.00092548f, 0.00089301f, 0.00131038f, 0.00127857f, 0.00219713f, 0.00100817f, 0.00054105f, 0.00368739f, 0.00047608f, 0.00102648f, 0.00094759f, 0.00069226f, 0.00999315f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00533241f, 0.00169359f, 0.00136609f, 0.00127915f, 0.00119152f, 0.00132844f, 0.00178697f, 0.00194579f, 0.00071553f, 0.01117956f, 0.00914460f, 0.00210897f, 0.00197461f, 0.00256159f, 0.00135781f, 0.00241601f, 0.00343452f, 0.00038538f, 0.00148001f, 0.02075171f, 0.0f, 0.0f, 0.0f},
  {0.00216889f,	0.00184720f, 0.00817702f, 0.01132359f, 0.00039639f, 0.00167504f, 0.00360769f, 0.00262066f, 0.00119141f, 0.00104814f, 0.00161121f, 0.00254914f, 0.00055605f, 0.00076850f, 0.00112832f, 0.00293312f, 0.00200316f,	0.00018511f, 0.00070483f, 0.00132262f, 0.01596521f, 0.0f, 0.0f},
  {0.00275660f,	0.00261199f, 0.00191481f, 0.00336792f, 0.00035029f, 0.00550866f, 0.01010847f, 0.00178646f, 0.00122834f, 0.00112545f, 0.00201308f, 0.00370364f, 0.00075826f, 0.00069031f, 0.00120831f, 0.00242526f, 0.00185643f, 0.00022204f, 0.00081370f, 0.00155771f, 0.00264137f, 0.01216585f, 0.0f},
  {1.0,	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
};

//DNA/RNA Alignment Models
static float nuc_initDistrib2Default[] = { 0.9615409374f, 0.0000004538f, 0.0000004538f, 0.0192291681f, 0.0192291681f };
static float nuc_gapOpen2Default[] = { 0.0082473317f, 0.0082473317f, 0.0107844425f, 0.0107844425f };
static float nuc_gapExtend2Default[] = { 0.3210460842f, 0.3210460842f, 0.3298229277f, 0.3298229277f };

static char nuc_alphabetDefault[] = "ACGUTN";
static float nuc_emitSingleDefault[6] = { 
  0.2270790040f, 0.2422080040f, 0.2839320004f, 0.2464679927f, 0.2464679927f, 0.0003124650f 
};

static float nuc_emitPairsDefault[6][6] = {
  { 0.1487240046f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f },
  { 0.0184142999f, 0.1583919972f, 0.0f, 0.0f, 0.0f, 0.0f },
  { 0.0361397006f, 0.0275536999f, 0.1979320049f, 0.0f, 0.0f, 0.0f },
  { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f, 0.0f, 0.0f },
  { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f, 0.1557479948f, 0.0f },
  { 0.0000375308f, 0.0000815823f, 0.0000824765f, 0.0000743985f, 0.0000743985f, 0.0000263252f }
};


HMM::HMM(char type)
{
	float *initDistribDef, *gapOpen, *gapExtend, *emitSingle, **emitPair=NULL;
	char *alphabet;
	short alpha_len;
	//typedef float prot_emit[22][22];
	
	_num_states = 5;
	_num_ins_states = 2;
	if (type=='P')
	{
		emitPair= new float*[23];
		for(int i = 0; i < 23; ++i)
		      emitPair[i] = new float[23];
		initDistribDef= &prot_initDistrib2Default[0];
		gapOpen = &prot_gapOpen2Default[0];
		gapExtend= &prot_gapExtend2Default[0];
		emitSingle =&prot_emitSingleDefault[0];
		for(int i=0; i<23; i++)
		  for(int j=0; j<23; j++)
		    emitPair[i][j]=prot_emitPairsDefault[i][j];
		alphabet = &prot_alphabetDefault[0];
		alpha_len =23;//sizeof(prot_alphabetDefault);
	}
	if (type=='N')
	{ 
		emitPair= new float*[6];
		for(int i = 0; i < 6; ++i)
		      emitPair[i] = new float[6];
		initDistribDef= &nuc_initDistrib2Default[0];
		gapOpen = &nuc_gapOpen2Default[0];
		gapExtend= &nuc_gapExtend2Default[0];
		emitSingle =&nuc_emitSingleDefault[0];
		for(int i=0; i<6; i++)
		  for(int j=0; j<6; j++)
		    emitPair[i][j]=nuc_emitPairsDefault[i][j];
		alphabet = &nuc_alphabetDefault[0];
		alpha_len =6;//sizeof(nuc_alphabetDefault);
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
HMM::calculate_insertion_probs(const BioTools::Seq::Alignment &aln, BioTools::Utils::Matrix<float> &ins_probs)
{
	size_t aln_len = aln.size();
	size_t i, j;

	//calculate new insertion scores
	ins_probs.resize(aln_len,1);
	for (j=0; j<aln_len; ++j)
		ins_probs[j][0]=0;
	std::vector<int> non_gap_counter(aln_len, 0);
	int totalSeq=0;
	size_t n_seqs = aln.n_seqs();
	for (i=0; i<n_seqs; ++i)
		totalSeq += aln[i].ocurrence();

	char c;
	for (i=0; i<n_seqs; ++i)
	{
		for (j=0; j<aln_len; ++j)
		{
			c=tolower(aln[i][j]);
			if ((c!='-') && (c!='.'))
			{
				non_gap_counter[j]+=aln[i].ocurrence();
				ins_probs[j][0] += _insProb[c][0]*aln[i].ocurrence();
			}
		}

	}
	for (j=0; j<aln_len; ++j)
	{
		ins_probs[j][0] /= non_gap_counter[j];
		//-->Normalization that doesn't seem to improve always the results
		/*if (non_gap_counter[j] != totalSeq)
			ins_probs[j][0] += log(1-(1.0*(non_gap_counter[j])/totalSeq));*/
	}
}



void
HMM::calculate_match_probs(const BioTools::Seq::Alignment &aln1, const BioTools::Seq::Alignment &aln2, Matrix<float> &match_probs)
{

	size_t aln_len1=aln1.size();
	size_t aln_len2=aln2.size();
//	Matrix<float> prof1(aln_len1, 26, 0);
//	Matrix<float> prof2(aln_len2, 26, 0);
	std::vector<std::unordered_map<short, int> > prof1(aln_len1);
	size_t i,j;
	std::vector<std::unordered_map<short, int> > prof2(aln_len2);
	std::unordered_map<short, int>::iterator it;

	size_t n_seq1=aln1.n_seqs();
	size_t n_seq2=aln2.n_seqs();



	vector<size_t> totaNonGap1(aln_len1);
	vector<size_t> totaNonGap2(aln_len2);
	size_t total1 = 0, total2 = 0;
	for (i=0; i<n_seq1; ++i)
		total1+=aln1[i].ocurrence();
	for (i=0; i<n_seq2; ++i)
		total2+=aln2[i].ocurrence();


	short c;
	std::vector<int> observed(26,0);
	size_t occ;
	for (i=0; i<n_seq1; ++i)
	{
		const Sequence &seq = aln1[i];
		occ = seq.ocurrence();
		for (j=0; j<aln_len1; ++j)
		{
			if ((seq[j] != '-') && (seq[j]!='.'))
			{
				totaNonGap1[j] += occ;
				c=tolower(seq[j]);
				if ((it =prof1[j].find(c)) != prof1[j].end())
					it->second+=occ;
				else
					prof1[j][c]=occ;
			}
		}
	}

	for (i=0; i<n_seq2; ++i)
	{
		const Sequence &seq = aln2[i];
		occ = seq.ocurrence();
		for (j=0; j<aln_len2; ++j)
			if ((seq[j] != '-') && (seq[j] != '.'))
			{
				totaNonGap2[j] += occ;
				c=tolower(seq[j]);
				if ((it =prof2[j].find(c)) != prof2[j].end())
					++it->second+=occ;
				else
					prof2[j][c]=occ;
			}
	}

	match_probs.resize(aln_len1, aln_len2);
	for (i=0; i<aln_len1; ++i)
	{
		for (j=0; j<aln_len2; ++j)
			match_probs[i][j] = 0;
	}

	std::unordered_map<short, int>::iterator it1,it2,it1_end,it2_end;
	double tmp, eqW;
	for (i=0; i<aln_len1; ++i)
	{
		it1_end=prof1[i].end();
		for (j=0; j<aln_len2; ++j)
		{
			it2_end=prof2[j].end();
			tmp=0;
			for (it1=prof1[i].begin(); it1!=it1_end; ++it1)
			{
				for (it2=prof2[j].begin(); it2!=it2_end; ++it2)
				{
					if(it1->first == it2->first){  eqW=it1->second * it2->second; }
					else{ eqW=1; }
					match_probs[i][j] += _matchProb[it1->first][it2->first] * it1->second * it2->second;
					tmp += it1->second * it2->second * eqW;
				}
			}
			match_probs[i][j] /= tmp;
			//match_probs[i][j] += log((1.0*totaNonGap1[i]/total1) * (1.0*totaNonGap2[j]/total2));  //-->Normalization that doesn't seem to improve always the results
		}
	}



}



}
}


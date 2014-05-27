

#include "aln_analysis.h"


namespace BioTools
{
namespace Seq
{

using namespace std;
using namespace BioTools::Utils;

// Analysis methods


/*
// old version
pair<size_t, size_t>
identity(const Alignment &aln)
{
	size_t n_seqs = aln.n_seqs();
	size_t aln_length = aln.size();
	size_t identical_pairs = 0;
	size_t total_pairs = 0;
	char c1,c2;
	size_t i,j,k;
	int chunk = 30;
	#pragma omp parallel shared(n_seqs) private(i,j,k,c1,c2) reduction(+:identical_pairs) reduction(+:total_pairs)
	{
	#pragma omp for schedule(dynamic,chunk) nowait

	for (i=0; i<n_seqs; ++i)
	{
		const Sequence &seq_1 = aln[i];

		for (j=i+1; j<n_seqs; ++j)
		{
			const Sequence &seq_2 = aln[j];
			for (k=0; k<aln_length; ++k)
			{
				if (((c1=tolower(seq_1[k])) != '-') && ((c2=tolower(seq_2[k])) != '-'))
				{
					++total_pairs;
					if (c1==c2)
						++identical_pairs;
				}
			}
		}
	}
	} //END OpenMP
	return pair<size_t, size_t>(identical_pairs, total_pairs);
}*/

pair<size_t, size_t>
msaIdentity(const Alignment &aln)
{
	size_t n_seqs = aln.n_seqs();
	size_t aln_length = aln.size();
	size_t identical_pairs = 0;
	size_t total_pairs = 0;
	char c1;
	size_t i,j,k;
	int chunk = 300;
	Matrix<size_t> counter(aln_length, 27, 0);
	#pragma omp parallel shared(n_seqs) private(i,j,k,c1) reduction(+:identical_pairs) reduction(+:total_pairs)
	{
		Matrix<size_t> tmp_counter(aln_length, 27, 0);
		#pragma omp for schedule(dynamic,chunk) nowait
		for (i=0; i<n_seqs; ++i)
		{
			const Sequence &seq = aln[i];
			for (k=0; k<aln_length; ++k)
			{
				if ((c1=seq[k]) != '-')
					++tmp_counter[k][tolower(c1)-'a'];
			}
			#pragma omp critical(dataupdate)
			{
				for (k=0; k<aln_length; ++k)
				{
					for (j=0; j<26; ++j)
					{
						counter[k][j] += tmp_counter[k][j];
						tmp_counter[k][j] = 0;
					}
				}
			}

		}
	} //END OpenMP

	size_t tmp_total_in_column;
	for (k=0; k<aln_length; ++k)
	{
		tmp_total_in_column=0;
		for (j=0; j<26; ++j)
		{
			identical_pairs += (counter[k][j] * (counter[k][j]-1)/2);
			tmp_total_in_column += counter[k][j];
		}
		total_pairs += (tmp_total_in_column * (tmp_total_in_column-1)/2);
	}
	return pair<size_t, size_t>(identical_pairs, total_pairs);
}


double
sum_of_pairs_score(const Alignment &aln, const Scoring_Matrix &matrix, double gop, double gep)
{
	size_t n_seqs = aln.n_seqs();
	size_t aln_length = aln.size();
	double score = 0;

	size_t i,j,k;
	char c1,c2;
	short gap_continue=0;
	for (i=0; i < n_seqs; ++i)
	{
		const Sequence &seq1 = aln[i];
		for (j=i+1; j < n_seqs; ++j)
		{
			const Sequence &seq2 = aln[j];
			for (k=0; k < aln_length; ++k)
			{
				c1 = seq1[k];
				c2 = seq2[k];
				if ((c1 != '-') && (c2 != '-'))
				{
					score += matrix[tolower(c1)-97][tolower(c2)-97];
					gap_continue=0;
				}
				else if ((c1 != '-') || (c2 != '-'))
				{
					if (((c1 == '-') && (gap_continue == 1)) || ((c2 == '-') && (gap_continue == 2)))
						score += gep;
					else
						gap_continue = 0;
					
					if (!gap_continue)
					{
						score += gep + gop;
						gap_continue = (c1=='-')?1:2;
					}
				}

			}

		}
	}

	return score;
}


/*
 * \brief Defines which columns to used depending on the gap_threshold
 * \param aln The alignment.
 * \param The gap threshold in percentage.
 * \return An array of the same length as aln. It contains 0 for columns having a number of gaps larger than the thresold, else 1.
 */
size_t *
use_columns(const Alignment &aln, double gap_threshold)
{
	gap_threshold/=100;
	size_t n_seqs = aln.n_seqs();
	size_t aln_length = aln.size();

	size_t i,j;
	size_t *n_gaps = new size_t[aln_length];
	for (j=0; j<aln_length; ++j)
		n_gaps[j]=0;
	for (i=0; i<n_seqs; ++i)
	{
		const Sequence &tmp_seq = aln[i];
		for (j=0; j<aln_length; ++j)
		{
			if (tmp_seq[j] == '-')
				++n_gaps[j];
		}
	}

	for (j=0; j<aln_length; ++j)
		n_gaps[j] = ((static_cast<double>(n_gaps[j])/static_cast<double>(n_seqs)) <= gap_threshold) ? 1 : 0;

	return n_gaps;
}


size_t *
seqs_in_test(const Alignment &ref, const Alignment &test, bool ignore_missing_seqs)
{
	size_t n_seqs = ref.n_seqs();
	size_t *ids = new size_t[n_seqs];

	size_t n_seqs_test = test.n_seqs();
	size_t i,j;
	string name;
	for (i=0; i<n_seqs; ++i)
	{
		name = ref[i].name();
		for (j=0; j <n_seqs_test; ++j)
		{
			if (test[j].name() == name)
			{
				ids[i] = j;
				break;
			}
		}
		if (j==n_seqs_test)
		{
			if (ignore_missing_seqs)
				ids[i] =j;
			else
			{
				delete[] ids;
				return NULL;
			}
		}

	}
	return ids;
}



/*
 * \brief Calculates the sum of pairs similarity of two alignments.
 * @param[in] ref The first alignment
 * @param[in] test The second alignment
 * @param[in] id Connecting the sequences of the reference alignment to the test alignment.
 * @param[in] use_columns The columns to use.
 * @param[out] similarity The similarity of the alignment. The first is the number of common pairs, the second the number of total aligned pairs.
 */
void
pair_sim(const Alignment &ref, const Alignment &test, size_t *id, size_t *use_columns, pair<double, double> &similarity)
{
	size_t i,j,k,l;
	size_t n_columns=ref.size();
	size_t n_seqs = ref.n_seqs();
	size_t test_size = test.size();
	int seq1_ref_pos, seq2_ref_pos, seq1_test_pos, seq2_test_pos;
	double total_pair = 0, same_pair =0;
	int chunk = 10;
	#pragma omp parallel shared(ref, test,chunk, use_columns, test_size,id) private(i,l,k,j, seq1_ref_pos,seq2_ref_pos ,seq1_test_pos ,seq2_test_pos)
	{
	#pragma omp for schedule(dynamic,chunk) reduction(+:same_pair) reduction(+:total_pair) nowait
	for (i=0; i<n_seqs; ++i)
	{
		if (id[i]==test.n_seqs())
			continue;
		const Sequence &seq1_ref = ref[i];
		const Sequence &seq1_test = test[id[i]];

		for (j=i+1; j<n_seqs; ++j)
		{
			if (id[j]==test.n_seqs())
				continue;
			const Sequence &seq2_ref = ref[j];
			const Sequence &seq2_test = test[id[j]];
			seq1_ref_pos = seq2_ref_pos = seq1_test_pos = seq2_test_pos = -1;
			l=0;
			for (k=0; k<n_columns; ++k)
			{
				if (seq1_ref[k] != '-')
					++seq1_ref_pos;
				if (seq2_ref[k] != '-')
					++seq2_ref_pos;
				if (!use_columns[k])
					continue;

				if ((seq1_ref[k] != '-') && (seq2_ref[k] != '-'))
				{
					++total_pair;
					while (l<test_size)
					{
						if (seq1_test[l] != '-')
							++seq1_test_pos;
						if (seq2_test[l] != '-')
							++seq2_test_pos;
						++l;
						if (seq1_test_pos==seq1_ref_pos)
							break;
					}
					if ((seq1_test[l-1] != '-') && (seq2_test[l-1] != '-') && (seq1_test_pos==seq1_ref_pos) && (seq2_test_pos==seq2_ref_pos) )
						++same_pair;
				}

			}

		}

	}
	}//end OpenMP
	similarity.first=same_pair;
	similarity.second=total_pair;
}


/*
 * \brief Calculates the column similarity of two alignments.
 * @param[in] ref The first alignment
 * @param[in] test The second alignment
 * @param[in] id Connecting the sequences of the reference alignment to the test alignment.
 * @param[in] use_columns The columns to use.
 * @param[out] similarity The similarity of the alignment. The first is the number of common columns, the second the number of total columns.
 */
void
column_sim(const Alignment &ref, const Alignment &test, size_t *id, size_t *use_columns, pair<double, double> &similarity)
{

	// Sequences are encoded by position values to make comparison easier.
	size_t n_seqs_ref=ref.n_seqs();
	size_t n_seqs_test=test.n_seqs();
	size_t n_cols_ref = ref.size();
	size_t n_cols_test = test.size();

	// get memory for the recoding
	int **ref_num = new int*[n_seqs_ref];
	int **test_num = new int*[n_seqs_ref];
	size_t i,j;
	for (i = 0; i<n_seqs_ref; ++i)
	{
		if (id[i] != n_seqs_test)
		{
			ref_num[i] = new int[n_cols_ref];
			test_num[i] = new int[n_cols_test];
		}
		else
		{
			ref_num[i] = NULL;
			test_num[i] = NULL;
		}
	}

	// recode sequences each character is recoded with the position in the gapless sequence
	// gaps are encoded with -1
	int pos;
	int *num;
	for (i=0; i<n_seqs_ref; ++i)
	{
		if (ref_num[i]!=NULL)
		{
			const Sequence &ref_seq= ref[i];
			num = ref_num[i];
			pos=-1;
			for (j=0; j<n_cols_ref; ++j)
				num[j] = (ref_seq[j] == '-')? -1 : ++pos;

			const Sequence &test_seq= test[id[i]];
			num = test_num[i];
			pos=-1;
			for (j=0; j<n_cols_test; ++j)
				num[j] = (test_seq[j] == '-')? -1 : ++pos;
		}
	}


	// scan the alignment for similar columns.
	double same_column=0, total_column=0;
	size_t test_col_pos = 0 ,k;
	int val;
	for (j=0; j<n_cols_ref; ++j)
	{
		if (!use_columns[j])
			continue;
		++total_column;
		for (i=0; i<n_seqs_ref; ++i)
			if ((ref_num[i] != NULL) && (ref_num[i][j] >= 0))
				break;
		if (i==n_seqs_ref)
			continue;
		val=ref_num[i][j];
		while ((test_num[i][test_col_pos] < val) && (test_col_pos < n_cols_test))
			++test_col_pos;
		for (k=0; k<n_seqs_ref; ++k)
		{
			if ((ref_num[k] != NULL) && (ref_num[k][j] != test_num[k][test_col_pos]))
				break;
		}
		if (k==n_seqs_ref)
			++same_column;
	}

	similarity.first=same_column;
	similarity.second=total_column;

	// free memory
	for (i=0; i<n_seqs_ref; ++i)
	{
		if (test_num[i] != NULL)
		{
			delete[] test_num[i];
			delete[] ref_num[i];
		}
	}
	delete[] test_num;
	delete[] ref_num;
}




pair<double, double>
sim(const Alignment &ref, const Alignment &test, double gap_limit, const string &mode, bool ignore_missing_seqs)
{
	pair<double, double> similarity;
	size_t *use_column = use_columns(ref, gap_limit);
	size_t *id = seqs_in_test(ref, test, ignore_missing_seqs);
	if (id==NULL)
		throw Alignment_excep("Sequence not found!");

	size_t n_seqs = ref.n_seqs();

	for (size_t i=0; i<n_seqs; ++i)
	{
		if (id[i] != test.n_seqs())
		{
			if (!seq_check(ref[i], test[id[i]]))
				throw Alignment_excep("Sequences differ!");
		}
	}
	if (mode == "pair")
		pair_sim(ref, test, id, use_column, similarity);
	else
		column_sim(ref, test, id, use_column, similarity);

	delete[] use_column;
	delete[] id;

	return similarity;
}



Sequence *
simple_consensus(const Alignment &aln, float threshold, char ambiguous)
{
	size_t n_seqs = aln.n_seqs();
	size_t aln_length = aln.size();
	Sequence *consensus= new Sequence(aln.file(), "consensus", aln_length);
	Matrix<size_t> counter = Matrix<size_t>(26,aln_length,0);
	size_t i,j;

	// counting occurrences of characters
	for (i=0; i<n_seqs; ++i)
	{
		const Sequence &seq=aln[i];
		for (j=0; j<aln_length; ++j)
			if ((seq[j]!='-') && (seq[j]!='.'))
//			{
				//cout << seq[j] << endl;
				++counter[tolower(seq[j])-97][j];
//			}
	}

	// find character surpassing threshold
	size_t min=ceil(threshold*n_seqs);
	char c;
	for (j=0; j<aln_length; ++j)
	{
		for (i=0; i<26; ++i)
		{
			if (counter[i][j] >= min)
				break;
		}
		c= (i==26) ? ambiguous : i+'A';
		consensus->append(c);
	}
	return consensus;
}




}
}

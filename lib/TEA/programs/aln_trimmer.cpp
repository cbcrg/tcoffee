
// C-header
#include <cstdlib>
#include <cstdio>

// C++ header
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// boost header
#include <boost/program_options.hpp>

//BioTools header
#include "../lib/Sequence/Alignment.h"
#include "../lib/Sequence/aln_analysis.h"



using namespace std;
namespace po = boost::program_options;
using namespace BioTools::Seq;
using namespace BioTools::Utils;


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


void
pair_sim(const Alignment &ref, const Alignment &test, size_t *id, size_t *use_columns, double *similarity)
{
	size_t i,j,k,l;
	size_t n_columns=ref.size();
	size_t n_seqs = ref.n_seqs();
	size_t test_size = test.size();
	int seq1_ref_pos, seq2_ref_pos, seq1_test_pos, seq2_test_pos;
	size_t total_pair = 0, same_pair =0;
	int chunk = 20;
	#pragma omp parallel shared(ref, test,chunk, use_columns, test_size,id) private(same_pair, total_pair, i,l,k,j, seq1_ref_pos,seq2_ref_pos ,seq1_test_pos ,seq2_test_pos)
	{
		#pragma omp for schedule(dynamic,chunk) nowait
		for (i=0; i<n_seqs; ++i)
		{
			if (id[i]==test.n_seqs())
				continue;
			const Sequence &seq1_ref = ref[i];
			const Sequence &seq1_test = test[id[i]];
			for (j=i+1; j<n_seqs; ++j)
			{
				same_pair=0;
				total_pair=0;
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
				if (total_pair)
				{
					#pragma omp critical (pair_sim)
					{
						similarity[i]+=(static_cast<double>(same_pair)/static_cast<double>(total_pair));
						similarity[j]+=(static_cast<double>(same_pair)/static_cast<double>(total_pair));
					}
				}
			}


		}
	}//end OpenMP
// 	for (size_t l=0; l<n_seqs; ++l)
// 	{
// 		printf(" %f", similarity[l]);
// 	}
// 	printf("\n");
}




double*
sim(const Alignment &ref, const Alignment &test)
{

	size_t *use_column = new size_t[ref.size()];
	size_t i;
	for (i=0; i<ref.size(); ++i)
		use_column[i] = 1;
	size_t *id = seqs_in_test(ref, test, 0);
	if (id==NULL)
		throw Alignment_excep("Sequence not found!");
	size_t n_seqs = ref.n_seqs();

	for (i=0; i<n_seqs; ++i)
	{
		if (id[i] != test.n_seqs())
		{
			if (!seq_check(ref[i], test[id[i]]))
				throw Alignment_excep("Sequences differ!");
		}
	}

	double *similarity = new double[ref.n_seqs()];
	for (i=0; i<ref.n_seqs(); ++i)
		similarity[i]=0;
	pair_sim(ref, test, id, use_column, similarity);


	delete[] use_column;
	delete[] id;

	return similarity;
}


void
write_seqs(string &aln_out_f, Alignment &aln)
{
	size_t n_seqs = aln.n_seqs();
	size_t n_cols = aln.size();
	size_t j;
	FILE *aln_F = my_fopen(aln_out_f.c_str(), "w");
	char c;
	for (size_t i=0; i<n_seqs; ++i)
	{
		const Sequence &seq = aln[i];
		if (seq.comment().empty())
			fprintf(aln_F, ">%s\n", seq.name().c_str());
		else
			fprintf(aln_F, ">%s %s\n", seq.name().c_str(), seq.comment().c_str());
		for (j=0; j<n_cols; ++j)
		{
			if ((c=seq[j])!='-')
				fprintf(aln_F, "%c", c);
		}
		fprintf(aln_F, "\n");
	}
	fclose(aln_F);
}


int compare2_my (const void * a, const void * b)
{
	return ( ((pair<double, size_t>*)a)->first - ((pair<double, size_t>*)b)->first );
}

int compare_my (const void * a, const void * b)
{
	return ( ((pair<double, size_t>*)b)->first - ((pair<double, size_t>*)a)->first );
}

void
disagree_trim(const Alignment *alns, size_t n_alns, double percentage, vector<size_t> &to_del)
{
	size_t i,j;
	size_t n_seqs=alns[0].n_seqs();
	pair<double, size_t> *similarity = new pair<double, size_t>[n_seqs];

	for (j=0; j<n_seqs; ++j)
	{
		similarity[j].first=0;
		similarity[j].second=j;
	}
	double *sim_tmp;
	size_t k;
	for (i=0; i<n_alns; ++i)
	{
		for (j=i+1; j<n_alns; ++j)
		{
			sim_tmp = sim(alns[i], alns[j]);
			for (k=0; k<n_seqs; ++k)
				similarity[k].first += sim_tmp[k]/(n_seqs-1);
			delete[] sim_tmp;
			sim_tmp = sim(alns[j], alns[i]);
			for (k=0; k<n_seqs; ++k)
				similarity[k].first += sim_tmp[k]/(n_seqs-1);
			delete[] sim_tmp;
		}
	}


	qsort(similarity, n_seqs, sizeof(pair<double, size_t>), compare2_my);
	double comp_sim = 0;
	for (k=0; k<n_seqs; ++k)
	{
		similarity[k].first =1-(similarity[k].first/(n_alns*(n_alns-1)));
		comp_sim += similarity[k].first;
	}


	double so_far = 0;
	k=0;
	while (so_far/comp_sim < percentage)
	{
		so_far+=similarity[k].first;
		to_del.push_back(similarity[k].second);
		++k;
	}
}


/**
 * \brief Identifies sequences which are introducing the most gaps.
 * \param aln The alignment
 * \param percentage The thresold to use
 * \param[out] The ids of the sequences to delete
 */
void
gap_trim(const Alignment &aln, double percentage, vector<size_t> &to_del)
{
	size_t aln_l = aln.size();
	size_t n_seqs = aln.n_seqs();
	size_t i,j;

	// counting gaps in each column
	size_t *gap_counter = new size_t[aln_l];
	for (j=0; j<aln_l; ++j)
		gap_counter[j] = 0;

	for (i=0; i<n_seqs; ++i)
	{
		const Sequence &seq = aln[i];
		for (j=0; j<aln_l; ++j)
		{
			if (seq[j]== '-')
				++gap_counter[j];
		}
	}

	pair<double, size_t> *seq_val = new pair<double, size_t>[n_seqs];
	for (i=0; i<n_seqs; ++i)
	{
		seq_val[i].first = 0; // gap responsibility of this sequence
		seq_val[i].second = i; //sequence id
	}

	//calculates for each sequence the responsibility of gaps
	for (i=0; i<n_seqs; ++i)
	{
		const Sequence &seq = aln[i];
		for (j=0; j<aln_l; ++j)
		{
			if ((seq[j]!='-') && (gap_counter[j] != n_seqs))
				seq_val[i].first+=(static_cast<double>(gap_counter[j])/static_cast<double>(n_seqs-gap_counter[j])); // (number of gaps)/(number of non-gaps)
		}
	}

	// calculate total value
	double total=0;
	for (i=0; i<n_seqs; ++i)
		total+=seq_val[i].first;

	// sort by responsibility
	qsort(seq_val, n_seqs, sizeof(pair<double, size_t>), compare_my);

	// delete sequences until threshold is reached
	double so_far = 0;
	i=0;
	while (so_far/total < percentage)
	{
		so_far+=seq_val[i].first;
		to_del.push_back(seq_val[i].second);
		++i;
	}


	delete[] gap_counter;
}




int
main(int argc, char *argv[])
{
	unsigned int n_threads;
	double threshold;
	vector<string> aln_fs;
	string out_f, del_out_f, keep_f="";
	po::options_description general("General options");
	general.add_options()
	("help,h", "Produces this help message")
	("in,i", po::value<vector<string> >(&aln_fs)->multitoken(), "Input file (fasta format)")
	("out,o", po::value<string>(&out_f), "The trimmed sequence set")
	("del_out,d", po::value<string>(&del_out_f), "The deleted sequences names")
	("threshold,t", po::value<double>(&threshold)->default_value(5), "Percentage threshold to use.")
	("n_threads,n", po::value<unsigned int>(&n_threads)->default_value(1), "Number of threads to use.")
	("keep,k", po::value<string>(&keep_f), "File in fasta format with sequences to keep anyway")
	;

	po::options_description all("aln_trimmer v1.0.\n\nAllowed options are displayed below.");
	all.add(general);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << all<< "\n";
		return EXIT_SUCCESS;
	}

	omp_set_num_threads(n_threads);
	threshold/=100;

	size_t n_alns=aln_fs.size();
	Alignment *alns = new Alignment[n_alns];
	for (size_t i=0; i<n_alns; ++i)
	{
		try
		{
			alns[i].read(aln_fs[i]);
		}
		catch (bad_alloc&)
		{
			fprintf(stderr, "ERROR: Not able to allocate enough memory for alignment %s.\n", aln_fs[i].c_str());
			exit(EXIT_FAILURE);
		}
		catch (Alignment_excep &aln_e)
		{
			fprintf(stderr, "ERROR! Format could not be identified %s\n", aln_fs[i].c_str());
			exit(EXIT_FAILURE);
		}
		catch (My_IO_Exception &e)
		{
			fprintf(stderr, "ERROR! Alignment file %s could not be opened: %s\n", aln_fs[i].c_str(), e.what());
			exit(EXIT_FAILURE);
		}
		alns[i].sort("name");
	}

	set<string> names2keep;
	if (!keep_f.empty())
	{
		Alignment keep_aln;
		try
		{
			keep_aln.read(keep_f);
		}
		catch (bad_alloc&)
		{
			fprintf(stderr, "ERROR: Not able to allocate enough memory for alignment %s.\n", keep_f.c_str());
			exit(EXIT_FAILURE);
		}
		catch (Alignment_excep &aln_e)
		{
			fprintf(stderr, "ERROR! Format could not be identified %s\n", keep_f.c_str());
			exit(EXIT_FAILURE);
		}
		catch (My_IO_Exception &e)
		{
			fprintf(stderr, "ERROR! Alignment file %s could not be opened: %s\n", keep_f.c_str(), e.what());
			exit(EXIT_FAILURE);
		}
		size_t n_keep = keep_aln.n_seqs();
		for (size_t i=0; i<n_keep; ++i)
			names2keep.insert(keep_aln[i].name());
	}

	vector<size_t> to_del;
	if (n_alns > 1)
		disagree_trim(alns, n_alns, threshold, to_del);
	else
		gap_trim(alns[0], threshold, to_del);


	// Delete keep names form the delete list

	set<string>::iterator it_end=names2keep.end();
	for (int i = to_del.size()-1; i>=0; --i)
	{
		if (names2keep.find(alns[0][to_del[i]].name()) != it_end)
			to_del.erase (to_del.begin()+i);
	}

	// Delete sequences from the alignment when necessary!
	FILE *del_out_F = NULL;
	if (del_out_f.empty())
		del_out_F = stderr;
	else
		del_out_F = fopen(del_out_f.c_str(), "w");
	if (!to_del.empty())
	{
		for (size_t i=0; i<to_del.size(); ++i)
			fprintf(del_out_F, ">%s\n", alns[0][to_del[i]].name().c_str());
		alns[0].delete_seqs(to_del);
	}
	if (!del_out_f.empty())
		fclose(del_out_F);

	try
	{
		write_seqs(out_f, alns[0]);
	}
	catch (My_IO_Exception &e)
	{
		fprintf(stderr, "ERROR! Alignment file %s could not be opened: %s\n", out_f.c_str(), e.what());
		exit(EXIT_FAILURE);
	}

	delete[] alns;

	return EXIT_SUCCESS;
}

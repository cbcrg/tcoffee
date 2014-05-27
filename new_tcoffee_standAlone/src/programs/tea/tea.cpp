#include <iostream>
#include "tea.h"
#include "Version.hpp"

using namespace std;

using namespace BioTools::Seq;
using namespace BioTools::Utils;
using namespace BioTools::HMM;
using namespace BioTools::Clustering;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
#define NDEBUG

bool verbose;
bool debug;


Alignment *
produce_aln_from_sequences(SequenceSet &seq_set, const Aln_opts &aln_opts, const Tree_opts tree_opts)
{
	// pairwise alignments to build library
	if (verbose)
		printf("Calculate pairwise alignments\n");

	size_t n_seqs=seq_set.n_seqs();

	Library lib(seq_set);
	Matrix<float> dist_mat(n_seqs, n_seqs);	
	char seq_type=toupper(aln_opts.seq_type);
	
	if (aln_opts.pair_aln_method=="hmm")
		all_hmm_pairs(seq_set, lib, dist_mat, seq_type);
	else if (aln_opts.pair_aln_method == "gotoh")
		all_gotoh_pairs(seq_set, lib, aln_opts.score_mat, aln_opts.gapopen, aln_opts.gapext);

	size_t i;
	// tree construction
	if (verbose)
		printf("Calculate/Read tree\n");

	//Matrix<float> dist_mat(n_seqs, n_seqs);

	Tree tree;
	if (!tree_opts.tree_i_f.empty())
	{
		tree.read_newick(tree_opts.tree_i_f);
		tree.match_seq_set(seq_set);
	}
	else
	{
		// calculate distance matrix if needed
		if (tree_opts.dist_mat_method=="kmer")
			make_dist_mat(seq_set, tree_opts.tree_k, tree_opts.alphabet, dist_mat);

		// calculate tree
		vector<std::string> names;
		for (i=0; i<n_seqs; ++i)
			names.push_back(seq_set[i].name());

		if (tree_opts.tree_method=="nj")
			tree.nj(dist_mat, names);
		else if (tree_opts.tree_method=="upgma")
			tree.upgma(dist_mat, names);
		else
		{
			fprintf(stderr, "Error! Tree method unknown\n");
			exit(EXIT_FAILURE);
		}
	}


	if (!tree_opts.tree_o_f.empty())
		tree.print_newick(tree_opts.tree_o_f);

	if (verbose)
		printf("Relaxing library\n");
	lib.relax();


	// calculate msa using scores from pairwise alignments
	std::map<Match, int> match_points;
	std::map<Match, int>::iterator it, it_end;

	vector<Alignment*> alns;
	alns.resize(n_seqs);
	for (i=0; i< n_seqs; ++i)
	{
		seq_set[i].aln_id(i);
		alns[i] = new Alignment(seq_set, i);
	}
//	printf("NSEQS: %i\n", n_seqs);
	//lib.print(seq_set, "library.txt");
	Alignment *aln = progressive_aln(tree.root() , lib, alns);
	
	return aln;
}



Alignment *
produce_aln_from_profiles(vector<Alignment*> &aln_set, const Aln_opts &aln_opts, const Tree_opts tree_opts)
{
	// pairwise alignments to build library
	if (verbose)
		printf("Calculate pairwise alignments\n");

	size_t n_alns=aln_set.size();

	Sequence *tmp_seq;
	SequenceSet consensus_set;
	size_t i;

	for (i=0; i<n_alns; ++i)
	{
		tmp_seq = BioTools::Seq::simple_consensus(*(aln_set[i]));
		tmp_seq->id(i);
		consensus_set.addSeq(tmp_seq);
	}

	Library lib(aln_set);
	Matrix<float> dist_mat(n_alns, n_alns);
	char seq_type=toupper(aln_opts.seq_type);
	
	if (aln_opts.pair_aln_method=="hmm")
		all_hmm_pairs(aln_set, lib, dist_mat, seq_type, 0, n_alns-1);
	//else if (aln_opts.pair_aln_method == "gotoh")
	//	all_gotoh_pairs(consensus_set, lib, aln_opts.score_mat, aln_opts.gapopen, aln_opts.gapext);


	// tree construction
	if (verbose)
		printf("Calculate/Read tree\n");

	//Matrix<float> dist_mat(n_seqs, n_seqs);

	Tree tree;
	if (!tree_opts.tree_i_f.empty())
	{
		tree.read_newick(tree_opts.tree_i_f);
		tree.match_seq_set(consensus_set);
	}
	else
	{
		// calculate distance matrix if needed
		if (tree_opts.dist_mat_method=="kmer")
			make_dist_mat(consensus_set, tree_opts.tree_k, tree_opts.alphabet, dist_mat);

		// calculate tree
		vector<std::string> names;
		for (i=0; i<n_alns; ++i)
			names.push_back(consensus_set[i].name());

		if (tree_opts.tree_method=="nj")
			tree.nj(dist_mat, names);
		else if (tree_opts.tree_method=="upgma")
			tree.upgma(dist_mat, names);
		else
		{
			fprintf(stderr, "Error! Tree method unknown\n");
			exit(EXIT_FAILURE);
		}
	}


	if (!tree_opts.tree_o_f.empty())
		tree.print_newick(tree_opts.tree_o_f);


	if (verbose)
		printf("Relaxing library\n");
	//lib.print("lib.txt");
	lib.relax();

	// calculate msa using scores from pairwise alignments
	std::map<Match, int> match_points;
	std::map<Match, int>::iterator it, it_end;

/*	vector<Alignment*> alns;
	alns.resize(n_seqs);
	for (i=0; i< n_seqs; ++i)
		alns[i] = new Alignment(seq_set, i);

	//lib.print(seq_set, "library.txt");*/
	Alignment *result_aln = progressive_aln(tree.root() , lib, aln_set);
	return result_aln;
}



/**
 * Calculates
 * @param[in] root The root of the tree to use for the alignment
 * @param[in] lib The library to use.
 * @param[in] alns The input alignments.
 * @return The final alignment.
 */
Alignment *
progressive_aln(TreeNode *root, Library &lib, vector<Alignment*> alns)
{
	Alignment *aln1=NULL, *aln2;
	stack<pair<TreeNode*, unsigned int> > to_do;
	to_do.push(pair<TreeNode*, int>(root, 0));
	std::map<Match, int> match_points;
	Gotoh_Int_array dp_mat;
	Gotoh_Trace_array trace_mat;


	vector<pair<unsigned int, unsigned int> > gap1, gap2;
	while (!to_do.empty())
	{
		pair<TreeNode*, unsigned int> &current = to_do.top();
		if (current.first->children.empty())
			to_do.pop();
		else if (current.second == current.first->children.size())
		{
			aln1 = alns[current.first->children[0]->id];
			aln2 = alns[current.first->children[1]->id];
			lib.get(*aln1 , *aln2, match_points);
			gotoh_chaining(aln1->size(), aln2->size(), match_points, dp_mat, trace_mat, 0, 0);
			gotoh_chain_trace(gap1, gap2, trace_mat);
			aln1->merge(gap1, *aln2, gap2);
			current.first->id = current.first->children[0]->id;
			to_do.pop();
		}
		else
		{
			to_do.push(pair<TreeNode*, int>(current.first->children[current.second], 0));
			++current.second;
		}
	}

	return aln1;
}


Alignment *
large_scale_aln(SequenceSet &seq_set, Aln_opts &aln_opts, Tree_opts &tree_opts)
{
	//Alignment *result_aln=new Alignment;

	// produce km-tree
	std::vector<Vec_double_ptr>* vec_set = seqset2vecs_kmer(seq_set, tree_opts.tree_k, tree_opts.alphabet);
	KM_node* root = hierarchical_kmeans(*vec_set, seq_set, tree_opts.cluster_size, "kkz", 0.001);

	// parse tree
	stack<pair<KM_node*, size_t> > to_do;
	to_do.push(pair<KM_node*, size_t>(root,0));
	KM_node *current_node;
	size_t child_id; //the child to process
	size_t j;
	int dcount = 0;
	map<int,Alignment*> alns;
	vector<int> ids;
	while (!to_do.empty())
	{
		current_node = to_do.top().first;
		child_id = to_do.top().second;
		if (current_node->children.empty())
		{       
			if(current_node->end-current_node->start > 1)
			{
				// produce sequence alignment
				SequenceSet tmp_seq_set;
				for (size_t i=current_node->start; i<current_node->end; ++i)
				{ 	
					//tmp_seq_set.share(seq_set, i);
					tmp_seq_set.share(seq_set, (*vec_set)[i]->id());
					
					/*const BioTools::Seq::Sequence &seq =seq_set[(*vec_set)[i]->id()];
					size_t len=seq.size();
					for (j=0; j<len; ++j)
					      cout<< seq[j];
					cout<< endl;*/
				}
				Alignment *seq_aln = produce_aln_from_sequences(tmp_seq_set, aln_opts, tree_opts);
				alns.insert(pair<int, Alignment*>(current_node->id,seq_aln));
			}
			else
			{
				// turn single sequence into alignment
				//Alignment *seq_aln = new Alignment(seq_set, current_node->start);
				Alignment *seq_aln = new Alignment(seq_set, (*vec_set)[current_node->start]->id());
				alns.insert(pair<int, Alignment*>(current_node->id,seq_aln));
			}
			to_do.pop();
		}
		else
		{
			if (child_id == current_node->children.size())
			{
				// produce profile profile alignment
				vector<Alignment *> tmp_aln_set;
				for (j=0; j<child_id; ++j)
				{

					tmp_aln_set.push_back(alns[current_node->children[j]->id]);
					if (debug)
					     tmp_aln_set[j]->write(std::to_string(static_cast<long long>(dcount))+"_"+std::to_string(static_cast<long long>(j)), "fasta");
					tmp_aln_set[j]->id(j);
					replace_char('-','.', *tmp_aln_set[j]);

				}
				alns.insert(pair<int, Alignment*>(current_node->id,produce_aln_from_profiles(tmp_aln_set, aln_opts, tree_opts)));
				to_do.pop();
				dcount++;
			}
			else
			{
				++to_do.top().second;
				to_do.push(pair<KM_node*, size_t>(current_node->children[child_id] ,0));
			}
		}
	}
	replace_char('.','-', *alns[0]);
	return alns[0];
}


/**
 * \brief Removes identical sequences from the set and inreases the occuerence pointer for the remaining one.
 * @param seq_set The sequence set to group.
 * @param seqGroup Stores the duplicates.
 */
void
group(SequenceSet &seq_set, map<string, SequenceSet> &seqGroup)
{
	map<string, Sequence*> identicals;
	for (int i=seq_set.n_seqs()-1; i>=0; --i)
	{
		//Reduction to 14 ("A,C,D,EQ,FY,G,H,IV,KR,LM,N,P,ST,W");
		
		string seq2="", seq1=seq_set[i].sequence();
		//int ssize=seq1.length(); 
		
		//Reduction to 10 ("AST,C,DEN,FY,G,H,ILMV,KQR,P,W");
		/*for(int j=0; j<ssize; j++){
		  char c=toupper(seq1[j]);
		  if(c=='S' || c=='T'){ seq2+='A';}
		  else if(c=='E' || c=='N'){ seq2+='D';}
		  else if(c=='Y'){ seq2+='F';}
		  else if(c=='L' || c=='M' || c=='V'){ seq2+='I';}
		  else if(c=='Q' || c=='R'){ seq2+='K';}
		  else{ seq2+=c;}
		}*/
		
		//Reduction to 6 ("AST,CP,DEHKNQR,FWY,G,ILMV");
		/*for(int j=0; j<ssize; j++){
		  char c=toupper(seq1[j]);
		  if(c=='S' || c=='T'){ seq2+='A';}
		   else if(c=='P'){ seq2+='C';}
		  else if(c=='E' || c=='H' || c=='K' || c=='N' || c=='Q' || c=='R'){ seq2+='D';}
		  else if(c=='Y' || c=='W'){ seq2+='F';}
		  else if(c=='L' || c=='M' || c=='V'){ seq2+='I';}
		  else{ seq2+=c;}
		}*/
		
		auto tmp=identicals.insert(pair<string, Sequence*>(seq_set[i].sequence(), &(seq_set[i])));
		//auto tmp=identicals.insert(pair<string, Sequence*>(seq2, &(seq_set[i])));
		if (tmp.second==false)
		{
			tmp.first->second->ocurrence(tmp.first->second->ocurrence()+1);
			seqGroup[tmp.first->second->name()].transfer(seq_set, i);
		}
	}

}

void
ungroup(Alignment &aln, map<string, SequenceSet> &seqGroup)
{
	map<string, SequenceSet>::iterator it, itEnd=seqGroup.end();
	size_t j;
	std::vector<std::pair<unsigned int, unsigned int> > gaps;
	size_t alnSize=aln.n_seqs();
	size_t alnLength = aln[0].size();
	unsigned int gapLength;
	unsigned int charPos;
	for (size_t i =0; i<alnSize; ++i)
	{
		if ((it=seqGroup.find(aln[i].name()))!=itEnd)
		{
			SequenceSet &set=it->second;
			gaps.clear();
			charPos=0;
			for (j=0; j<alnLength; ++j)
			{
				if (aln[i][j] == '-')
				{
					gapLength=1;
					while ((++j<alnLength) && (aln[i][j] == '-'))
						++gapLength;
					gaps.push_back(pair<unsigned int, unsigned int>(charPos, gapLength));

				}
				++charPos;
			}
			for (j=0; j<set.n_seqs(); ++j)
				set[j].insert_gaps(gaps);
			aln.transfer(set);
		}
	}
}

int
main(int argc, char *argv[])
{
	string seq_f;
	string out_f;
	string out_format;

	Tree_opts tree_opts;
	Aln_opts aln_opts;
	bool noSeqSort;
	vector<string> profile_fs;

	int n_threads;
	po::options_description general("General options");
	general.add_options()
		("help,h", "Produces this help message")
		("in,i", po::value<string>(&seq_f), "Input file (fasta format)")
		("profiles,p", po::value<vector<string> >(&profile_fs)->multitoken(), "Profile files")
		("out,o", po::value<string>(&out_f), "Output file")
		("format,f", po::value<string>(&out_format)->default_value("fasta"), "The output format to use. Options: fasta, fasta_seq, clustalw, msf, phylip_s, phylip_i")
		("n_threads,n", po::value<int>(&n_threads)->default_value(1), "Number of threads to use")
		("seq_type,s", po::value<char>(&aln_opts.seq_type)->default_value('A'), "Type of input sequences. Options: auto detect (A), protein (P), nucleotide (N)")
		("verbose,v", po::value<bool>(&verbose)->default_value(false)->zero_tokens(),"Print information")
		("debug,d", po::value<bool>(&debug)->default_value(false)->zero_tokens(),"Print the profiles created on every step")
	;

	po::options_description tree_options("Guide tree options");
	tree_options.add_options()

		//("tree_in", po::value<string>(&tree_opts.tree_i_f), "A file with the tree in newick format (not supported yet)")
		("tree_out", po::value<string>(&tree_opts.tree_o_f), "Filename to write the tree into (beta)")
		("tree_method,t", po::value<string>(&tree_opts.tree_method)->default_value("nj"), "The tree construction method to use. Options: nj, upgma")
		("alphabet,a", po::value<string>(&tree_opts.alphabet)->default_value("full"), "The alphabet to use")
		("dist_mat,", po::value<string>(&tree_opts.dist_mat_method)->default_value("pairs"), "Distance matrix method")
		("tree_k_mer", po::value<unsigned short>(&tree_opts.tree_k)->default_value(2), "k-mer size to use for tree construction")
		("cluster_size", po::value<unsigned short>(&tree_opts.cluster_size)->default_value(100), "The maximal size of a cluster before subclustering");
	;

	po::options_description pw_options("Pairwise alignment options");
	pw_options.add_options()
		("pair_aln_method,m", po::value<string>(&aln_opts.pair_aln_method)->default_value("hmm"), "The pairwise method to use")
		("gapopen,g", po::value<int>(&aln_opts.gapopen)->default_value(-12), "Gap opening costs")
		("gapext,e", po::value<int>(&aln_opts.gapext)->default_value(-2), "Gap extension costs")
		("matrix,m", po::value<string>(&aln_opts.matrix_name)->default_value("BLOSUM62"), "Score matrix to use")
	;

	bool noIdCollapse;
	po::options_description hidden_options("The secrets...");
	hidden_options.add_options()
		("noIdCollapse", po::value<bool>(&noIdCollapse)->default_value(true)->zero_tokens(), "no collapse copy seqs")
		("noInputSort", po::value<bool>(&noSeqSort)->default_value(true)->zero_tokens(), "no sorting of the input sequences")
	;
	po::options_description all("tea (" + version + ")\n\nAllowed options are displayed below.");
	all.add(general).add(pw_options).add(tree_options);
	po::options_description allall;
	allall.add(general).add(pw_options).add(tree_options).add(hidden_options);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(allall).run(), vm);
	po::notify(vm);

	// Change alphabet and kmer in case of nucleotide seqs


	if (vm.count("help")) {
		cout << all<< "\n";
		return EXIT_SUCCESS;
	}

	// No input files given
	if ((seq_f.empty()) && (profile_fs.empty()))
	{
		fprintf(stderr, "Error! No input file given.\n");
		return EXIT_FAILURE;
	}

	// Profiles and sequence files given (not supported)
	if ((!seq_f.empty()) && (!profile_fs.empty()))
	{
		fprintf(stderr, "Error! Sequence files and alignment files cannot be used together.\n");
		return EXIT_FAILURE;
	}

	omp_set_num_threads(n_threads);
	char *bioTools_data= getenv("BIOTOOLS_DATA");

	if (aln_opts.pair_aln_method != "hmm")
	{
		fs::path matrix_n(aln_opts.matrix_name);
		fs::file_status mat_stat = status(matrix_n);
		if ((!aln_opts.matrix_name.empty()) && ((!exists(mat_stat)) || is_directory(mat_stat)))
		{
			fs::path matrix_path;
			string mat_path;
			if (bioTools_data == NULL)
			{
				matrix_path = fs::path(string(getenv("HOME")));
				matrix_path /= fs::path(".BioTools++");
			}
			else
				matrix_path = fs::path(string(bioTools_data));

			matrix_path /= fs::path("data");
			matrix_path /= matrix_n;
			aln_opts.matrix_name = matrix_path.string();
		}

		//Scoring_Matrix score_mat;
		try
		{
			aln_opts.score_mat.read(aln_opts.matrix_name);
		}
		catch (const My_IO_Exception& ex)
		{
			fprintf(stderr, "Error! Could not open file %s: %s\n", aln_opts.matrix_name.c_str(), ex.what());
			exit(EXIT_FAILURE);
		}
	}

	size_t i;
	size_t n_seqs;
	SequenceSet seq_set;
	Alignment *result_aln;
	vector<size_t> input_order;
	if (!seq_f.empty())
	{
		// read the sequence file
		try
		{
			seq_set.read(seq_f);
		}
		catch (const My_IO_Exception& ex)
		{
			fprintf(stderr, "Error! Could not open input file %s: %s\n", seq_f.c_str(), ex.what());
			exit(EXIT_FAILURE);
		}
	
		n_seqs = seq_set.n_seqs();

		if (!noSeqSort)
		{
			input_order.resize(n_seqs);
			seq_set.sort("seq");
			for (size_t i=0; i<n_seqs; ++i)
			{
				input_order[i] = seq_set[i].id();
				seq_set[i].id(i);
			}
		}
		if (aln_opts.seq_type=='A')
			aln_opts.seq_type=identify_seq_type(seq_set[0]);
		if(aln_opts.seq_type!='P')
		{
		    tree_opts.alphabet="nfull";
		    tree_opts.tree_k=4;
		}


		if (n_seqs > tree_opts.cluster_size)
		{
			// identical sequnces are grouped
			map<string, SequenceSet> seqGroups;
			if (!noIdCollapse)
				group(seq_set, seqGroups);
			tree_opts.cluster_size=min<size_t>(tree_opts.cluster_size, seq_set.n_seqs());
			result_aln = large_scale_aln(seq_set, aln_opts, tree_opts);
			if (!noIdCollapse)
				ungroup(*result_aln, seqGroups);
		}
		else
			result_aln = produce_aln_from_sequences(seq_set, aln_opts, tree_opts);
	}
	else
	{
		vector<Alignment*> aln_set;
		n_seqs = profile_fs.size();
		//aln_set.resize(n_seqs,NULL);
		for (i=0; i< n_seqs; ++i)
		{
			aln_set.push_back(new Alignment());
			aln_set[i]->read(profile_fs[i]);
			replace_char('-','.', *aln_set[i]);
			aln_set[i]->id(i);
		}
		if (aln_opts.seq_type=='A')
			aln_opts.seq_type=identify_seq_type((*aln_set[0])[0]);
		if(aln_opts.seq_type!='P')
		{
		    tree_opts.alphabet="nfull";
		    tree_opts.tree_k=4;
		}
		result_aln = produce_aln_from_profiles(aln_set, aln_opts, tree_opts);
		replace_char('.','-', *result_aln);
	}

	// write alignment
	/*if (input_order.size() != 0)
	{
		result_aln->sort("input");
		for (size_t i=0; i<input_order.size(); ++i)
			(*result_aln)[i].id(input_order[i]);
		result_aln->sort("input");
	}*/
	result_aln->write(out_f, out_format);
	return EXIT_SUCCESS;
}

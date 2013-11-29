/*
 * msaa.cpp
 *
 *  Created on: Oct 9, 2011
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 *
 *
 * This file is part of BioTools.
 *
 * BioTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * BioTools is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

// C header
#include <cstdlib>
#include <cstdio>

// C++ header
#include <string>
#include <vector>
#include <algorithm>

// boost header
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// other headers
#include <omp.h>

// BioTools headers
#include "../lib/Sequence/Alignment.h"
#include "../lib/Sequence/aln_analysis.h"
#include "../lib/utils/filesystem.h"
#include "../lib/utils/ScoringMatrix.h"
#include "../lib/utils/Matrix.h"
#include "../lib/utils/graph_algorithms.h"

using namespace std;
using namespace BioTools::Seq;
using namespace BioTools::Utils;

namespace po = boost::program_options;
namespace fs = boost::filesystem;


typedef struct
{
	double column_trim;
	bool to_upper;
	bool to_lower;
	vector<string> seq_trim;
	vector<size_t> id_to_delete;
	vector<string> seq_to_delete;
} Mod_options;


typedef struct
{
	size_t *aln_length;
	double *id;
	size_t *identical;
	size_t *total;
	double *s_o_p;
	double gop;
	double gep;
	Scoring_Matrix *matrix;
} Analysis_options;


void
get_files_in_dir(const string & dir_path,
				 vector<fs::path> &files,
				 const vector<string> &endings,
				 bool recursive)
{
	if ( !fs::exists( dir_path ) )
		return ;
	vector<fs::path> subdirs;
	subdirs.push_back(dir_path);
	fs::path extension;
	unsigned int n_extension = endings.size();
	unsigned int i;
	while (!subdirs.empty())
	{
		fs::directory_iterator end_itr; // default construction yields past-the-end
		fs::directory_iterator itr( subdirs.back() );
		subdirs.pop_back();
		for (; itr != end_itr; ++itr )
		{
			if ( is_directory(itr->status()) )
			{
				if (recursive)
					subdirs.push_back(itr->path());
			}
			else
			{
				if (endings.empty())
				{
					files.push_back(itr->path());
				}
				else
				{
					if (itr->path().has_extension())
					{
						extension=itr->path().extension();
						for (i = 0; i < n_extension; ++i)
						{
							if (extension == endings[i])
							{
								files.push_back(itr->path());
								break;
							}
						}
					}
				}
			}
		}

	}
}

void
read_aln(const string &aln_f, Alignment &aln, FILE *log_F, bool check)
{
	aln.clear();
	try
	{
		aln.read(aln_f, check);
	}
	catch (bad_alloc&)
	{
	  fprintf(log_F, "ERROR: Not able to allocate enough memory for alignment %s.\n", aln_f.c_str());
	  exit(EXIT_FAILURE);
	}
	catch (Alignment_excep &aln_e)
	{
		fprintf(log_F, "ERROR in %s! %s\n", aln_f.c_str(), aln_e.what());
		exit(EXIT_FAILURE);
	}
	catch (My_IO_Exception &e)
	{
		fprintf(log_F, "ERROR! Alignment file %s could not be opened: %s\n", aln_f.c_str(), e.what());
		exit(EXIT_FAILURE);
	}
}

void
read_aln(const string &aln_f, Alignment &aln, const vector<string> &names, FILE *log_F, bool check)
{
	aln.clear();
	try
	{
		aln.read(aln_f, names, check);
	}
	catch (bad_alloc&)
	{
	  fprintf(log_F, "ERROR: Not able to allocate enough memory for alignment %s.\n", aln_f.c_str());
	  exit(EXIT_FAILURE);
	}
	catch (Alignment_excep &aln_e)
	{
		fprintf(log_F, "ERROR in %s! %s\n", aln_f.c_str(), aln_e.what());
		exit(EXIT_FAILURE);
	}
	catch (My_IO_Exception &e)
	{
		fprintf(log_F, "ERROR! Alignment file %s could not be opened: %s\n", aln_f.c_str(), e.what());
		exit(EXIT_FAILURE);
	}
}


bool greater_eq(float x, float y){
	return x >= y;
}


bool less_eq(float x, float y){
	return x <= y;
}


void seq_trim(Alignment &aln, vector<string> options)
{
	size_t n_seqs=aln.n_seqs();
	size_t i,j,k;
	pair<size_t, size_t> (*seq_compare_function)(const Sequence &seq1, const Sequence &seq2) = NULL;
	bool (*compare_function)(float x, float y) = NULL;
	string trim_option;
	float threshold;
	vector<size_t> seqs_to_del;
	aln.sort("seq");
	//aln.write("arg","fasta");
	pair<size_t, size_t> tmp;
	for (i=0; i<options.size(); ++i)
	{
		seqs_to_del.clear();
		Matrix<bool> mat(n_seqs, n_seqs, false);
		trim_option=options[i];
		threshold = atof(options[++i].c_str());
		if ((trim_option== "min_cov") || (trim_option== "max_cov"))
			seq_compare_function=coverage;
		else
			seq_compare_function=id;

		if ((trim_option== "min_sim") || (trim_option== "min_cov"))
			compare_function= greater_eq;
		else
			compare_function= less_eq;


		for (j=0; j<n_seqs; ++j)
		{
			for (k=j+1; k<n_seqs; ++k)
			{
				tmp=seq_compare_function(aln[j], aln[k]);
				mat[k][j]=mat[j][k]=compare_function(tmp.first*100.0/tmp.second, threshold);
			}
		}

		set<size_t> result;
		max_clique(mat, result);
		set<size_t>::iterator it, it_end=result.end();
		for (it=result.begin(); it != it_end; ++it)
			seqs_to_del.push_back(*it);
		aln.keep_seqs(seqs_to_del);
	}
	aln.sort("input");
}



void modify_aln(Alignment &aln, const Mod_options &options)
{
	if (!options.seq_trim.empty())
		seq_trim(aln, options.seq_trim);
	if (options.column_trim>=0)
		aln.column_trim(options.column_trim);
	if (options.to_upper)
		aln.to_upper();
	if (options.to_lower)
		aln.to_lower();
}


void
run_analysis(const Alignment &aln, Analysis_options &ana_opts, bool avg_only, FILE *analysis_F)
{
	double tmp_d;
	if (!avg_only)
		fprintf(analysis_F, "%s", aln.file().c_str());
	if (ana_opts.aln_length != NULL)
	{
		size_t tmp = aln.size();
		if (!avg_only)
			fprintf(analysis_F, "\t%li", tmp);
		*ana_opts.aln_length += tmp;
	}
	if (ana_opts.id != NULL)
	{
		pair<size_t, size_t> id = identity(aln);
		tmp_d = 100*exp(log(id.first)-log(id.second));
		if (!avg_only)
			fprintf(analysis_F, "\t%.1f", tmp_d);
		*ana_opts.id += tmp_d;
		*ana_opts.identical += id.first;
		*ana_opts.total += id.second;
	}
	if (ana_opts.s_o_p != NULL)
	{
		double tmp = sum_of_pairs_score(aln, *ana_opts.matrix, ana_opts.gop, ana_opts.gep);
		if (!avg_only)
			fprintf(analysis_F,"\t%.1f", tmp);
		*ana_opts.s_o_p += tmp;
	}
	if (!avg_only)
		fprintf(analysis_F,"\n");
}


void
aln_compare(const string &ref_aln_f, const vector<fs::path> &test_fs, double gap_limit, const string &mode, bool ignore_missing_seqs, FILE *analysis_F, FILE *log_F, bool check)
{
	Alignment ref_aln, test_aln;
	read_aln(ref_aln_f, ref_aln, log_F, check);

	size_t n_alns = test_fs.size();
	pair<double, double> similarity;
	for (size_t i=0; i<n_alns; ++i)
	{
		read_aln(test_fs[i].string(), test_aln, log_F, check);
		similarity = sim(ref_aln,test_aln, gap_limit, mode, ignore_missing_seqs);
		fprintf(analysis_F, "%s %.1f\n", test_fs[i].string().c_str(), 100*similarity.first/similarity.second);
	}
}


int
main(int argc, char *argv[])
{
	vector<string> input_list;
	vector<string> extract_names;
	string out_f, out2_f, out_format, analysis_f;
	vector<string> extensions;
	string new_suffix;
	bool score_aln, avg_only, no_avg, no_header, aln_length, del_gaps;
	string matrix_name, compare_mode, ref_aln;
	double gap_limit;
	bool compute_identities, write_f = false;
	int n_threads;
	Mod_options mod_options;
	mod_options.column_trim=-1;
	Analysis_options ana_opts;
	string log_f;
	bool ignore_missing_seqs, check;

	po::options_description general("General options");
	general.add_options()
		("help,h", "Produces this help message")
		("alignments,a", po::value<vector<string> >(&input_list)->multitoken(), "Alignment file(s) or directory(s)")
		("extensions,e", po::value<vector<string> >(&extensions)->multitoken(), "File extension to consider when given a directory. If none given all files will be used.")
		("output,o", po::value<string>(&out_f), "The file or directory to write the output to")
		("format,f", po::value<string>(&out_format)->default_value("FASTA"), "The output format to use")
		("suffix,s", po::value<string>(&new_suffix), "Replaces extension with a new one. In case no extension exists it will be added.")
		("analysis,A", po::value<string>(&analysis_f), "File for analysis output")
		("num_threads,n", po::value<int>(&n_threads)->default_value(1), "Number of threads to use")
		("logging,l", po::value<string>(&log_f), "Log errors in this file")
		("no-check", po::value<bool>(&check)->default_value(false)->zero_tokens(), "No checking of the alignments will be done");
	;

	po::options_description analysis_opts("Alignment analysis");
	analysis_opts.add_options()
		("average_only,v", po::value<bool> (&avg_only)->default_value(false)->zero_tokens(), "Only the average is printed")
		("no-average", po::value<bool> (&no_avg)->default_value(false)->zero_tokens(), "Do not print an average")
		("no-header", po::value<bool>(&no_header)->default_value(false)->zero_tokens(), "Do not print the header")
		("identity,i", po::value<bool> (&compute_identities)->default_value(false)->zero_tokens(), "Computes the average identity of the alignments")
		("score,S", po::value<bool>(&score_aln)->default_value(false)->zero_tokens(), "Compute the Sum-of-Pairs score")
		("length,L", po::value<bool>(&aln_length)->default_value(false)->zero_tokens(), "Computes the length of the alignment")
		("matrix,m", po::value<string>(&matrix_name)->default_value("BLOSUM62"), "Score matrix to use!")
		("gop", po::value<double>(&ana_opts.gop)->default_value(-11.0), "Gap opening costs")
		("gep", po::value<double>(&ana_opts.gep)->default_value(-1.0), "Gap extension costs")
	;

	po::options_description mani_opts("Alignment manipulation");
	mani_opts.add_options()
		("upper", po::value<bool> (&mod_options.to_upper)->default_value(false)->zero_tokens(), "Turns all sequences to uppercase")
		("lower", po::value<bool> (&mod_options.to_lower)->default_value(false)->zero_tokens(), "Turns all sequences to lowercase")
		("extract,E", po::value<vector<string> >(&extract_names)->multitoken(), "Sequences to extract")
		("column_trim,c", po::value<double >(&mod_options.column_trim), "Deletes the columns which contain more gaps then the given percentage")
		("delete_gaps,D", po::value<bool>(&del_gaps)->default_value(false)->zero_tokens(), "Turns back the alignment into sequences by deleteing all gap characters.")
		("seq_trim,t", po::value<vector<string> >(&mod_options.seq_trim)->multitoken(), "trim_option threshold. The trimming option can be one of the following min_cov, max_cov, min_sim, max_sim")
	;

	po::options_description aln_compare_opts("Alignment comparison");
	aln_compare_opts.add_options()
		("ref_aln,r", po::value<string> (&ref_aln), "The reference alignment")
		("gap_limit,G", po::value<double> (&gap_limit)->default_value(100), "Columns with less then the given percentages of gaps are considered")
		("compare_mode,C", po::value<string> (&compare_mode)->default_value("pair"), "Compare mode: pair/column")
		("ignore_mis_seqs,q", po::value<bool>(&ignore_missing_seqs)->zero_tokens()->default_value(false), "Ignores missing sequences in found in the reference but not in the test alignment")
	;

	po::options_description all("MSAA - Multiple sequence alignment analyzer v1.0  Copyright (C) 2012  Carsten Kemena\nThis program comes with ABSOLUTELY NO WARRANTY;\n\nAllowed options are displayed below. If using alignment comparison options all other options are ignored. For more information please see the Manual.");
	all.add(general).add(analysis_opts).add(mani_opts).add(aln_compare_opts);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << all<< "\n";
		return EXIT_SUCCESS;
	}

	omp_set_num_threads(n_threads);
	check = !check;
	if (del_gaps)
		out_format="fasta_seq";
	
	FILE *analysis_F;
	FILE *log_F = NULL;
	if (log_f.empty())
		log_F = stderr;
	else
	{
		try
		{
			log_F = my_fopen(log_f, "w");
		}
		catch (My_IO_Exception &e)
		{
			fprintf(log_F, "ERROR! Logfile %s could not be opened: %s\n", analysis_f.c_str(), e.what());
			exit(EXIT_FAILURE);
		}
	}
	if (analysis_f.empty())
		analysis_F = stdout;
	else
	{
		try
		{
			analysis_F = my_fopen(analysis_f, "w");
		}
		catch (My_IO_Exception &e)
		{
			fprintf(log_F, "ERROR! Analysis output file %s could not be opened: %s\n", analysis_f.c_str(), e.what());
			exit(EXIT_FAILURE);
		}
	}
	char *bioTools_data= getenv("BIOTOOLS_DATA");

	vector<fs::path> aln_fs;
	for (unsigned int i=0; i<input_list.size(); ++i)
	{
		if (fs::is_directory(input_list[i]))
			get_files_in_dir(input_list[i], aln_fs, extensions, false);
		else
			aln_fs.push_back(input_list[i]);
	}

	if (aln_fs.empty())
	{
		fprintf(stderr, "ERROR! No input files have been given/found! Please check the arguments. You can use -h/--help to get more information.\n");
		return EXIT_FAILURE;
	}
	sort(aln_fs.begin(),aln_fs.end());

	if (!ref_aln.empty())
	{
		aln_compare(ref_aln, aln_fs, gap_limit, compare_mode, ignore_missing_seqs, analysis_F, log_F, check);
		return EXIT_SUCCESS;
	}

	fs::path matrix_n(matrix_name);
	fs::file_status mat_stat = status(matrix_n);
	if ((!matrix_name.empty()) && ((!exists(mat_stat)) || is_directory(mat_stat)))
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
		matrix_name = matrix_path.string();
	}

	if (score_aln)
	{
		ana_opts.matrix = new Scoring_Matrix;
		try
		{
			ana_opts.matrix->read(matrix_name);
		}
		catch (My_IO_Exception &e)
		{
			fprintf(log_F, "ERROR! Matrix file %s could not be opened: %s\n", matrix_name.c_str(), e.what());
			exit(EXIT_FAILURE);
		}
	}
	Alignment aln;
	string aln_out_f;
	size_t pos;
	bool use_new_dir = false;
	bool change_suffix = false;
	bool analysis = false;

	fs::path new_out;

	if (!out_f.empty())
	{
		write_f = true;
		new_out = fs::path(out_f);
		if (fs::is_directory(new_out))
			use_new_dir = true;
	}

	if(!new_suffix.empty())
	{
		write_f =true;
		change_suffix = true;
	}

	if ((aln_fs.size() >1)&& (!use_new_dir) && (!out_f.empty()))
	{
		fprintf(log_F, "ERROR! If more than one file is given output has to be more than one file!\n");
		exit(EXIT_FAILURE);
	}

	//  analysis
	string header = "NAME";
	if (aln_length)
	{
		analysis = true;
		header+="\taln_length";
		ana_opts.aln_length = new size_t;
		*ana_opts.aln_length = 0;
	}
	else
		ana_opts.aln_length = NULL;

	if (compute_identities)
	{
		analysis = true;
		header +="\tid";
		ana_opts.id = new double;
		*ana_opts.id = 0;
		ana_opts.identical = new size_t;
		*ana_opts.identical = 0;
		ana_opts.total = new size_t;
		*ana_opts.total = 0;
	}
	else
		ana_opts.id = NULL;

	if (score_aln)
	{
		analysis = true;
		header +="\t S-o-P";
		ana_opts.s_o_p = new double;
		*ana_opts.s_o_p = 0;
	}
	else
		ana_opts.s_o_p = NULL;


	//check if any modification option has been chosen
	bool modify = ((mod_options.to_upper) || (mod_options.to_lower) || (!extract_names.empty()) || (mod_options.column_trim !=-1) || (!mod_options.seq_trim.empty()));


	if ((analysis) && (!no_header))
		fprintf(analysis_F, "%s\n", header.c_str());

	// if extract names is a file replace with names in file
	if (extract_names.size()==1)
	{
		if (fs::is_regular_file(fs::path(extract_names[0])))
		{
			SequenceSet set;
			set.read(extract_names[0]);
			size_t n_names=set.n_seqs();
			extract_names.clear();
			extract_names.reserve(n_names);
			for (size_t i = 0; i<n_names; ++i)
				extract_names.push_back(set[i].name());
		}
	}


	// read and apply modifications do analysis for each alignment
	size_t aln_processed = 0;
	for (unsigned int i=0; i<aln_fs.size(); ++i)
	{
		if (!extract_names.empty())
			read_aln(aln_fs[i].string(), aln, extract_names, log_F, check);
		else
			read_aln(aln_fs[i].string(), aln, log_F, check);

		if (use_new_dir)
		{
			aln_out_f = out_f;
			aln_out_f += "/" + aln_fs[i].filename().string();
		}
		else if (!out_f.empty())
			aln_out_f = out_f;
		else
			aln_out_f = aln_fs[i].string();
		if (change_suffix)
		{
			pos =aln_out_f.rfind(".");
			if (pos != string::npos)
				aln_out_f.replace(pos+1, aln_out_f.size()-pos, new_suffix);
			else
				aln_out_f.append("."+change_suffix);
		}

		if (modify)
			modify_aln(aln, mod_options);

		if (analysis)
		{
			run_analysis(aln, ana_opts, avg_only, analysis_F);
			++aln_processed;
		}
		if (write_f)
		{
			try
			{
				aln.write(aln_out_f, out_format);
			}
			catch (My_IO_Exception &aln_e)
			{
				fprintf(log_F, "ERROR! Alignment could not be written to file %s: %s\n", aln.file().c_str(), aln_e.what());
				continue;
			}
			catch (Alignment_excep &aln_e)
			{
				fprintf(log_F, "ERROR! %s\n", aln_e.what());
				exit(EXIT_FAILURE);
			}
		}
	}


	if ((analysis) && (!no_avg))
	{
		fprintf(analysis_F, "AVG");
		if (ana_opts.aln_length != NULL)
			fprintf(analysis_F, "\t%li", *ana_opts.aln_length / aln_processed);
		if (ana_opts.id != NULL)
		{
			fprintf(analysis_F, "\t%.1f", *ana_opts.id / aln_processed);
			fprintf(analysis_F, "/%.1f", 100*exp(log(*ana_opts.identical)-log(*ana_opts.total)));
		}
		if (ana_opts.s_o_p != NULL)
			fprintf(analysis_F, "\t%.1f", *ana_opts.s_o_p / aln_processed);
		fprintf(analysis_F, "\n");
		if (aln_processed != aln_fs.size())
			fprintf(log_F, "%li of %li alignments processed!\n", aln_processed, aln_fs.size());
	}
	if (!log_f.empty())
		fclose(log_F);

	return EXIT_SUCCESS;

}


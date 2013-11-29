/*
 * setAnalyzer.cpp
 *
 *  Created on: Jun 14, 2012
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

#include "setAnalyser.h"

using namespace std;
using namespace BioTools::Seq;
using namespace BioTools::Utils;

namespace po = boost::program_options;
namespace fs = boost::filesystem;



typedef struct
{
	size_t *aln_length;
	double *seq_length;
	double *id;
	size_t *num_seqs;
	size_t *identical;
	size_t *total;
	size_t *composition;
	map<string, size_t> genome_len;
} Analysis_options;


size_t genome_size(const string &genome_f)
{
	size_t genome_length=0;
	FILE *genome_F =my_fopen(genome_f, "r");
	const int LINE_LENGTH = 1001;
	char line[LINE_LENGTH];
	char *pos;
	while (fgets(line, LINE_LENGTH, genome_F) != NULL)
	{
		if (line[0] != '>')
		{
			pos=&line[0];
			while (*pos != '\0')
			{
				if ((*pos != '\n') && (*pos != '\0'))
					++genome_length;
				++pos;
			}
		}
	}
	fclose(genome_F);

	return genome_length;
}



void
run_analysis(const Alignment &aln, Analysis_options &ana_opts, bool detailed, FILE *out_F)
{
	double tmp_d;
	size_t tmp;
	if (detailed)
		fprintf(out_F, "%s", aln.file().c_str());
	if (ana_opts.num_seqs != NULL)
	{
		tmp = aln.n_seqs();
		if (detailed)
			fprintf(out_F, "\t%li", tmp);
		*ana_opts.num_seqs += tmp;
	}
	if (ana_opts.aln_length != NULL)
	{
		tmp = aln.size();
		if (detailed)
			fprintf(out_F, "\t%li", tmp);
		*ana_opts.aln_length += tmp;
	}
	if (!ana_opts.genome_len.empty())
	{

	}
	if (ana_opts.seq_length != NULL)
	{
		tmp_d = aln.avg_size();
		if (detailed)
			fprintf(out_F, "\t%f", tmp_d);
		*ana_opts.seq_length += tmp_d;
	}
	if (ana_opts.id != NULL)
	{
		pair<size_t, size_t> id = identity(aln);
		tmp_d = 100*exp(log(id.first)-log(id.second));
		if (detailed)
			fprintf(out_F, "\t%.1f", tmp_d);
		*ana_opts.id += tmp_d;
		*ana_opts.identical += id.first;
		*ana_opts.total += id.second;
	}
	if (detailed)
		fprintf(out_F,"\n");
	if (ana_opts.composition)
		composition(aln, ana_opts.composition);
	if (!ana_opts.genome_len.empty())
		coverage(aln, ana_opts.genome_len);
}


void
run_analysis(const SequenceSet &seqSet, Analysis_options &ana_opts, bool detailed, FILE *out_F)
{
	if (detailed)
		fprintf(out_F, "%s", seqSet.file().c_str());
	size_t tmp;
	double tmp_d;
	if (ana_opts.num_seqs != NULL)
	{
		tmp = seqSet.n_seqs();
		if (detailed)
			fprintf(out_F, "\t%li", tmp);
		*ana_opts.num_seqs += tmp;
	}
	if (ana_opts.seq_length != NULL)
	{
		tmp_d = seqSet.avg_size();
		if (detailed)
			fprintf(out_F, "\t%f", tmp_d);
		*ana_opts.seq_length += tmp_d;
	}
	if (detailed)
		fprintf(out_F,"\n");
	if (ana_opts.composition)
		composition(seqSet, ana_opts.composition);
	if (!ana_opts.genome_len.empty())
		coverage(seqSet, ana_opts.genome_len);
}



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
				string name = itr->path().string();
				if (name[name.size()-1]=='~')
					continue;
				if (endings.empty())
				{
					files.push_back(name);
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
								files.push_back(name);
								break;
							}
						}
					}
				}
			}
		}

	}
}



int
main(int argc, char *argv[])
{
	vector<string> seqs_input_list, aln_input_list, genome_list;
	vector<string> extract_names;
	string out_f;
	vector<string> extensions;
	bool details, aln_length, seq_length;
	bool compute_identities;
	bool histogramm, num_seqs;
	bool extended, multiple;
	string diagramm_format;
	int n_threads;
	Analysis_options ana_opts;
	string log_f;
	bool composition;

	po::options_description general("General options");
	general.add_options()
		("help,h", "Produces this help message")
		("alignments,a", po::value<vector<string> >(&aln_input_list)->multitoken(), "Alignment file(s) or directory(s)")
		("sequences,s", po::value<vector<string> >(&seqs_input_list)->multitoken(), "Sequence file(s) or directory(s)")
		("genomes,g", po::value<vector<string> > (&genome_list)->multitoken(), "A list of files containing the genomes to calculate the coverage")
		("extensions,E", po::value<vector<string> >(&extensions)->multitoken(), "File extension to consider when given a directory. If none given all files will be used.")
		("output,o", po::value<string>(&out_f), "The prefix of the output files")
		("num_threads,n", po::value<int>(&n_threads), "Number of threads to use.")
	;

	po::options_description analysis_opts("Analysis options");
	analysis_opts.add_options()
		("details,d", po::value<bool> (&details)->default_value(false)->zero_tokens(), "A analysis is printed for each alignment")
		("identity,i", po::value<bool> (&compute_identities)->default_value(false)->zero_tokens(), "Computes the average identity of the alignments")
		("aln_length,l", po::value<bool>(&aln_length)->default_value(false)->zero_tokens(), "Computes the length of the alignment")
		("seq_length,S", po::value<bool>(&seq_length)->default_value(false)->zero_tokens(), "Computes the length of the sequences")
		("num_seqs,u", po::value<bool>(&num_seqs)->default_value(false)->zero_tokens(), "Number of sequences")
		("composition,c", po::value<bool>(&composition)->default_value(false)->zero_tokens(), "Calculates the composition of the sequences")
		("extended,x", po::value<bool>(&extended)->default_value(false)->zero_tokens(), "Use of extendend alphabet")
	;

	po::options_description graph_opts("Diagramm options");
	graph_opts.add_options()
		("diagramm,D", po::value<bool>(&histogramm)->default_value(false)->zero_tokens(), "Produces a graphic representation of the data")
		("multiple,m", po::value<bool>(&multiple)->default_value(false)->zero_tokens(), "All graphics will be written into the same field")
		("dia_format,f", po::value<string>(&diagramm_format)->default_value("pdf"), "The format to be used for the diagram")
	;

	po::options_description all("SetAnalyser v1.0.\n\nAllowed options are displayed below. For more information please see the Manual.");
	all.add(general).add(analysis_opts).add(graph_opts);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		cout << all<< "\n";
		return 1;
	}


	bool is_aln;
	vector<string> *in_list;
	if ((!seqs_input_list.empty()) && (!aln_input_list.empty()))
	{
		printf("ERROR! Alignments and Sequences cannot be analyzsed a the same time!");
		return EXIT_FAILURE;
	} else {
		is_aln = seqs_input_list.empty();
		if (is_aln)
			in_list=&aln_input_list;
		else
			in_list=&seqs_input_list;
	}

	char *script_path=getenv("BioTools++_DATA");
	if (script_path==NULL)
		script_path=getenv("HOME");
	string biotools_data=script_path;
	biotools_data.append("/.BioTools++/R/");

	vector<fs::path> input_fs;
	for (unsigned int i=0; i<in_list->size(); ++i)
	{
		if (fs::is_directory((*in_list)[i]))
			get_files_in_dir((*in_list)[i], input_fs, extensions, false);
		else
			input_fs.push_back((*in_list)[i]);
	}

	// get genome sizes if genome files given
	map<string, size_t> genome_sizes;
	size_t n_genomes=genome_list.size();
	for (size_t i=0; i<n_genomes; ++i)
	{
		genome_sizes[genome_list[i]] = genome_size(genome_list[i]);
		ana_opts.genome_len[genome_list[i]]=0;
	}


	//  analysis
	string header = "NAME";
//	if (input_fs.size()==1)
//		header+="\tseq_size";

	bool print_something = false;
	if (num_seqs)
	{
		print_something = true;
		header+="\tnum_seqs";
		ana_opts.num_seqs = new size_t;
		*ana_opts.num_seqs = 0;
	}
	else
		ana_opts.num_seqs = NULL;


	if (aln_length)
	{
		print_something = true;
		header+="\taln_length";
		ana_opts.aln_length = new size_t;
		*ana_opts.aln_length = 0;
	}
	else
		ana_opts.aln_length = NULL;

	if (seq_length)
	{
		print_something = true;
		header+="\tseq_length";
		ana_opts.seq_length = new double;
		*ana_opts.seq_length = 0;
	}
	else
		ana_opts.seq_length = NULL;

	if (compute_identities)
	{
		print_something = true;
		header +="\tidentity";
		ana_opts.id = new double;
		*ana_opts.id = 0;
		ana_opts.identical = new size_t;
		*ana_opts.identical = 0;
		ana_opts.total = new size_t;
		*ana_opts.total = 0;
	}
	else
		ana_opts.id = NULL;

	if (composition)
	{
		ana_opts.composition=new size_t[127];
		for (int i=0; i<127; ++i)
			ana_opts.composition[i]=0;
	}
	else
		ana_opts.composition=NULL;
	FILE *out_F = my_fopen(out_f, "w");
	if (print_something)
		fprintf(out_F, "==Summary==\n%s\n", header.c_str());

	FILE *log_F=stderr;

	// read and apply modifications do analysis for each alignment
	size_t aln_processed = 0;
	Alignment aln;

	char seq_type='x';
	if (!is_aln)
	{
		SequenceSet seqs;
		for (unsigned int i=0; i<input_fs.size(); ++i)
		{
			seqs.clear();
			try
			{
				seqs.read(input_fs[i].string());
				seq_type = seqs.seq_type();
			}
			catch (bad_alloc&)
			{
				fprintf(log_F, "ERROR: Not able to allocate enough memory.\n");
				continue;
			}
			catch (Alignment_excep &aln_e)
			{
				fprintf(log_F, "ERROR! Format could not be identified %s\n", input_fs[i].string().c_str());
				continue;
			}
			catch (My_IO_Exception &e)
			{
				fprintf(log_F, "ERROR! Sequence file %s could not be opened: %s\n", input_fs[i].string().c_str(), e.what());
				continue;
			}

			run_analysis(seqs, ana_opts, details, out_F);
			++aln_processed;
		}
	}
	else
	{
		for (unsigned int i=0; i<input_fs.size(); ++i)
		{
			aln.clear();
			try
			{
				aln.read(input_fs[i].string(), extract_names);
				seq_type = aln.seq_type();
			}
			catch (bad_alloc&)
			{
				fprintf(log_F, "ERROR: Not able to allocate enough memory.\n");
				continue;
			}
			catch (Alignment_excep &aln_e)
			{
				fprintf(log_F, "ERROR! Format could not be identified %s\n", input_fs[i].string().c_str());
				continue;
			}
			catch (My_IO_Exception &e)
			{
				fprintf(log_F, "ERROR! Alignment file %s could not be opened: %s\n", input_fs[i].string().c_str(), e.what());
				continue;
			}

			run_analysis(aln, ana_opts, details, out_F);
			++aln_processed;
		}
	}
	if (histogramm)
	{
		fclose(out_F);
		string r_call="R --no-restore --no-save --args ";
		string graph;
		if (multiple)
			graph="multiple";
		else
			graph="single";
		r_call.append(out_f + " " + out_f + " " + graph + " " +diagramm_format);
		r_call.append(" < " + biotools_data + "histogram.r >/dev/null 2>/dev/null");
		system(r_call.c_str());
		out_F = my_fopen(out_f, "a");
	}

	if (print_something)
		fprintf(out_F, "AVG");
	if (ana_opts.aln_length != NULL)
		fprintf(out_F, "\t%li", *ana_opts.aln_length / aln_processed);
	if (ana_opts.seq_length != NULL)
		fprintf(out_F, "\t%f", *ana_opts.seq_length / aln_processed);
	if (ana_opts.id != NULL)
	{
		fprintf(out_F, "\t%.1f", *ana_opts.id / aln_processed);
//		fprintf(out_F, "/%.1f", 100*exp(log(*ana_opts.identical)-log(*ana_opts.total)));
	}
	if (print_something)
		fprintf(out_F, "\n\n");

	if (composition)
	{
		size_t i;
		vector<int> char_use;
		if (seq_type=='P')
		{
			for (i=65; i<91; ++i)
			{
				if ((i != 66) && (i!=74) && (i!=88) && (i!=79) && (i!=90) && (i!=85))
					char_use.push_back(i);
			}
			if (extended)
			{
				char_use.push_back(66);
				char_use.push_back(90);
				char_use.push_back(88);
			}
		}
		else
		{
			char_use.push_back('A');
			char_use.push_back('C');
			char_use.push_back('G');
			if (char_use['U']!=0)
				char_use.push_back('U');
			else
				char_use.push_back('T');
			if (extended)
			{
				char_use.push_back('R');
				char_use.push_back('Y');
				char_use.push_back('M');
				char_use.push_back('K');
				char_use.push_back('W');
				char_use.push_back('S');
				char_use.push_back('B');
				char_use.push_back('D');
				char_use.push_back('H');
				char_use.push_back('V');
				char_use.push_back('N');
			}
		}
		size_t total=0;
		size_t n_chars = char_use.size();
		for (i=0; i<n_chars; ++i)
			total+=ana_opts.composition[i];
		fprintf(out_F, "==Composition==\nalphabet");
		for (i=0; i<n_chars; ++i)
			fprintf(out_F, " %c", char_use[i]);
		fprintf(out_F, "\ncount ");
		for (i=0; i<n_chars; ++i)
			fprintf(out_F, " %li", ana_opts.composition[char_use[i]]);
		fprintf(out_F, "\n\n");

	}

	if (!genome_sizes.empty())
	{
		map<string, size_t>::iterator it, it_end = genome_sizes.end();
		fprintf(out_F, "==Coverage==\n");
		for (it=genome_sizes.begin(); it!=it_end; ++it)
		{
			fprintf(out_F, "%s %f\n", it->first.c_str(), ana_opts.genome_len[it->first]*100.0/it->second);
		}
	}


	fclose(out_F);

	//printf("%li of %li alignments processed!\n", aln_processed, input_fs.size());

	if (aln_length)
		delete ana_opts.aln_length;

	if (seq_length)
		delete ana_opts.seq_length;

	if (compute_identities)
	{
		delete ana_opts.id;
		delete ana_opts.identical;
		delete ana_opts.total;
	}

	if (num_seqs)
		delete ana_opts.num_seqs;

	if (composition)
		delete[] ana_opts.composition;

	if (!log_f.empty())
		fclose(log_F);

	return EXIT_SUCCESS;
}






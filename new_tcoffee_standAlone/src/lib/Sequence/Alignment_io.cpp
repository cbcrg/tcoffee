/*
 * Alignment_io.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: ck
 */


#include "Alignment.h"
using namespace std;
using namespace BioTools::Utils;


namespace BioTools
{
namespace Seq
{


bool
Alignment::_len_check()
{
	size_t n_seqs = this->n_seqs();
	for (size_t i=0; i<n_seqs; ++i)
	{
		if ((*this)[i].size() != _aln_length)
			return false;
	}
	return true;
}

void
Alignment::read(const std::string &seq_f, const vector<string> &seq_names, bool check, short format)
{
	SequenceSet::read(seq_f, seq_names, check, format);
	_aln_length=(*this)[0].size();
	//if (!seq_names.empty())
	_delete_gap_columns();

	if (check)
	{
		if (!_len_check())
			throw Alignment_excep("Sequences have not the same length");
	}
}

void
Alignment::write(const string &aln_f, const string format) const
{
	string format_lc;
	format_lc.reserve(format.size());
	for (size_t i= 0; i<format.size(); ++i)
		format_lc.push_back(tolower(format[i]));

	FILE *aln_F;
	if ((format_lc == "fasta") || (format_lc == "fasta_seq") || (format_lc == "clustalw") || (format_lc == "msf") || (format_lc == "phylip_s") || (format_lc == "phylip_i"))
	{
		if (aln_f.empty())
			aln_F=stdout;
		else
			aln_F = my_fopen(aln_f.c_str(), "w");
	}
	else
		throw(Alignment_excep("Unknown alignment format!"));

	if (format_lc == "fasta")
		_write_fasta(aln_F);
	else if (format_lc == "fasta_seq")
		_write_fasta_seq(aln_F);
	else if (format_lc == "clustalw")
		_write_clustalw(aln_F);
	else if (format_lc == "msf")
		_write_msf(aln_F);
	else if (format_lc == "phylip_s")
		_write_phylip_sequential(aln_F);
	else if (format_lc == "phylip_i")
		_write_phylip_interleaved(aln_F);

	if (!aln_f.empty())
		fclose(aln_F);
}



void
Alignment::_write_fasta(FILE *aln_F, unsigned int line_break) const
{
	unsigned int num_seqs = this->n_seqs();
	size_t current_length;
	const char *tmp_seq;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		const Sequence &seq = (*this)[i];
		if (seq.comment().empty())
			fprintf(aln_F, ">%s\n", seq.name().c_str());
		else
			fprintf(aln_F, ">%s %s\n", seq.name().c_str(), seq.comment().c_str());
		current_length = 0;
		tmp_seq = seq.sequence().c_str();
		while (current_length < _aln_length)
		{
			fprintf(aln_F, "%.*s\n", line_break, &tmp_seq[current_length]);
			current_length += line_break;
		}
	}
}

void
Alignment::_write_fasta_seq(FILE *aln_F, unsigned int line_break) const
{
	unsigned int num_seqs = this->n_seqs();
	size_t current_length, seq_length;
	string tmp_s;
	const char *tmp_seq;
	size_t j, pos;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		const Sequence &seq = (*this)[i];
		if (seq.comment().empty())
			fprintf(aln_F, ">%s\n", seq.name().c_str());
		else
			fprintf(aln_F, ">%s %s\n", seq.name().c_str(), seq.comment().c_str());
		current_length = 0;
		seq_length = seq.size();
		tmp_s=seq.sequence();
		pos=0; 
		for (j=0; j<seq_length; ++j)
		{
			if (tmp_s[j] != '-')
				tmp_s[pos++]=tmp_s[j];
		}
		tmp_s.resize(pos);
		seq_length = tmp_s.size();
		tmp_seq = tmp_s.c_str();
		
		while (current_length < seq_length)
		{
			fprintf(aln_F, "%.*s\n", line_break, &tmp_seq[current_length]);
			current_length += line_break;
		}
	}
}


void
Alignment::_write_clustalw(FILE *aln_F, unsigned int line_break) const
{
	vector<unsigned int> vec_col[26];
	vector<unsigned int> vec_dot[26];

	if(this->seq_type()=='P')
	{
		vec_col[0].push_back(0);	// A
		vec_col[1].push_back(13);	// B
		vec_col[2].push_back(11);	// C
		vec_col[3].push_back(3);	// D
		vec_col[4].push_back(1);	// E
		vec_col[4].push_back(3);
		vec_col[5].push_back(6);	// F
		vec_col[5].push_back(8);
		vec_col[6].push_back(9);	// G
		vec_col[7].push_back(2);	// H
		vec_col[7].push_back(4);
		vec_col[7].push_back(7);
		vec_col[8].push_back(5);	// I
		vec_col[8].push_back(6);
		vec_col[9].push_back(13);	// J
		vec_col[10].push_back(1);	// K
		vec_col[10].push_back(2);
		vec_col[10].push_back(4);
		vec_col[11].push_back(5);	// L
		vec_col[11].push_back(6);
		vec_col[12].push_back(5);	// M
		vec_col[12].push_back(6);
		vec_col[13].push_back(1);	// N
		vec_col[13].push_back(2);
		vec_col[13].push_back(3);
		vec_col[14].push_back(13);	// O
		vec_col[15].push_back(10);	// P
		vec_col[16].push_back(1);	// Q
		vec_col[16].push_back(2);
		vec_col[16].push_back(3);
		vec_col[16].push_back(4);
		vec_col[17].push_back(4);	// R
		vec_col[18].push_back(0);	// S
		vec_col[19].push_back(0);	// T
		vec_col[20].push_back(5);	// U
		vec_col[21].push_back(13);	// V
		vec_col[22].push_back(8);	// W
		vec_col[23].push_back(13);	// X
		vec_col[24].push_back(7);	// Y
		vec_col[24].push_back(8);	// Y
		vec_col[25].push_back(11);	// Z


		vec_dot[0].push_back(0);	// A
		vec_dot[0].push_back(1);
		vec_dot[0].push_back(2);
		vec_dot[0].push_back(4);
		vec_dot[1].push_back(12);	// B
		vec_dot[2].push_back(0);	// C
		vec_dot[3].push_back(5);	// D
		vec_dot[3].push_back(6);
		vec_dot[3].push_back(7);
		vec_dot[4].push_back(6);	// E
		vec_dot[4].push_back(7);
		vec_dot[4].push_back(8);
		vec_dot[5].push_back(9);	// F
		vec_dot[5].push_back(10);
		vec_dot[6].push_back(2);	// G
		vec_dot[6].push_back(5);
		vec_dot[7].push_back(7);	// H
		vec_dot[7].push_back(8);
		vec_dot[7].push_back(10);
		vec_dot[8].push_back(9);	// I
		vec_dot[9].push_back(12);	// J
		vec_dot[10].push_back(3);	// K
		vec_dot[10].push_back(6);
		vec_dot[10].push_back(7);
		vec_dot[10].push_back(8);
		vec_dot[11].push_back(9);	// L
		vec_dot[12].push_back(9);	// M
		vec_dot[13].push_back(3);	// N
		vec_dot[13].push_back(5);
		vec_dot[13].push_back(6);
		vec_dot[13].push_back(7);
		vec_dot[13].push_back(8);
		vec_dot[14].push_back(12);	// O
		vec_dot[15].push_back(4);	// P
		vec_dot[16].push_back(6);	// Q
		vec_dot[16].push_back(7);
		vec_dot[16].push_back(8);
		vec_dot[17].push_back(8);	// R
		vec_dot[18].push_back(0);	// S
		vec_dot[18].push_back(2);
		vec_dot[18].push_back(3);
		vec_dot[18].push_back(4);
		vec_dot[18].push_back(5);
		vec_dot[18].push_back(6);
		vec_dot[19].push_back(1);	// T
		vec_dot[19].push_back(3);
		vec_dot[19].push_back(4);
		vec_dot[20].push_back(12);	// U
		vec_dot[21].push_back(1);	// V
		vec_dot[21].push_back(9);
		vec_dot[22].push_back(12);	// W
		vec_dot[23].push_back(12);	// X
		vec_dot[24].push_back(10);	// Y
		vec_dot[25].push_back(12);	// Z
	}

	fprintf(aln_F, "CLUSTAL W (1.82) multiple sequence alignment\n\n\n");
	unsigned int aln_length = size();
	unsigned int num_seqs = this->n_seqs();
	unsigned int pos = 0;
	unsigned int end = 0, i, j;
	const char *tmp_seq;
	Fixed_multi_array<unsigned int> col_groups(line_break, 14);
	Fixed_multi_array<unsigned int> dot_groups(line_break, 13);
	short c;
	size_t x;
	size_t *positions = new size_t[num_seqs];
	for (i=0; i<num_seqs; ++i)
		positions[i] = 0;
	char *id = new char[line_break];
	while (pos < aln_length)
	{
		end += line_break;
		if (end > aln_length)
		{
			line_break = aln_length-pos;
			end=aln_length;
		}
		dot_groups.set_all(0);
		col_groups.set_all(0);

		for (i = 0; i < num_seqs; ++i)
		{
			tmp_seq = (*this)[i].sequence().c_str()+pos;
			if (i==0)
			{
				for (j = 0; j < line_break; ++j)
					id[j] = tmp_seq[j];
			}
			else
			{
				for (j = 0; j < line_break; ++j)
					if (id[j] != tmp_seq[j])
						id[j] = '-';
			}

			for (j = 0; j < line_break; ++j)
			{
				c = tolower(tmp_seq[j])-97;
				if (c >= 0)
				{
					for (x=0; x < vec_col[c].size(); ++x)
						++col_groups[j][vec_col[c][x]];
					for (x=0; x < vec_dot[c].size(); ++x)
						++dot_groups[j][vec_dot[c][x]];
				}
			}
			fprintf(aln_F, "%-10.10s      ", (*this)[i].name().c_str());//, tmp_seq, end");
			for (x=0; x<line_break; ++x)
			{
				if (tmp_seq[x] != '-')
					++positions[i];
				fprintf(aln_F, "%c", tmp_seq[x]);
			}
			fprintf(aln_F, " %li\n", positions[i]);
		}
		fprintf(aln_F, "                ");

		bool printed = false;
		for (j = 0; j < line_break; ++j)
		{
			if (id[j] != '-')
			{
				fprintf(aln_F, "*");
				continue;
			}
			printed=false;
			for (i=0; i < 14; ++i)
			{
				if (col_groups[j][i] == num_seqs)
				{
					fprintf(aln_F, ":");
					printed = true;
					break;
				}
			}
			if (!printed)
			{
				for (i=0; i<13; ++i)
				{
					if (dot_groups[j][i] == num_seqs)
					{
						fprintf(aln_F, ".");
						printed = true;
						break;
					}
				}
			}
			if (!printed)
				fprintf(aln_F, " ");
		}
		fprintf(aln_F, "\n\n");
		pos += line_break;
	}
	delete[] id;
	delete[] positions;
}


void
Alignment::_write_msf(FILE *aln_F) const
{
	unsigned int aln_length = this->size();
	unsigned int num_seqs = this->n_seqs();
	fprintf(aln_F, "PileUp\n\n MSF: %i Type: %c Check:  0  ..\n\n", aln_length, this->seq_type());

	for (unsigned int i = 0; i < num_seqs; ++i)
		fprintf(aln_F, " Name: %s Len: %i, Check: 0, Weight: 1\n", (*this)[i].name().c_str(), aln_length);

	fprintf(aln_F, "\n//\n\n");

	unsigned int pos = 0;
	unsigned int num = 0;
	unsigned int current_length = 0;
	const char *tmp_seq = NULL;

	unsigned int maxNameLength=0;
	for (unsigned int i = 0; i < num_seqs; ++i)
	{
		if ((*this)[i].name().size()>maxNameLength)
			maxNameLength=(*this)[i].name().size();
	}

	while (pos < aln_length)
	{
		for (unsigned int i = 0; i < num_seqs; ++i)
		{
			tmp_seq = (*this)[i].sequence().c_str()+pos;
			num = 0;
			current_length = pos;
			fprintf(aln_F, "%-*s  ", maxNameLength, (*this)[i].name().c_str());
			while ((current_length < aln_length) && (num < 5))
			{
				fprintf(aln_F, " %.10s", tmp_seq);
				tmp_seq += 10;
				++num;
				current_length += 10;
			}
			fprintf(aln_F ,"\n");
		}
		fprintf(aln_F, "\n");
		pos = current_length;
	}
}


void
Alignment::_write_phylip_sequential(FILE *aln_F) const
{
	size_t num_seqs = this->n_seqs();
	fprintf(aln_F, "%li %li\n", num_seqs, _aln_length);
	const char *tmp_seq;
	size_t j;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		tmp_seq = (*this)[i].sequence().c_str();
		fprintf(aln_F, "%-10.10s", (*this)[i].name().c_str());
		j=0;
		while (j < _aln_length)
		{
			fprintf(aln_F, "%.10s", &tmp_seq[j]);
			j+=10;
		}
		fprintf(aln_F, "\n");
	}
}

void
Alignment::_write_phylip_interleaved(FILE *aln_F) const
{
	size_t num_seqs = this->n_seqs();
	fprintf(aln_F, "%li %li\n", num_seqs, _aln_length);
	const char *tmp_seq;
	size_t j=0;
	size_t line_break = 60;
	size_t write_length =line_break-10 < _aln_length ? line_break-10 : _aln_length;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		tmp_seq = (*this)[i].sequence().c_str();
		fprintf(aln_F, "%-10.10s ", (*this)[i].name().c_str());
		j=0;
		while (j < write_length)
		{
			fprintf(aln_F, "%.10s ", &tmp_seq[j]);
			j+=10;
		}
		fprintf(aln_F, "\n");
	}
	fprintf(aln_F, "\n");

	write_length = line_break / 10;
	size_t k, pos;

	while (j < _aln_length)
	{
		for (size_t i = 0; i < num_seqs; ++i)
		{
			pos = j;
			k=0;
			tmp_seq = (*this)[i].sequence().c_str();
			while ((k < write_length) && (pos < _aln_length))
			{
				fprintf(aln_F, "%.10s ", &tmp_seq[pos]);
				pos+=10;
				++k;
			}
			fprintf(aln_F, "\n");
		}
		j+=write_length * 10;
		fprintf(aln_F, "\n");
	}
}


}
}

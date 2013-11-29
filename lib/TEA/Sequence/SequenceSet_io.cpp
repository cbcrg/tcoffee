/*
 * SequenceSet_io.cpp
 *
 *  Created on: Jun 26, 2012
 *      Author: Carsten Kemena
 *
 *   Copyright 2012 Carsten Kemena
 *
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


#include "SequenceSet.h"
#include "SequenceSet_io.h"



namespace BioTools
{
namespace Seq
{

using namespace std;
using namespace BioTools::Utils;
// Read Functions

void
SequenceSet::_gap_replace()
{
	size_t num_seqs = this->n_seqs();
	size_t n_cols;
	size_t j;

	for (size_t i=0; i<num_seqs; ++i)
	{
		Sequence &tmp_seq = *_seqs[i];
		n_cols = tmp_seq.size();
		for (j=0; j<n_cols; ++j)
		{
			if (!isalpha(tmp_seq[j]))
					tmp_seq[j]='-';
		}
	}
}

short
SequenceSet::_identify_aln_format(FILE *aln_F)
{
	const unsigned int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	char *tmp, *pos;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '>')
		{
			fgets(line, LINE_LENGTH, aln_F);
			if (line[0] == '>')
				return 6; // AMPS
			else
				return 1; // Fasta
		}
		if (strstr(line, "CLUSTAL ") != NULL)
			return 2; // Clustal
		if (strstr(line, "MSF:") != NULL)
			return 3; // MFS
		if (strstr(line, "# STOCKHOLM 1.0") != NULL)
			return 4; // Stockholm
		if (strstr(line, "ENTRY") != NULL )
			return 5; // Codata

		//check Phylip
		StrTok strTok(line);
		tmp = strTok.next(" \n\t");
		if (tmp == NULL)
			continue;
		pos=tmp;
		while (*pos != '\0')
			if (!isdigit(*(pos++)))
				continue;
		tmp = strTok.next(" \n\t");
		if (tmp == NULL)
			continue;
		pos=tmp;
		while (*pos != '\0')
			if (!isdigit(*(pos++)))
				continue;
		tmp = strTok.next(" \n\t");
		if (tmp == NULL)
			return 7; // Phylip
	}
	return 0;
}


void
SequenceSet::read(const std::string &seq_f, const vector<string> &seq_names, bool check, short format)
{
	FILE *aln_F = my_fopen(seq_f.c_str(), "r");

	_file = seq_f;
	if (format <0)
		format = _identify_aln_format(aln_F);
	map<string, short> extract_only;
	size_t n_extract_seqs = seq_names.size();
	for (unsigned int i = 0; i < n_extract_seqs; ++i)
		extract_only.insert(pair<string,short>(seq_names[i],1));

	fseek ( aln_F , 0 , SEEK_SET );
	switch (format)
	{
		case 1:
			_read_fasta_f(aln_F, extract_only);
			break;
		case 2:
			_read_clustalw_f(aln_F, extract_only);
			break;
		case 3:
			_read_msf_f(aln_F, extract_only);
			break;
		case 4:
			_read_stockholm_f(aln_F, extract_only);
			break;
		case 5:
			_read_codata_f(aln_F, extract_only);
			break;
		case 6:
			_read_amps_f(aln_F, extract_only);
			break;
		case 7:
			_read_phylip_f(aln_F, extract_only);
			break;
		default:
			throw Alignment_excep("Format could not be identified");
	}
	fclose(aln_F);

	if ((check) && (!check_set(*this)))
		throw Alignment_excep("Sequence contains bad character");
	this->seq_type(identify_seq_type(*(this->_seqs)[0]));
}


void
SequenceSet::_read_fasta_f(FILE *aln_F, const map<string, short> &seq_names)
{
	const unsigned int LINE_LENGTH = 500;
	char line[LINE_LENGTH];
	Sequence *tmp_seq = NULL;
	char *comment = NULL, *name = NULL;
	bool read_sequence = true;
	size_t id=0;
	size_t seq_length=0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '>')
		{
			if (this->n_seqs())
				seq_length = tmp_seq->size();
			StrTok tokenizer(&line[1]);
			name = tokenizer.next(" \n");
			comment = tokenizer.next("\n");
			if (seq_names.empty() || (seq_names.count(name)>0))
			{
				read_sequence = true;
				if (comment == NULL)
					tmp_seq = new Sequence(name, "", seq_length, id++);
				else
					tmp_seq = new Sequence(name, comment, seq_length, id++);
				_seqs.push_back(Seq_ptr(tmp_seq));
			}
			else
				read_sequence = false;
		}
		else if (line[0] == '/')
			break;
		else
		{
			if (read_sequence)
			{
				tmp_seq->append(line);
				if ((*tmp_seq)[tmp_seq->size()-1] == '\n')
					tmp_seq->resize(tmp_seq->size()-1);
			}
		}
	}
}


void
SequenceSet::_read_clustalw_f(FILE *aln_F, const map<string, short> &seq_names)
{
	const unsigned int LINE_LENGTH = 200;
	char line[LINE_LENGTH];
	fgets(line, LINE_LENGTH, aln_F);
	long t_pos = 0;
	size_t id = 0;
	while ((fgets(line, LINE_LENGTH, aln_F) != NULL) && (line[0] == '\n'))
		t_pos = ftell(aln_F);
	fseek(aln_F, t_pos, SEEK_SET);

	//read first block in format
	Sequence *seq_p = NULL;
	char *name = NULL, *tmp_seq = NULL;
	vector<short> read_sequence;
	StrTok tokenizer;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if ((line[0] == ' ') || (line[0] == '\n'))
			break;
		else
		{
			tokenizer.set(line);
			name = tokenizer.next(" \n");
			if (seq_names.empty() || seq_names.count(name) > 0)
			{
				read_sequence.push_back(1);
				tmp_seq = tokenizer.next(" \n");
				seq_p = new Sequence(name, "", tmp_seq, id++);
				_seqs.push_back(Seq_ptr(seq_p));
			}
			else
				read_sequence.push_back(0);
		}
	}

	// Read rest of the alignment
	unsigned int counter = 0;
	unsigned int seq_id = 0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if ((line[0] == ' ') || (line[0] == '\n'))
		{
			seq_id = 0;
			counter = 0;
		}
		else
		{
			if (read_sequence[counter])
			{
				tokenizer.set(line);
				tokenizer.next(" \n");
				tmp_seq = tokenizer.next(" \n");
				_seqs[seq_id++]->append(tmp_seq);
			}
			++counter;
		}
	}
}


void
SequenceSet::_read_msf_f(FILE *aln_F, const map<string, short> &seq_names)
{
	char *msf;
	const unsigned int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	while ((fgets(line, LINE_LENGTH, aln_F) != NULL) && ((msf=strstr(line, "MSF:")) == NULL));
	msf +=4;
	size_t seq_length=atoi(msf);
	char *seq_name;
	Sequence *tmp_seq;
	StrTok tokenizer;
	vector<short> read_sequence;
	size_t id = 0;

	//read comment block
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if ((line[0] == '/') && (line[1] == '/'))
			break;
		if (((seq_name=strstr(line, "Name:")) != NULL) || ((seq_name=strstr(line, "NAME:")) != NULL))
		{
			tokenizer.set(seq_name+5);
			seq_name = tokenizer.next(" \n");
			if (seq_names.empty() || (seq_names.count(seq_name)>0))
			{
				read_sequence.push_back(1);
				tmp_seq = new Sequence(seq_name, "", seq_length, id++);
				_seqs.push_back(Seq_ptr(tmp_seq));

			}
			else
				read_sequence.push_back(0);
		}
	}

	//read alignment
	char *pos;
	char c;
	unsigned int counter = 0;
	unsigned int seq_id = 0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '\n')
		{
			counter = 0;
			seq_id = 0;
		}
		else
		{
			if (read_sequence[counter]==1)
			{
				pos = line;
				while (*(++pos) != ' ');
				while (1)
				{
					while ((c=*pos) != '\0')
					{
						if ((c=='-') || (c=='.') || (isalpha(c)))
							_seqs[seq_id]->append(c);
						++pos;
					}
					if (*(pos-1)=='\n')
						break;
					else
					{
						fgets(line, LINE_LENGTH, aln_F);
						pos = &line[0];
					}
				}
				++seq_id;
			}
			else
			{
				while (1)
				{
					pos = line;
					while (*pos != '\0')
						++pos;
					if (*(pos-1)=='\n')
						break;
					fgets(line, LINE_LENGTH, aln_F);
				}
			}
			++counter;
		}
	}
	_gap_replace();
}


void
SequenceSet::_read_stockholm_f(FILE *aln_F, const map<string, short> &seq_names)
{
	const unsigned int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	char *part1_p = NULL, *part2_p =NULL;
	StrTok tokenizer;
	Sequence *seq_p = NULL;
	char *end=&line[LINE_LENGTH-1];
	bool use = false;
	bool append = false;
	size_t seq_id = 0;
	vector<short> read_seq;
	*end=8;
	size_t id = 0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '\n')
		{
			append = true;
			seq_id = 0;
		}
		else if (line[0] == '#')
		{
			while ((*end != 8) && (*end != '\n'))
			{
				*end = 8;
				fgets(line, LINE_LENGTH, aln_F);
			}
		}
		else if (line[0] != '/')
		{
			tokenizer.set(line);
			part1_p = tokenizer.next(" \n");
			part2_p = tokenizer.next(" \n");
			if (part2_p != NULL)
			{
				if (seq_names.empty() || (seq_names.count(part1_p)))
				{
					use=true;
					if (append)
					{
						seq_p = &(*_seqs[seq_id]);
						seq_p->append(part2_p);
					} else {
						seq_p = new Sequence(part1_p, "", part2_p, id++);
						_seqs.push_back(Seq_ptr(seq_p));
					}
					++seq_id;
				}
				else
					use=false;
			}
			else
			{
				if (use)
					seq_p->append(part1_p);
			}
		} else
			use=false;
		*end=8;
	}

	//Turn "." into "-"
	_gap_replace();
}


void
SequenceSet::_read_codata_f(FILE *aln_F, const map<string, short> &seq_names)
{
	const unsigned int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	char *name=NULL;
	size_t seq_length = 0;
	StrTok tokenizer;
	Sequence *seq_p = NULL;
	size_t id = 0;
	char c;
	unsigned int pos1 = 0, pos2 = 1;
	bool extract_seq = false;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (!strncmp(line, "ENTRY", 5))
		{
			tokenizer.set(&line[6]);

			name = tokenizer.next(" \n");
		}
		else if (!strncmp(line, "SEQUENCE", 8))
		{
			if (seq_names.empty() || seq_names.count(name))
			{
				seq_p = new Sequence(name, "", seq_length, id++);
				_seqs.push_back(Seq_ptr(seq_p));
				fgets(line, LINE_LENGTH, aln_F);
				extract_seq = true;
			}
			else
				extract_seq = false;
		}
		else if (strncmp(line, "///", 3) == 0)
		{
			if (extract_seq)
				seq_length = seq_p->size();
		}
		else
		{
			if (extract_seq)
			{
				pos1 = 0;
				pos2 = 0;
				while ((c=line[pos1++]) != '\0')
				{
					if (isalpha(c) || (c == '-'))
						line[pos2++] = c;
				}
				line[pos2] = '\0';
				seq_p->append(line);
			}
		}
	}
}


void
SequenceSet::_read_amps_f(FILE *aln_F, const map<string, short> &seq_names)
{
	const unsigned int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	vector<short> extract_seq;
	StrTok tokenizer;
	Sequence *seq_p;
	size_t id = 0;
	char *name = NULL;
	unsigned int num_seqs = 0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '*')
			break;
		else if (line[0] == '>')
		{
			++num_seqs;
			tokenizer.set(&line[1]);
			name = tokenizer.next("\n");
			if (seq_names.empty() || seq_names.count(name))
			{
				seq_p = new Sequence(name, "", 0, id++);
				_seqs.push_back(Seq_ptr(seq_p));
				extract_seq.push_back(1);
			}
			else
				extract_seq.push_back(0);
		}
	}
	char c;
	char *seq_line = new char[num_seqs+2];
	unsigned int n_seqs2 = num_seqs+2;
	unsigned int i;
	unsigned int seq_id = 0;
	while (fgets(seq_line, n_seqs2, aln_F) != NULL)
	{
		if (seq_line[0] == '*')
			break;
		else
		{
			seq_id = 0;
			for (i=0; i<num_seqs; ++i)
			{
				if (extract_seq[i])
				{
					if ((c = seq_line[i]) == ' ')
						_seqs[seq_id]->append('-');
					else
						_seqs[seq_id]->append(c);
					++seq_id;
				}
			}
		}
	}
	delete[] seq_line;
}


// reads interleaved and sequential Phylip format
void
SequenceSet::_read_phylip_f(FILE *aln_F, const std::map<std::string, short> &seq_names)
{
	const unsigned int LINE_LENGTH = 500;
	char line[LINE_LENGTH];
	vector<int> use_seq;
	size_t num_seqs, seq_length;
	fgets(line, LINE_LENGTH, aln_F);
	sscanf(line, "%lu %lu", &num_seqs, &seq_length);
	char name[11];
	name[10] = '\0';
	char *pos;
	Sequence *seq_p;
	size_t seq_num = num_seqs;
	size_t id = 0;
	char c;
	size_t line_num = 0;
	int x =-1;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '\n')
			continue;

		//Read sequence name
		if (line_num < num_seqs)
		{
			strncpy(name, line, 10);
			pos=name;
			while ((*pos != '\0') && (*pos != ' '))
				++pos;
			*pos = '\0';
			if ((seq_names.empty()) || (seq_names.count(name)))
			{
				use_seq.push_back(++x);
				seq_p = new Sequence(name, "", seq_length, id++);
				_seqs.push_back(Seq_ptr(seq_p));
			}
			else
			{
				use_seq.push_back(-1);
			}
			pos=&line[10];
			++line_num;
		}
		else
			pos=&line[0];

		//read_sequence
		++seq_num;
		if (seq_num >= num_seqs)
			seq_num=0;
		if (use_seq[seq_num] >= 0)
		{

			seq_p=&(*_seqs[use_seq[seq_num]]);
			while (1)
			{
				while ((c=*pos) != '\0')
				{
					if ((isalpha(c)) || (c=='-'))
						seq_p->append(c);
					++pos;
				}

				if (*(pos-1) != '\n')
				{
					if (fgets(line, LINE_LENGTH, aln_F) == NULL)
						break;
					pos=&line[0];
				}
				else
					break;
			}
		}
		else
		{
			while (line[strlen(line)-1] != '\n')
			{
				if (fgets(line, LINE_LENGTH, aln_F) == NULL)
					break;
			}
		}

	}
}



// Write functions

void
SequenceSet::write(const string &aln_f, const string format) const
{
	string format_lc;
	format_lc.reserve(format.size());
	for (size_t i= 0; i<format.size(); ++i)
		format_lc.push_back(tolower(format[i]));

	FILE *aln_F;
	if ((format_lc == "fasta"))
		aln_F = my_fopen(aln_f.c_str(), "w");
	else
		throw(Alignment_excep("Unknown alignment format!"));

	if (format_lc == "fasta")
		_write_fasta(aln_F);

	fclose(aln_F);
}


void
SequenceSet::_write_fasta(FILE *aln_F, unsigned int line_break) const
{
	unsigned int num_seqs = this->n_seqs();
	size_t current_length, seq_length;
	string tmp_seq;
// 	char symbol='\0';

	for (unsigned int i = 0; i < num_seqs; ++i)
	{
		if (_seqs[i]->comment().empty())
			fprintf(aln_F, ">%s\n", _seqs[i]->name().c_str());
		else
			fprintf(aln_F, ">%s %s\n", _seqs[i]->name().c_str(), _seqs[i]->comment().c_str());
		current_length = 0;
		tmp_seq = _seqs[i]->sequence();
		seq_length=_seqs[i]->size();
		while (current_length < seq_length)
		{
			fprintf(aln_F, "%.*s\n", line_break, &tmp_seq.c_str()[current_length]);
			current_length += line_break;
		}

	}
}




} // end namespace Sequence
} // end namespace Biotools








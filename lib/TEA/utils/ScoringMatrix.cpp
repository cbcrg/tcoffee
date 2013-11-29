/*
 * ScoringMatrix.cpp
 *
 *  Created on: Dec 4, 2011
 *      Author: ck
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

#include "ScoringMatrix.h"

namespace BioTools {
namespace Utils {

using namespace std;


Scoring_Matrix::Scoring_Matrix():_size(27),_matrix_name(""),_matrix(new double*[_size])
{
	_matrix=new double*[_size];
	for (int i = 0; i < _size;++i)
		_matrix[i] = new double[_size];
}

Scoring_Matrix::~Scoring_Matrix() {
	for (int i = 0; i < _size;++i)
		delete[] _matrix[i];
	delete[] _matrix;
}

Scoring_Matrix::Scoring_Matrix(const Scoring_Matrix &mat):_size(mat._size),_matrix_name(mat._matrix_name),_matrix(new double*[_size])
{
	for (int i = 0; i < _size;++i)
		_matrix[i] = new double[_size];
	for (int i = 0; i < _size;++i)
		memcpy(_matrix[i], mat._matrix[i], sizeof(double)*_size );
}

Scoring_Matrix&
Scoring_Matrix:: operator=(const Scoring_Matrix &mat)
{
	int i;
	for (i=0; i<_size; ++i)
		delete[] _matrix;
	delete _matrix;
	_matrix = new double*[_size];
	for (i=0; i<_size; ++i)
		_matrix[i] = new double[_size];
	for (i=0; i<_size; ++i)
		memcpy(_matrix[i], mat._matrix[i], sizeof(double)*_size );
	return *this;
}

void
Scoring_Matrix::read(const std::string &matrix_f)
{
	FILE *matrix_F = my_fopen(matrix_f.c_str(), "r");
	const unsigned int LINE_LENGTH = 500;
	char line[LINE_LENGTH];
	StrTok tokenizer;
	char *tmp;
	int order[256];
	int pos,step;
	while (fgets(line, LINE_LENGTH, matrix_F) != NULL)
	{
		if ((line[0] == '#') || (line[0] == '\n'))
			continue;
		if (line[0] == ' ')
		{
			tokenizer.set(line);
			pos = 0;
			while ((tmp = tokenizer.next(" \n\t")) != NULL)
			{
				step = std::tolower(tmp[0])-97;
				if (step < 0)
					order[pos++]=26;
				else
					order[pos++]=step;

			}
			continue;
		}
		tokenizer.set(line);
		tmp=tokenizer.next(" \n\t");
		pos = std::tolower(tmp[0])-97;
		if (pos<0)
			pos=26;
		step=0;
		while ((tmp = tokenizer.next(" \n\t")) != NULL)
			_matrix[pos][order[step++]] = atof(tmp);
	}
	fclose(matrix_F);



}


} /* namespace Utils */
} /* namespace BioTools */

/*
 * ScoringMatrix.h
 *
 *  Created on: Dec 4, 2011
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

#ifndef SCORINGMATRIX_H_
#define SCORINGMATRIX_H_


// C headers
#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <cstring>

// C++ headers
#include <string>

// My headers
#include "filesystem.h"
#include "utils.h"


namespace BioTools {
namespace Utils {


class Scoring_Matrix {
private:
	int _size;
	std::string _matrix_name;
	double **_matrix;

public:
	Scoring_Matrix();
	virtual ~Scoring_Matrix();
	Scoring_Matrix(const Scoring_Matrix &mat);

	Scoring_Matrix& operator=(const Scoring_Matrix &mat);


	double *operator[](unsigned int index)
	{
		return _matrix[index];
	}
	const double *operator[](unsigned int index) const
	{
		return _matrix[index];
	}



	void read(const std::string &matrix_f);



};

} /* namespace Utils */
} /* namespace BioTools */

#endif /* SCORINGMATRIX_H_ */

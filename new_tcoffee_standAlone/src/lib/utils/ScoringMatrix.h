

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

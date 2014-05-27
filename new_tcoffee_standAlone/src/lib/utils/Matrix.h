


#ifndef MATRIX_H_
#define MATRIX_H_


// C header
#include <cstdlib>


// C++ header
#include <vector>

namespace BioTools {
namespace Utils {

/**
 * \brief A simple class to produce 2 dimensional matrices.
 */
template <typename DataType>
class Matrix {


private:
	std::vector<std::vector<DataType> > _matrix;
public:
	Matrix()
	{}

	/**
	 * \brief Constructor setting the size.
	 * @param dim1 Size of the first dimension.
	 * @param dim2 Size of the second dimension.
	 */
	Matrix(size_t dim1, size_t dim2);

	/**
	 * \brief Constructor allowing initialization of fields.
	 * @param dim1 Size of the first dimension.
	 * @param dim2 Size of the second dimension.
	 * @param init Value to initalize the fields with.
	 */
	Matrix(size_t dim1, size_t dim2, DataType init);

	/**
	 * \brief Destructor
	 */
	virtual ~Matrix()
	{

	}

	/**
	 * \brief Access operator
	 * @param index The index to acess.
	 * @return Reference to the field.
	 */
	std::vector<DataType> &operator[](unsigned int index)
	{
		return _matrix[index];
	}

	/**
	 * \brief Access operator
	 * @param index The index to acess.
	 * @return Reference to the field.
	 */
	const std::vector<DataType> &operator[](unsigned int index) const
	{
		return _matrix[index];
	}


	/**
	 * \brief Resizes the matrix to the new dimensions.
	 * @param dim_1 The first dimension
	 * @param dim_2 The second diemension
	 */
	void resize(size_t dim_1, size_t dim_2)
	{
		_matrix.resize(dim_1);
		for (size_t i=0; i<dim_1; ++i)
			_matrix[i].resize(dim_2);
	}

	/**
	 * \brief Returns the size of the first dimension.
	 * @return The size of the first dimension.
	 */
	size_t dim1() const
	{
		return _matrix.size();
	}

	/**
	 * \brief Returns the size of the second dimension.
	 * @return The size of the second dimension.
	 */
	size_t dim2() const
	{
		return _matrix[0].size();
	}

	void
	fill(const DataType &value)
	{
		size_t dim_1 = _matrix.size();
		size_t dim_2= (dim_1>0) ? _matrix[0].size() : 0;
		size_t j;
		for(size_t i=0; i<dim_1; ++i)
		{
			std::vector<DataType> &vec=_matrix[i];
			for (j=0; j<dim_2; ++j)
				vec[j] = value;
		}
	}
};


template <typename DataType>
Matrix<DataType>::Matrix(size_t dim_1, size_t dim_2, DataType init)
{
	_matrix.resize(dim_1);
	for (size_t i=0; i<dim_1; ++i)
		_matrix[i].resize(dim_2, init);
}


template <typename DataType>
Matrix<DataType>::Matrix(size_t dim_1, size_t dim_2)
{
	_matrix.resize(dim_1);
	for (size_t i=0; i<dim_1; ++i)
		_matrix[i].resize(dim_2);
}

} /* namespace Utils */
} /* namespace BioTools */

#endif /* MATRIX_H_ */

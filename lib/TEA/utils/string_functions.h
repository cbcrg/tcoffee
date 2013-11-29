/*
 * string_functions.h
 *
 *  Created on: Oct 10, 2011
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 */

#ifndef STRING_FUNCTIONS_H_
#define STRING_FUNCTIONS_H_


namespace BioTools {
namespace Utils {


template<typename DataType>
class Fixed_multi_array
{
private:
	size_t _dim1;
	size_t _dim2;
	DataType **_data_p;

public:
	Fixed_multi_array(size_t dim1, size_t dim2);
	~Fixed_multi_array();


	// Operators
	DataType* &operator[](size_t index)
	{
		return _data_p[index];
	}

	const DataType* &operator[](size_t index) const
	{
		return _data_p[index];
	}

	void set_all(const DataType &value);

};


template<typename DataType>
Fixed_multi_array<DataType>::Fixed_multi_array(size_t dim1, size_t dim2):_dim1(dim1),_dim2(dim2)
{
	_data_p = new DataType*[_dim1];
	for (size_t i = 0; i < _dim1; ++i)
		_data_p[i] = new DataType[dim2];
}

template<typename DataType>
Fixed_multi_array<DataType>::~Fixed_multi_array()
{
	for (size_t i = 0; i < _dim1; ++i)
		delete[] _data_p[i];
	delete[] _data_p;
}


template<typename DataType>
void
Fixed_multi_array<DataType>::set_all(const DataType &value)
{
	size_t j;
	DataType *data_p = NULL;
	for (size_t i = 0; i < _dim1; ++i)
	{
		data_p = _data_p[i];
		for (j = 0; j < _dim2; ++j)
			data_p[j] = value;
	}
}



} /* namespace Utils */
} /* namespace BioTools */




#endif /* STRING_FUNCTIONS_H_ */

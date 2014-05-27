

#ifndef VECTOR_H_
#define VECTOR_H_


/*! \file Vector.h
    \brief Contains the Vector class which can be used in the different Clustering algorithms
*/

// C header
#include <cstdlib>
#include <cmath>


// C++ Header
#include <vector>
#include <algorithm>

namespace BioTools
{
namespace Clustering
{


/**
 * \brief A vector class which made for the different clustering algorithms.
 * \todo changed internal vector to boost vector
 */
template<typename DataType>
class Vector {
private:
	std::vector<DataType> _data;
	size_t _vec_num;
	size_t _assignment;
	size_t _seqsize;
	size_t _seqdiv;

public:
	/**@{*/
	//! \name Constructors & Destructors
	/**
	 * \brief Standard Vector constructor
	 */
	Vector();

	/**
	 * \brief Copy constructor.
	 * \param vec The vector.
	 */
	Vector(const Vector &vec);

	/**
	 * \brief Constructor
	 * \param len The length of the vector
	 */
	Vector(size_t len);

	/**
	 * \brief Constructor
	 * \param len The length of the vector.
	 * \param value The value each field should be given.
	 * \param vec_num The vector id.
	 */
	Vector(size_t len, DataType value, size_t vec_num);
	/**@}*/

	// Operators
	/**@{*/
	//! \name Operators
	DataType &operator[](unsigned int index)
	{
		return _data[index];
	}

	const DataType &operator[](unsigned int index) const
	{
		return _data[index];
	}
	/**@}*/

	/**
	 * \brief Sets the cluster the vector belongs to.
	 * \param value The Cluster id.
	 */
	void assignment(size_t value)
	{
		_assignment = value;
	}

	/**
	 * \brief Returns the cluster membership.
	 * \return cluster membership.
	 */
	size_t assignment() const
	{
		return _assignment;
	}

	/**
	 * \brief Sets the vector size of a seq.
	 * \param value The sum of all contents of a vector.
	 */
	void seqsize(size_t value)
	{
		_seqsize = value;
	}

	/**
	 * \brief Returns the vector size of a seq.
	 * \return Sum of all vector cells.
	 */
	size_t seqsize() const
	{
		return _seqsize;
	}
	
	/**
	 * \brief Sets the number of different cells filled in a vector which reflects seq diversity.
	 * \param value Sum of different cells filled in a vector.
	 */
	void seqdiv(size_t value)
	{
		_seqdiv = value;
	}

	/**
	 * \brief Returns the Sum of different cells filled in a vector which reflects seq diversity.
	 * \return Sum of different cells filled in a vector.
	 */
	size_t seqdiv() const
	{
		return _seqdiv;
	}
	
	/**
	 * \brief Sets the id.
	 * \param value The new id.
	 */
	void id(size_t value)
	{
		_vec_num = value;
	}

	/**
	 * \brief Returns the id of this vector.
	 * \return the id of the vector.
	 */
	size_t id() const
	{
		return _vec_num;
	}

	virtual ~Vector();

	/**
	 * \brief Returns the dimension of the vector.
	 * \return The dimension.
	 */
	size_t size() const
	{
		return _data.size();
	}

	void resize(size_t new_size)
	{
		_data.resize(new_size);
	}
};


template<typename DataType>
Vector<DataType>::Vector() {
}

template<typename DataType>
Vector<DataType>::Vector(size_t len) {
	_data.resize(len);
}

template<typename DataType>
Vector<DataType>::Vector(size_t len, DataType value, size_t vec_num) : _data(len,value)
{
	_vec_num=vec_num;
	for (size_t i = 0; i < len; ++i)
		_data[i] = 0;

}


template<typename DataType>
Vector<DataType>::Vector(const Vector &vec)
{
	size_t len = vec.size();
	_data.resize(len);
	_vec_num = vec.id();
	_assignment = vec.assignment();
	_seqsize = vec.seqsize();
	_seqdiv = vec.seqdiv();
	
	for (size_t i = 0; i < len; ++i)
		_data[i] = vec[i];
}



template<typename DataType>
Vector<DataType>::~Vector() {
	// TODO Auto-generated destructor stub
}

/**
 * \brief Calculates the sum of all vector cells, which corresponds to seqence_length-1.
 * \param vec1 The vector.
 * \return The sum of all vector cells.
 */
template<typename DataType>
DataType
vec_sum(const Vector<DataType> &vec1 )
{
	size_t vec_len = vec1.size();
	unsigned int sum = 0;
	for (size_t i=0; i<vec_len; ++i)
		sum += vec1[i];
	return sum;
}

/**
 * \brief Calculates the sum of different vector cells that are filled, which corresponds to seqence diversity.
 * \param vec1 The vector.
 * \return The sum of non-zero vector cells.
 */
template<typename DataType>
DataType
vec_div(const Vector<DataType> &vec1 )
{
	size_t vec_len = vec1.size();
	unsigned int sum = 0;
	for (size_t i=0; i<vec_len; ++i){
	  if( vec1[i] != 0)
		sum++;
	}
	return sum;
}

/**
 * \brief Calculates the sum of different vector cells that are filled, which corresponds to seqence diversity.
 * \param vec1 The vector.
 * \return The sum of non-zero vector cells.
 */
template<typename DataType>
DataType
vec_div2(const Vector<DataType> &vec1 )
{
	size_t vec_len = vec1.size();
	unsigned int sum = 0;
	for (size_t i=0; i<vec_len; ++i){
	  if( vec1[i] != 0)
		sum+=i;
	}
	return sum;
}


//Distance functions between two vectors.


/**
 * \brief Calculates the squared distance between two vectors.
 * \param vec1 The first vector.
 * \param vec2 The second vector.
 * \return The squared distance.
 */
template<typename DataType>
DataType
sq_dist(const Vector<DataType> &vec1, const Vector<DataType> &vec2)
{
	size_t vec_len = vec1.size();
	double dist = 0;
	double tmp;
	for (size_t i=0; i<vec_len; ++i)
	{
		tmp = vec1[i]-vec2[i];
		dist += (tmp * tmp);
	}
	return dist;
}

/**
 * \brief Calculates the euclidean distance between two vectors.
 * \param vec1 The first vector.
 * \param vec2 The second vector.
 * \return The euclidean distance.
 */
template<typename DataType>
double
euclidean_dist(const Vector<DataType> &vec1, const Vector<DataType> &vec2)
{
	size_t vec_len = vec1.size();
	double dist = 0;
	double tmp;
	for (size_t i=0; i<vec_len; ++i)
	{
		tmp = vec1[i]-vec2[i];
		dist += (tmp * tmp);
	}
	return sqrt(dist);
}


/**
 * \brief Calculates the angle between two vectors.
 * \param vec1 The first vector.
 * \param vec2 The second vector.
 * \return The angle between the two values.
 */
template<typename DataType>
double
angle_dist(const Vector<DataType> &vec1, const Vector<DataType> &vec2)
{
	size_t vec_len = vec1.size();
	double dist = 0;
	for (size_t i=0; i<vec_len; ++i)
	{
		dist += vec1[i]*vec2[i];
	}
	return acos(dist);
}

/**
 * \brief Calculates eucleidan norm of a vector.
 * \param vec The vector.
 * \return The euclidean norm.
 *
 */
template<typename DataType>
double
euclidean_norm(const Vector<DataType> &vec)
{
	size_t vec_len = vec.size();
	double dist = 0;
	for (size_t i=0; i<vec_len; ++i)
	{
		dist += vec[i]*vec[i];
	}
	return sqrt(dist);
}


/**
 * \brief Calculates eucleidan norm of a vector.
 * \param vec The vector.
 * \return The euclidean norm.
 *
 */
template<typename DataType>
double
muscle_dist(const Vector<DataType> &vec1, const Vector<DataType> &vec2)
{
	size_t vec_len = vec1.size();
	size_t common=0;
	size_t l1=0, l2=0;
	for (size_t i=0; i<vec_len; ++i)
	{
		common += std::min(vec1[i],vec2[i]);
		l1+=vec1[i];
		l2+=vec2[i];
	}
	return 1-(common/std::min(l1,l2));
}


} // namespace Clustering
} // namespace BioTools

#endif /* VECTOR_H_ */

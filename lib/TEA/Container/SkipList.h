/*
 * SkipList.h
 *
 *  Created on: May 1, 2012
 *      Author: ck
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
#ifndef SKIPLIST_H_
#define SKIPLIST_H_

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <cstdio>

/*! \file SkipList.h
    \brief Contains a container Class.
*/


namespace BioTools {
namespace Container {


// SkipNode
template <typename KeyType, typename ValueType>
struct SkipNode {
	KeyType key;
    ValueType value;
    SkipNode<KeyType, ValueType> **next; // array of pointers

    SkipNode(const KeyType &k, const ValueType &v, int level):key(k),value(v)
    {
        next = new SkipNode<KeyType, ValueType>*[level+1];
        memset(next, 0, sizeof(SkipNode<KeyType, ValueType>*)*(level+1)); //Initializing pointers to NULL
    }

    ~SkipNode()
    {
        delete[] next;
    }
};



// SkipList

/**
 * \brief A class to store key value pairs.
 *
 */
template <typename KeyType, typename ValueType>
class SkipList {

private:
	const int _MAX_LEVEL;
	const float _PROBABILITY;
	size_t _size;
	int _current_max;
	SkipNode<KeyType, ValueType> **_head;

	// Private methods
	int random_level() const {
		float rand_val = static_cast<float>(rand())/RAND_MAX;
	    int lvl = (int)(log(rand_val)/log(1.-_PROBABILITY));
	    return lvl < _MAX_LEVEL ? lvl : _MAX_LEVEL-1;
	}

	SkipNode<KeyType, ValueType>* private_find(const KeyType &key) const;
	SkipNode<KeyType, ValueType>* private_find_smaller(const KeyType &key) const;


	// Iterator
	class SkipList_iterator
	{
	private:
		SkipNode<KeyType, ValueType> *_node;

	public:
		SkipList_iterator(SkipNode<KeyType, ValueType> *node):_node(node)
		{}

		SkipList_iterator& operator++()
		{
			_node = _node->next[0];
			return *this;
		}

		SkipList_iterator operator++(int)
		{
			SkipList::SkipList_iterator temp = *this;
			_node = _node->next[0];
			return temp;
		}

		SkipNode<KeyType, ValueType>& operator*()
		{
			return *_node;
		}

		SkipNode<KeyType, ValueType>* operator->()
		{
			return _node;
		}

		bool operator!=(const SkipList_iterator& rhs) const
		{
			return _node != rhs._node;
		}

		bool operator==(const SkipList_iterator& rhs) const
		{
			return _node == rhs._node;
		}
	};

public:
    typedef SkipList_iterator iterator;
    typedef const SkipList_iterator const_iterator;

	SkipList(int max_level, float probability) : _MAX_LEVEL(max_level), _PROBABILITY(probability), _size(0), _current_max(0)
	{
		_head = new SkipNode<KeyType, ValueType>*[_MAX_LEVEL];
		memset(_head, 0, sizeof(SkipNode<KeyType, ValueType>*)*(_MAX_LEVEL));
	}

	virtual ~SkipList();

	/**
	 * \brief Searches for an element in the SkipList.
	 * \param key The key value to search for.
	 * \return Iterator to the found element. If key is not found an iterator to the SkipList::end is returned.
	 */
    iterator find(const KeyType &key)
    {
    	return iterator(private_find(key));
    }


	/**
	 * \brief Searches for an element in the SkipList.
	 * \param key The key value to search for.
	 * \return Const iterator to the found element. If key is not found an iterator to the SkipList::end is returned.
	 */
    const_iterator find(const KeyType &key) const
    {
    	return const_iterator(private_find(key));
    }

	/**
	 * \brief Searches for an element in the SkipList.
	 * \param key The key value to search for.
	 * \return Iterator to the found element. If key is not found an iterator to the SkipList::end is returned.
	 */
    iterator find_smaller(const KeyType &key)
    {
    	return iterator(private_find_smaller(key));
    }


	/**
	 * \brief Inserts an element into the SkipList.
	 * \param key The key value.
	 * \param value The value.
	 */
    void insert(const KeyType &key, const ValueType &value);

	/**
	 * \brief Erases the element with a given key.
	 * \param key The key of the element to erase.
	 *
	 */
    void erase(const KeyType &key);

    /**
     * \brief Returns the number of elements in this object.
     * \return The number of elements.
     */
    size_t size() const
    {
    	return _size;
    }

    /**
     * \brief Checks if the container is empty.
     * \return True if container is empty, else false
     */
    bool empty() const
    {
    	return _size==0;
    }

	// Iterator
	/**@{*/
	//! \name Iterators

    /**
     * \brief Returns an iterator to the first element in the SkipList.
     */
    iterator begin() { return iterator(_head[0]); }
    /**
	 * \brief Returns a const iterator to the first element in the SkipList.
	 */
    const_iterator begin() const { return iterator(_head[0]); }
    /**
	 * \brief Returns an iterator to the element behind the last one in the SkipList.
	 */
    iterator end() { return iterator(0); }
    /**
	 * \brief Returns a const iterator to the element behind the last one in the SkipList.
	 */
    const_iterator end() const { return iterator(0); }
    /**@}*/
};


template <typename KeyType, typename ValueType>
SkipList<KeyType, ValueType>::~SkipList()
{
	SkipNode<KeyType, ValueType> *tmp;
	while (_head[0] != NULL)
	{
		tmp = _head[0]->next[0];
		delete _head[0];
		_head[0]=tmp;
	}
    delete[] _head;
}



template <typename KeyType, typename ValueType>
void SkipList<KeyType, ValueType>::insert(const KeyType &key, const ValueType &value)
{
	SkipNode<KeyType, ValueType> ***tmp = &_head;
	SkipNode<KeyType, ValueType> ****update = new SkipNode<KeyType, ValueType> ***[_MAX_LEVEL];

	// Find elements, remember path
	int i;
	for (i=_MAX_LEVEL-1; i >= 0; --i)
	{
		while (((*tmp)[i] != NULL) && ((*tmp)[i]->key< key))
			tmp = &((*tmp)[i]->next);
		update[i] = tmp;
	}

	// Update pointers and add new value
	if (((*tmp)[0] == NULL) || ((*tmp)[0]->key != key))
	{
		int lvl = random_level();
		if (lvl>_current_max)
			_current_max=lvl;
		SkipNode<KeyType, ValueType> *tmp_node= new SkipNode<KeyType, ValueType>(key, value, lvl);
		for (i=0; i<=lvl; ++i)
		{
			tmp_node->next[i]=(*(update[i]))[i];
			(*(update[i]))[i] = tmp_node;
		}
	}
	else
		(*tmp)[0]->value=value;
	delete[] update;
	++_size;
}


template <typename KeyType, typename ValueType>
void SkipList<KeyType, ValueType>::erase(const KeyType &key)
{
	SkipNode<KeyType, ValueType> ***tmp = &_head;
	SkipNode<KeyType, ValueType> ****update = new SkipNode<KeyType, ValueType> ***[_MAX_LEVEL];

	// Find elements, remember path
	for (int i=_current_max; i >= 0; --i)
	{
		while (((*tmp)[i] != NULL) && ((*tmp)[i]->key< key))
			tmp = &((*tmp)[i]->next);
		update[i] = tmp;
	}

	// Update pointers and delete element
	SkipNode<KeyType, ValueType>* p=(*tmp)[0];
	if (((*tmp)[0] != NULL) & ((*tmp)[0]->key == key))
	{
		for (int i = 0; i <= _current_max; ++i)
		{
			if ((*(update[i]))[i] == (*tmp)[0])
				(*(update[i]))[i] = (*tmp)[0]->next[i];
			else
				break;
		}
		delete p;
		while ((_current_max > 0) && (_head[_current_max] == NULL))
		    --_current_max;
	}
	--_size;
	delete[] update;

}


template <typename KeyType, typename ValueType>
SkipNode<KeyType, ValueType>*
SkipList<KeyType, ValueType>::private_find(const KeyType &key) const
{
	SkipNode<KeyType, ValueType> ** const *tmp = &_head;
	for (int i = _current_max; i>=0; --i)
	{
		while (((*tmp)[i] != NULL) && ((*tmp)[i]->key < key))
			tmp = &(*tmp)[i]->next;
	}
	if (((*tmp)[0] != NULL) && ((*tmp)[0]->key==key))
		return (*tmp)[0];
	else
		return NULL;
}


template <typename KeyType, typename ValueType>
SkipNode<KeyType, ValueType>*
SkipList<KeyType, ValueType>::private_find_smaller(const KeyType &key) const
{
	if ((empty()) || (_head[0]->key >= key))
		return NULL;
	SkipNode<KeyType, ValueType> ** const *tmp = &_head;
	for (int i = _current_max; i>=0; --i)
	{
		while (((*tmp)[i]->next[i] != NULL) && ((*tmp)[i]->next[i]->key < key))
			tmp = &(*tmp)[i]->next;
	}
	if (((*tmp)[0] != NULL) && ((*tmp)[0]->key<key))
		return (*tmp)[0];
	else
		return NULL;
}



} /* namespace Container */
} /* namespace BioTools */
#endif /* SKIPLIST_H_ */

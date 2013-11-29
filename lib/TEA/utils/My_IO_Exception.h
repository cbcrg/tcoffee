/*
 * My_IO_Exception.h
 *
 *  Created on: Feb 23, 2012
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

#ifndef MY_IO_EXCEPTION_H_
#define MY_IO_EXCEPTION_H_


class My_IO_Exception : public std::exception
{
public:
private:
	const char *_message;

public:
	virtual const char* what() const throw()
	{
		return _message;
	}

	My_IO_Exception(const char *message):_message(message)
	{}

	virtual ~My_IO_Exception() throw()
	{}
};




#endif /* MY_IO_EXCEPTION_H_ */

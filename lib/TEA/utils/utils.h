/*
 * utils.h
 *
 *  Created on: Jul 28, 2012
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

/*! \file utils.h
    \brief Several small classes and functions which do not fit anywhere else.
*/

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <cstring>
#include <cctype>



namespace BioTools {
namespace Utils {



/**
 * \brief A string tokenizer class which in difference to standard C is thread save.
 */
class StrTok
{
private:
	char *_str;

public:
	/**
	 * \brief Constructor
	 * \param str The string to use.
	 */
	StrTok(char *str);

	/**
	 * \brief Constructor.
	 */
	StrTok();

	/**
	 * \brief Set a new string.
	 * \param str The string.
	 */
	void set(char *str)
	{
		_str=str;
	}

	/**
	 * \brief Gets the next part of the string.
	 * \param sep The seperator to use.
	 */
	char* next(const char *sep);
};



/**@{*/
//! \name String functions

/**
 * \brief Converts all characters to lowercase.
 * \param[in,out] str The string to convert
 */
void str_upper(std::string &str);

/**
 * \brief Converts all characters to uppercase.
 * \param[in,out] str The string to convert
 */
void str_lower(std::string &str);

/**@}*/

/**
 * \brief Returns a reduced alphabet
 * \param alphabet The alphabet to use.
 * \return The clustering of the alphabet.
 */
std::string get_alphabet(std::string alphabet);

/**
 * \brief Encodes an alphabet into numbers
 * \param string  The alphabet to use.
 * \return The encoding of an alphabet.
 */
short *encode(std::string alphabet);

}
}

#endif /* UTILS_H_ */

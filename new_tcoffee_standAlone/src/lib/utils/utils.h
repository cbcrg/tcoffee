

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

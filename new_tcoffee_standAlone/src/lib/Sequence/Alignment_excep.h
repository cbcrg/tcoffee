/*
 * Alignment_excep.h
 *
 *  Created on: Feb 23, 2012
 *      Author: ck
 */

#ifndef ALIGNMENT_EXCEP_H_
#define ALIGNMENT_EXCEP_H_



namespace BioTools {
namespace Seq {



/**
 * \brief A class for throwing exceptions
 */
class Alignment_excep: public std::exception
{
private:
	const char *_message;

public:
	virtual const char* what() const throw()
	{
		return _message;
	}

	Alignment_excep(const char *message):_message(message)
	{}


	virtual ~Alignment_excep() throw()
	{}

};


}
}

#endif /* ALIGNMENT_EXCEP_H_ */

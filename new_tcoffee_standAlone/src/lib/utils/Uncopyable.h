

#ifndef UNCOPYABLE_H_
#define UNCOPYABLE_H_

namespace BioTools {
namespace Utils {

/*! \file Uncopyable.h
    \brief Contains a single class which is not copyable.
*/


/**
 * \brief A class to which is not copyable.
 *
 * Any class derived from this class cannot be copied. The copy and the acess operator are declared
 * but not implemented thus any copy/assignment operation will fail on compilation.
 *
 */
class Uncopyable {

private:
	Uncopyable(const Uncopyable&);
	Uncopyable& operator=(const Uncopyable&);

public:
	Uncopyable()
	{}

	virtual ~Uncopyable()
	{}
};

} /* namespace Utils */
} /* namespace BioTools */
#endif /* UNCOPYABLE_H_ */

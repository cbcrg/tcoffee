/*
 * Uncopyable.h
 *
 *  Created on: May 26, 2012
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

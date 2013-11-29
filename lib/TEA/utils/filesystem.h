/*
 * filesystem.h
 *
 *  Created on: Oct 9, 2011
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

#ifndef FILESYSTEM_H_
#define FILESYSTEM_H_


// C header
#include <cstdio>
#include <cstring>
#include <cerrno>

// C++ header
#include <string>
#include <vector>
#include <exception>

// BioTools header
#include "./My_IO_Exception.h"

/*! \file filesystem.h
    \brief Contains some filesystem functions.
*/




namespace BioTools {
namespace Utils {



inline FILE *
my_fopen(std::string name_f, std::string mode)
{
	FILE *name_F = fopen(name_f.c_str(), mode.c_str());
	if (name_F == NULL)
	{
		throw My_IO_Exception(strerror(errno));
	}
	else
		return name_F;
}

} /* namespace Utils */
} /* namespace BioTools */

//void
//get_files_in_dir(const std::string &directory, std::vector<std::string> &files, std::vector<std::string> ending, bool recursive);

#endif /* FILESYSTEM_H_ */



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

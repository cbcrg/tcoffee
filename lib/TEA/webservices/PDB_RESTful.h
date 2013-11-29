/*
 * PDB_RESTful.h
 *
 *  Created on: Feb 18, 2012
 *      Author: ck
 */

#ifndef PDB_RESTFUL_H_
#define PDB_RESTFUL_H_

#include <string>
#include <curl/curl.h>

namespace BioTools {
namespace Web {

class PDB_RESTful {
private:


public:
	PDB_RESTful();
	virtual ~PDB_RESTful();

	void get_pdb_file(const std::string &pdb_id, const std::string &out_f);
};



} /* namespace Web */
} /* namespace BioTools */
#endif /* PDB_RESTFUL_H_ */

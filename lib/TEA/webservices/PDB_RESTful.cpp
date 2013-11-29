/*
 * PDB_RESTful.cpp
 *
 *  Created on: Feb 18, 2012
 *      Author: ck
 */

#include "PDB_RESTful.h"

namespace BioTools {
namespace Web {

using namespace std;

PDB_RESTful::PDB_RESTful() {
	// TODO Auto-generated constructor stub

}

PDB_RESTful::~PDB_RESTful() {
	// TODO Auto-generated destructor stub
}


void PDB_RESTful::get_pdb_file(const string &pdb_id, const string &out_f)
{

	FILE *out_F = fopen(out_f.c_str(), "w");
	CURL* easyhandle = curl_easy_init();
	string url("http://www.pdb.org/pdb/files/");
	url.append(pdb_id);
	url.resize(url.size()-2);
	url.append(".pdb");
	printf("%s\n", url.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_URL, url.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  out_F);
	curl_easy_perform(easyhandle);
	curl_easy_cleanup(easyhandle);
	fclose(out_F);
}


} /* namespace Web */
} /* namespace BioTools */

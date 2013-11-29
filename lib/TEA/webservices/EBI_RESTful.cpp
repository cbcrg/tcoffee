/*
 * EBI_RESTful.cpp
 *
 *  Created on: Feb 18, 2012
 *      Author: ck
 */

#include "EBI_RESTful.h"

namespace BioTools {
namespace Web {

using namespace std;
using namespace BioTools::Seq;

std::string EBI_RESTful::_email="";
void EBI_RESTful::lib_start(const string &email)
{

	unsigned int pos = email.find('@');
	unsigned int email_len = email.size();
	_email.reserve(email_len+2);
	_email.append(email, 0, pos);
	_email.append("%40");
	_email.append(email, pos+1, email_len-(pos+1));
	curl_global_init(CURL_GLOBAL_ALL);
}



void writefunc(void *ptr, size_t size, size_t nmemb, string *s)
{
	s->assign((char*)ptr);
}


short
EBI_RESTful::get_job_status(const string &job_id)
{
	CURL *easyhandle = curl_easy_init();
	string status;
	string url("http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/");
	url.append(job_id);
	curl_easy_setopt(easyhandle, CURLOPT_URL, url.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_WRITEFUNCTION, writefunc);
	curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  &status);
	curl_easy_perform(easyhandle);
	curl_easy_cleanup(easyhandle);
	if (status == "RUNNING")
		return 0;
	if (status == "FINISHED")
		return 1;
	else
		return -1;
}


void
EBI_RESTful::get_job_result(const string &job_id,
						 const string &out_f)
{
	short stat;
	FILE *out_F = fopen(out_f.c_str(), "w");
	while (!(stat = get_job_status(job_id)))
	{
		sleep(5);
	}
	if (stat == -1)
	{
		fclose(out_F);
		return;
	}

	CURL* easyhandle = curl_easy_init();
	string status;
	string url("http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/");
	url.append(job_id);
	url.append("/out");
// 		printf("%s\n", job_id.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_URL, url.c_str());
// 	curl_easy_setopt(easyhandle, CURLOPT_WRITEFUNCTION, writefunc);
	curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  out_F);
	curl_easy_perform(easyhandle);
	curl_easy_cleanup(easyhandle);
	fclose(out_F);
}


void
EBI_RESTful::submit_blast_job(const std::string &program,
						   const std::string &type,
						   const std::string &database,
						   double ecut,
						   const Sequence &seq,
						   std::string &job_id)
{
	char exp[10];
	sprintf(exp, "%.1f",ecut);
	string request;
	request.reserve(seq.size()+200);
	request.append("sequence=%3E");
	request.append(seq.name());
	request.append("%0A");
	request.append(seq.sequence());
	request.append("&database=");
	request.append(database);
	request.append("&stype=");
	request.append(type);
	request.append("&exp=");
	request.append(exp);
	request.append("&align=7");
	request.append("&program=");
	request.append(program);
	request.append("&email=");
	request.append(_email);
	printf("%s\n", request.c_str());
	CURL* easyhandle = curl_easy_init();

	if(easyhandle) {
		curl_easy_setopt(easyhandle, CURLOPT_POSTFIELDS, request.c_str());
		curl_easy_setopt(easyhandle, CURLOPT_WRITEFUNCTION, writefunc);
		curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  &job_id);
		curl_easy_setopt(easyhandle, CURLOPT_URL, "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/");
		curl_easy_perform(easyhandle);
		curl_easy_cleanup(easyhandle);
	}
}



} /* namespace Web */
} /* namespace BioTools */

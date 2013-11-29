/*
 * EBI_RESTful.h
 *
 *  Created on: Feb 18, 2012
 *      Author: ck
 */

#ifndef EBI_RESTFUL_H_
#define EBI_RESTFUL_H_

#include<cstdio>
#include<string>

#include <curl/curl.h>

namespace BioTools {
namespace Web {



/** \file EBI_RESTful.h
 * \brief Contains function to use functinality of the EBI RESTful service via C++.
 *
 * \warning This functions should not be used in a high number at the same time (e.g. on a cluster) as too many connections at the same time may block you or your whole insitute from acessing
 * this servers. For more information please see the <a href="http://www.ebi.ac.uk/Tools/webservices/"> EBI webservice</a> webpage.
 */

/**
 * \brief A class to access the EBI_RESTful service.
 *
 * This class needs the curl library available from <a href="http://curl.haxx.se/libcurl/">http://curl.haxx.se/libcurl/</a> which is already included in many LINUX distributions.
 * \warning This functions should not be used in a high number at the same time (e.g. on a cluster) as too many connections at the same time may block you or your whole insitute from acessing
 * this servers. For more information please see the <a href="http://www.ebi.ac.uk/Tools/webservices/"> EBI webservice</a> webpage.
 *
 */
class EBI_RESTful {

private:
	static std::string _email;

public:

	/**@{*/
	//! \name Class start/end

	// Static Methods
	/**
	 * \brief This function has to be used before any of the other functions can be used.
	 * \param email The email address needs to be given if you want to use the EBI RESTful service.
	 */
	void static lib_start(const std::string &email);

	/**
	 * \brief This function should be called after all jobs have been finished as the curl lib neads to be cleaned up.
	 *
	 * It just calles the curl_global_cleanup() function.
	 */
	void static lib_stop()
	{
		//printf("CLEANING UP!\n");
		curl_global_cleanup();
	}
	/**@}*/

	/**@{*/
	//! \name Sending and retreiving jobs.

	// Sending retreiving jobs.
	/**
	 * \brief Submits a job to the EBI server.
	 *
	 * \param[in] program The blast program to use.
	 * \param[in] type The type of sequence.
	 * \param[in] database The database to blast against.
	 * \param[in] ecut The evalue threshold to use.
	 * \param[in] seq The sequence to blast.
	 * \param[out] job_id The job_id given by the EBI.
	 */
	void static submit_blast_job(const std::string &program, const std::string &type, const std::string &database, double ecut, const BioTools::Seq::Sequence &seq, std::string &job_id);

	/**
	 * \brief Asks the EBI webserver for the jobstatus of a job.
	 *
	 * \param job_id The job_id to check for.
	 * \return The job status: 0 job is still running, 1 job is finished, -1 else
	 */
	short static get_job_status(const std::string &job_id);

	/**
	 * \brief Gets the job results.
	 *
	 * \param job_id The job_id.
	 * \param out_f The file to write the string to.
	 */
	void static get_job_result(const std::string &job_id, const std::string &out_f);
	/**@}*/

};


} /* namespace Web */
} /* namespace BioTools */
#endif /* EBI_RESTFUL_H_ */

#include <stdlib.h>

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <sched.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

static int entry_len;

/**
 * \file util_constraints_list.c
 * All utilities for the Constraint_list.
 */


int compare_constraint_list_entry ( const void*vx, const void*vy);


/*********************************************************************/
/*                                                                   */
/*                        Post Process Constraint_list               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

/*********************************************************************/
/*                                                                   */
/*                         PRODUCE IN LIST                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list* make_test_lib (Constraint_list *CL);
 

Constraint_list *fork_line_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode,Job_TC *job, int nproc,FILE *local_stderr);
Constraint_list *fork_cell_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode,Job_TC *job, int nproc,FILE *local_stderr);
Constraint_list *fork_subset_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, int nproc, FILE *local_stderr);
int job2first_seq(Job_TC *job);

Constraint_list *produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode)
{
	Job_TC *job=NULL;
	FILE *local_stderr;
	int njob;
	int nproc;

	store_seq_type (S->type);
	if ( CL==NULL)CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?vtmpfile():NULL, NULL);
	local_stderr=(CL->local_stderr!=NULL)?CL->local_stderr:stderr;

	CL->local_stderr=vfopen("/dev/null", "w");
	job=queue2heap(method2job_list ( method,S,weight, CL->lib_list,CL->DM, CL));
	njob=queue2n(job)+1;

	nproc=get_nproc();

	if (strstr ( CL->multi_thread, "jobcells"))return fork_cell_produce_list (CL, S, method, weight, mem_mode,job,nproc,local_stderr);
	else if (strstr ( CL->multi_thread, "joblines"))return fork_line_produce_list (CL, S, method, weight, mem_mode,job, nproc,local_stderr);
	else if (strstr ( CL->multi_thread, "jobs"))return fork_subset_produce_list (CL, S, method, weight, mem_mode,job, nproc,local_stderr); //Recommended default
	else return fork_subset_produce_list (CL, S, method, weight, mem_mode,job,nproc,local_stderr); //Recommended default
}
int job2first_seq(Job_TC *job)
{
	int *seqlist;
	int r;

	if (!job) return -1;
	else if ( !job->param)return -1;
	else if ( !(job->param)->seq_c) return -1;
	seqlist=string2num_list ((job->param)->seq_c);
	if (seqlist[0]<2)r=-1;
	else r=seqlist[2];
	vfree (seqlist);
	return r;
}

Constraint_list *fork_subset_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, int nproc, FILE *local_stderr)
{
	//forks lines of the matrix
	int a;
	Job_TC *heap,*end,*start, ***jl;
	TC_method *M;
	char **pid_tmpfile;
	int   *pid_list;
	int pid, npid, njob, max_nproc;
	int max=max_n_pid();
	
	max_nproc=nproc;
	max_nproc*=2;
	/*OUT_MODE:
	A->  alignment provided in a file
	L->  list      provided in a file
	aln-> alignment computed and provided in a file
	list->list      computed and provided in a file
	*/

	if ( job->jobid==-1)
	{
		M=(job->param)->TCM;
		fprintf (local_stderr, "\n\tMethod %s: No Suitable Sequences [Type: %s]\n", method,M->seq_type);
		return CL;
	}

	job=queue2heap (job);
	heap=job;
	njob=queue2n (job);

	M=(job->param)->TCM;
	if (M)M->PW_CL=method2pw_cl ( M, CL);
	pid_tmpfile=(char**)vcalloc (MAX(njob,nproc)+1, sizeof (char*));
	pid_list   =(int *)vcalloc (max, sizeof (int *));

	fprintf ( local_stderr, "\n\tMulti Core Mode (Compute): %d processor(s) [subset]\n", nproc);

	jl=split_job_list(job,nproc);
	a=npid=0;
	while (jl[a])
	{
	  start=job=jl[a][0];
	  end=jl[a][1];
	  pid_tmpfile[a]=vtmpnam(NULL);
	  pid=vvfork(NULL);
	  
	  if (pid==0)//child process
	    {
	      int done, todo;
	      
	      freeze_constraint_list (CL);//record the current state, so as not to dump everything
	      initiate_vtmpnam(NULL);
	      vfclose(vfopen (pid_tmpfile[a],"w"));
	      
	      todo=0;
	      while (job!=end){todo++;job=job->c;}
	      job=start;
	      
	      done=0;
	      while (job!=end)
		{
		  if (a==0)output_completion ( local_stderr,done,todo,1, "Submit   Job");
		  job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",CL,submit_lib_job, retrieve_lib_job, CL->multi_thread );
		  
		  job=submit_job (job);
		  retrieve_job (job);
		  
		  job=job->c;
		  done++;
		}
	      dump_constraint_list (CL, pid_tmpfile[a], "a");
	      unfreeze_constraint_list (CL);
	      myexit (EXIT_SUCCESS);
	    }
	  else
	    {
	      pid_list[pid]=npid;
	      //set_pid(pid);
	      npid++;
	      a++;
	    }
	}
	
	//wait for all processes to finish
	for (a=0; a<npid; a++)
	{
		pid=vwait(NULL);

		CL=undump_constraint_list (CL, pid_tmpfile[pid_list[pid]]);
		remove(pid_tmpfile[pid_list[pid]]);
	}

	vfree (pid_list);
	vfree (pid_tmpfile);

	job=heap;
	while (job)	job=delete_job (job);
	CL->local_stderr=local_stderr;
	free_queue  (heap);
	return CL;
}

    
Constraint_list *fork_line_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, int nproc,FILE *local_stderr)
{
	//forks lines of the matrix
	int a;
	Job_TC *heap;
	TC_method *M;

	char **pid_tmpfile;
	int   *pid_list;
	int pid,npid, njob;
	int max_nproc, submited;
	int cseq, seq, nlines;
	int max=max_n_pid();
	max_nproc=nproc;
	max_nproc*=2;
	/*OUT_MODE:
	A->  alignment provided in a file
	L->  list      provided in a file
	aln-> alignment computed and provided in a file
	list->list      computed and provided in a file
	*/



 
	if ( job->jobid==-1)
	{
		M=(job->param)->TCM;
		fprintf (local_stderr, "\n\tMethod %s: No Suitable Sequences [Type: %s]\n", method,M->seq_type);
		return CL;
	}


	job=queue2heap (job);

	heap=job;
	M=(job->param)->TCM;
	if (M)M->PW_CL=method2pw_cl ( M, CL);


	/* Cf. parse method for possible out_mode flags*/

	njob=queue2n(job)+1;
	pid_tmpfile=(char**)vcalloc (njob, sizeof (char*));

	pid_list   =(int *)vcalloc (max, sizeof (int *));
	fprintf ( local_stderr, "\n\tMulti Core Mode (Compute): %d processor(s) [jobline]\n", nproc);


	//count the number of lines
	cseq=-1;
	nlines=0;
	while (job)
	{
		nlines++;
		seq=job2first_seq(job);
		if ( seq!=cseq)
		{
			cseq=seq;
			while (job && cseq==job2first_seq(job))job=job->c;
		}
	}
	job=heap;

	npid=submited=0;
	cseq=-1;
	while (job)
	{
		seq=job2first_seq(job);
		if ( seq!=cseq)
		{
			cseq=seq;
			pid_tmpfile[npid]=vtmpnam(NULL);

			pid=vvfork(NULL);

			if (pid==0)//Child Process
			{
				initiate_vtmpnam(NULL);
				while (job && cseq==job2first_seq(job))
				{
					job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",CL,submit_lib_job, retrieve_lib_job, CL->multi_thread );
					job=submit_job (job);
					retrieve_job (job);
					dump_constraint_list ((job->io)->CL,pid_tmpfile[npid], "a");
					job=job->c;
				}
				myexit (EXIT_SUCCESS);
			}
			else //parent process
			{

				pid_list[pid]=npid;
				//set_pid (pid);
				npid++;
				submited++;
				fprintf(stderr, "%i %i\n", submited, max_nproc);
				if (submited>max_nproc)
				{
					//wait for nproc
					for (a=0; a<nproc; a++)
					{
						int index;
						local_stderr=output_completion ( local_stderr,npid,nlines,1, "Processed   Job");
						pid=vwait(NULL);
						index=pid_list[pid];
						CL=undump_constraint_list (CL,pid_tmpfile[index]);
						remove (pid_tmpfile[index]);
						submited--;
					}
				}
			}
		}
		else
		{
			job=job->c;
		}
	}

	for (a=0; a<submited; a++)
	{
		int index;
		local_stderr=output_completion ( local_stderr,npid-(submited-a),nlines,1, "Processed   Job");
		pid=vwait(NULL);
		index=pid_list[pid];
		CL=undump_constraint_list(CL,pid_tmpfile[index]);
		remove (pid_tmpfile[index]);
	}

	vfree (pid_list);
	vfree (pid_tmpfile);

	CL->local_stderr=local_stderr;

	free_queue  (heap);

	return CL;
}

Constraint_list *fork_cell_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job,int nproc, FILE *local_stderr)
{
	//forks cells of the matrix
	int a;
	Job_TC *heap;
	TC_method *M;

	int *pid_list;
	char **pid_tmpfile;
	int pid,npid, njob;
	int max_nproc;
	int submited;
	int max=max_n_pid();
	max_nproc=nproc;

	/*OUT_MODE:
	A->  alignment provided in a file
	L->  list      provided in a file
	aln-> alignment computed and provided in a file
	list->list      computed and provided in a file
	*/


	if ( job->jobid==-1)
	{
		M=(job->param)->TCM;
		fprintf (local_stderr, "\n\tMethod %s: No Suitable Sequences [Type: %s]\n", method,M->seq_type);
		return CL;
	}


	job=queue2heap (job);

	heap=job;
	M=(job->param)->TCM;
	if (M)M->PW_CL=method2pw_cl ( M, CL);


	/* Cf. parse method for possible out_mode flags*/

	njob=queue2n(job)+1;
	pid_tmpfile=(char**)vcalloc (njob, sizeof (char*));
	pid_list   =(int *)vcalloc  (max, sizeof (int *));

	fprintf ( local_stderr, "\n\tMulti Core Mode (compute): %d processor(s):\n", nproc);
	npid=0;
	submited=0;
	while (job)
	{
		job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",CL,submit_lib_job, retrieve_lib_job, CL->multi_thread );
		pid_tmpfile[npid]=vtmpnam(NULL);
		pid=vvfork (NULL);
		if ( pid==0)
		{
			initiate_vtmpnam (NULL);
			job=submit_job (job);
			retrieve_job (job);
			dump_constraint_list ((job->io)->CL,pid_tmpfile[npid], "w");
			myexit (EXIT_SUCCESS);
		}
		else
		{
			job=job->c;
			pid_list[pid]=npid;
			//set_pid(pid);
			npid++;
			submited++;

			if (submited>max_nproc)
			{
				for (a=0; a<nproc; a++)
				{
					int index;
					local_stderr=output_completion ( local_stderr,npid,njob,1, "Processed   Job");
					pid=vwait(NULL);
					index=pid_list[pid];
					CL=undump_constraint_list (CL,pid_tmpfile[index]);
					remove (pid_tmpfile[index]);
					submited--;
				}
			}
		}
	}

	for (a=0; a<submited; a++)
	{
		int index;
		local_stderr=output_completion ( local_stderr,npid-a, npid,1, "Processed   Job");
		pid=vwait(NULL);
		index=pid_list[pid];
		CL=undump_constraint_list (CL,pid_tmpfile[index]);
		remove (pid_tmpfile[index]);
	}
	vfree (pid_list);
	vfree (pid_tmpfile);

	CL->local_stderr=local_stderr;

	free_queue  (heap);

	return CL;
}

Job_TC *retrieve_lib_job ( Job_TC *job)
{
	Job_param_TC *p;
	Job_io_TC *io;
	TC_method *M;

	p=job->param;
	io=job->io;
	M=(job->param)->TCM;


	if ( job->status==EXIT_SUCCESS)
	{
		if ( !M) return job;
		else if (strm2(M->out_mode, "aln", "A"))
		{
			io->CL=read_constraint_list (io->CL,io->out,"aln","disk",M->weight);
		}
		else if (strm2(M->out_mode, "lib","L"))
		{

			io->CL=read_constraint_list (io->CL,io->out,"lib","disk",M->weight);

		}
		return job;
	}
	else
		return NULL;
}

int add_method_output2method_log (char *l,char *command,Alignment *A, Constraint_list *CL, char *io_file)
{
	static int header, set;
	static char *file, *log_file;

	if ( set && log_file==NULL && l==NULL) return 0;
	if (!set )
	  {

	    log_file=get_string_variable ("method_log");if (log_file && strm (log_file, "no"))log_file=NULL; set=1;
	  }

	if (!file)file=vtmpnam (NULL);

	if ( l);
	else if (!l && log_file) l=log_file;
	else return 0;


	if (!header){printf_file ( l, "w", "# CONC_MSF_FORMAT_01\n");header=1;}
	if (command)printf_file (l,  "a", "%s\n#----------------------------------------------\n#%s\n", TC_REC_SEPARATOR,command);


	if ( A)
	{

		io_file=file;
		output_fasta_aln (io_file, A);
	}
	else if (CL)
	{
		io_file=file;
		vfclose (save_constraint_list ( CL, 0, CL->ne,io_file, NULL, "ascii",CL->S));
	}
	else
		file_cat (io_file,l);


	return 1;
}



Job_TC *submit_lib_job ( Job_TC *job)
{
	Job_param_TC *p;
	Job_io_TC *io;
	TC_method *M;
	static char *tdir=NULL;
	static char *cdir=NULL;
	

	
	if (!tdir)tdir=get_tmp_4_tcoffee();
	if (!cdir)cdir=get_pwd(NULL);
	  
	tdir=get_tmp_4_tcoffee();
	cdir=get_pwd(NULL);
	
	p=job->param;
	io=job->io;
	M=(job->param)->TCM;
	add_method_output2method_log (NULL, p->aln_c, NULL, NULL, NULL);
	if ( getenv4debug ("DEBUG_LIBRARY"))fprintf ( stderr, "\n[DEBUG_LIBRARY:produce_list] Instruction: %s\n", p->aln_c);

	if ( !M)
	  {
	    return job;
	  }
	else if (strm4 (M->out_mode,"A", "L", "aln", "lib"))
	  {
	    char *com;
	    static int do_flip;
	    int flipped=0;

	    if (!do_flip)
	      {
		do_flip=get_int_variable ("flip");
		if (!do_flip)do_flip=-1;
	      }

	    seq_list2in_file ( M, (io->CL)->S, p->seq_c, io->in);
	    if (do_flip!=-1)
	      {
		if ((rand()%100)<do_flip)
		  {
		    invert_seq_file (io->in);
		    flipped=1;
		  }
	      }

	    com=(char*)vcalloc ( strlen (p->aln_c)+100, sizeof (char));
	    sprintf (com, "%s", p->aln_c);
	    substitute (com, "//", "/");
	    substitute (com, tdir, "./");
	    chdir (tdir);
	    printf_system ("%s ::IGNORE_FAILURE::", com);
	    vfree (com);
	    chdir (cdir);
	    add_method_output2method_log (NULL,NULL, NULL, NULL, io->out);

	    if (!evaluate_sys_call_io (io->out,p->aln_c, "") || (strm (M->out_mode, "aln") && !(is_aln (io->out) || is_seq(io->out))) )
	      {
		job->status=EXIT_FAILURE;
		return job;
	      }
	    if (flipped==1)invert_aln_file (io->out);
	  }
	else if ( strm2 (M->out_mode, "fA", "fL"))
	  {

	    
	    io->CL= seq2list(job);
	    if (!io->CL)
	      {
		add_warning (stderr, "FAILED TO EXECUTE:%s [SERIOUS:%s]", p->aln_c, PROGRAM);
		job->status=EXIT_FAILURE;
	      }
	  }
	else
	  {
	    myexit(fprintf_error ( stderr, "\nERROR: Unknown out_mode=%s for method[FATAL:%s]\n", M->out_mode, M->executable));
	  }

	return job;
}


      

Job_TC* method2job_list ( char *method_name,Sequence *S, char *weight, char *lib_list, Distance_matrix *DM, Constraint_list *CL)
{
  int preset_method;
  static char *fname, *bufS, *bufA;
  char *in,*out;
  TC_method *method;
  char aln_mode[100];
  char out_mode[100];
  Job_TC *job;
  int hijack_P_jobs=1;
  int debugio=atoigetenv ("DEBUG_METHOD_4_TCOFFEE");

  /*A method can be:
    1- a pre computed alignment out_mode=A
    2- a precomputed Library    out_mode=L
    3- a method producing an alignment out_mode=aln
    4- a method producing an alignment out_mode=list
    5- a function producing an alignment out_mode=faln
    6- a function producing a library    out_mode=flist
  */

  if ( !fname)
    {
      fname=(char*)vcalloc ( 1000, sizeof (char));
      bufS=(char*)vcalloc ( S->nseq*10, sizeof (char));
    }

  /*Make sure that fname is a method*/

 
  sprintf(fname, "%s", method_name);

  if ( fname[0]=='A' || fname[0]=='L')
    {
      
      method=method_file2TC_method("no_method");
      sprintf ( method->out_mode, "%c", fname[0]);

      if (!strm (weight, "default"))sprintf ( method->weight, "%s", weight);

      return print_lib_job(NULL,"param->out=%s param->TCM=%p",fname+1, method);
    }
  else if ( fname[0]=='M' && is_in_pre_set_method_list (fname+1))
    {
      
      preset_method=1;
      fname++;
    }
  else if ( is_in_pre_set_method_list (fname))
    {
      
      preset_method=1;
    }
  else
    {
      char buf[1000];
      
      if ( check_file_exists ( fname));
      else if (fname[0]=='M' && check_file_exists(fname+1));
      else
	{
	  sprintf ( buf, "%s/%s", get_methods_4_tcoffee(), fname);
	  if( check_file_exists(buf)){sprintf ( fname, "%s", buf);}
	  else
	    {
	      myexit (fprintf_error ( stderr, "%s is not a valid method", fname));
	    }
	}
    }

  
  method=method_file2TC_method(fname);
  job=print_lib_job (NULL, "param->TCM=%p", method);
  job->jobid=-1;

  if (!strm (weight, "default"))sprintf ( method->weight, "%s", weight);

  sprintf ( aln_mode, "%s", method->aln_mode);
  sprintf ( out_mode, "%s", method->out_mode);


  if (lib_list && lib_list[0])
    {
      static char **lines, **list=NULL;
      int a,i, x, n, nl;



      if ( lines) free_char (lines, -1);


      if ( strstr (lib_list, "prune"))
	{
	  lines=file2lines (list2prune_list (S,DM->similarity_matrix));
	}
      else
	{
	  lines=file2lines (lib_list);
	}

      nl=atoi (lines[0]);
      for (a=1; a<nl; a++)
	{

	  if (list) free_char (list, -1);
	  if (isblanc (lines[a]))continue;

	  list=string2list (lines[a]);n=atoi(list[1]);
	  if ( n> 2 && strm (aln_mode, "pairwise"))continue;
	  if ( n==2 && strm (aln_mode, "multiple"))continue;

	  for (i=2; i<n+2; i++)
	    {
	      if ( is_number (list[i]));
	      else if ((x=name_is_in_list (list[i], S->name, S->nseq, 100))!=-1)sprintf(list[i], "%d", x);
	      else
		{
		  add_warning ( stderr, "%s is not part of the sequence dataset", list[i]);
		  continue;
		}
	    }
	  sprintf ( bufS, "%s", list[1]);
	  for ( i=2; i<n+2; i++) {strcat (bufS, " ");strcat ( bufS, list[i]);}


	  bufA=make_aln_command (method, in=vtmpnam(NULL),out=vtmpnam(NULL));

	  if (!debugio && strrchr(bufA, '>')==NULL)strcat (bufA,TO_NULL_DEVICE);
	  if ( check_seq_type ( method, bufS, S))
	    {
	      job->c=print_lib_job (NULL, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ", method, fname, bufA, bufS, in, out);
	      job=queue_cat (job, job->c);
	    }
	  vfree (bufA);

	}
    }
  else if ( strcmp (aln_mode, "multiple")==0)
    {
      int d;
      char buf[10000];

      sprintf (bufS, "%d",S->nseq);
      for (d=0; d< S->nseq; d++)
	{
	  sprintf ( buf," %d",d);
	  strcat ( bufS, buf);
	}

      bufA=make_aln_command (method, in=vtmpnam(NULL),out=vtmpnam(NULL));
      if (!debugio && strrchr(bufA, '>')==NULL)strcat (bufA,TO_NULL_DEVICE);

      if ( check_seq_type ( method, bufS, S))
	{
	  job->c=print_lib_job (NULL, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ", method, fname, bufA, bufS, in, out, S->template_file);

	  job=queue_cat (job, job->c);

	}
      vfree (bufA);
    }
  else if ( strstr(aln_mode, "o2a"))
    {
      int x, n;
      static char *tmpf;
      int byte=CL->o2a_byte;
      int max=0;
      FILE *fp;
     
      for (x=0; x<(CL->S)->nseq; x++)max+=(CL->master[x]);

      if (CL->o2a_byte>=max)byte=(max/get_nproc())+1;

      if (!tmpf)tmpf=vtmpnam (NULL);

      fp=vfopen (tmpf, "w");
      for (n=0,x=0; x<(CL->S)->nseq; x++)
	{
	  if (CL->master[x]){fprintf (fp, "%d ", x);n++;}
	  if (n==byte || (n && x==(CL->S)->nseq-1))
	    {
	      vfclose (fp);
	      sprintf (bufS, "%d %s", n,file2string (tmpf));n=0;
	      bufA=make_aln_command (method, in=vtmpnam(NULL),out=vtmpnam(NULL));
	      if (!debugio && strrchr(bufA, '>')==NULL)strcat (bufA,TO_NULL_DEVICE);
	      if ( check_seq_type ( method, bufS, S))
		{
		  job->c=print_lib_job (NULL, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ", method, fname, bufA, bufS, in, out, S->template_file);
		  job=queue_cat (job, job->c);
		}
	      vfree (bufA);
	      fp=vfopen (tmpf, "w");
	    }
	}
      vfclose (fp);
    }


  else if ( strstr(aln_mode, "pairwise"))
    {
      int do_mirror, do_self, x, y, id;
      int **set;

      do_mirror=(strstr(aln_mode, "m_"))?1:0;
      do_self=(strstr(aln_mode, "s_"))?1:0;
      set=declare_int (S->nseq, S->nseq);
      if (method->minid!=0 || method->maxid!=100)
	DM=CL->DM=cl2distance_matrix ( CL,NOALN,NULL,NULL,1);

      for (x=0; x< S->nseq; x++)
	{
	  if (!CL->master[x])continue;
	  for ( y=0; y< S->nseq; y++)
	    {
	      if (!do_mirror && set[y][x])continue;
	      set[x][y]=1;

	      id=(DM && DM->similarity_matrix)?DM->similarity_matrix[x][y]:-1;

	      if ( x==y && !do_self);
	      else if ( id!=-1 && !is_in_range(id,method->minid, method->maxid));
	      else
		{
		  sprintf (bufS, "2 %d %d",x,y);
		  bufA=make_aln_command (method,in=vtmpnam(NULL),out=vtmpnam (NULL));

		  if (!debugio && strrchr(bufA, '>')==NULL)strcat (bufA, TO_NULL_DEVICE);
		  if (check_seq_type (method, bufS, S))
		    {
		      job->c=print_lib_job (job->c, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ",method,fname,bufA, bufS, in, out);
		      job=queue_cat (job, job->c);
		    }
		  else if ( method->seq_type[0]=='P' && hijack_P_jobs)
		    {
		      //Hijack _P_ jobs without enough templates
		      static TC_method *proba_pairM;

		      add_information(stderr, "Method %s cannot be applied to [%s vs %s], proba_pair will be used instead",(method->executable2)?method->executable2:method->executable, (CL->S)->name[x], (CL->S)->name [y]);
		      if (!proba_pairM)
			{
			  proba_pairM=method_file2TC_method(method_name2method_file ("proba_pair"));
			  proba_pairM->PW_CL=method2pw_cl(proba_pairM, CL);
			}
		      job->c=print_lib_job (job->c, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ",proba_pairM,fname,bufA, bufS, in, out);
		      job=queue_cat (job, job->c);
		    }

		  vfree (bufA);
		}
	    }
	}
      free_int (set, -1);
    }

  return job;
}

int check_seq_type (TC_method *M, char *list,Sequence *S)
{
	char t1, t2;
	int s1, s2, n1, nseq, ntype, i;
	int *slist;
	Template *T1, *T2;


	slist=string2num_list (list);

	nseq=slist[1];
	ntype=strlen (M->seq_type);
	t1=M->seq_type[0];
	t2=M->seq_type[1];
	n1=0;

	/*Profiles and Sequences MUST NOT be distinguished so that sequences and profiles can easily be aligned*/
	if ( tolower(t1)=='r')t1='S';
	if ( tolower(t2)=='r')t2='S';


	if ( strm ( M->aln_mode, "pairwise") && nseq>2)n1=0;
	else if (strm ( M->aln_mode, "multiple" ) && ntype>1)n1=0;
	else if (ntype==1)
	{

		for (n1=0, i=0; i<nseq; i++)
		{
			s1=slist[i+2];
			T1=S->T[s1];

			n1+=(strchr (T1->seq_type,t1) || check_profile_seq_type (S, s1, t1))?1:0;
		}
		n1=(n1==nseq)?1:0;
	}
	else if (ntype==2)
	{
		int s1_has_t1;
		int s1_has_t2;
		int s2_has_t1;
		int s2_has_t2;

		s1=slist[2];
		s2=slist[3];
		T1=(S->T[s1]);
		T2=(S->T[s2]);

		s1_has_t1=(strchr ( T1->seq_type, t1) || check_profile_seq_type (S, s1, t1))?1:0;
		s1_has_t2=(strchr ( T1->seq_type, t2) || check_profile_seq_type (S, s1, t2))?1:0;
		s2_has_t1=(strchr ( T2->seq_type, t1) || check_profile_seq_type (S, s2, t1))?1:0;
		s2_has_t2=(strchr ( T2->seq_type, t2) || check_profile_seq_type (S, s2, t2))?1:0;
		n1=((s1_has_t1 && s2_has_t2) || (s1_has_t2 && s2_has_t1))?1:0;
	}

	vfree (slist);
	return n1;
}

int check_profile_seq_type (Sequence *S, int i, char t)
{
	Alignment *A;
	Template *T;
	int a;

	/*returns 1 if the sequence S is associated with a profile containing the right sequence*/
	A=seq2R_template_profile (S, i);
	if (A==NULL || A->S==NULL) return 0;
	for ( a=0; a< A->nseq; a++)
	{

		T=(A->S)->T[a];
		if ( T && strchr( T->seq_type,t))return 1;

	}
	return 0;
}


char **method_list2method4dna_list ( char **list, int n)
{
	int a;
	static char *buf;

	if ( !buf)buf=(char*)vcalloc ( 1000, sizeof (char));

	if ( !list || n==0)return list;
	buf=(char*)vcalloc ( 1000, sizeof (char));

	for ( a=0; a< n; a++)
	{

		sprintf ( buf,"%s",list[a]);

		if ( strm ( list[a], "4dna"));
		else
		{

			char **para ;
			int b;

			para=string2list2 (list[a], "@");
			sprintf ( buf, "%s4dna", para[1]);
			if (method_name2method_file (buf) || method_name2method_file(buf+1))
			{
				sprintf ( list[a],"%s", buf);
				for (b=2; b< atoi (para[0]); b++)
				{
					strcat (list[a], "@");
					strcat (list[a], para[b]);
				}
			}

			free_char (para, -1);
		}
	}
	return list;
}

int is_in_pre_set_method_list ( char *method)
{
	char *new_name;




	new_name=method_name2method_file (method);

	if ( !new_name) return 0;
	else
	{

		sprintf ( method, "%s", new_name);
		return 1;
	}
}
char *** display_method_names (char *mode, FILE *fp)
{
	char ***list, ***l2;
	int a, ml1=0, ml2=0, ml3=0;

	list=produce_method_file (NULL);
	l2=(char***)declare_arrayN(3,sizeof (char), 1000, 10, 100);

	fprintf ( fp, "\n#######   Compiling the list of available methods ... (will take a few seconds)\n");
	a=0;

	while (list[a])
	{


		sprintf (l2[a][0], "%s", method_file_tag2value (list[a][1],"ADDRESS"));
		sprintf (l2[a][1], "%s", method_file_tag2value (list[a][1],"PROGRAM"));
		sprintf (l2[a][2], "%s", method_file_tag2value (list[a][1],"ALN_MODE"));
		sprintf (l2[a][3], "%s", method_file_tag2value (list[a][1],"SEQ_TYPE"));
		ml1=MAX((strlen (list[a][0])),ml1);
		ml2=MAX((strlen (l2  [a][0])),ml2);
		ml3=MAX((strlen (l2  [a][1])),ml3);
		l2[a][4][0]= check_program_is_installed (l2[a][1],NULL,NULL,l2[a][0],NO_REPORT);
		if (strm (l2[a][0], "built_in")){sprintf (l2[a][5], "built_in");}
		else if (l2[a][4][0])l2[a][5]=pg2path (l2[a][1]);

		a++;
	}
	fprintf ( fp, "\n#######   Methods For which an Interface is available in T-Coffee\n");
	fprintf ( fp, "You must install the packages yourself when required (use the provided address)\n");
	fprintf ( fp, "Contact us if you need an extra method to be added [%s]\n", EMAIL);


	fprintf ( fp, "\n****** Pairwise Sequence Alignment Methods:\n");
	fprintf ( fp, "--------------------------------------------\n");
	a=0;
	while (list[a])
	{
		if ( strm (l2[a][2], "pairwise") && !strstr (l2[a][3], "P"))
			fprintf ( fp, "%-*s %-*s [pg: %*s is %s Installed][%s]\n", ml1,list[a][0],ml2, l2[a][0], ml3,l2[a][1], (l2[a][4][0]==0)?"NOT":"",(l2[a][5])?l2[a][5]:"");
		a++;
	}

	fprintf ( fp, "\n****** Pairwise Structural Alignment Methods:\n");
	fprintf ( fp, "--------------------------------------------\n");
	a=0;
	while (list[a])
	{
		if ( strm (l2[a][2], "pairwise") && strstr (l2[a][3], "P"))
			fprintf ( fp, "%-*s %-*s [pg: %*s is %s Installed][%s]\n", ml1,list[a][0],ml2, l2[a][0], ml3,l2[a][1], (l2[a][4][0]==0)?"NOT":"",(l2[a][5])?l2[a][5]:"");
		a++;
	}
	fprintf ( fp, "\n****** Multiple Sequence Alignment Methods:\n");
	fprintf ( fp, "--------------------------------------------\n");
	a=0;
	while (list[a])
	{
		if ( strm (l2[a][2], "multiple"))
			fprintf ( fp, "%-*s %-*s [pg: %*s is %s Installed][%s]\n", ml1,list[a][0],ml2, l2[a][0], ml3,l2[a][1], (l2[a][4][0]==0)?"NOT":"",(l2[a][5])?l2[a][5]:"");
		a++;
	}
	fprintf ( fp, "\n#######   Prediction Methods available to generate Templates\n");
	fprintf ( fp, "-------------------------------------------------------------\n");
	a=0;
	while (list[a])
	{
		if ( strm (l2[a][2], "predict"))
			fprintf ( fp, "%-*s %-*s [pg: %*s is %s Installed][%s]\n", ml1,list[a][0],ml2, l2[a][0], ml3,l2[a][1], (l2[a][4][0]==0)?"NOT":"",(l2[a][5])?l2[a][5]:"");
		a++;
	}
	fprintf ( fp, "\n\n\nAll these Methods are supported by T-Coffee, but you HAVE to install them yourself [use the provided address]\n\n");
	fprintf ( fp, "\nThese methods were selected because they are freeware opensource, easy to install and well supported");
	fprintf ( fp, "\nContact us if you need an extra method to be added [%s]\n", EMAIL);

	return l2;
}

char* method_name2method_file (char *method)
{
	char *fname=NULL;
	char ***mlist, *p;
	char address[100];
	char program[100];
	int a;

	if ( check_file_exists (method) || (toupper(method[0])==method[0] && check_file_exists (method+1)))return NULL;

	if ( (p=strstr (method, "@"))!=NULL && !strstr (method, "em@"))p[0]='\0';
	mlist=produce_method_file (method);
	


	a=0;
	while (mlist[a])
	{
		if ( lstrstr (method,mlist[a][0])){fname=mlist[a][1];break;}
		else {a++;}
	}
	if (p)p[0]='@';
	if ( fname==NULL) return NULL;
	else
	{
		sprintf (address, "%s", method_file_tag2value (fname,"ADDRESS"));
		sprintf (program, "%s", method_file_tag2value (fname,"PROGRAM"));
		check_program_is_installed (program,NULL,NULL,address,INSTALL_OR_DIE);
		if ( (method=strstr (method, "EP@"))!=NULL)
		{
			int a;
			char **list;
			FILE *fp;
			list=string2list2 ( method, "@");
			fp=vfopen (fname, "a");
			for ( a=2; a<atoi (list[0]); a+=2)
			{
				fprintf (fp, "%s %s\n", list[a], list[a+1]);

			}

			vfclose (fp);
			free_char (list, -1);
		}
		return fname;
	}
}
char *method_file_tag2value (char *method, char *tag)
{
	FILE *fp;
	char c, *p;
	static char *buf, *line;
	if (!buf) buf=(char*)vcalloc ( 1000, sizeof (char));
	else buf[0]='\0';

	if (!line)
	{
		line=(char*)vcalloc (LONG_STRING+1, sizeof ( char));
	}

	fp=vfopen (method, "r");
	while ( (c=fgetc (fp))!=EOF)
	{
		ungetc ( c, fp);
		fgets ( line,LONG_STRING, fp);

		if ( (line && (line[0]=='*' || line[0]=='#' || line[0] == '$'|| line[0]==' ' || line[0]=='\0' )));
		else if ( (p=strstr (line, tag    )))

		{
			sscanf (p+strlen (tag), " %s"   , buf);
			vfclose (fp);
			return buf;
		}
	}
	vfclose ( fp);
	return buf;
}
struct TC_method * method_file2TC_method( char *method)
{
	TC_method *m;
	static char *subcommand;
	static char *line;
	char *p;
	FILE *fp;
	int c;
	//<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PARAM>
	if (!line)
	{
		line=(char*)vcalloc (LONG_STRING+1, sizeof ( char));
		subcommand=(char*)vcalloc ( LONG_STRING, sizeof (char));
	}

	m=(TC_method*)vcalloc ( 1, sizeof (TC_method));

	/*set default parameter values*/
	m->gop=m->gep=UNDEFINED;
	sprintf (m->seq_type, "S");
	sprintf (m->weight, "sim");
	m->minid=0;
	m->maxid=100;

	fp=vfopen (method, "r");
	
	while ( (c=fgetc (fp))!=EOF)
	{
		ungetc ( c, fp);
		fgets ( line,LONG_STRING, fp);
		

		line=substitute (line, "\n", " ");
		line=substitute (line, "%s", " ");
		line=substitute (line, "%e", "=");
		line=substitute (line, "%m", "-");

		if ( (line && (line[0]=='*' || line[0]=='#' || line[0] == '$'|| line[0]==' ' || line[0]=='\0' )))subcommand[0]='\0';
		//Parse PARAM, PARM1 and PARAM2 first because they may contain keywords
		else if ( (p=strstr (line, "PARAM1"     )))
		{
			sprintf (subcommand, " %s ", p+6);
			strcat ( m->param1, subcommand);
		}
		else if ( (p=strstr (line, "PARAM2"     )))
		{
			sprintf (subcommand, " %s ", p+6);
			strcat ( m->param2, subcommand);
		}
		else if ( (p=strstr (line, "PARAM"     )))
		{
			sprintf (subcommand, " %s ", p+5);
			strcat ( m->param, subcommand);
		}else if ( (p=strstr (line, "EXECUTABLE2" )))
		{
			sscanf (p, "EXECUTABLE2 %s", m->executable2);
		}
		else if ( (p=strstr (line, "EXECUTABLE" )))
		{
			sscanf (p, "EXECUTABLE %s", m->executable);
		}
		else if ( (p=strstr (line, "METHOD"     ))) sscanf (p, "METHOD %s"   , m->method);
		else if ( (p=strstr (line, "ALN_MODE"   ))) sscanf (p, "ALN_MODE %s"  , m->aln_mode);
		else if ( (p=strstr (line, "PRFMODE"    ))) sscanf (p, "PRFMODE %s"   , m->prfmode);
		else if ( (p=strstr (line, "IN_FLAG2"   ))) sscanf (p, "IN_FLAG2 %s"   , m->in_flag2);
		else if ( (p=strstr (line, "IN_FLAG"    ))) sscanf (p, "IN_FLAG %s"   , m->in_flag);
		else if ( (p=strstr (line, "OUT_FLAG"   ))) sscanf (p, "OUT_FLAG %s"  , m->out_flag);
		else if ( (p=strstr (line, "OUT_MODE"   ))) sscanf (p, "OUT_MODE %s"  , m->out_mode);
		else if ( (p=strstr (line, "SEQ_TYPE"   ))) sscanf (p, "SEQ_TYPE %s"  , m->seq_type);
		else if ( (p=strstr (line, "WEIGHT"     ))) sscanf (p, "WEIGHT %s"  , m->weight);
		else if ( (p=strstr (line, "MATRIX"     ))) sscanf (p, "MATRIX %s"  , m->matrix);
		else if ( (p=strstr (line, "GOP"        ))) sscanf (p, "GOP %d"  , &m->gop);
		else if ( (p=strstr (line, "GEP"        ))) sscanf (p, "GEP %d"  , &m->gep);
		else if ( (p=strstr (line, "MAXID"      ))) sscanf (p, "MAXID %d"  , &m->maxid);
		else if ( (p=strstr (line, "MINID"      ))) sscanf (p, "MINID %d"  , &m->minid);
		else if ( (p=strstr (line, "EXTEND_SEQ"      ))) sscanf (p, "EXTEND_SEQ %d"  , &m->extend_seq);
		else if ( (p=strstr (line, "REVERSE_SEQ"      ))) sscanf (p, "REVERSE_SEQ %d"  , &m->reverse_seq);

		

	}
	
	
	vfclose ( fp);
	if ( !m-> executable && ! m->executable2)return NULL;
	return m;
}

int TC_method2method_file( struct TC_method*m,char *fname )
{
	FILE *fp;
	if ( !m) return 0;
	fp=vfopen ( fname, "w");
	if ( m->executable[0])fprintf (fp, "EXECUTABLE %s\n", m->executable);
	if (m->in_flag[0])fprintf (fp, "IN_FLAG %s\n", m->in_flag);
	if (m->out_flag[0])fprintf (fp, "OUT_FLAG %s\n", m->out_flag);
	if (m->out_mode[0])fprintf (fp, "OUT_MODE %s\n", m->out_mode);
	if (m->aln_mode[0])fprintf (fp, "ALN_MODE %s\n", m->aln_mode);
	if (m->seq_type)fprintf (fp, "SEQ_TYPE %s\n", m->seq_type);
	if (m->weight[0])fprintf (fp, "WEIGHT %s\n", m->weight);
	if (m->matrix[0])fprintf (fp, "MATRIX %s\n", m->matrix);
	if (m->gop!=UNDEFINED)fprintf (fp, "GOP %d\n", m->gop);
	if (m->gep!=UNDEFINED)fprintf (fp, "GEP %d\n", m->gep);
	if (m->minid!=0  )fprintf (fp, "MINID %d\n", m->minid);
	if (m->maxid!=100)fprintf (fp, "MAXID %d\n", m->maxid);
	if (m->param[0])fprintf (fp, "PARAM %s\n", m->param);
	if (m->param1[0])fprintf (fp, "PARAM1 %s\n", m->param1);
	if (m->param2[0])fprintf (fp, "PARAM1 %s\n", m->param2);
	if (m->in_flag2[0])fprintf (fp, "IN_FLAG2 %s\n", m->in_flag2);
	if (m->prfmode[0])fprintf (fp, "PRFMODE %s\n", m->prfmode);
	if (m->method[0])fprintf (fp, "METHOD %s\n", m->method);
	vfclose ( fp);
	return 1;
}

char *make_aln_command(TC_method *m, char *seq, char *aln)
{
	char *command;
	char buf[1000];

	//      sprintf ( buf, "%s %s %s%s %s%s %s", m->executable, m->param1, m->in_flag, seq,m->param2, m->out_flag,aln, m->param);

	sprintf ( buf, "%s %s %s%s %s %s%s %s", m->executable, m->param1, m->in_flag, seq,m->param2, m->out_flag,aln, m->param);
	command=(char*)vcalloc ( strlen (buf)+100, sizeof (char));
	sprintf ( command, "%s", buf);
	

	command=substitute (command, "&bnsp", " ");
	command=substitute (command, "no_name", "");

	return command;
}





/*********************************************************************/
/*                                                                   */
/*                         WRITE IN LIST                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *  unfreeze_constraint_list (Constraint_list *CL)
{
	free_int (CL->freeze, -1);
	CL->freeze=NULL;
}
Constraint_list *  freeze_constraint_list (Constraint_list *CL)
{
	int a, b, d=0;
	Sequence *S=CL->S;
	int **freeze=declare_int2 (S->nseq, S->len, 1);

	for ( a=0; a<S->nseq; a++)
	  {
	    b=1;
	    while (CL->residue_index[a][b])
	      {
		  d+=CL->residue_index[a][b][0];
		  freeze[a][b]=CL->residue_index[a][b][0];
		  b++;
	      }
	  }
	CL->freeze=freeze;
	return CL;
}
Constraint_list *  empty_constraint_list (Constraint_list *CL)
{
  //reset all the indexes
  int a, b;

  if ( !CL || ! CL->residue_index) return CL;
  for (a=0; a<(CL->S)->nseq; a++)
    {
      b=0;
      while (CL->residue_index[a][b])
	{
	  CL->residue_index[a][b][0]=1;
	  b++;
	}
    }
  CL->ne=0;
  return CL;
}

Constraint_list *  undump_constraint_list (Constraint_list *CL, char *file)
{

	int *entry, b, c, e, tot;
	FILE *fp;
	
	if (!CL || !CL->residue_index)return CL;
	entry=(int*)vcalloc ( CL->entry_len+1, sizeof (int));
	
	
	fp=vfopen (file, "rb");
	while ((b=fread (entry, sizeof(int), CL->entry_len, fp))==CL->entry_len)
	  {
	    CL=add_entry2list2 (entry, CL);
	  }

	vfree(entry);
	vfclose (fp);
	remove(file);
	return CL;
}

int safe_dump_constraint_list (Constraint_list *CL,char *file, char *mode, Sequence *RS)
{
	Sequence *S;
	int **cache=NULL;
	FILE *fp;

	int *entry =(int*)vcalloc (CL->entry_len+1, sizeof (int));
	int b,c,s1, r1, s2, r2;
	int d=0;
	char bmode[10];
	if (!CL || !CL->S || !CL->residue_index || CL->ne==0) return 0;
	S=CL->S;
	if (RS)cache=fix_seq_seq (S, RS);

	if      (mode[0]=='a')sprintf (bmode, "ab");
	else if (mode[0]=='w')sprintf (bmode, "wb");
	else sprintf (bmode, "%s", mode);

	fp=vfopen (file, bmode);
	for (s1=0; s1<S->nseq; s1++)
	  {
	    if (cache && cache[s1][0]==-1)continue;
	    for (r1=1; r1<=S->len[s1]; r1++)
	      {
		entry[SEQ1]=(cache)?cache[s1][0]:s1;
		entry[R1]=(cache)?cache[s1][r1]:r1;
		if (entry[R1]<=0)continue;
		
		b=(CL->freeze)?CL->freeze[s1][r1]:1;
		for (;b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
		  {
		    s2=CL->residue_index[s1][r1][b+SEQ2];
		    r2=CL->residue_index[s1][r1][b+R2];
		    
		    entry[SEQ2]=(cache)?cache[s2][0]:s2;
		    if ( entry[SEQ2]==-1)continue;
		    else
		      {
			entry[R2]=(cache)?cache[s2][r2]:r2;
			if (entry[R2]<=0)continue;
			else
			  {
			    d++;
			    entry[WE]=CL->residue_index[s1][r1][b+WE];
			    entry[CONS]=CL->residue_index[s1][r1][b+CONS];
			    entry[MISC]=CL->residue_index[s1][r1][b+MISC];
			    fwrite(entry, sizeof (int),CL->entry_len,fp);
			  }
		      }
		  }
	      }
	  }
	vfclose (fp);
	
	free_int (cache, -1);
	vfree (entry);
	return d;
}


int dump_constraint_list (Constraint_list *CL, char *file, char *mode)
{
	return safe_dump_constraint_list (CL, file, mode, NULL);
}


FILE* display_constraint_list (Constraint_list *CL, FILE *fp, char *tag)
{
	Sequence *S=CL->S;
	int b, s1, r1, s2, r2, w2, n;

	n=0;
	for (s1=0; s1<S->nseq; s1++)
	  {
	    fprintf (fp, "SEQUENCE %d\n", s1);

	    r1=1;
	    while (CL->residue_index[s1][r1])
	      {
		for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
		  {
		    s2=CL->residue_index[s1][r1][b+SEQ2];
		    r2=CL->residue_index[s1][r1][b+R2];
		    w2=CL->residue_index[s1][r1][b+WE];
		    fprintf ( fp, "\t%sS1:%5d - R1:%5d S2:%5d R2:%5d W:%5d\n", tag,s1,r1, s2, r2,w2);
		  }
		r1++;
	      }

	  }
	return fp;
}



/*********************************************************************/
/*                                                                   */
/*                         LIST EXTENTION                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

/**
 * Determine maximal constraints for normalisation.
 *
 * Sets Constaint_list::max_value to the largest edge weight occuring in the complete Constraint_list,
 * and Constraint_list::max_ext_value to the largest sum of outgoing edges of all residues.
 */
Constraint_list * evaluate_constraint_list_reference ( Constraint_list *CL)
{
	int a, b, s1, r1;
	int **max_res;

	if ( CL->M)
	  {
	    CL->max_value=CL->max_ext_value=20;

	  }
	else
	  {
	    Sequence *S=CL->S;


	    CL->max_value=CL->max_ext_value=0;
	    max_res=(int**)vcalloc ( (CL->S)->nseq, sizeof (int*));

	    for ( a=0; a< (CL->S)->nseq; a++)
	      {
		max_res[a]=(int*)vcalloc ( strlen ((CL->S)->seq[a])+1, sizeof (int));
	      }


	    for (s1=0; s1<S->nseq; s1++)
	      {
		for ( r1=1; r1<=S->len[s1]; r1++)
		  {

		    for (a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
		      {
			int s2=CL->residue_index[s1][r1][a+SEQ2];
			int r2=CL->residue_index[s1][r1][a+R2];
			int w2=CL->residue_index[s1][r1][a+WE];


			if ( w2==UNDEFINED || ( (CL->moca) && (CL->moca)->forbiden_residues  && ((CL->moca)->forbiden_residues[s1][r1]==UNDEFINED || (CL->moca)->forbiden_residues[s2][r2]==UNDEFINED)));
			else
			  {
			    max_res[s1][r1]+=w2;
			    max_res[s2][r2]+=w2;
			    CL->max_value=MAX(w2, CL->max_value);
			  }
		      }
		  }
	      }

	    for ( a=0; a< (CL->S)->nseq; a++)
	      for ( b=1; b<=(CL->S)->len[a]; b++)
		{
		  CL->max_ext_value=MAX(max_res[a][b],CL->max_ext_value);
		}

	    free_int (max_res,-1);
	    CL->max_ext_value=MAX(1,CL->max_ext_value);
	  }

	if (CL->normalise)
	  {
	    CL->nomatch=(CL->nomatch*CL->normalise)/CL->max_ext_value;
	  }

	return CL;
}

/*********************************************************************/
/*                                                                   */
/*                         ENTRY MANIPULATION                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list * add_list_entry2list (Constraint_list *CL, int n_para, ...)
{
	int a;
	int *entry;
	int field, val;
	va_list ap;

	if (n_para>LIST_N_FIELDS)
	{
		myexit (fprintf_error ( stderr, "Too Many Fields in List"));
	}

	va_start (ap,n_para);
	entry=(int*)vcalloc (CL->entry_len+1, sizeof (int));

	for ( a=0; a<n_para; a++)
	{
		field=va_arg(ap, int);
		val  =va_arg(ap, CLIST_TYPE);
		entry[field]=val;
	}
	va_end (ap);
	add_entry2list(entry, CL);
	vfree(entry);
	return CL;
}

int *extract_entry (Constraint_list *CL)
{
	static int s=0;
	static int r=0;
	static int l=1;
	static int *entry;

	if (!entry)entry=(int*)vcalloc (100, sizeof (int));
	if (!CL){s=0;r=1;l=1; return NULL;}
	

	for (; s<(CL->S)->nseq; s++)
	  {
	    while ((CL->residue_index[s][r]))
	      {
		for (; l<CL->residue_index[s][r][0];)
		  {

		    entry[SEQ1]=s;
		    entry[R1]=r;
		    entry[SEQ2]=CL->residue_index[s][r][l+SEQ2];
		    entry[R2]=  CL->residue_index[s][r][l+R2];
		    entry[WE]=  CL->residue_index[s][r][l+WE];
		    entry[CONS]=CL->residue_index[s][r][l+CONS];
		    entry[MISC]=CL->residue_index[s][r][l+MISC];
		    entry[INDEX]=l;
		    l+=ICHUNK;
		    return entry;
		  }
		l=1;
		r++;
	      }
	    r=0;
	  }
	s=0;
	r=0;
	l=1;

	return NULL;
}
#ifdef FAFAFA
int next_entry (Constraint_list *CL, int *s, int*r, int *l);
int *extract_entry (Constraint_list *CL)
{
	static int s=0;
	static int r=1;
	static int l=1;
	static int *entry;
	int v;

	if (!entry)entry=vcalloc (100, sizeof (int));
	if (!CL){s=0;r=1;l=1; return NULL;}


	while ((v=next_entry (CL, &s, &r, &l)))
	  {
	    if (v==1)

		{

			entry[SEQ1]=s;
			entry[R1]=r;
			entry[SEQ2]=CL->residue_index[s][r][l-ICHUNK+SEQ2];
			entry[R2]=  CL->residue_index[s][r][l-ICHUNK+R2];
			entry[WE]=  CL->residue_index[s][r][l-ICHUNK+WE];
			entry[CONS]=CL->residue_index[s][r][l-ICHUNK+CONS];
			entry[MISC]=CL->residue_index[s][r][l-ICHUNK+MISC];
			entry[INDEX]=l-ICHUNK;
			return entry;
		}
	  }
	s=0; r=1; l=1;
	return NULL;
}
int next_entry (Constraint_list *CL, int *s, int*r,int *l)
{
	Sequence *S=CL->S;

	if (s[0]>=S->nseq)return 0;
	else if (!CL->residue_index[s[0]][r[0]])
	{
		r[0]=1;
		l[0]=1;
		s[0]++;
		return -1;
	}
	else if ( l[0]>=CL->residue_index[s[0]][r[0]][0])
	{
		r[0]++;
		l[0]=1;
		return -1;
	}
	else
	{

		l[0]+=ICHUNK;
		return 1;
	}
}
#endif

int CLisCompacted (Constraint_list *CL, char *t)
{
	int s1, r1, s2, r2,ps2,pr2,b,c;
	Sequence *S=CL->S;

	for (s1=0; s1<S->nseq; s1++)
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][b+SEQ2];
				r2=CL->residue_index[s1][r1][b+R2];

				if (b>1)
				{
					if (s2==ps2 && r2==pr2)
					{

					  fprintf (stderr,"%s -- NOT COMPACTED",t);
					  exit(0);
					}
					else if ( s2<=ps2 && r2<pr2)
					{
					  fprintf (stderr, "%s -- NOT SORTED",t);
					  fprintf ( stderr, "\n");
					  for (c=1; c<CL->residue_index[s1][r1][0]; c++)
					    {
					      fprintf (stderr, "%5d ",CL->residue_index[s1][r1][c]);
					    }
					  
					  exit(0);
					}
				}
				ps2=s2;
				pr2=r2;
			}
		}
		return 1;
}

int checkCL (Constraint_list *CL, char *t)
{
	int fail=0;
	int a, b;

	for (a=0; a<(CL->S)->nseq; a++)
		for (b=1; b<=(CL->S)->len[a]; b++)
		  {
		    HERE ("%s %4d ==> %d",t,CL->residue_index[a][b][0], CL->residue_index[a][b][0]%ICHUNK);
		  }
	
	
	if (!CL){HERE ("CL is Not declared");fail=1;}
		else if ( !CL->S){HERE ("S is not declared in CL");fail=1;}
		else if ( !CL->residue_index){HERE ("residue index is not declared");fail=1;}
		else if (read_array_size_new (CL->residue_index)!=(CL->S)->nseq){HERE ("CL not well declared (S)");fail=1;}
		else
		{
			for (a=0; a<(CL->S)->nseq; a++)
			{
				int s=(CL->S)->len[a]+1;
				int j=read_array_size_new (CL->residue_index[a]);
				if (s!=j)fail=1;

			}
			if (fail)
			{
				for (a=0; a<(CL->S)->nseq; a++)
				{
					int s=(CL->S)->len[a]+1;
					int j=read_array_size_new (CL->residue_index[a]);
					if (s!=j)fail=1;
					HERE ("\t %s %d %d",(CL->S)->name[a], s-1,j);
				}
			}

			if ( fail)
			{
			  fprintf ( stderr,"******** CHECKED CL: %s ******", t);
				exit (0);
			}
		}
		return 1;
}

Constraint_list *add_entry2list( CLIST_TYPE *entry, Constraint_list *CL)
{
  //adds an entry and its mirror to the list
  //if INDEX is set the entry replaces the entry with a similar index
  //otherwise the entry (and its mirror) are added
  int s1=entry[SEQ1];
  int s2=entry[SEQ2];
  int r1=entry[R1];
  int r2=entry[R2];

  if (entry[INDEX])return add_entry2list2(entry, CL);

  entry[SEQ1]=s2;
  entry[SEQ2]=s1;
  entry[R1]=r2;
  entry[R2]=r1;
  
  add_entry2list2 (entry, CL);

  entry[SEQ1]=s1;
  entry[SEQ2]=s2;
  entry[R1]=r1;
  entry[R2]=r2;
  add_entry2list2 (entry, CL);
  return CL;
}


int get_entry_index (int *entry, Constraint_list *CL);
Constraint_list *insert_entry (CLIST_TYPE *entry, Constraint_list *CL, int i);
Constraint_list *update_entry (CLIST_TYPE *entry, Constraint_list *CL, int i);
Constraint_list *reset_entry (CLIST_TYPE *entry, Constraint_list *CL, int i);
Constraint_list *remove_entry (CLIST_TYPE *entry, Constraint_list *CL, int i);
Constraint_list *add_entry2list2 (CLIST_TYPE *entry, Constraint_list *CL)
{
  //adds an entry to the list
  int s2=entry[SEQ2];
  int r2=entry[R2];
  int s1=entry[SEQ1];
  int r1=entry[R1];
  int i;
  
  if (r1>(CL->S)->len[s1])
    {
      myexit (fprintf_error ( stderr, "Library out of bounds: %s::%d [%d]* vs %s::%d", (CL->S)->name[s1], r1,(CL->S)->len[s1], (CL->S)->name[s2], r2));
    }
  else if (r2>(CL->S)->len[s2])
    {
      myexit (fprintf_error ( stderr, "Library out of bounds: %s::%d vs %s::%d [%d]*", (CL->S)->name[s1], r1, (CL->S)->name[s2], r2,(CL->S)->len[s2] ));
    }

  
  i=CL->residue_index[entry[SEQ1]][entry[R1]][0];
  if (entry[INDEX]){return reset_entry (entry, CL, entry[INDEX]); }
  else if (i==1){CL=insert_entry (entry, CL,1);}
  else
    {
      i=get_entry_index(entry,CL);
      if (i<0)insert_entry (entry,CL,-i);
      else update_entry (entry, CL, i);
    }
  
  return CL;
}


int get_entry_index_serial (int *entry, Constraint_list *CL);
int get_entry_index_dico (int *entry, Constraint_list *CL);
int get_entry_index (int *entry, Constraint_list *CL)
{
  return get_entry_index_dico (entry, CL);
}
int get_entry_index_serial (int *entry, Constraint_list *CL)
{
  //return the index of an entry
  //positive value: the entry exists on position i
  //negative value: the entry must be created on poistion -i
  int s1=entry[SEQ1];
  int r1=entry[R1];
  int s2=entry[SEQ2];
  int r2=entry[R2];
  int a;
  int *r= CL->residue_index[s1][r1];

  static int tot;
  static int pr;
  if (r[0]==1)return -1;//corresponding entry undeclared->must be inserted
  else
    {
      for (a=1; a<r[0]; a+=ICHUNK)
	{
	  tot++;
	  pr++;
	  if (r[a+SEQ2]==s2 && r[a+R2]==r2)return a;
	  else if (r[a+SEQ2]==s2 && r[a+R2]>r2)return -a;
	  else if (r[a+SEQ2]>s2)return -a;
	}
    }
  return -a;
}

int get_entry_index_dico (int *entry, Constraint_list *CL)
{
  //return the index of an entry
  //positive value: the entry exists on position i
  //negative value: the entry must be created on poistion -i
  int s1=entry[SEQ1];
  int r1=entry[R1];
  int s2=entry[SEQ2];
  int r2=entry[R2];
  int dir;
  int delta;
  int i;
  int *r= CL->residue_index[s1][r1];
  int ps2, pr2, ns2,nr2,p,n;
  static int tot;
  static int pr;

  if (r[0]==1)return -1;//corresponding entry undeclared->must be inserted
  dir=1;
  i=1;
  delta=((r[0]-1)/ICHUNK);
  delta=MAX(1,(delta/2));

  while (1==1)
    {
      pr++;
      tot++;
      i+=(delta*dir*ICHUNK);
      i=MAX(i,1);
      i=MIN(i,(r[0]));
      if (i<r[0] && s2==r[i+SEQ2] && r2==r[i+R2])return i;
      else
	{
	  p=1; n=1;
	  if (i>1)
	    {
	      ps2=r[i-ICHUNK+SEQ2];
	      pr2=r[i-ICHUNK+R2];
	      p=(s2>ps2 || (s2==ps2 && r2>pr2))?1:0;
	    }
	  if (i<r[0])
	    {
	      ns2=r[i+SEQ2];
	      nr2=r[i+R2];
	      n=(s2<ns2 || (s2==ns2 && r2<nr2))?1:0;
	    }

	  if ( p && n) return -i;
	  else if (p)
	    {
	      if (dir==-1)delta=MAX(1,(delta/2));
	      dir=1;
	    }
	  else if (n)
	    {
	      if ( dir==1)delta=MAX(1,(delta/2));
	      dir=-1;
	    }

	}
    }
}

Constraint_list *remove_entry (CLIST_TYPE *entry, Constraint_list *CL, int i)
{
	//insert entry right after i;
	static char *buf;
	static int bsize;
	int *r=CL->residue_index[entry[SEQ1]][entry[R1]];
	int s;

	i+=ICHUNK;
	s=r[0]-i;
	if (bsize<s){vfree(buf); bsize=s;buf=(char*)vcalloc (bsize,sizeof (int));}
	memcpy(buf,r+i, s*sizeof (int));
	memcpy(r+i-ICHUNK, buf, s*sizeof (int));
	r[0]-=ICHUNK;
	CL->residue_index[entry[SEQ1]][entry[R1]]=(int*)vrealloc(r,r[0]*sizeof (int));
	CL->ne--;
	return CL;
}

Constraint_list *insert_entry (CLIST_TYPE *entry, Constraint_list *CL, int i)
{
	//insert entry right after i;
	static char *buf;
	static int bsize;
	int *r=CL->residue_index[entry[SEQ1]][entry[R1]];
	int s;
	//inserts a new entry between i-1 and i

	
	r=(int*)vrealloc (r, (r[0]+ICHUNK)*sizeof (int));
	CL->residue_index[entry[SEQ1]][entry[R1]]=r;
	CL->ne++;
	
	
	s=r[0]-i;
	if (bsize<s){vfree(buf); bsize=s;buf=(char*)vcalloc (bsize,sizeof (int));}
	if (s)memcpy(buf,r+i, s*sizeof (int));
	if (s)memcpy(r+i+ICHUNK, buf, s*sizeof (int));
	memcpy (r+i, entry, ICHUNK*sizeof (int));
	r[0]+=ICHUNK;
	return CL;
}

Constraint_list *reset_entry (CLIST_TYPE *entry, Constraint_list *CL, int i)
{
  int *r=CL->residue_index[entry[SEQ1]][entry[R1]];
  memcpy (r+i, entry, ICHUNK*sizeof (int));
  return CL;
}
Constraint_list *update_entry (CLIST_TYPE *entry, Constraint_list *CL, int i)
{
  int s1=entry[SEQ1];
  int r1=entry[R1];

  CL->residue_index[s1][r1][i+SEQ2]=entry[SEQ2];
  CL->residue_index[s1][r1][i+R2]  =entry[R2];
  CL->residue_index[s1][r1][i+WE]  =MAX(entry[WE],CL->residue_index[s1][r1][i+WE]);
  CL->residue_index[s1][r1][i+CONS]+=entry[CONS];
  CL->residue_index[s1][r1][i+MISC]=entry[MISC];
  return CL;
}

/*********************************************************************/
/*                                                                   */
/*                         SEARCH IN LIST (ARRAY AND FILE)           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

int compare_residues_between_clists (int s1, int r1,Constraint_list *CL1, Constraint_list *CL2,int *id, int *tot)
{
  static int *lu;
  static int nseq;
  int n, a, s2, r2;
  //Note: residues: 1->N
  //Note Sequences: 0->N-1 !
  if ((CL1->S)->nseq !=nseq)
    {
      if (lu) vfree (lu);
      lu=(int*)vcalloc ((CL1->S)->nseq, sizeof (int));
      nseq=(CL1->S)->nseq;
    }
  
  
  n=CL1->residue_index[s1][r1][0];
  for (a=1; a<n; a+=ICHUNK)
    {
      s2=CL1->residue_index[s1][r1][a+SEQ2];
      r2=CL1->residue_index[s1][r1][a+R2];
      lu[s2]=r2;
      tot[0]++;
    }
  n=CL2->residue_index[s1][r1][0];
  for (a=1; a<n; a+=ICHUNK)
    {
      s2=CL2->residue_index[s1][r1][a+SEQ2];
      r2=CL2->residue_index[s1][r1][a+R2];
      if (lu[s2]==r2)id[0]++;
    }

  n=CL1->residue_index[s1][r1][0];
  for (a=1; a<n; a+=ICHUNK)
    {
      s2=CL1->residue_index[s1][r1][a+SEQ2];
      r2=CL1->residue_index[s1][r1][a+R2];
      lu[s2]=0;
    }
  return id[0];
}
	
  

CLIST_TYPE *main_search_in_list_constraint ( int *key,int *p,int k_len,Constraint_list *CL)
{
	static CLIST_TYPE *l=NULL;
	int a, s1, s2, r1, r2, ni;

	if (!l)l=(int*)vcalloc (CL->entry_len+1, sizeof (int));
	for (a=0; a<CL->entry_len; a++)l[a]=key[a];

	l[INDEX]=1;

	s1=key[SEQ1];
	r1=key[R1];
	s2=key[SEQ2];
	r2=key[R2];
	ni=CL->residue_index[s1][r1][0];
	for (a=1; a<ni; a+=ICHUNK)
	{
		if (CL->residue_index[s1][r1][a+SEQ2]==s2 && CL->residue_index[s1][r1][a+R2]==r2)
		{
			l[SEQ1]=s1;
			l[R1]=r1;
			l[SEQ2]=s2;
			l[R2]=r2;
			l[WE]=CL->residue_index[s1][r1][a+WE];
			l[INDEX]=a;
			return l;
		}
	}
	return NULL;
}


/*********************************************************************/
/*                                                                   */
/*                                                                   */
/*      LIST SORTING                                                 */
/*                                                                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *sort_constraint_list_inv (Constraint_list *CL, int start, int len)
{
	return sort_constraint_list (CL, start, len);
}

Constraint_list *invert_constraint_list (Constraint_list *CL, int start,int len)
{
	return CL;
}

Constraint_list *sort_constraint_list(Constraint_list *CL, int start, int len)
{
	int s1, r1;
	Sequence *S=CL->S;
	
	for (s1=0; s1<S->nseq; s1++)
	{
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			int nE=(CL->residue_index[s1][r1][0]-1)/ICHUNK;
			int sizE=sizeof(int)*ICHUNK;
			int *start=CL->residue_index[s1][r1]+1;
			qsort ((void*)start,nE,sizE,compare_constraint_list_entry);
		}
	}
	return CL;
}

int compare_constraint_list_entry ( const void*vx, const void*vy)
{
	int a;
	const int *x=(int*)vx, *y=(int*)vy;
	for (a=0; a<ICHUNK; a++)
	{
		if (x[a]<y[a])return -1;
		else if (x[a]>y[a]) return 1;
	}
	return 0;
}
/*********************************************************************/
/*                                                                   */
/*                         LIST PARSING                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list*  fork_read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source,int nproc);

/**
 * Read several types of input files and build the library.
 *
 * Decides on how many processors to use and redirects to ::fork_read_n_constraint_list,
 * so please follow th elibnk to continue reading.
 */

Constraint_list*CL2simCL (Constraint_list *CL)
{
  int s1, s2,r1, r2, i,b;
  int *entry1, *entry2;
  Sequence *S;
  
  //CL MUST be symetrical
  //add missing symetrical entries to CL
  
  if (!CL || !CL->residue_index)return CL;
  S=CL->S;
  if (!S)return CL;
  
  entry1=(int*)vcalloc ( CL->entry_len+1, sizeof (int));
  entry2=(int*)vcalloc ( CL->entry_len+1, sizeof (int));
  
  for (s1=0; s1<S->nseq; s1++)
    {
      for (r1=1; r1<=S->len[s1]; r1++)
	{
	  entry1[SEQ1]=entry2[SEQ2]=s1;
	  entry1[R1]  =entry2[R2]  =r1;
	  if (entry1[R1]<=0)continue;
	  
	  for (b=1;b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
	    {
	      s2=CL->residue_index[s1][r1][b+SEQ2];
	      r2=CL->residue_index[s1][r1][b+R2];
	      
	      entry1[SEQ2]=entry2[SEQ1]=s2;
	      if ( entry1[SEQ2]==-1)continue;
	      else
		{
		  entry1[R2]=entry2[R1]=r2;
		  if (entry1[R2]<=0)continue;
		  else
		    {
		      entry1[WE]  =entry2[WE]   =CL->residue_index[s1][r1][b+WE];
		      entry1[CONS]=entry2[CONS] =CL->residue_index[s1][r1][b+CONS];
		      entry1[MISC]=entry2[MISC] =CL->residue_index[s1][r1][b+MISC];
		      CL=add_entry2list2 (entry1, CL);
		      CL=add_entry2list2 (entry2, CL);
		    }
		}
	    }
	}
    }
  vfree (entry1); vfree(entry2);
  return CL;
}
Constraint_list*reload_constraint_list (Constraint_list *CL)
{
  
  char *tmp=vtmpnam (NULL);
  if (!CL)return CL;
  CL->ne=0;
  CL->residue_index=declare_residue_index(CL->S);
  CL=undump_constraint_list (CL,tmp);
  return CL;
}

Constraint_list* read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source)
{

  
  if (1==1 && !strstr (CL->multi_thread, "methods"))//deprecated
    return fork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode,type,local_stderr, CL, seq_source, 1);
  else
    return fork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode,type,local_stderr, CL, seq_source, get_nproc());
}


/**
 * Distributes the reading process for each input file.
 *
 * This function forks the process (using ::vvfork) of reading several input files into
 * several subroutines. Each of these subroutines will eventually call ::read_constraint_list,
 * which is the central reading method and recommended reading!
 *

 * \param[in] fname       The names of the input files
 * \param[in] n_list      Number of input files
 * \param[in,out] CL      Points to the global ::Constraint_list object
 * \param[in] nproc       Number of cores to use, per default the output of ::get_nproc
 *
 */
Constraint_list* fork_read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL,char *seq_source, int nproc)
{
	int a, b;
	Sequence *S;
	char **tmp_list;
	int*proclist;
	int  ns;
	int max=max_n_pid();

	proclist=(int*)vcalloc (max, sizeof (int));
	tmp_list=(char**)vcalloc (n_list+1, sizeof (char*));
	for (a=0; a<n_list; a++)tmp_list[a]=vtmpnam(NULL);

	if (!(CL->S) && (S=read_seq_in_n_list (fname, n_list,type,seq_source))==NULL)
	{
		myexit (fprintf_error ( stderr, "NO SEQUENCE WAS SPECIFIED"));
	}
	else if (CL->S==NULL)
	{
		CL->S=S;
	}



	/*CHECK IF THERE IS A MATRIX AND GET RID OF OTHER METHODS*/
	for (b=0, a=0; a< n_list; a++)if (is_matrix(fname[a]) ||is_matrix(fname[a]+1) )b=a+1;

	if ( b)
	{
		if ( b==1);
		else sprintf ( fname[0], "%s", fname[b-1]);
		n_list=1;

		return read_constraint_list (CL, fname[0], in_mode, mem_mode,weight_mode);
// 		return fork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode, type,local_stderr, CL, seq_source,1);

	}

	if (!CL)CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?tmpfile():NULL, NULL);

	if (CL->ne)
	{
		dump_constraint_list(CL,tmp_list[n_list], "w");
		CL->ne=0;
	}

	CL->local_stderr=local_stderr;
	fprintf ( local_stderr, "\n\tMulti Core Mode (read): %d processor(s):\n", nproc);
	for (ns=0,a=0; a< n_list; a++)
	{
		int pid;
		ns++;
		pid=vvfork (NULL);
		if ( pid==0)
		{
			int in;
			initiate_vtmpnam (NULL);
			CL->local_stderr=vfopen("/dev/null", "w");
			in=CL->ne;
			CL=read_constraint_list (CL, fname[a], in_mode, mem_mode,weight_mode);
			if (CL->ne>in)dump_constraint_list(CL,tmp_list[a], "w");
			myexit (EXIT_SUCCESS);
		}
		else
		{
			//set_pid (pid);
			fprintf ( local_stderr, "\n\t--- Process Method/Library/Aln %s", fname[a], ns);
			proclist[pid]=a;
			if (ns>=nproc)
			{
				b=proclist[vwait(NULL)];
				fprintf (local_stderr, "\n\txxx Retrieved %s",fname[a]);
				if (tmp_list[b] && check_file_exists (tmp_list[b]))
				{
					CL=undump_constraint_list(CL,tmp_list[b]);
				}
				ns--;
			}
		}
	}

	while (ns)
	{
		int pid2;
		pid2=vwait(NULL);
		a=proclist[pid2];
		fprintf (local_stderr, "\n\txxx Retrieved %s",fname[a]);
		if (tmp_list[a] && check_file_exists (tmp_list[a]))
		{
			CL=undump_constraint_list (CL,tmp_list[a]);
		}
		ns--;
	}
	fprintf ( local_stderr, "\n\n\tAll Methods Retrieved\n");

	if (tmp_list[n_list] && check_file_exists (tmp_list[n_list]))
	{
		CL=undump_constraint_list(CL,tmp_list[n_list]);
	}

	CL->local_stderr=local_stderr;

	vfree (proclist);
	vfree (tmp_list);
	return CL;
}


/**
 * Central method for reading different types of files.
 *
 * \todo This MUST be documented! Explain what types of input are possible and so on.
 * \todo It should also contain odds that occur when reading several inputs, like keeping the maximum edge weight.
 *
 * \param[in,out] CL         points to the global ::Constrain_list object.
 * \param[in]     in_fname   Name of the input file to be read.
 */
Constraint_list* read_constraint_list(Constraint_list *CL,char *in_fname,char *in_mode, char *mem_mode,char *weight_mode)
{
	Sequence *SL=NULL, *TS=NULL;
	int a;
	Constraint_list *SUBCL=NULL;
	char *read_mode;
	char *fname;

	fname=in_fname;
	read_mode=(char*)vcalloc ( STRING, sizeof (char));

	if ( is_lib_list (in_fname))sprintf ( read_mode, "lib_list");
	else if ( in_mode)sprintf (read_mode, "%s", in_mode);
	else if ( fname[0]=='A'){sprintf ( read_mode, "aln");fname++;}
	else if ( fname[0]=='L'){sprintf ( read_mode, "lib");fname++;}
	else if ( fname[0]=='M'){sprintf ( read_mode, "method");fname++;}
	else if ( fname[0]=='S'){sprintf ( read_mode, "sequence");return CL;}
	else if ( fname[0]=='P'){sprintf ( read_mode, "pdb")     ;return CL;}
	else if ( fname[0]=='R'){sprintf ( read_mode, "profile") ;return CL;}
	else if ( fname[0]=='X'){sprintf ( read_mode, "matrix");++fname;}
	else if ( fname[0]=='W'){sprintf ( read_mode, "structure");fname++;}
	else
	{
		fprintf ( stderr, "\nERROR: The descriptor %s could not be identified as a file or a method.[FATAL]\nIf it is a method file please indicate it with M%s\n", fname, fname);
		myexit (EXIT_SUCCESS);
	}

	fprintf (CL->local_stderr, "\n\t%s [%s]\n", fname, read_mode);


	if ( strm (read_mode, "lib_list"))
	{
		int n, a;
		char **l;
		l=read_lib_list (fname, &n);
		for ( a=0; a<n; a++)
			CL=read_constraint_list (CL,l[a],in_mode, mem_mode,weight_mode);
		free_char (l,-1);
	}
	else if (strm(read_mode, "binary"))
	{
		fprintf ( stderr, "\nERROR: Library %s: binary mode is not any more supported [FATAL:%s]\n", fname,PROGRAM);
		myexit (EXIT_FAILURE);
	}
	else if ( strm (fname, "make_test_lib"))
	{
		CL=make_test_lib (CL);
	}
	else if ( strm2 (read_mode,"ascii","lib"))
	{

		SUBCL=read_constraint_list_file(NULL, fname);

	}
	else if (strm(read_mode, "method"))
	{
		CL=produce_list ( CL, CL->S, fname,weight_mode,mem_mode);
	}
	else if (strm(read_mode, "matrix"))
	{
	  CL->residue_index=NULL;
	  CL->extend_jit=0;
	  CL->M=read_matrice ( fname);
	}
	else if ( strm ( read_mode, "structure"))
	{
		if ( CL->ne>0)
		{
			fprintf ( stderr, "\nERROR: Wstructure must come before Mmethod or Aaln [FATAL:%s]",PROGRAM);
			myexit (EXIT_FAILURE);
		}

		if ( !(CL->STRUC_LIST))
		{
			CL->STRUC_LIST=declare_sequence (1,1,10000);
			(CL->STRUC_LIST)->nseq=0;
		}
		SL=CL->STRUC_LIST;

		if ( check_file_exists(fname))
		{
			TS=main_read_seq ( fname);
			for (a=0; a<TS->nseq; a++)sprintf (SL->name[SL->nseq++], "%s", TS->name[a]);
			free_sequence (TS, TS->nseq);
		}
		else
		{
			sprintf (SL->name[SL->nseq++], "%s", fname);
		}
	}
	else if (strm (read_mode, "aln"))
	{
		CL=aln_file2constraint_list ( fname,CL,weight_mode);
	}
	else
	{
		SUBCL=read_constraint_list_file(SUBCL, fname);
	}


	if (SUBCL)
	{
		CL=merge_constraint_list    (SUBCL, CL, "default");
		free_constraint_list_full (SUBCL);
	}

	return CL;
}



Sequence *precompute_blast_db (Sequence *S, char **ml, int n)
{
	int a;
	for (a=0; a<n; a++)
	{
		if ( strstr (ml[a], "blast"))
		{
			S->blastdb=vtmpnam(NULL);
			seq2blastdb(S->blastdb,S);
			return S;
		}
		a++;
	}
	S->blastdbS=S;
	return S;
}


#define is_seq_source(Symbol,Mode,SeqMode)            (Symbol==Mode && (SeqMode==NULL || strm (SeqMode, "ANY") || (SeqMode[0]!='_' && strchr (SeqMode,Symbol)) || (SeqMode[0]=='_' && !strchr (SeqMode,Symbol))))


/**
 * Reads sequence files and alignments or can extract the sequence from a library.
 *
 * This is basically a decision function that calls the appropriate read function depending on the type of the input file.
 * \sa ::main_read_aln to read Alignment files or sequence files
 * \sa ::read_seq_in_list to read from a library
 * \sa several auxiliary functions like ::aln2seq, ::seq2unique_name_seq and ::merge_seq
 */
Sequence * read_seq_in_n_list(char **fname, int n, char *type, char *SeqMode)
{
	int nseq=0;
	int a, b;
	Alignment *A;
	char **sequences=NULL;
	char **seq_name=NULL;
	Genomic_info *genome_co = NULL;
	Sequence *S=NULL;
	Sequence *S1;
	char mode;

	

	/*THE TYPE OF EACH FILE MUST BE INDICATED*/
	/*SeqMode indicates the type of file that can be used as sequence sources*/
	/*
	ANY: any mode
	SL: only sequences from Libraries and Sequences
	_A: anything BUT sequences from A(lignments)
	*/

	if ( n==0)
	  {
	    myexit (fprintf_error ( stderr, "NO in FILE"));
	  }
	else
	  {
	    for ( a=0; a< n ; a++)
	      {
		static char *buf;
		char *lname;
		if (buf)vfree (buf);
		
		//prevent silly bug when the file Sname also exists...		    
		if (fname[a][0]=='S' && file_exists(NULL, fname[a]+1))
		  {
		    buf=(char*)vcalloc (strlen (fname[a])+1, sizeof(char));
		    sprintf (buf, "%s", fname[a]);
		  }
		else
		  buf=name2type_name(fname[a]);
		
		mode=buf[0];lname=buf+1;
		
		if (is_seq_source ('A', mode, SeqMode))
		  {
		    
		    
		    A=main_read_aln (lname,NULL);
		    
		    S1=aln2seq(A);
		    
		    S1=seq2unique_name_seq (S1);
		    
		    if ((S=merge_seq ( S1, S))==NULL){fprintf ( stderr, "\nERROR: Sequence Error in %s [FATAL:%s]\n",lname, PROGRAM); myexit(EXIT_FAILURE);}
		    free_aln (A);
		    free_sequence (S1, S1->nseq);
		    
		  }
		else if ( is_seq_source ('R', mode, SeqMode))
		  {
		    
		    S=add_prf2seq (lname, S);
		    
		  }
		else if (is_seq_source ('P', mode, SeqMode))
		  {
		    int i;
		    
		    S1=get_pdb_sequence (lname);
		    if (S1==NULL)
		      {
			add_warning ( stderr, "Could not use PDB: %s", lname);
		      }
		    else
		      {
			if ((S=merge_seq ( S1, S))==NULL){fprintf ( stderr, "\nERROR: Sequence Error in %s [FATAL:%s]\n",lname, PROGRAM); myexit(EXIT_FAILURE);}
			i=name_is_in_list (S1->name[0], S->name, S->nseq, 100);
			(S->T[i])->P=fill_P_template (S->name[i], lname, S);
		      }
		    free_sequence (S1, S1->nseq);
		  }
		else if ( mode=='M');
		else if ( mode=='X');
		else if ( mode=='W');
		
		else if (is_seq_source ('S', mode, SeqMode))
		  {
		    /*1 Try with my routines (read t_coffee and MSF)*/
		    
		    if ( (A=main_read_aln ( lname, NULL))!=NULL)
		      {
			
			
			S1=aln2seq(A);
			free_aln(A);
		      }
		    else
		      {
			S1=main_read_seq (lname);
		      }
		    
		    for ( b=0; b< S1->nseq; b++)ungap(S1->seq[b]);
		    S1=seq2unique_name_seq (S1);
		    
		    
		    if ((S=merge_seq ( S1, S))==NULL){fprintf ( stderr, "\nSequence Error in %s [FATAL:%s]\n",lname,PROGRAM); myexit(EXIT_FAILURE);}
		    
		    free_sequence (S1,S1->nseq);
		  }
		else if (is_seq_source ('L', mode, SeqMode))
		  {
		    S1=read_seq_in_list(lname);
		    
		    if ((S=merge_seq( S1, S))==NULL)
		      {
			fprintf ( stderr, "\nSequence Error in %s [FATAL:%s]\n",lname,PROGRAM);
			myexit(EXIT_FAILURE);
		      }
		    free_sequence(S1, S1->nseq);
		  }
		
		else if ( !strchr ( "ALSMXPRWG", mode))
		  {
		    myexit (fprintf_error ( stderr, "%s is neither a file nor a method [FATAL:%s]\n", lname));
		  }
	      }
	   
	    S=remove_empty_sequence (S);
	    
	    
	    if ( type && type[0] )sprintf ( S->type, "%s", type);
	    else S=get_sequence_type (S);
	    
	    if ( strm (S->type, "PROTEIN_DNA"))
	      {
		for ( a=0; a< S->nseq; a++)
		  {
		    if (strm ( get_string_type ( S->seq[a]), "DNA") ||strm ( get_string_type ( S->seq[a]), "RNA")  );
		    else if ( strm ( get_string_type ( S->seq[a]), "PROTEIN"))
		      {
			S->seq[a]=thread_aa_seq_on_dna_seq (S->seq[a]);
			S->len[a]=strlen (S->seq[a]);
			S->max_len=MAX(S->max_len, S->len[a]);
		      }
		  }
	      }
	    
	    for (a=0; a<S->nseq; a++)
	      {
		int l=S->len[a]=strlen (S->seq[a]);
		S->max_len=(S->max_len>l)?S->max_len:l;
	      }
	    return S;
	  }
		
	return NULL;
}

int read_cpu_in_list ( char *fname)
{
	FILE *fp;
	int c;
	int cpu=0;

	fp=vfopen ( fname, "r");
	while ( (c=fgetc(fp))!='#');
	while ( (c=fgetc(fp))!='C' && c!=EOF);
	if ( c=='C')fscanf( fp, "PU %d\n", &cpu);
	vfclose ( fp);
	return cpu;
}



/*

char * expand_constraint_list_file ( char *file)
{
  char *new_file;
  FILE *IN, *OUT;
  int a, b, c, n;
  char **list;
  static char *buf;
  int line=0;
  if ((!token_is_in_file (file,"+BLOCK+")) && (!token_is_in_file (file,"++")))
	  return file;

  new_file=vtmpnam (NULL);
  IN=vfopen ( file,"r");
  OUT=vfopen (new_file, "w");

  while ( (c=fgetc (IN))!=EOF)
    {
      ungetc (c, IN);
      buf=vfgets (buf, IN);
      line++;
	  if (( !strstr (buf, "+BLOCK+")) && ( !strstr (buf, "++")))
			fprintf (OUT, "%s", buf);
      else
	{
	  list=string2list (buf);
	  n=atoi (list[2]);
	  for (a=0; a< n; a++)
	    {
	      if (atoi(list[3])<0 || atoi(list[4])<0)
		myexit (fprintf_error (stderr,"Illegal Coordinates: File %s, Line %d : %s", file,line, buf));
	      fprintf ( OUT, "%5d %5d ",atoi(list[3])+a, atoi(list[4])+a);
	      for (b=5; b<atoi(list[0]); b++)
		{
		  if (atoi(list[b])<0){HERE ("*** %s ****", buf); exit (0);}
		  fprintf ( OUT, "%s ", list[b]);
		}
	      fprintf (OUT, "\n");
	    }
	  free_char (list, -1);
	}
    }
  vfclose (IN); vfclose (OUT);
  return new_file;
}*/


Constraint_list * make_test_lib(Constraint_list *CL)
{
	int a, b, c, l1;
	Sequence *S;
	int *entry;

	entry=(int*)vcalloc (CL->entry_len+1, sizeof (int));
	HERE ("cncn: making artificial lib");//Keep this warning alive
	S=CL->S;

	fprintf ( stderr, "\n");
	for (a=0; a< S->nseq-1; a++)
	{
		fprintf ( stderr, "[%d-%d]", a, CL->ne);
		for ( b=a+1; b< S->nseq; b++)
		{

			l1=MIN(S->len[a], S->len[b]);
			l1=MIN(10, l1);
			for (c=0; c<l1; c++)
			{
				entry[SEQ1]=a;
				entry[SEQ1]=b;
				entry[R1]=c+1;
				entry[R2]=c+1;
				entry[WE]=100;
				add_entry2list (entry, CL);
			}
		}
	}
	return CL;
}



Constraint_list * old_read_constraint_list_file(Constraint_list *CL, char *fname)
{
	Sequence *NS;
	int **index, *entry;
	int c,line=0,s1, s2, r1, r2, misc, cons,x,we;
	FILE *fp;
	char *buf=NULL;
	int error=0;
	int keepNS=0;
	int start=0;
	
	NS=read_seq_in_n_list (&fname, 1,NULL, NULL);
	
	if (!CL)
	  {
	    CL=declare_constraint_list_simple(NS);
	    keepNS=1;
	  }
	

	index=fix_seq_seq (NS, CL->S);
	entry=(int*)vcalloc ( CL->entry_len+1, sizeof (int));
	
	

	fp=vfopen(fname,"r");
	while ((c=fgetc(fp))!='#' && c!=EOF){line+=(c=='\n')?1:0;}


	ungetc (c, fp);
	while ((buf=vfgets ( buf, fp))!=NULL)
	{
		line++;
		if (buf[0]=='!')
		  {
		    if (strstr (buf, "!CMT:"))
		      {
			if ((CL->comment && ! strstr (CL->comment, buf))|| !CL->comment)
			  CL->comment=vcat(CL->comment, buf);
		      }
		    else continue;
		  }
		else if (buf[0]!='#' && !start)continue;
		else if (buf[0]=='#')
		  {
		    start=1;
		    sscanf ( buf, "#%d %d", &s1, &s2);
		    s1--; s2--;
		    if (s1>NS->nseq || s2>NS->nseq)error=1;
		  }
		else
		{
			cons=misc=r1=r2=we=0;
			x=sscanf (buf, "%d %d %d %d %d",&r1,&r2,&we,&cons,&misc);

			if (r1>NS->len[s1] || r2>NS->len[s2] || x<3)error=1;
			else
			{
				int is1, is2, ir1, ir2;

				is1=entry[SEQ1]=index[s1][0];//0-N-1
				is2=entry[SEQ2]=index[s2][0];//0-N-1
				ir1=entry[R1]=index[s1][r1];//1-N
				ir2=entry[R2]=index[s2][r2];//1-N
				entry[WE]=we;//0-1000
				entry[CONS]=cons;
				entry[MISC]=misc;
				if (ir1>(CL->S)->len[is1] || ir2>(CL->S)->len[is2])
				{

					myexit(fprintf_error (stderr, "%s::%d -- %s::%d ---- %s::%d %s::%d", NS->name[s1],r1, NS->name[s2],r2, (CL->S)->name[is1],ir1, (CL->S)->name[is2],ir2));
				}
				if (entry[SEQ1]>-1 && entry[SEQ2]>-1 && entry[R1]>0 && entry[R2]>0 && entry[WE]>0)
					add_entry2list (entry, CL);
			}
		}
		if (error)
			printf_exit (EXIT_FAILURE,stderr,"Parsing Error [L:%d F:%s S:%s]",line,fname,buf);
	}

	vfclose (fp);

	if (!keepNS)free_sequence (NS,-1);
	vfree (entry);
	vfree (buf);
	free_int (index, -1);
	return CL;
}

Constraint_list *read_constraint_list_file(Constraint_list *CL, char *fname)
{
  Sequence *NS;
  int **index, *entry;
  int c,line=0,s1, s2, r1, r2, misc, cons,x,we;
  FILE *fp;
  char *buf=NULL;
  int error=0;
  int keepNS=0;
  unsigned int i;
  char arg[10];
  unsigned int length;
  int start=0;
  int a;
  int d=0;
  NS=read_seq_in_n_list (&fname, 1,NULL, NULL);
  
  if (!CL)
    {
      CL=declare_constraint_list_simple(NS);
      keepNS=1;
    }
  
  index=fix_seq_seq (NS, CL->S);
  
  

  entry=(int*)vcalloc ( CL->entry_len+1, sizeof (int));
  
  
  fp=vfopen(fname,"r");
  //while ((c=fgetc(fp))!='#' && c!=EOF){line+=(c=='\n')?1:0;}
  //ungetc (c, fp);
  
  while ((buf=vfgets ( buf, fp))!=NULL)
    {
      line++;
      if (buf[0]=='!')
	{
	  if (strstr (buf, "!CMT:"))
	    {
	      if ((CL->comment && ! strstr (CL->comment, buf))|| !CL->comment)
		CL->comment=vcat(CL->comment, buf);
	    }
	  else continue;
	}	
      else if (buf[0]!='#' && !start)continue;
      else if (buf[0]=='#')
	{
	  start=1;
	  sscanf ( buf, "#%d %d", &s1, &s2);
	  s1--; s2--;
	  if (s1>NS->nseq || s2>NS->nseq)error=1;
	}
      else
	{
	  cons=misc=r1=r2=we=0;
	  if (buf[0] != '+')
	    {
	      x=sscanf (buf, "%d %d %d %d %d",&r1,&r2,&we,&cons,&misc);
	      length = 1;
	    }
	  else
	    {
	      x=sscanf (buf, "%s %d %d %d %d %d",arg, &length, &r1,&r2,&we,&cons,&misc);
	    }
	  
	  for (i = 0; i < length; ++i)
	    {
	      if (r1>NS->len[s1] || r2>NS->len[s2] || x<3)error=1;
	      else
		{
		  int is1, is2, ir1, ir2;
		  int ax, bx, cx, dx, ex, fx;
		  
		  is1=entry[SEQ1]=index[s1][0];//0-N-1
		  is2=entry[SEQ2]=index[s2][0];//0-N-1
		  ir1=entry[R1]=index[s1][r1];//1-N
		  ir2=entry[R2]=index[s2][r2];//1-N
		  entry[WE]=we;//0-1000
		  entry[CONS]=cons;
		  entry[MISC]=misc;
		  if(ir1>0){;}
		  if (ir2>0){;}
		  if (is1>0){;}
		  if (is2>0){;}
		  if((CL->S)->len[is1]){;}
		  if((CL->S)->len[is2]){;}
		  
		  if (ir1>(CL->S)->len[is1] || ir2>(CL->S)->len[is2])
		    {
		      
		      myexit(fprintf_error (stderr, "%s::%d -- %s::%d ---- %s::%d %s::%d", NS->name[s1],r1, NS->name[s2],r2, (CL->S)->name[is1],ir1, (CL->S)->name[is2],ir2));
		    }
		  if (entry[SEQ1]>-1 && entry[SEQ2]>-1 && entry[R1]>0 && entry[R2]>0 && entry[WE]>0)
		    {
		      int ine=CL->ne;
		      int id=d;
		      add_entry2list (entry, CL);
		      d++;
		    }
		}
	      
	      ++r1;
	      ++r2;
	    }
	  
	}
      if (error)
	{
	  printf_exit (EXIT_FAILURE,stderr,"Parsing Error [Line:%d F:%s S:%s][r1=%d r2=%d L1=%d L2=%d NP:%d]\[%s]\n[%s]",line,fname,buf, r1, r2,NS->len[s1],r2>NS->len[s2], x, NS->seq[s1], NS->seq[s2] );
	}
    }
  vfclose (fp);
  
  if (!keepNS)free_sequence (NS,-1);
  vfree (entry);
  vfree (buf);
  free_int (index, -1);
  return CL;
}
Sequence *constraint_list2seq (char *lname)
{
  return read_seq_in_list(lname);
}

Sequence * read_seq_in_list (char *fname)
{
  char *name=NULL;
  char *seq=NULL;
  char *buf=NULL;
  FILE *fp1,*fp2;
  char *tmp=vtmpnam(NULL);
  int a,l;
  int nseq;

  fp1=vfopen (fname, "r");
  fp2=vfopen (tmp, "w");
  
  fp1=skip_commentary_line_in_file ('!', fp1);
  fscanf (fp1, "%d\n", &nseq);
  for (a=0; a<nseq; a++)
    {
      buf=vfgets(buf, fp1);
      name=csprintf (name, "%s", buf);
      seq =csprintf (seq , "%s", buf);
      sscanf (buf, "%s %d %s\n", name, &l,seq);
      fprintf ( fp2, ">%s\n%s\n",name, seq);
    }
  vfclose (fp1); vfclose (fp2);
  vfree (name); vfree(seq); vfree(buf);
  return get_fasta_sequence (tmp, NULL);
}

/*********************************************************************/
/*                                                                   */
/*                        EXTENDED LIST OUTPUT                       */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

FILE * save_extended_constraint_list      (  Constraint_list *CL, char *mode, FILE *fp)
{
	int a, b;

	if (!CL || !CL->S || CL->residue_index)return fp;
	if ( !fp)fp=stdout;

	fp=save_list_header (fp,CL);

	for ( a=0; a< (CL->S)->nseq; a++)
	{
		for ( b=a; b<(CL->S)->nseq; b++)
		{

			if ( a==b && !CL->do_self)continue;
			fp=save_extended_constraint_list_pair(CL, mode, (CL->S)->name[a], (CL->S)->name[b], fp);
		}
	}
	fprintf (fp, "! SEQ_1_TO_N\n");
	return fp;
}


FILE * save_extended_constraint_list_pair (  Constraint_list *CL, char *mode, char* seq1, char * seq2,FILE *fp)
{
	int a, b, t;
	int s1, s2, score;
	char *p;


	if ((p=strstr (mode, "THR"))!=NULL)t=atoi(p+3);
	else t=0;

	s1=name_is_in_list (seq1,(CL->S)->name, (CL->S)->nseq, 100);
	s2=name_is_in_list (seq2,(CL->S)->name, (CL->S)->nseq, 100);

	if ( s1==-1)
	{
		fprintf ( stderr, "Output Error: %s is not a sequence [FATAL:%s]\n", seq1, PROGRAM);
		crash ("");
	}
	if ( s2==-1)
	{
		fprintf ( stderr, "Output Error: %s is not a sequence [FATAL:%s]\n", seq2, PROGRAM);
		crash ("");
	}

	if ( strstr (mode, "pair"))fprintf (fp, "# 1 2\n");
	else if ( strstr (mode, "lib"))fprintf (fp, "# %d %d\n", s1+1, s2+1);

	for ( a=0; a<(CL->S)->len[s1]; a++)
	{

		for ( b=0; b<(CL->S)->len[s2]; b++)
		{
			if ( a>=b && s1==s2)continue;
			if ( strstr (mode, "pc"))score=residue_pair_extended_list_pc (CL, s1,a+1, s2, b+1);
			else if ( strstr (mode, "raw"))score=residue_pair_extended_list_raw (CL, s1,a+1, s2, b+1);
			else
				score=CL->evaluate_residue_pair (CL, s1,a+1, s2, b+1);

			if (score<=t) continue;
			fprintf (fp, "%5d %5d %5d \n", a+1, b+1, score);

		}
	}
	return fp;
}


/*********************************************************************/
/*                                                                   */
/*                         LIST OUTPUT                               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
#ifdef MMMMMMMM
FILE *save_extended_constraint_list ( Constraint_list *CL,Sequence *S, char *fname)
{
	int a, b, c, d;
	int *tr, *ns;
	int **pos0, **l_s;
	int epsilon=0;
	Alignment *A;
	FILE *fp;

	fp=vfopen (fname, "w");
	fp=save_sub_list_header(fp, S->nseq, S->name, CL);

	tr=vcalloc (S->nseq+1, sizeof (int));
	for ( b=0,a=0; a< S->nseq; a++)
	{
		int i;
		if ( (i=name_is_in_list(S->name[a],(CL->S)->name,(CL->S)->nseq, 100))==-1)
		{
			printf_exit (EXIT_FAILURE, stderr, "\nERROR: Sequence %s is not part of the sequence dataset [FATAL:%s]", S->name[a], PROGRAM);

		}
		else
		{
			tr[a]=i;
		}
	}

	A=declare_aln (S);
	pos0=vcalloc ( S->nseq, sizeof (int*));
	for (a=0; a<S->nseq; a++)
	{
		int l;
		l=strlen (S->seq[a]);
		A->seq_al[a]=S->seq[a];
		pos0[a]=vcalloc (l+1, sizeof (int));
		for (b=0; b<l; b++)pos[a][b]=b+1;
	}
	l_s=declare_int (2,2);
	ns=vcalloc ( 2, sizeof (int));


	for ( a=0; a< S->nseq-1; a++)
		for ( b=a+1; b<S->nseq; b++)
		{
			int pos_i, pos_j, s;
			l_s[0]=tr[a];l_s[1]=tr[b];
			for ( pos_i=0; pos_i< S->len[a]; pos_i++)
				for (pos_j=0; pos_j<S->len[b]; pos_j++)
				{
					s=(CL->get_dp_cost) ( A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],pos_j-1, CL);
					if (s>epsilon)fprintf (fp, "%d %d %d", i, j, s);
				}
		}
		return fp;
}
#endif
int save_contact_constraint_list (Constraint_list *CL, char *name)
{
  FILE *fp;
  Sequence *S;
  int s1, s2;
  

  S=CL->S;
  fp=vfopen (name, "w");
  fprintf ( fp, "! TC_LIB_FORMAT_01\n");
  if (!CL->comment)fprintf ( fp, "!CMT: Intra sequence contacts\n");
  else fprintf ( fp, "%s", CL->comment);
  
  fprintf (fp, "%d\n", S->nseq);
  for (s1=0; s1<S->nseq; s1++)fprintf ( fp, "%s %d %s\n",S->name[s1], S->len[s1], S->seq[s1]);
  for (s1=0; s1<S->nseq; s1++)
    {
      for (s2=s1; s2<S->nseq; s2++)
	{
	  
	  int r1, r2, we;
	  int **contact;
	  int head=0;
	  
	  
	  if (s1==s2)contact=declare_int (S->len[s1]+1,S->len[s1]+1);
	  
	  for (r1=1; r1<=S->len[s1]; r1++)
	    {
	      int b;
	      for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
		{
		  int rs2;
		  rs2=CL->residue_index[s1][r1][b+SEQ2];
		  r2 =CL->residue_index[s1][r1][b+R2];
		  we =CL->residue_index[s1][r1][b+WE];
		  if (rs2==s2)
		    {
		      if (!head)
			{
			  fprintf (fp,"#%d %d\n",s1+1,s2+1);
			  head=1;
			}
		      if (s1==s2 && !contact[r1][r2])
			{
			  if (r1<r2)fprintf (fp, "%d %d %d\n", r1, r2, we);
			  else fprintf (fp, "%d %d %d\n", r2, r1, we);
			  contact[r1][r2]=contact[r2][r1]=1;
			}
		      else fprintf (fp, "%d %d %d\n", r1, r2, we);
		    }
		}
	    }
	  if (s1==s2)free_int (contact,-1);
	}
    }
  fprintf ( fp, "! SEQ_1_TO_N\n");
  vfclose (fp);
  return 1;
}
int save_contact_constraint_list_old (Constraint_list *CL, char *name)
{
  FILE *fp;
  Sequence *S1;
  int a,b;
  

  S1=CL->S;
  fp=vfopen (name, "w");
  fprintf ( fp, "! TC_LIB_FORMAT_01\n");
  if (!CL->comment)fprintf ( fp, "!CMT: Intra sequence contacts\n");
  else fprintf ( fp, "%s", CL->comment);
  
  fprintf (fp, "%d\n", S1->nseq);
  for (a=0; a<S1->nseq; a++)fprintf ( fp, "%s %d %s\n",S1->name[a], S1->len[a], S1->seq[a]);
  for (a=0; a<S1->nseq; a++)
    {
      int r1, r2, we;
      int **contact;
      fprintf (fp,"#%d %d\n",a+1,a+1); 
      contact=declare_int (S1->len[a]+1,S1->len[a]+1);
      
      for (r1=1; r1<=S1->len[a]; r1++)
	{
	  for (b=1; b<CL->residue_index[a][r1][0]; b+=ICHUNK)
	    {
	      r2=CL->residue_index[a][r1][b+R2];
	      we=CL->residue_index[a][r1][b+WE];
	      if (!contact[r1][r2])
		{
		  if (r1<r2)fprintf (fp, "%d %d %d\n", r1, r2, we);
		  else fprintf (fp, "%d %d %d\n", r2, r1, we);
		}
	      contact[r1][r2]=contact[r2][r1]=1;
	    }
	}
      free_int (contact,-1);
    }

  fprintf ( fp, "! SEQ_1_TO_N\n");
  vfclose (fp);
  return 1;
}
FILE * save_constraint_list ( Constraint_list *CL,int start, int len, char *fname, FILE *fp,char *mode, Sequence *S)
{
	int a, b;
	static int* translation;


	if ( fp==NULL)
	{
		if ( translation!=NULL)vfree(translation);
		translation=(int*)vcalloc ( (CL->S)->nseq+1, sizeof (int));
		for ( b=0,a=0; a< (CL->S)->nseq; a++)
		{
			if ( name_is_in_list((CL->S)->name[a],S->name,S->nseq, 100)==-1)
			{
				(CL->S)->len[a]=-1;
				translation [a]=-1;
			}
			else
			{
				translation[a]=b++;
			}
		}

	}
	if (strm2(mode, "lib","ascii"))
	{
		if ( fp==NULL)fp=vfopen ( fname, "w");
		fp=save_list_header (fp,CL);
		fp=save_constraint_list_ascii(fp, CL, 0, CL->ne, translation);
	}

	else
	{

		myexit(fprintf (stderr,"\nUNKOWN MODE FOR OUTPUT: %s", mode));
	}
	return fp;
}

FILE * save_sub_list_header ( FILE *OUT, int n, char **name, Constraint_list *CL)
{
	int a,b;
	int nseq=0;



	for ( a=0; a<(CL->S)->nseq; a++)
		for ( b=0; b<n; b++)
			if (strm (name[b] ,  (CL->S)->name[a]))
				nseq+=((CL->S)->len[a]!=-1);

			fprintf ( OUT, "! TC_LIB_FORMAT_01\n%d\n",nseq);
		for ( a=0; a<n; a++)
			for ( b=0; b<(CL->S)->nseq; b++)
				if (strm (name[a] ,  (CL->S)->name[b]))
					if ((CL->S)->len[b]!=-1) fprintf ( OUT, "%s %d %s\n", (CL->S)->name[b], (CL->S)->len[b],(CL->S)->seq[b]);

					return OUT;
}

FILE * save_list_header ( FILE *OUT,Constraint_list *CL)
{
	int a;
	int nseq=0;

	for ( a=0; a<(CL->S)->nseq; a++)nseq+=((CL->S)->len[a]!=-1);


	fprintf ( OUT, "! TC_LIB_FORMAT_01\n%d\n",nseq);
	for ( a=0; a<(CL->S)->nseq; a++)
		if ((CL->S)->len[a]!=-1)
		{
			fprintf ( OUT, "%s %d %s\n", (CL->S)->name[a], (CL->S)->len[a],(CL->S)->seq[a]);

		}
		return OUT;
}

FILE *save_list_footer (FILE *OUT,Constraint_list *CL)
{
	if ( CL->cpu)fprintf (OUT, "! CPU %d\n",get_time());
	fprintf (OUT, "! SEQ_1_TO_N\n");
	return OUT;
}
FILE * save_constraint_list_ascii ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation)
{
	int a, b, s1, s2, r1, r2;
	int d1=0,d2=0;
	if (len==start && CL->cpu!=-1)
	{
		fprintf (OUT, "! CPU %d\n",get_time());
		return OUT;
	}
	else
	{
		Sequence *S=CL->S;
		int ***cacheI=(int***)vcalloc (S->nseq, sizeof (int**));
		int **cacheR=(int**)vcalloc (S->nseq, sizeof (int*));
		int *tot=(int*)vcalloc ( S->nseq, sizeof (int));
		int *max=(int*)vcalloc ( S->nseq, sizeof (int));
		for (a=0;a<S->nseq; a++)
		  {
		    cacheI[a]=(int**)vcalloc (S->len[a], sizeof (int*));
		    cacheR[a]=(int*)vcalloc (S->len[a], sizeof (int));
		  }
		for (a=0; a<S->nseq; a++)max[a]=S->len[a];

		for (s1=0; s1<S->nseq; s1++)
		  {
		  for (r1=1; r1<=S->len[s1]; r1++)
		    {
		      for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{

			  s2=CL->residue_index[s1][r1][b+SEQ2];
			  if (tot[s2]>=max[s2])
			    {
			      max[s2]+=100;
			      cacheI[s2]=(int**)vrealloc (cacheI[s2], max[s2]*sizeof (int*));
			      cacheR[s2]=(int*)vrealloc (cacheR[s2], max[s2]*sizeof (int*));
			    }

			  cacheI[s2][tot[s2]]=CL->residue_index[s1][r1]+b;
			  cacheR[s2][tot[s2]]=r1;
			  tot[s2]++;
			  d1++;
			}
		    }

		  for (s2=0;s2<S->nseq; s2++)
		    {
		      int x1=translation[s1];
		      int x2=translation[s2];
		      if (tot[s2] && x1!=-1 && x2!=-1 && x1<x2)
			{
			  fprintf ( OUT, "#%d %d\n", x1+1, x2+1);
			  for ( r2=0; r2<tot[s2];r2++)
			    {
			      int *v=cacheI[s2][r2];
			      fprintf (OUT, "%5d %5d %5d %5d %5d\n",cacheR[s2][r2], v[R2], v[WE], v[CONS], v[MISC]);
			      d2++;
			    }
			}
		      tot[s2]=0;
		    }
		  }
		for (a=0; a<S->nseq; a++)
		  {
		    vfree (cacheI[a]);
		    vfree (cacheR[a]);
		  }
		vfree (cacheI); vfree (cacheR);
		vfree (tot);
		vfree (max);
	}

	return save_list_footer (OUT, CL);
	
}


/*********************************************************************/
/*                                                                   */
/*                         LIST CONVERTION                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

/**
 * Apply Consistency.
 *
 * Similar to ::relax_constraint_list, but here even edges that had zero weight before (i.e. that did not exist)
 * will be reweighted.
 * \param[in,out] CL The global Constraint_list object
 */

Constraint_list * extend_constraint_list ( Constraint_list *CL)
{
	Sequence *S;
	int **cache;
	int tot=0;
	char *tmp;
	FILE *fp;
	int a,e,b,c;
	int s1, r1, s2, r2,w2, s3, r3, w3;

	static int *entry;
	if (!CL || !CL->residue_index || !CL->S)return CL;
	S=CL->S;

	tmp=vtmpnam (NULL);
	fp=vfopen (tmp, "w");
	cache=declare_int (S->nseq, S->max_len+1);
	if (!entry)entry=(int*)vcalloc ( CL->entry_len+1, sizeof (int));
	for (s1=0; s1< S->nseq; s1++)
	{
		for ( r1=1; r1<=S->len[s1]; r1++)
		{
			cache[s1][r1]=1;

			for ( a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][a+SEQ2];
				r2=CL->residue_index[s1][r1][a+R2];
				w2=CL->residue_index[s1][r1][a+WE];
				cache[s2][r2]=w2;
			}

			for ( a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][a+SEQ2];
				r2=CL->residue_index[s1][r1][a+R2];
				w2=CL->residue_index[s1][r1][a+WE];

				for (b=1; b<CL->residue_index[s2][r2][0]; b+=ICHUNK)
				{
					s3=CL->residue_index[s2][r2][b+SEQ2];
					r3=CL->residue_index[s2][r2][b+R2];
					w3=CL->residue_index[s2][r2][b+WE];

					if (!cache[s3][r3] && s3>s1)
					{
						entry[SEQ1]=s1;
						entry[R1]=r1;
						entry[SEQ2]=s3;
						entry[R2]=r3;
						entry[WE]=MIN(w2, w3);
						entry [CONS]=1;
						tot++;c++;
						for (e=0; e<CL->entry_len; e++)fprintf ( fp, "%d ", entry[e]);

						entry[SEQ1]=s3;
						entry[R1]=r3;
						entry[SEQ2]=s1;
						entry[R2]=r1;
						for (e=0; e<CL->entry_len; e++)fprintf ( fp, "%d ", entry[e]);
					}
				}
			}
			cache[s1][r1]=0;
			for ( a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][a+SEQ2];
				r2=CL->residue_index[s1][r1][a+R2];
				cache[s2][r2]=0;
			}
		}
	}
	vfclose (fp);
	CL=undump_constraint_list(CL, tmp);
	free_int (cache, -1);
	return CL;
}

Constraint_list * fork_relax_constraint_list (Constraint_list *CL, int njobs,int nproc);


/**
 * Distributes the relaxation, see ::fork_relax_constraint_list.
 */
Constraint_list * relax_constraint_list (Constraint_list *CL)
{
  int njobs,nproc;


  if (!CL || !CL->S || !CL->residue_index) return CL;
  if (strstr ( CL->multi_thread, "relax")){njobs=nproc=get_nproc();}
  else {njobs=nproc=1;}

  return fork_relax_constraint_list (CL, njobs,nproc);
}

/**
 * First central step of Consistency evaluation.
 *
 * Here we evaluate the consistency of edges in the Constraint_list::residue_index.
 * This function distributes this process on several cores using ::vvfork.
 * On each core, the alogrithm looks like this:
 * \code
 * for (A,resA) in [set of Sequences] { // Iterate over Sequences A and all residues in this sequence
 *
 *     n_resA = len(Constraint_list(A,resA)) // number edges/constraint going out from this residue in sequence A
 *
 *     for (B, resB, w_resA_resB) in Constraint_list(A,resA) { // Iterate over all matches that can be seen from this residue (usually between 0 and 4 per external sequence ehen using proba_pair)
 *
 *         edges_from_A [B, resB] = w_resA_resB // keep in mind outgoing edges from sequence A.
 *         }
 *
 *     for (B, resB, w_resA_resB) in Constraint_list(A,resA) {
 *
 *         n_resB = len(Constraint_list(B,resB)) // number edges/constraint going out from resB in sequence B
 *
 *         for (C, resC, w_resB_resC) in Constraint_list(B,resB) {
 *
 *             if( C == A && resC == resA ) { // looking back at the same edge
 *
 *                 score += w_resA_resB
 *                 }
 *             else if ( isset edges_from_A[] ) { // found a 3-circle!
 *
 *                 score += MIN (w_resA_resB, w_resB_resC) // How to weight this consistency? Here we take the minimum.
 *                 }
 *             }
 *
 *         save (A,resA)--(B,resB) = score / MIN(n_resA, n_resB)      // Now we have a new score for edge resA -- resB. Save it
 *
 *         }
 *     unset edges_from_A
 *     }
 *
 * // when all jobs are finished
 * Update CL using the new scores.
 * \endcode
 *
 * During the relaxation step, the new scores are written to a file and updated not before all loops have finished.
 * It is important not to update edges while the computation is still going on - this would lead to unexpected results.
 * \note \b Normalisation: The new score, which is a sum of existing scores, is divided by the minimum of
 * 		\c n_resA nad n_resB, which are the numbers of outgoing edges (to all sequences).
 *
 * Finally, the Constrint_list is filtered with a hardcoded threshold of 10, see ::filter_constraint_list.
 * \attention The relaxation differs from the extension in the fact that only previously existing edges are reweighted.
 *            Non-existing edges (i.e. edges with zero weight) cannot be created, although they might be highly consistent.
 * \sa ::extend_constraint_list
 */

Constraint_list * fork_relax_constraint_list (Constraint_list *CL, int njobs,int nproc)
    {
    
  int a, s1, s2, r1, r2,j;

  int thr=10;
  FILE *fp;
  char **pid_tmpfile;
  int sjobs;
  int **sl;
  Sequence *S;
  int in;

  static int **hasch;
  static int max_len;

  int t_s1, t_s2, t_r1, t_r2,x;
  double score;
  int normalisation_mode=0;

  
  if (!CL || !CL->residue_index)return CL;
  fprintf ( CL->local_stderr, "\nLibrary Relaxation: Multi_proc [%d]\n ", nproc);

  
  
  if ( !hasch || max_len!=(CL->S)->max_len)
    {
      max_len=(CL->S)->max_len;
      if ( hasch) free_int ( hasch, -1);
      hasch=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
    }

  S=CL->S;
  in=CL->ne;

  sl=n2splits (njobs, (CL->S)->nseq);
  pid_tmpfile=(char**)vcalloc (njobs, sizeof (char*));

  for (sjobs=0,j=0;j<njobs; j++)
    {
      pid_tmpfile[j]=vtmpnam(NULL);

      // Child process
      if (vvfork (NULL)==0)
	{
	  int norm,norm1,norm2;

	  initiate_vtmpnam(NULL);
	  fp=vfopen (pid_tmpfile[j], "w");

	  for (s1=sl[j][0]; s1<sl[j][1]; s1++)
	    {
	      if (j==0)output_completion (CL->local_stderr,s1,sl[0][1],1, "Relax Library");
	      for (r1=1; r1<S->len[s1]; r1++)
		{
		  norm1=0;
		  for (x=1; x< CL->residue_index[s1][r1][0]; x+=ICHUNK)
		    {
		      t_s1=CL->residue_index[s1][r1][x+SEQ2];
		      t_r1=CL->residue_index[s1][r1][x+R2];
		      hasch[t_s1][t_r1]=CL->residue_index[s1][r1][x+WE];
		      norm1++;

		    }
		  for ( a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
		    {
		      score=0;
		      norm2=0;
		      s2=CL->residue_index[s1][r1][a+SEQ2];
		      r2=CL->residue_index[s1][r1][a+R2];

		      for (x=1; x< CL->residue_index[s2][r2][0]; x+=ICHUNK)
			{
			  t_s2=CL->residue_index[s2][r2][x+SEQ2];
			  t_r2=CL->residue_index[s2][r2][x+R2];
			  if (t_s2==s1 && t_r2==r1)score+=(float)CL->residue_index[s2][r2][x+WE];
			  else if (hasch[t_s2][t_r2])
			    {
			      score+=MIN((((float)hasch[t_s2][t_r2])),(((float)CL->residue_index[s2][r2][x+WE])));
			    }
			  norm2++;
			}
		      if (normalisation_mode==0)
			norm=MIN(norm1,norm2);
		      else if (normalisation_mode==1)
			norm=MAX(norm1,norm2);
		      else if (normalisation_mode==2)
			norm=(CL->S)->nseq;
			      
		      score=((norm)?score/norm:0);
		      fprintf (fp, "%d ",(int)(score));
		    }
		  for (x=1; x< CL->residue_index[s1][r1][0]; x+=ICHUNK)
		    {
		      t_s1=CL->residue_index[s1][r1][x+SEQ2];
		      t_r1=CL->residue_index[s1][r1][x+R2];
		      hasch[t_s1][t_r1]=0;
		    }
		}
	    }
	  vfclose (fp);
	  myexit (EXIT_SUCCESS);
	}
      // Parent process
      else
	{
	  sjobs++;
	}
    }

//Constraint_list * fork_relax_constraint_list (Constraint_list *CL, int njobs,int nproc)
//    {
//  int a, s1, s2, r1, r2,j;
//
//  /** Sascha: No filtering please */
//  int thr=0;
//  FILE *fp;
//  char **pid_tmpfile;
//  int sjobs;
//  int **sl;
//  Sequence *S;
//  int in;
//
//  static int **hasch;
//  static int max_len;
//
//  int t_s1, t_s2, t_r1, t_r2,x;
//  double score;
//  int np;
//
////   HERE ("%d", nproc);
//  if (!CL || !CL->residue_index)return CL;
//  fprintf ( CL->local_stderr, "\nLibrary Relaxation: Multi_proc [%d]\n ", nproc);
//
//  if ( !hasch || max_len!=(CL->S)->max_len)
//    {
//      max_len=(CL->S)->max_len;
//      if ( hasch) free_int ( hasch, -1);
//      hasch=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
//    }
//
//  S=CL->S;
//  in=CL->ne;
//
//  sl=n2splits (njobs, (CL->S)->nseq);
//  pid_tmpfile=vcalloc (njobs, sizeof (char*));
//
//  for (sjobs=0,j=0;j<njobs; j++)
//    {
//      pid_tmpfile[j]=vtmpnam(NULL);
//      if (vvfork (NULL)==0)
//	{
//	  int norm,norm1,norm2;
//
//	  initiate_vtmpnam(NULL);
//	  fp=vfopen (pid_tmpfile[j], "w");
//
//	  for (s1=sl[j][0]; s1<sl[j][1]; s1++)
//	    {
//	      if (j==0)output_completion (CL->local_stderr,s1,sl[0][1],1, "Relax Library");
//	      for (r1=1; r1<S->len[s1]; r1++)
//		{
//		  norm1=0;
//		  for (x=1; x< CL->residue_index[s1][r1][0]; x+=ICHUNK)
//		    {
//		      t_s1=CL->residue_index[s1][r1][x+SEQ2];
//		      t_r1=CL->residue_index[s1][r1][x+R2];
//		      hasch[t_s1][t_r1]=CL->residue_index[s1][r1][x+WE];
//		      /** Sascha: Sum up all weights / probabilities going out from here. */
//		      norm1 += hasch[t_s1][t_r1];
//
//		    }
//		  for ( a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
//		    {
//		      score=0;
//		      norm2=0;
//		      s2=CL->residue_index[s1][r1][a+SEQ2];
//		      r2=CL->residue_index[s1][r1][a+R2];
//
//		      for (x=1; x< CL->residue_index[s2][r2][0]; x+=ICHUNK)
//			{
//			  t_s2=CL->residue_index[s2][r2][x+SEQ2];
//			  t_r2=CL->residue_index[s2][r2][x+R2];
//			  if (t_s2==s1 && t_r2==r1) score+= (float) CL->residue_index[s2][r2][x+WE];
//			  else if (hasch[t_s2][t_r2])
//			    {
//			      /** Sascha: Product of probabilities! See ProbCons paper */
//			      score+= ((float)hasch[t_s2][t_r2]) * ((float)CL->residue_index[s2][r2][x+WE]);
//			    }
//			  /** Sascha: Sum up all weights / probabilities going out from here. */
//			  norm2 += CL->residue_index[s2][r2][x+WE];
//			}
//		      /** Sascha: Normalized by the product of posterior probabilities */
//		      //norm=(float)norm1*norm2;
//		      /** Sascha: Take only nominator, the marginal probability of Aligning A with B through C */
//		      norm = 1000;
//		      score=((score)?score/norm:0);
//		      fprintf (fp, "%d ",(int)(score));
//		    }
//		  for (x=1; x< CL->residue_index[s1][r1][0]; x+=ICHUNK)
//		    {
//		      t_s1=CL->residue_index[s1][r1][x+SEQ2];
//		      t_r1=CL->residue_index[s1][r1][x+R2];
//		      hasch[t_s1][t_r1]=0;
//		    }
//		}
//	    }
//	  vfclose (fp);
//	  myexit (EXIT_SUCCESS);
//	}
//      else
//	{
//	  sjobs++;
//	}
//    }



  while (sjobs>=0){vwait(NULL); sjobs--;}//wait for all jobs to complete
  for (j=0; j<njobs; j++)
    {
      fp=vfopen (pid_tmpfile[j], "r");
      for (s1=sl[j][0]; s1<sl[j][1]; s1++)
	for (r1=1; r1<S->len[s1]; r1++)
	  for ( a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
	    {
	      if (!(fscanf ( fp, "%d ", &CL->residue_index[s1][r1][a+WE])))
		{
		  printf_exit (EXIT_FAILURE,stderr, "Could not complete relaxation cycle");
		}
	    }
      vfclose (fp);
      remove (pid_tmpfile[j]);
    }

  CL=filter_constraint_list (CL,WE,thr);
  fprintf ( CL->local_stderr, "\nRelaxation Summary: [%d]--->[%d]\n", in,CL->ne);
  vfree (pid_tmpfile);
  free_int (sl, -1);
  return CL;
}



// relax constraint list for gene prediction
Constraint_list * expand_constraint_list_4gp (Constraint_list *CL, int T)
{
	int s1, s2, r1, r2, a, b, c, d, w;
	int *entry;
	Sequence *S;
	

	entry=(int*)vcalloc (CL->entry_len+1, sizeof (int));

	S=CL->S;
	for (a=0; a<S->nseq; a++)//loop sequences
	{
		for (b=1; b<=S->len[a];b++)
		{
			for (c=1; c<CL->residue_index[a][b][0];c+=ICHUNK)
			{


				s2=a;
				r2=b;
				s1=CL->residue_index[a][b][c+SEQ2];;
				r1=CL->residue_index[a][b][c+R2];
				w=residue_pair_extended_list_4gp (CL,s1, r1,s2, r2);
				if (w>T)
				{
					entry[SEQ2]=s2;
					entry[SEQ1]=s1;
					entry[R1]=r2;
					entry[R2]=r1;
					entry[WE]=w;
					add_entry2list (entry, CL);
				}
				for (d=c+3; d<CL->residue_index[a][b][0]; d+=ICHUNK)
				{
					s2=CL->residue_index[a][b][d];
					r2=CL->residue_index[a][b][d+1];
					w=residue_pair_extended_list_4gp (CL,s1, r1,s2, r2);
					if (w>T)
					{
						entry[SEQ2]=s2;
						entry[SEQ1]=s1;
						entry[R1]=r2;
						entry[R2]=r1;
						entry[WE]=w;
					}
				}
			}
		}
	}
	vfree (entry);
	return CL;
}


Constraint_list * nfork_relax_constraint_list_4gp (Constraint_list *CL);
Constraint_list * fork_relax_constraint_list_4gp (Constraint_list *CL);
Constraint_list * relax_constraint_list_4gp (Constraint_list *CL)
{
	if ( get_nproc()==1)return  nfork_relax_constraint_list_4gp (CL);
	else if (strstr ( CL->multi_thread, "relax"))return fork_relax_constraint_list_4gp (CL);
	else return nfork_relax_constraint_list_4gp (CL);
}

Constraint_list * fork_relax_constraint_list_4gp (Constraint_list *CL)
{
	return nfork_relax_constraint_list_4gp(CL);
}
Constraint_list * nfork_relax_constraint_list_4gp (Constraint_list *CL)
{
	int s, r, a, t_s, t_r;
	int th=0;
	static int *entry;
	char *tmp;
	Sequence *S=CL->S;
	FILE *fp;
	
	tmp=vtmpnam(NULL);
	entry=(int*)vcalloc (CL->entry_len, sizeof (int));

	if (!CL || !CL->residue_index)return CL;

	fp=vfopen (tmp, "w");
	for (s=0; s< S->nseq; s++)
	{
		for ( r=1; r<=S->len[s]; r++)
		{
			entry[SEQ1]=s;
			entry[R1]=r;

			for ( a=1; a<CL->residue_index[s][r][0]; a+=ICHUNK)
			{

				t_s=CL->residue_index[s][r][a+SEQ2];
				t_r=CL->residue_index[s][r][a+WE];
				fprintf (fp, "%d ",residue_pair_extended_list_pc (CL,s, r,t_s, t_r));
			}
		}
	}
	vfclose (fp);
	fp=vfopen (tmp, "r");

	for (s=0; s< S->nseq; s++)
	{
		for ( r=1; r<=S->len[s]; r++)
		{
			entry[SEQ1]=s;
			entry[R1]=r;

			for ( a=1; a<CL->residue_index[s][r][0]; a+=ICHUNK)
			{
				fscanf (fp, "%d ", &CL->residue_index[s][r][a+WE]);
			}
		}
	}
	vfree (entry);
	vfclose (fp);
	return filter_constraint_list (CL,WE,th);
}

int constraint_list2max  ( Constraint_list *CL)
{
	int max=0;
	int b, s1, r1, s2 ,r2 ,w2;
	Sequence *S=CL->S;
	
	for (s1=0; s1<S->nseq; s1++)
	{
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][b+SEQ2];
				r2=CL->residue_index[s1][r1][b+R2];
				w2=CL->residue_index[s1][r1][b+WE];
				max=MAX(w2,max);
			}
		}
	}
	return max;
}
int constraint_list2ne  ( Constraint_list *CL)
{
	int max=0;
	int s1,r1;
	Sequence *S;
	
	if (!CL || !CL->residue_index || !CL->S)return 0;

	S=CL->S;
	for (s1=0; s1<S->nseq; s1++)
	{
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			max+=CL->residue_index[s1][r1][0]-1;
		}
	}
	return max/ICHUNK;
}



float constraint_list2connectivity (Constraint_list *CL)
{
	float **mat;
	float tot=0;
	float ntot=0;
	int s1, s2,r1,b;
	Sequence *S=CL->S;

	mat=declare_float (S->nseq, S->nseq);
	for (s1=0; s1<S->nseq; s1++)
		for ( r1=1;r1<=S->len[s1]; r1++)
		{
			for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{
				mat[s1][CL->residue_index[s1][r1][b+SEQ2]]++;
			}
		}
		for (s1=0; s1<S->nseq; s1++)
			for ( s2=0; s2<S->nseq; s2++)
			{
				if ( s1==s2)continue;
				mat[s1][s2]*=(float)2;
				mat[s1][s2]/=(float)(S->len[s1]+S->len[s2]);
				tot+=mat[s1][s2];
				ntot++;
			}
			free_float (mat, -1);
		tot=(float)tot/(float)(ntot);
		return tot;
}




int constraint_list2avg ( Constraint_list *CL)
{
	int max=0, n=0;
	int b, s1, r1, s2, r2, w2;
	Sequence *S=CL->S;
	
	for (s1=0; s1<S->nseq; s1++)
	{
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][b+SEQ2];
				r2=CL->residue_index[s1][r1][b+R2];
				w2=CL->residue_index[s1][r1][b+WE];
				max+=w2;
				n++;
			}
		}
	}
	return (n==0)?0:max/n;
}

Constraint_list * constraint_list2sub_constraint_list (Constraint_list *CL, Sequence *SMALL)
{
  char *tmp_in=vtmpnam(NULL);
  char *tmp_out=vtmpnam(NULL);
  char *buf=NULL;
  FILE *out;
  FILE *in;
  Sequence *S;
  int a,ns,i;
  int *lu;
  int do_print=0;
  
  S=CL->S;
  lu=(int*)vcalloc (S->nseq, sizeof(int));
  
  for (ns=0,a=0; a<SMALL->nseq; a++)
    {
      if ((i=name_is_in_list (SMALL->name[a], S->name, S->nseq,100))!=-1)lu[i]=++ns;
    }
  
 
  out=vfopen (tmp_out, "w");
  fprintf (out,"%d\n", ns);
  for (a=0; a<S->nseq; a++)
    {
      if (lu[a])fprintf (out, "%s %d %s\n",S->name[a],S->len[a], S->seq[a]);
    }
  
  vfclose (save_constraint_list ( CL, 0, CL->ne,tmp_in, NULL, "ascii",CL->S));
  in=vfopen (tmp_in,"r");

  do_print=0;
  while ((buf=vfgets (buf,in))!=NULL)
    {
      int s1, s2;
      if (buf[0]=='#')
	{
	  sscanf(buf, "#%d %d", &s1, &s2);
	  if (lu[s1-1] && lu[s2-1])
	    {
	      do_print=1;
	      fprintf (out, "#%d %d\n", lu[s1-1], lu[s2-1]);
	    }
	  else 
	    do_print=0;
	}
      else if (buf[0]=='!' || do_print==1)fprintf (out, "%s", buf);
    }
  vfclose (out);
  vfclose (in);
  return read_constraint_list_file (NULL,tmp_out);
}
  

/**
 * Deletes unimportant edges from the Constraint_list.
 *
 * This function reduces the size of th eConstraint_list, which is useful to speed
 * up calculations such as the extension step, by kicking out all constraints/edges
 * with a weight smaller or equal to a given threshold.
 * Currently there is a hardcoded value of 10 in ::fork_relax_constraint_list, but it
 * can also be invokled before the relaxation by setting the \c -filter_lib parameter.
 *
 * \param[in,out] CL    The global Constraint_list object.
 * \param[in]     field Internal, constant number to access the CHUNKS of Constraint_list::residue_index correctly.
 * \param[in]     T     threshold: Edges with weight smaller or equal to T will be removed.
 */
Constraint_list * filter_constraint_list (Constraint_list *CL,int field, int T)
{
	int s1,r1,b,c,d;
	Sequence *S;

	if (!CL || !CL->residue_index || !CL->S)return CL;
	S=CL->S;
	CL->ne=0;
	for (s1=0; s1<S->nseq; s1++)
	{
		for ( r1=1; r1<=S->len[s1]; r1++)
		{
			int *r=CL->residue_index[s1][r1];
			for (d=1,b=1; b<r[0]; b+=ICHUNK)
			{
				if (r[b+field]>T)
				{
					if (d!=b)for (c=0; c<ICHUNK; c++){r[d+c]=r[b+c];}
					d+=ICHUNK;
					CL->ne++;
				}
			}
			r[0]=d;
			CL->residue_index[s1][r1]=(int*)vrealloc (r, d*sizeof(int));;
		}
	}

	return CL;
}


Constraint_list * shrink_constraint_list (Constraint_list *CL)
{
	int a, b, n, tot;
	Constraint_list *CL2;
	Alignment *A, *B;
	int *ns, **ls;

	if (!CL || !CL->S || !CL->residue_index) return CL;

	ns=(int*)vcalloc (2, sizeof (int));
	ls=declare_int ((CL->S)->nseq, 2);

	A=seq2aln (CL->S,NULL, RM_GAP);
	B=seq2aln (CL->S,NULL, RM_GAP);
	CL2=declare_constraint_list (CL->S,NULL, NULL, 0,NULL, NULL);
	n=(CL->S)->nseq;
	tot=((n*n)-n)/2;
	fprintf ( CL->local_stderr, "\n\n\tSHRINK Constraint List [%d element(s)]", CL->ne);
	for (n=0,a=0; a<(CL->S)->nseq-1; a++)
		for (b=a+1; b<(CL->S)->nseq; b++, n++)
		{
			output_completion (CL->local_stderr,n, tot, 100, "slow_pair");
			ns[0]=ns[1]=1;
			ls[0][0]=a;
			ls[1][0]=b;
			ungap (A->seq_al[a]);
			ungap (A->seq_al[b]);
			linked_pair_wise (A, ns, ls, CL);
			B->seq_al[0]=A->seq_al[a];
			B->seq_al[1]=A->seq_al[b];
			sprintf (B->name[0], "%s", A->name[a]);
			sprintf (B->name[1], "%s", A->name[b]);
			B->nseq=2;
			B->len_aln=strlen (B->seq_al[0]);
			CL2=aln2constraint_list (B, CL2, "sim");
		}
		CL->ne=CL2->ne;
		return CL;
}


Constraint_list *aln_file2constraint_list (char *alname, Constraint_list *CL,char *weight_mode)
{
	Alignment *A;
	A=main_read_aln ( alname, NULL);
	CL=aln2constraint_list (A, CL, weight_mode);
	free_aln (A);
	return CL;
}

int *seqpair2weight (int s1, int s2, Alignment *A,Constraint_list *CL, char *weight_mode, int *weight)
{
	int *col;
	int a,c, ref_weight;


	if ( !weight)weight=(int*)vcalloc (MAX(2,A->len_aln), sizeof (int));

	weight[0]=FORBIDEN;
	if ( weight_mode==NULL || strcmp (weight_mode, "no")==0 || is_number (weight_mode))
	  {

	    if (is_number (weight_mode))ref_weight=atoi(weight_mode);
	    else ref_weight=1;
	    weight[1]=ref_weight;

	  }
	else if ( strstr ( weight_mode, "cons"))
	  {
	    ref_weight=weight[1]=1000;
	  }
	else if ( strstr ( weight_mode, "OW"))
	  {
	    int ow;
	    sscanf ( weight_mode, "OW%d", &ow);
	    weight[1]=ow*get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-", NULL);

	  }
	else if ( strstr ( weight_mode, "len"))
	  {
	    weight[1]=A->len_aln;
	  }
	else if ( strstr (weight_mode, "winsim"))
	  {
	    weight=get_seq_winsim ( A->seq_al[s1], A->seq_al[s2], "-", weight_mode+6, weight);
	  }
	else if ( strstr ( weight_mode, "sim") || strstr (weight_mode, "default"))
	  {
	    char *sim_mode;
	    if ( strstr(weight_mode, "sim"))sim_mode=strstr(weight_mode, "sim")+3;
	    else sim_mode=NULL;
	    ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-", sim_mode);
	    if (ref_weight == 0)
	      ref_weight = 1;
	    weight[1]=ref_weight;

	  }
	else if ( strstr ( weight_mode, "subset"))
	  {
	    ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-",NULL);
	    weight[1]=ref_weight;
	  }
	else if (  strstr ( weight_mode, "cdna"))
	  {
		ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-", weight_mode+4);
		col=(int*)vcalloc ( A->len_aln+1, sizeof (int));
		if (A->cdna_cache)
			for ( a=0; a<=A->len_aln; a++)col[a]=A->cdna_cache[0][a];
			else
				for ( a=0; a<=A->len_aln; a++)col[a]=1;
				for ( c=0; c< A->len_aln; c++)weight[c]=ref_weight*col[c];
				vfree (col);
	  }

	else if ( strstr (weight_mode, "overaln"))
	  {
	    ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-","idmat");
	    //weight=pw_aln2clean_aln_weight (A->seq_al[s1], A->seq_al[s2], ref_weight,0, 0, 0, 0, NULL);
	    printf_exit (EXIT_FAILURE, stderr,"ERROR: mode overaln not currently supported [FATAL:%s]", PROGRAM);
	  }
	else
	  {
	    fprintf ( stderr, "\nERROR: Weight Mode %s is unknown [FATAL:%s]", weight_mode, PROGRAM);
	    crash ("");
	  }
	return weight;
}


Constraint_list *aln2constraint_list_full    (Alignment *A, Constraint_list *CL,char *in_weight_mode)
  {
  return aln2constraint_list_generic    (A,CL, in_weight_mode, 0);
  }
Constraint_list *aln2constraint_list    (Alignment *A, Constraint_list *CL,char *in_weight_mode)
  {
  return aln2constraint_list_generic    (A,CL, in_weight_mode, 1);
  }
Constraint_list *aln2constraint_list_generic    (Alignment *A, Constraint_list *CL,char *in_weight_mode, int top)
{
	int a, b, c;
	int *weight=NULL;
	int s1, s2;
	int fixed_nres1, fixed_nres2;
	int do_pdb=0;
	int pdb_weight=0;
	int set_misc;
	char *alp=NULL, *p;
	char weight_mode [100];
	int *cache;

	int *entry;
	int **fixed;
	entry=(int*)vcalloc (CL->entry_len+1, sizeof (int));


	//MSA are now read as one to all (+ extra bits if needed) libraries
	//if ( atoigetenv ("TOP4TC")){HERE ("TOP=1");top=1;}

	sprintf ( weight_mode , "%s", (!in_weight_mode || strm (in_weight_mode, "default"))?"sim":in_weight_mode);

	if ( !A)return CL;
	if (strstr (weight_mode, "extend"))
	  {
	    extend_seqaln(NULL, A);
	    extend_seqaln(CL->S,NULL);
	    if (CL->S!=A->S)extend_seqaln(A->S,NULL);
	  }

	fixed=fix_aln_seq_new (A, (CL->S));

	if ( !CL)
	  {
	    Sequence *S;
	    S=aln2seq (A);
	    CL=declare_constraint_list (S,NULL, NULL, 0,NULL, NULL);
	    CL->S=S;
	  }

	cache=(int*)vcalloc (A->len_aln, sizeof (int));

	if ( (p=strstr (weight_mode, "_subset_")))
	  {
	    alp=strchr (weight_mode, '_')+1;
	    p[0]='\0';
	  }

	for ( a=0; a<A->nseq-1; a++)
	  {
	    if ((s1=fixed[a][0])==-1)continue;

	    for (set_misc=0,b=a+1; b< A->nseq; b++)
	      {

		int use_pair;
		int nres1=0;
		int nres2=0;

		if ((s2=fixed[b][0])==-1)continue;
		weight=seqpair2weight (a, b, A, CL, weight_mode, weight);
		for (c=0; c< A->len_aln; c++)
		  {
		    int isgap1, isgap2;
		    isgap1=is_gap(A->seq_al[a][c]);
		    isgap2=is_gap(A->seq_al[b][c]);
		    nres1+=!isgap1;
		    nres2+=!isgap2;

		    if (cache[c]==-1 && top)continue;
		    if (!isgap1)cache[c]=1;
		    if (cache[c] && b==A->nseq-1)cache[c]=-1;

		    use_pair=1;
		    use_pair=use_pair && !is_gap(A->seq_al[a][c]);
		    use_pair=use_pair && !is_gap(A->seq_al[b][c]);
		    use_pair=use_pair && A->seq_al[b][c]!=UNDEFINED_RESIDUE;
		    use_pair=use_pair && A->seq_al[a][c]!=UNDEFINED_RESIDUE;
		    use_pair=use_pair && !(do_pdb && pdb_weight==0);
		    use_pair=use_pair && ((weight[0]==FORBIDEN)?weight[1]:weight[c]);

		    if (alp)use_pair=use_pair && is_in_set (A->seq_al[b][c], alp) && is_in_set (A->seq_al[a][c], alp);

		    if (use_pair)
		      {
			if ((fixed_nres1=fixed[a][nres1])<0)continue;
			if ((fixed_nres2=fixed[b][nres2])<0)continue;

			entry[SEQ1]=s1;
			entry[SEQ2]=s2;
			entry[R1]=fixed_nres1;
			entry[R2]=fixed_nres2;
			entry[CONS]=1;
			//lw=(do_pdb)?(NORM_F/MAXID)*pdb_weight:(((weight[0]==FORBIDEN)?weight[1]:weight[c]));
			//lw=(lw<0.1)?10:lw;
			//entry[WE]=(NORM_F/MAXID)*lw;
			if (do_pdb)entry[WE]=(NORM_F/MAXID)*pdb_weight;
			else entry[WE]=(NORM_F/MAXID)*((weight[0]==FORBIDEN)?weight[1]:weight[c]);
			add_entry2list (entry, CL);
		      }
		  }
	      }
	  }
	vfree (entry);
	vfree (cache);
	vfree (weight);
	free_int (fixed, -1);

	if (strstr (weight_mode, "extend"))
	  {
	    unextend_seqaln(NULL, A);
	    unextend_seqaln(CL->S,NULL);
	    if (A->S!=CL->S)unextend_seqaln(A->S,NULL);



	  }
	if (A->A)
	  {
	    return aln2constraint_list (A->A, CL, weight_mode);
	  }
	else
	  {
	    return CL;
	  }
}






int **list2residue_total_weight ( Constraint_list *CL)
{
	/*Returns
	tot_weight[nseq][maxlen]
	where each residue is associated with the total of its weights in CL
	####IMPORTANT

	-the numbering of the residues  goes from 1 to L:
	-the numbering of the sequences goes from 0 to N-1:
	*/

	int **tot_weight;
	int a,b,s1,r1,s2,r2,w2;
	Sequence *S=CL->S;

	if ( !CL || !CL->S || (CL->S)->nseq==0)return NULL;
	S=CL->S;
	tot_weight=declare_int2 (S->nseq, S->len, 1);
	for (s1=0; s1<S->nseq; s1++)
	{
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][b+SEQ2];
				r2=CL->residue_index[s1][r1][b+R2];
				w2=CL->residue_index[s1][r1][b+WE];
				tot_weight[s1][r1]+=w2;
				tot_weight[s2][r2]+=w2;
			}
		}
	}
	for (a=0; a<(CL->S)->nseq; a++)
		for (b=0; b<(CL->S)->len[a]; b++)
			tot_weight[a][b]/=(CL->S)->nseq;

		return tot_weight;
}



/*********************************************************************/
/*                                                                   */
/*                         LIST FUNCTIONS                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

Constraint_list *merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode)
{

  int **cache=NULL;
  int s1,r1,s2,r2,a,b;
  int *entry;

  entry=(int*)vcalloc (ICHUNK+10, sizeof (int));

  for (a=0;a<ICHUNK+10; a++)entry[a]=0;
  if (SL->S!= ML->S)cache=fix_seq_seq((SL->S),(ML->S));

  if (!ML || !SL)return ML;
  for (s1=0; s1<(SL->S)->nseq; s1++)
    {
      if (cache && cache[s1][0]==-1)continue;
      for (r1=1; r1<=((SL)->S)->len[s1]; r1++)
	{
	  entry[SEQ1]=(cache)?cache[s1][0]:s1;
	  entry[R1]=(cache)?cache[s1][r1]:r1;
	  if (entry[R1]<=0)continue;
	  b=(SL->freeze)?SL->freeze[s1][r1]:1;
	  for (;b<SL->residue_index[s1][r1][0]; b+=ICHUNK)
	    {
	      s2=SL->residue_index[s1][r1][b+SEQ2];
	      r2=SL->residue_index[s1][r1][b+R2];
	      entry[SEQ2]=(cache)?cache[s2][0]:s2;
	      if ( entry[SEQ2]==-1)continue;
	      else
		{
		  entry[R2]=(cache)?cache[s2][r2]:r2;
		  if (entry[R2]<=0)continue;
		  else
		    {
		      entry[WE]=SL->residue_index[s1][r1][b+WE];
		      entry[CONS]=SL->residue_index[s1][r1][b+CONS];
		      entry[MISC]=SL->residue_index[s1][r1][b+MISC];
		    }
		  add_entry2list2(entry, ML);
		}
	    }
	}
    }
  vfree (entry);
  if ( cache)free_int (cache, -1);
  return ML;
}




int ** seq2defined_residues ( Sequence *S, Constraint_list *CL)
{
	int **seq_count;
	int s1, r1, a;

	seq_count=declare_int (S->nseq, S->max_len+1);
	for (s1=0; s1<S->nseq; s1++)
	{
		for ( r1=1; r1<=S->len[s1]; r1++)
		{
			for (a=1; a<CL->residue_index[s1][r1][0]; a+=ICHUNK)
			{
				int s2=CL->residue_index[s1][r1][a+SEQ2];
				int r2=CL->residue_index[s1][r1][a+R2];
				seq_count[s1][r1]++;
				seq_count[s2][r2]++;
			}
		}
	}
	return seq_count;
}

int ** aln2defined_residues ( Alignment *A, Constraint_list *CL)
{

	int **seq_count;
	int **aln_count;
	int **pos;
	int a,ra, b;

	pos=aln2pos_simple(A, A->nseq);
	seq_count=seq2defined_residues(CL->S, CL);
	aln_count=declare_int (A->nseq, A->len_aln);
	for (a=0; a< A->nseq; a++)
	{
		ra=name_is_in_list(A->name[a], (CL->S)->name, (CL->S)->nseq, 100);
		if ( ra==-1) continue;
		for ( b=0; b<A->len_aln; b++)
			if (pos[a][b]>0 && seq_count[ra][pos[a][b]]>0)aln_count[a][b]=1;
	}

	free_int (seq_count, -1);
	free_int (pos,-1);
	return aln_count;
}

/*********************************************************************/
/*                                                                   */
/*                        DEBUG CONSTRAINT_LIST                       */
/*                                                                   */
/*                                                                   */
/*********************************************************************/


/*********************************************************************/
/*                                                                   */
/*                        PRUNE CONSTRAINT_LIST                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/


char * list2prune_list ( Sequence *S, int **sm)
{
	int a, b, c;
	int **mat, *used, *keep;
	int nk=0, n=0;
	int ns=4;
	char *file;
	FILE *fp;

	n=S->nseq;

	if (get_string_variable("prune_lib_mode"))
		ns=atoi(get_string_variable("prune_lib_mode"));


	if (ns==0)ns=n;
	else if (ns<0)ns=-(n*ns)/100;
	else if (ns>=n)ns=n;





	keep=(int*)vcalloc (n, sizeof (int));
	used=(int*)vcalloc (n, sizeof (int));
	mat=declare_int (n, n);
	file=vtmpnam (NULL);

	//1-Identify the seed sequence: the one on average the further away from the rest
	for (a=0; a<n;a++)
	{
		mat[a][0]=a;
		for (b=0; b<n; b++)
		{
			mat[a][1]+=sm[a][b];
		}
		sort_int_inv (mat, 2, 1, 0, n-1);
	}

	keep[nk++]=mat[0][0];
	used[mat[0][0]]=1;

	for (a=1; a<ns; a++)
	{
		for (c=0; c<n; c++)
		{
			mat[c][0]=c;
			if (used[c]){mat[c][1]=n*1000;continue;}
			for (b=0; b<nk; b++)
			{
				mat[c][1]+=sm[c][keep[b]];
			}
		}
		sort_int (mat, 2, 1, 0,n-1);
		keep[nk++]=mat[0][0];
		used[mat[0][0]]=1;
	}

	fp=vfopen (file, "w");
	for ( a=0; a<nk; a++)
	{
		for (b=0; b<S->nseq; b++)
			if (keep[a]!=b)
			{
				fprintf ( fp, "\n2 %d %d", keep[a], b);
			}
	}
	vfclose (fp);
	vfree (keep); vfree (used);free_int (mat, -1);

	return file;
}




/*********************************************************************/
/*                                                                   */
/*                        WEIGHT CONSTRAINT_LIST                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/


/**
 * Specify a weight for each sequence in the Constraint_list.
 *
 * If \c seq_weight is "t_coffee", weights (see ::Weights object) are computed according to
 * ::compute_t_coffee_weight and given to the Constraint_list::W reference.
 * If \c seq_weight is a file, this function will read the weight file (::read_seq_weight).
 * In every other case, all sequence weights are set to 1 and the weight mode Weight::mode is
 * set to "no_seq_weight".
 *
 * Unless the third case occurs, so whenever you want the sequences to be weighted,
 * the function ::re_weight_constraint_list takes care of it.
 *
 * \param[in] seq_weight can be "t_coffee" or the name of a weight file.
 * \param[in,out] CL points to the global Constraint_list object.
 */
Constraint_list *weight_constraint_list(Constraint_list * CL, char *seq_weight)

{
	Weights *W;

	if ( CL->ne==0){return CL;}
	else if ( strm(seq_weight, "t_coffee")) {W=compute_t_coffee_weight(CL);}
	else if (check_file_exists (seq_weight))
	{
		W=read_seq_weight ((CL->S)->name, (CL->S)->nseq, seq_weight);
	}
	else
	  {
	    int a;
	    W=declare_weights((CL->S)->nseq);
	    sprintf ( W->mode, "no_seq_weight");
	    for ( a=0; a<(CL->S)->nseq; a++)
	      {
		sprintf ( W->seq_name[a], "%s", (CL->S)->name[a]);
		W->SEQ_W[a]=1;
	      }
	    CL->W=W;
	    return CL;
	  }

	CL=re_weight_constraint_list (CL,W);

	CL->W=W;


	return CL;


}



Weights* compute_t_coffee_weight(Constraint_list * CL)
{
	int a, b;
	float p, d;
	Weights *W;
	int nseq;

	nseq=(CL->S)->nseq;
	W=declare_weights(nseq);
	sprintf ( W->mode, "t_coffee");
	for ( a=0; a< nseq; a++)
	{
		sprintf ( W->seq_name[a], "%s", (CL->S)->name[a]);
		W->SEQ_W[a]=1;
	}

	for (a=0; a< (CL->S)->nseq-1; a++)
		for ( b=a+1; b< (CL->S)->nseq; b++)
		{
			if ( b==a){d=1;}
			else if ( !(CL->S)->len[b] || !(CL->S)->len[a])d=1;
			else
			{
				d=((float)(CL->DM)->similarity_matrix[a][b]/MAXID)*10;
			}
			p=pow(d,3);

			W->SEQ_W[a]+=p;
			W->SEQ_W[b]+=p;

		}

		for ( p=0,b=0; b< (CL->S)->nseq; b++)
		{
			if ((CL->S)->len[b]==0)W->SEQ_W[b]=0;
			else W->SEQ_W[b]=2/W->SEQ_W[b];
			p+=W->SEQ_W[b];
		}
		for ( b=0; b< (CL->S)->nseq; b++)
		{
			W->SEQ_W[b]=W->SEQ_W[b]*((float)W->nseq/p);
		}


		return W;
}

Constraint_list *set_weight4constraint_list(Constraint_list * CL,int w)
{
	int a, s, r;
	Sequence *S;

	if ( !CL || !CL->S || !CL->residue_index)return CL;
	S=CL->S;
	for (s=0; s<S->nseq; s++)
	{
		for ( r=1; r<=S->len[s]; r++)
		{
			for (a=1; a<CL->residue_index[s][r][0]; a+=ICHUNK)
			{
				CL->residue_index[s][r][a+WE]*=w;
			}
		}
	}
	CL=evaluate_constraint_list_reference (CL);
	return CL;
}
Constraint_list *re_weight_constraint_list(Constraint_list * CL,Weights *W)
{
	float *weight;
	int s, r, a;
	Sequence *S=CL->S;

	weight=W->SEQ_W;

	for (s=0; s<S->nseq; s++)
	{
		for ( r=1; r<=S->len[s]; r++)
		{
			for (a=1; a<CL->residue_index[s][r][0]; a+=ICHUNK)
			{
				int t_s=CL->residue_index[s][r][a];
				CL->residue_index[s][r][a+WE]*=MIN(weight[s], weight[t_s]);
			}
		}
	}
	CL=evaluate_constraint_list_reference (CL);
	return CL;
}

Distance_matrix* cl2distance_matrix    (Constraint_list *CL, Alignment *A, char *in_mode, char *in_sim_mode, int print)
{
	char mode[100];
	char sim_mode [100];

	if ( !CL)return NULL;
	sprintf ( mode, "%s", (CL && in_mode==NULL)?CL->distance_matrix_mode:in_mode);
	sprintf ( sim_mode, "%s", (CL && in_sim_mode==NULL)?CL->distance_matrix_sim_mode:in_sim_mode);

	if ( !CL->DM ||!strm ((CL->DM)->mode, mode) || !strm ((CL->DM)->sim_mode, sim_mode) || A )
	  {
	    return seq2distance_matrix (CL, A, mode, sim_mode, print);
	  }
	else
	  {

	    return CL->DM;
	  }
}


Distance_matrix *seq2distance_matrix (Constraint_list *CL, Alignment *A,char *mode, char *sim_mode, int print)
{
	/*Compute the distance matrix associated with the Constraint List and the sequences*/
	/*Computation only occurs if the similiraty matrix is undefined : CL->similarity_matrix*/
	/*Undefine  CL->similarity_matrix to force computation*/

	int a, b;
	Alignment *B;
	Constraint_list *NCL;
	float score=0;
	int *ns;
	int **l_s;
	float id;
	int max_name=0;
	int  id_score;
	static float **g_matrix;
	float ref=0;
	int n_coor=0;
	Distance_matrix *DM;
	int **sim_table=NULL;

	//mode: computation mode
	//sim_mode: mode for computing the similarity

	//Composite modes

	if (strm (mode, "ktup2"))
	  {
	    B=seq2aln ( CL->S, NULL, 1);
	    B=very_fast_aln (B, B->nseq,NULL);
	    sprintf ( CL->distance_matrix_mode, "aln");
	    DM=cl2distance_matrix (CL, B, NULL, NULL, 1);
	    sprintf ( CL->distance_matrix_mode, "ktup2");
	    sprintf ( DM->mode, "%s", mode);
	    sprintf ( DM->sim_mode, "%s", sim_mode);
	    free_aln (B);
	    return DM;
	  }

	if ( !CL) return NULL;
	else
	  {
	    for ( max_name=0,a=0; a<  (CL->S)->nseq; a++)max_name=MAX(strlen ((CL->S)->name[a]), max_name);


	    if ( CL->DM)DM=CL->DM;
	    else
	      {
		DM=(Distance_matrix*)vcalloc ( 1, sizeof (Distance_matrix));
		DM->nseq=(CL->S)->nseq;
		DM->similarity_matrix=declare_int ( (CL->S)->nseq, (CL->S)->nseq);
		DM->distance_matrix  =declare_int ( (CL->S)->nseq, (CL->S)->nseq);
		DM->score_similarity_matrix=declare_int ( (CL->S)->nseq, (CL->S)->nseq);
		}

	    sprintf ( DM->mode, "%s", mode);
	    sprintf ( DM->sim_mode, "%s", sim_mode);

	    NCL=duplicate_constraint_list_soft (CL);
	    NCL->pw_parameters_set=1;

	    if (!A)
	      {
		if ( CL->tree_aln)B=CL->tree_aln;
		else B=seq2aln ( NCL->S, NULL, 1);
	      }
	    else
	      {
		B=copy_aln (A, NULL);
		B=reorder_aln (B, (CL->S)->name, (CL->S)->nseq);
	      }

	    if ( strm (mode, "very_fast"))
	      {
		sprintf ( NCL->dp_mode, "very_fast_pair_wise");
		NCL->evaluate_residue_pair=evaluate_matrix_score;
		if ( strm ((CL->S)->type, "DNA") ||strm ((CL->S)->type, "RNA")  )
		  {
				NCL->M=read_matrice ("idmat");
				NCL->gop=-10;
				NCL->gep=-1;
				CL->ktup=6;
		  }
		else
		  {
		    NCL->M=read_matrice ("blosum62mt");
		    NCL->gop=get_avg_matrix_mm (NCL->M, AA_ALPHABET)*10;
		    NCL->gep=-1;
		    CL->ktup=2;
		  }
		NCL->use_fragments=1;
		CL->diagonal_threshold=6;
	      }

	    else if ( strm (mode, "ktup"))
	      {

		NCL->ktup=6;
		sim_table=ktup_dist_mat((CL->S)->seq,(CL->S)->nseq,NCL->ktup, (CL->S)->type);
	      }


	    else if (strm (mode, "aln"))
	      {

		sim_table=aln2sim_mat (A, sim_mode);
	      }
	    else if ( strm (mode, "fast") || strm ("idscore", mode))
	      {
		sprintf ( NCL->dp_mode, "myers_miller_pair_wise");
		NCL->evaluate_residue_pair=evaluate_matrix_score;
		if ( strm ((CL->S)->type, "DNA") || strm ((CL->S)->type, "RNA"))
		  {
		    NCL->M=read_matrice ("idmat");
		    NCL->gop=-10;
		    NCL->gep=-1;
		  }
		else
		  {
		    NCL->M=read_matrice ("blosum62mt");
		    NCL->gop=get_avg_matrix_mm (NCL->M, AA_ALPHABET)*10;
		    NCL->gep=-1;
		  }
	      }

	    else if ( strm (mode, "cscore"))
	      {
		if (!CL || !CL->residue_index || CL->ne==0)
		  return seq2distance_matrix (CL, A,"idscore",sim_mode, print);
	      }
	    else if ( strm (mode, "geometric") );
	    else if (strm (mode, "slow"));
	    else if (strm (mode, "clustalw"));
	    else if (strm (mode, "no"))
	      print=1;
	    else if (strm (mode, "random"))
	      print=1;
	    else
	      {
			fprintf ( stderr, "\nError: %s is an unknown distance_matrix_mode [FATAL:%s]", mode,PROGRAM);
			crash ("");
	      }

		//Special Geometric Mode
	    if ( strm (NCL->distance_matrix_mode, "geometric"))
	      {
		free_arrayN(g_matrix, 2);
		g_matrix=declare_float ((CL->S)->nseq, 3);
		n_coor=MIN(3,((CL->S)->nseq));

		for ( a=0; a<(CL->S)->nseq; a++)
			{
			  for (b=0; b<n_coor; b++)
			    {
			      B=align_two_sequences ((CL->S)->seq[a], (CL->S)->seq[b], "pam250mt", -10, -1, "fasta_pair_wise");
			      g_matrix[a][b]=get_seq_sim ( B->seq_al[0], B->seq_al[1], "-", NULL);
			      free_aln(B);B=NULL;
			    }
			}
		ref=(float)sqrt((double)(10000*n_coor));
	      }


	    ns=(int*)vcalloc ( 2, sizeof(int));
	    l_s=declare_int ( 2, 1);
	    ns[0]=ns[1]=1;
	    l_s[0][0]=0;
	    l_s[1][0]=1;

	    if (CL->local_stderr && print>0)fprintf ( (CL->local_stderr), "\nCOMPUTE PAIRWISE SIMILARITY [dp_mode: %s] [distance_matrix_mode: %s][Similarity Measure: %s] \n", (NCL->dp_mode)?NCL->dp_mode:"NO_DP",mode, sim_mode);

	    for (a=0; a< (CL->S)->nseq; a++)
	      {
		if (CL->local_stderr && print>0)fprintf ( (CL->local_stderr), "\r\tSeq: %5d %20s -- [%3d %%]",a, (CL->S)->name[a], (int)((a*100)/(CL->S)->nseq));
		for ( b=a; b< (CL->S)->nseq; b++)
		  {
		    if ( b==a){DM->similarity_matrix[a][b]=MAXID;}
		    else
		      {
			l_s[0][0]=a;
			l_s[1][0]=b;
			if ( !strm(mode, "ktup2") && ! strm (mode, "geometric"))
			  {
			    ungap ( B->seq_al[a]);
			    ungap ( B->seq_al[b]);
			  }

			if ( strm (mode, "slow"))
			  {

			    B->score_aln=pair_wise (B, ns, l_s,NCL);

			    id=get_seq_sim ( B->seq_al[a], B->seq_al[b], "-", sim_mode);
			    if ( CL->residue_index)
			      {
				score=(int)(((float)B->score_aln)/(B->len_aln*SCORE_K));
				score=(int)(CL->normalise)?((score*MAXID)/(CL->normalise)):(score);
			      }
				else if ( CL->M)score=id;


			    if ( score>MAXID)score=(CL->residue_index)?sub_aln2sub_aln_score (B, CL, CL->evaluate_mode, ns, l_s):id;

			  }
			else if ( strm2 (mode,"fast", "very_fast"))
			  {
			    B->score_aln=pair_wise (B, ns, l_s,NCL);
			    id=get_seq_sim ( B->seq_al[a], B->seq_al[b], "-", sim_mode);
			    score=(int)(id)*SCORE_K;
			  }
			else if ( strm (mode, "cscore"))
			  {
			    ungap ( B->seq_al[a]);
				ungap ( B->seq_al[b]);
				score=(int)linked_pair_wise (B, ns, l_s, NCL);


				score/=(B->len_aln*SCORE_K);
				id=score/SCORE_K;
			      }
			    else if ( strm (mode, "idscore"))
			      {
				score=id=idscore_pairseq (B->seq_al[a], B->seq_al[b], NCL->gop, NCL->gep, NCL->M, sim_mode);
			      }
			    else if (strm (mode, "ktup"))
			      {
				id=sim_table[a][b];
				score=id*SCORE_K;

			      }
			    else if (strm (mode, "aln"))
			      {
				score=id=sim_table[a][b];
				score*=SCORE_K;
			      }

			    else if ( strm (mode, "geometric"))
			      {
				id=get_geometric_distance (g_matrix,n_coor, a, b, "euclidian");
				id=MAXID*(1-((id/ref)));
				score=(int)(id)*SCORE_K;
			      }
			    else if ( strm (mode, "no"))
			      {
				id=100;
				score=id*SCORE_K;
			      }
			    else if ( strm (mode, "random"))
			      {
				id=rand()%100;
				score=id*SCORE_K;
			      }
			    else
			      {
				id=B->score_aln=pair_wise (B, ns, l_s,NCL);
				score=id*SCORE_K;
			      }
			    /*Sim mat*/
			    DM->similarity_matrix[a][b]=DM->similarity_matrix[b][a]=(int)(id);
			    /*Dist mat*/
			    DM->distance_matrix[a][b]=DM->distance_matrix[b][a]=MAXID-(int)(id);
			    /*Score mat*/


			    DM->score_similarity_matrix[a][b]=DM->score_similarity_matrix[b][a]=(int)score;
			    id_score=id;
			    if (CL->local_stderr && print>1) fprintf (CL->local_stderr, "\n\t%-*s %-*s identity=%3d%% score=%3d [%3d %%]", max_name,(CL->S)->name[a], max_name,(CL->S)->name[b], id_score, (int)score,(int)((a*100)/(CL->S)->nseq));
			  }
		      }
		  }
		vfree (ns);
		free_int(l_s, -1);

	}


	if (CL->local_stderr) fprintf (CL->local_stderr, "\n");
	free_constraint_list (NCL);



	if (!CL->tree_aln)
	{
		free_aln (B);
	}

	free_int (sim_table, -1);



	return DM;
}
/*********************************************************************/
/*                                                                   */
/*                        RNA FUNCTIONS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char * seq2rna_lib ( Sequence *S, char *name)
{
	int a;
	FILE *fp;


	if (!name)name=vtmpnam (NULL);
	fp=vfopen (name, "w");
	for ( a=0; a<S->nseq; a++)
	{

		fprintf (fp, "%s\n", rna_struc2rna_lib(S->name[a], S->seq[a], NULL));
	}
	vfclose (fp);

	return name;
}

Constraint_list *read_contact_lib ( Sequence *S, char *fname, Constraint_list *R)
{
  char **list;
  int n=0,a;
  if (!S && !R)S=constraint_list2seq(fname);
  if (!R)R=declare_constraint_list ( S,NULL, NULL, 0,NULL, NULL);
  
  if  (fname && check_file_exists (fname) && is_lib(fname))return read_constraint_list_file (R, fname);
  else if (fname && check_file_exists (fname) && is_lib(fname)) list=read_lib_list ( fname, &n);
  else
    {
      X_template *F;
      list=(char**)vcalloc (S->nseq, sizeof (char*));
      for ( a=0; a<S->nseq; a++)
	{
	  if ((F=seq_has_template (S, a, "_F_")))
	    {
	      list[n++]=F->template_file;
	    }
	}
    }
  
  for (a=0; a< n; a++)
    {
      if (list[a])R=read_constraint_list_file (R, list[a]);
    }
  if (fname && fname[0]) save_contact_constraint_list (R, fname);
  return R;
}

Constraint_list * rna_lib_extension ( Constraint_list *CL, Constraint_list *R)
{
	CLIST_TYPE  *entry=NULL;
	int b, c, d, e, n1, n2, s1, r1, s2, r2, w2;
	int *list1;
	int *list2;
	Sequence *S=CL->S;
	static char *tmp;
	FILE *fp;

	list1=(int*)vcalloc ( 100, sizeof (int));
	list2=(int*)vcalloc ( 100, sizeof (int));

	if (!tmp) tmp=vtmpnam (NULL);
	fp=vfopen (tmp, "wb");
	entry=(int*)vcalloc ( 100, sizeof (int));
	for (s1=0; s1<S->nseq; s1++)
	{
		for (r1=1; r1<=S->len[s1]; r1++)
		{
			entry[SEQ1]=s1;
			for (c=1; c<CL->residue_index[s1][r1][0]; c+=ICHUNK)
			{
				s2=CL->residue_index[s1][r1][c+SEQ2];
				r2=CL->residue_index[s1][r1][c+R2];
				w2=CL->residue_index[s1][r1][c+WE];
				entry[SEQ2]=s2;
				entry[WE]=w2;

				n1=n2=0;
				list1[n1++]=r1;
				for (b=1; b<R->residue_index[s1][r1][0]; b+=ICHUNK)
				{
					list1[n1++]=R->residue_index[s1][r1][b+R2];
				}
				list2[n2++]=r2;

				for (b=1; b<R->residue_index[s2][r2][0]; b+=ICHUNK)
				{
					list2[n2++]=R->residue_index[s2][r2][b+R2];
				}

				for (b=1; b<n1; b++)
					for (d=1; d<n2; d++)
					{
						entry[R1]=list1[b];
						entry[R2]=list2[d];
						fwrite(entry, sizeof (int),CL->entry_len,fp);
						//for (e=0; e<CL->entry_len; e++) fprintf (fp, "%d ", entry[e]);
					}
			}
		}
	}

	vfclose (fp);
	return undump_constraint_list (CL,tmp);
}

char *** produce_method_file ( char *method)
{
	static char ***list;
	int n=0, a;
	FILE *fp;

	
	
	if (!list)

	  {
	    list=(char***)declare_arrayN(3, sizeof (char),1000,2,vtmpnam_size()+1);
	  }

	/*
	sprintf (list[n][0], "t_coffee");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE t_coffee\n");
	fprintf ( fp, "ADDRESS    built_in");
	fprintf ( fp, "ALN_MODE   any\n");
	fprintf ( fp, "OUT_MODE   L\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);

	vfclose (fp);}
	*/
	/*Space holder method to analyze very large dATASETS*/
	sprintf (list[n][0], "test_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE test_pair\n");
	fprintf ( fp, "DOC Fast alignmnents on the best diagonals\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "fast_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE fast_pair\n");
	fprintf ( fp, "DOC Fast alignmnents on the best diagonals\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if (strm (retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX dna_idmat\n");
		fprintf ( fp, "GOP -50\n");
		fprintf ( fp, "GEP -1\n");
	}
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "exon3_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE exon3_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   G\n");
	fprintf ( fp, "MATRIX     exon2mt\n");
	fprintf ( fp, "GOP        0\n");
	fprintf ( fp, "GEP        -1\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "exon2_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE exon2_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   G\n");
	fprintf ( fp, "MATRIX     exon2mt\n");
	fprintf ( fp, "GOP        -10\n");
	fprintf ( fp, "GEP        -1\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "exon_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE exon_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   G\n");
	fprintf ( fp, "MATRIX     exon2mt\n");
	fprintf ( fp, "GOP        -10\n");
	fprintf ( fp, "GEP        -1\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}


	sprintf (list[n][0], "blastr_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE slow_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "EXTEND_SEQ 1\n");
	fprintf ( fp, "MATRIX     blosumR\n");
	fprintf ( fp, "GOP        -20\n");
	fprintf ( fp, "GEP        0\n");
	fprintf ( fp, "WEIGHT     extend_sim\n");

	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "promo_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE slow_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "EXTEND_SEQ 1\n");
	fprintf ( fp, "MATRIX     promoter_tf1\n");
	fprintf ( fp, "GOP        -60\n");
	fprintf ( fp, "GEP        -1\n");
	fprintf ( fp, "WEIGHT     extend_sim\n");

	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "clean_slow_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC regular dynamic Programming\n");
	fprintf ( fp, "EXECUTABLE slow_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "WEIGHT     clean\n");
	if ( strm ( retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX   dna_idmat\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -1\n");
	}
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "slow_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC regular dynamic Programming\n");
	fprintf ( fp, "EXECUTABLE slow_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if ( strm ( retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX   dna_idmat\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -1\n");
	}
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	//hash_pair START
	sprintf (list[n][0], "hash_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC regular dynamic Programming\n");
	fprintf ( fp, "EXECUTABLE hash_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if ( strm ( retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX   dna_idmat\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -1\n");
	}
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}
	//hash_pair nsihed





	sprintf (list[n][0], "biphasic_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC bi-phasic dynamic Programming\n");
	fprintf ( fp, "EXECUTABLE biphasic_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if ( strm ( retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX   dna_idmat\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -1\n");
	}
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}
	


	sprintf (list[n][0], "proba_prfpair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE profile_pair\n");
	fprintf ( fp, "EXECUTABLE2 t_coffee\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -profile1=\n");
	fprintf ( fp, "IN_FLAG2   -profile2=\n");
	fprintf ( fp, "PARAM      -method=proba_pair\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   R\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	fprintf ( fp, "SUPPORTED  NO");
	vfclose (fp);}
	



	sprintf (list[n][0], "proba_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Probabilistic pairwise alignment\n");
	fprintf ( fp, "EXECUTABLE proba_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if ( strm ( retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "GOP   %d\n",CODE4DNA);//code for DNA
	}
	else
	{
		fprintf ( fp, "GOP   %d\n",CODE4PROTEINS);//Code for Proteins
	}
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}
	

	sprintf (list[n][0], "best_pair4prot");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Combination of the best template methods\n");
	fprintf ( fp, "EXECUTABLE best_pair4prot\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	if (strm (retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		printf_exit (EXIT_FAILURE, stderr, "\nERROR: The mode best_pair4prot is only suited for Proteins [FATAL:%s]\n", PROGRAM);
	}
	vfclose (fp);}


	sprintf (list[n][0], "best_pair4rna");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Combination of the best template methods\n");
	fprintf ( fp, "EXECUTABLE best_pair4rna\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	// 	  if (strm(retrieve_seq_type(), "RNA"))
	// 	  {
		// 		  printf_exit (EXIT_FAILURE, stderr, "\nERROR: The mode best_pair4rna is only suited for RNA [FATAL:%s]\n", PROGRAM);
		// 	  }
	vfclose (fp);}


	//Llaign ID PAIR
	sprintf (list[n][0], "lalign_id_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC local alignment reporting the N best pairwise local alignments\n");
	fprintf ( fp, "EXECUTABLE lalign_id_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	if ( strm (retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX     dna_idmat\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -1\n");
	}
	else
	{
		fprintf ( fp, "MATRIX     blosum50mt\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -4\n");
	}
	fprintf ( fp, "MAXID  100\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}



	//Llaign RS_S PAIR for Mocca
	sprintf (list[n][0], "lalign_rs_s_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC local alignment reporting the N best pairwise local alignments\n");
	fprintf ( fp, "EXECUTABLE lalign_id_pair\n");
	fprintf ( fp, "ALN_MODE   s_pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	if ( strm (retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))
	{
		fprintf ( fp, "MATRIX     dna_idmat\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -1\n");
	}
	else
	{
		fprintf ( fp, "MATRIX     blosum50mt\n");
		fprintf ( fp, "GOP   -10\n");
		fprintf ( fp, "GEP   -4\n");
	}
	fprintf ( fp, "MAXID  100\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}







	sprintf (list[n][0], "align_pdbpair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE align_pdb_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "lalign_pdbpair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE lalign_pdb_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "ktup_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE ktup_msa\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   s\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "seq_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE seq_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   s\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "extern_pdbpair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE tc_P_generic_method.pl\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   A\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "externprofile_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE tc_R_generic_method\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   A\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   R\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "thread_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	/*This runs thread_pair@EP@EXECUTABLE2@threpg*/
	fprintf ( fp, "EXECUTABLE thread_pair\n");
	fprintf ( fp, "EXECUTABLE2 t_coffee\n");
	fprintf ( fp, "ALN_MODE    pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile1=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   Ps\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "fugue_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	/*This runs thread_pair@EP@EXECUTABLE2@threpg*/
	fprintf ( fp, "EXECUTABLE thread_pair\n");
	fprintf ( fp, "EXECUTABLE2 fugueali\n");
	fprintf ( fp, "ALN_MODE    pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile1=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   Ps\n");
	fprintf ( fp, "ADDRESS    %s\n", FUGUE_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", FUGUE_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "pdb_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE pdb_pair\n");
	fprintf ( fp, "EXECUTABLE2 t_coffee\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}



	sprintf (list[n][0], "hh_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE profile_pair\n");
	fprintf ( fp, "EXECUTABLE2 hhalign\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -profile1=\n");
	fprintf ( fp, "IN_FLAG2   -profile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   R\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	fprintf ( fp, "SUPPORTED  NO");
	vfclose (fp);}

	sprintf (list[n][0], "co_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE profile_pair\n");
	fprintf ( fp, "EXECUTABLE2 clustalo\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -profile1=\n");
	fprintf ( fp, "IN_FLAG2   -profile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   R\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	fprintf ( fp, "SUPPORTED  NO");
	vfclose (fp);}

	sprintf (list[n][0], "cwprofile_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE profile_pair\n");
	fprintf ( fp, "EXECUTABLE2 clustalw2\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -profile1=\n");
	fprintf ( fp, "IN_FLAG2   -profile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   R\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	fprintf ( fp, "SUPPORTED  NO");
	vfclose (fp);}


	//Switch to TM_align if SAP is not installed
	//Intercept sap
	if (method && strm (method, "sap_pair") && !check_program_is_installed (SAP_4_TCOFFEE,NULL,NULL,SAP_ADDRESS,INSTALL))
	{
		static int issued;
		if (!issued)
		{
		  add_warning (stderr, "SAP is not installed: TMalign will be used instead");
			issued=1;
		}

		sprintf (list[n][0], "sap_pair");
		sprintf (list[n][1], "%s", vtmpnam(NULL));
		n++;

		fp=vfopen (list[n-1][1], "w");
		fprintf ( fp, "DOC: TM-Align: pairwise structural aligner [%s]\n", TMALIGN_ADDRESS);
		fprintf ( fp, "EXECUTABLE pdb_pair\n");
		fprintf ( fp, "EXECUTABLE2 TMalign\n" );
		fprintf ( fp, "ALN_MODE   pairwise\n");
		fprintf ( fp, "OUT_MODE   fL\n");
		fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
		fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
		fprintf ( fp, "OUT_FLAG   -outfile=\n");
		fprintf ( fp, "ADDRESS    %s\n", TMALIGN_ADDRESS);
		fprintf ( fp, "PROGRAM    %s\n", TMALIGN_4_TCOFFEE);
		fprintf ( fp, "SEQ_TYPE   P\n");
		vfclose (fp);
	}
	else
	{

		sprintf (list[n][0], "sap_pair");
		sprintf (list[n][1], "%s", vtmpnam(NULL));
		n++;
		if (method==NULL || strm (method, list[n-1][0]))
		{
			fp=vfopen (list[n-1][1], "w");
			fprintf ( fp, "DOC: sap: pairwise structural aligner [%s]\n", SAP_ADDRESS);
			fprintf ( fp, "EXECUTABLE sap_pair\n");
			fprintf ( fp, "ALN_MODE   pairwise\n");
			fprintf ( fp, "OUT_MODE   fL\n");
			fprintf ( fp, "IN_FLAG    no_name\n");
			fprintf ( fp, "OUT_FLAG   no_name\n");
			fprintf ( fp, "WEIGHT     100\n");
			fprintf ( fp, "SEQ_TYPE   P\n");
			fprintf ( fp, "ADDRESS    %s\n", SAP_ADDRESS);
			fprintf ( fp, "PROGRAM    %s\n", SAP_4_TCOFFEE);
			vfclose (fp);
		}
	}
	sprintf (list[n][0], "sara_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;
	if (method==NULL || strm (method, list[n-1][0])){
		fp=vfopen (list[n-1][1], "w");
		fprintf ( fp, "DOC: SARA: pairwise structural RNA aligner [%s]\n", ADDRESS_BUILT_IN);
		fprintf ( fp, "EXECUTABLE rnapdb_pair\n");
		fprintf ( fp, "EXECUTABLE2 runsara.py\n" );
		fprintf ( fp, "ALN_MODE   pairwise\n");
		fprintf ( fp, "OUT_MODE   fL\n");
		fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
		fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
		fprintf ( fp, "OUT_FLAG   -outfile=\n");
		fprintf ( fp, "SEQ_TYPE   P\n");
		fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
		fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
		vfclose (fp);
	}

	sprintf (list[n][0], "daliweb_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: Dalilite: pairwise structural aligner [%s]\n", DALILITEc_ADDRESS);
	fprintf ( fp, "EXECUTABLE pdbid_pair\n");
	fprintf ( fp, "EXECUTABLE2 daliweb\n" );
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", DALILITEc_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DALILITEc_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "dali_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: Dalilite: pairwise structural aligner [%s]\n", DALILITEc_ADDRESS);
	fprintf ( fp, "EXECUTABLE pdb_pair\n");
	fprintf ( fp, "EXECUTABLE2 DaliLite\n" );
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", DALILITEc_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DALILITEc_4_TCOFFEE);
	vfclose (fp);}


	sprintf (list[n][0], "mustang_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: Mustang: pairwise structural aligner [%s]\n", MUSTANG_ADDRESS);
	fprintf ( fp, "EXECUTABLE pdb_pair\n");
	fprintf ( fp, "EXECUTABLE2 mustang\n" );
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   P\n");
	fprintf ( fp, "ADDRESS    %s\n", MUSTANG_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MUSTANG_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "TMalign_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: TM-Align: pairwise structural aligner [%s]\n", TMALIGN_ADDRESS);
	fprintf ( fp, "EXECUTABLE pdb_pair\n");
	fprintf ( fp, "EXECUTABLE2 TMalign\n" );
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -pdbfile1=\n");
	fprintf ( fp, "IN_FLAG2   -pdbfile2=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "ADDRESS    %s\n", TMALIGN_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", TMALIGN_4_TCOFFEE);
	fprintf ( fp, "SEQ_TYPE   P\n");
	vfclose (fp);}

	sprintf (list[n][0], "cdna_fast_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE cdna_fast_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "cdna_cfast_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE cdna_cfast_pair\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    no_name\n");
	fprintf ( fp, "OUT_FLAG   no_name\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
	fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
	vfclose (fp);}

	sprintf (list[n][0], "blastp_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: BLAST multiple Aligner [%s]\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "EXECUTABLE seq_msa\n");
	fprintf ( fp, "EXECUTABLE2 blastp\n" );
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", NCBIBLAST_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "blastp_o2a");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: BLAST multiple Aligner [%s]\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "EXECUTABLE seq_msa\n");
	fprintf ( fp, "EXECUTABLE2 blastp\n" );
	fprintf ( fp, "ALN_MODE   o2a\n");
	fprintf ( fp, "OUT_MODE   fL\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", NCBIBLAST_4_TCOFFEE);
	vfclose (fp);}

	/*deprecated: Should now run clustalo via dynamic.pl	*/
    sprintf (list[n][0], "old_clustalo_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: clustalo [%s]\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "EXECUTABLE clustalo\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -i&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -o&bnsp\n");
    fprintf ( fp, "PARAM      --force --threads=1 \n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", CLUSTALO_4_TCOFFEE);
    vfclose (fp);}

    /*deprecated: Should now run clustalo via dynamic.pl	*/
    sprintf (list[n][0], "old_clustalo_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: clustalo [%s]\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "EXECUTABLE clustalo\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -i&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -o&bnsp\n");
    fprintf ( fp, "PARAM      --force --threads=1 &bnsp2>/dev/null\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", CLUSTALO_4_TCOFFEE);
    vfclose (fp);}

    // The following methods are T-Coffee modes made available via thscript dynamic.pl
    // This script can select a method based on the number of sequences -dynamic_config
    // It is also the easiest way to use a T-Coffee mode as child aligner when running the regressive algo
    // The script also supports all the standard T-Coffee modes
    sprintf (list[n][0], "dynamic_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Will use a method depending on the provides number of sequences");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}
    
    sprintf (list[n][0], "3dcoffee_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs 3d coffee from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method 3dcoffee_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}
    
    sprintf (list[n][0], "expresso_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs expresso from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method expresso_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}

    sprintf (list[n][0], "accurate_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Runs t_coffee -mode accurate from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method accurate_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}
    
    sprintf (list[n][0], "psicoffee_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs psicoffee from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method psicoffee_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}


    //FAMSA runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "famsa_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs famsa from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method famsa_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", FAMSA_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", FAMSA_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "famsa_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs fmasa from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method famsa_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", FAMSA_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", FAMSA_4_TCOFFEE);
    vfclose (fp);}
    
    
    
    //CLUSTALO runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "clustalo_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs clustalo from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method clustalo_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", CLUSTALO_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "clustalo_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs clustalo from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method clustalo_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", CLUSTALO_4_TCOFFEE);
    vfclose (fp);}

    

    
    //MAFFT runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "mafft_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafft_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "mafft_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafft_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}

    
    //MAFFT runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "mafftginsi_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftginsi_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "mafftginsi_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftginsi_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}


    //MAFFT runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "mafftfftns1_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftfftns1_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "mafftfftns1_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftfftns1_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}

    //MAFFT runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "mafftfftnsi_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftfftnsi_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "mafftfftnsi_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftfftnsi_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}

     //MAFFT runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "mafftnwnsi_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftnwnsi_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "mafftnwnsi_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftfftnwnsi_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}


     //MAFFT runs through dynamic.pl
    //Guide Tree mode can be selected
    sprintf (list[n][0], "mafftsparsecore_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftsparsecore_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}
   
    sprintf (list[n][0], "mafftnwnsi_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pair\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method mafftffsparsecore_pair\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
    vfclose (fp);}



    sprintf (list[n][0], "mafftsparsecore_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE mafft-sparsecore.rb\n");
	fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -i\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM1 -C \"--globalpair --maxiterate 100 --anysymbol\" \n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "mafftsparsecore_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE mafft-sparsecore.rb\n");
	fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -i\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM1 -C \"--globalpair --maxiterate 100 --anysymbol\" \n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
	vfclose (fp);}

	 //MAFFT runs through dynamic.pl
	//Guide Tree mode can be selected
	sprintf (list[n][0], "mafftlinsi_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	  fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
	  fprintf ( fp, "EXECUTABLE dynamic.pl\n");
	  fprintf ( fp, "ALN_MODE   multiple\n");
	  fprintf ( fp, "OUT_MODE   aln\n");
	  fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
	  fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
	  fprintf ( fp, "PARAM      -method mafftlinsi_msa\n");
	  fprintf ( fp, "SEQ_TYPE   S\n");
	  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
	  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
	  vfclose (fp);}
	
	sprintf (list[n][0], "mafftlinsi_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	  fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
	  fprintf ( fp, "EXECUTABLE dynamic.pl\n");
	  fprintf ( fp, "ALN_MODE   pair\n");
	  fprintf ( fp, "OUT_MODE   aln\n");
	  fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
	  fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
	  fprintf ( fp, "PARAM      -method mafftlinsi_pair\n");
	  fprintf ( fp, "SEQ_TYPE   S\n");
	  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
	  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
	  vfclose (fp);}

	//MAFFT runs through dynamic.pl
	//Guide Tree mode can be selected
	sprintf (list[n][0], "maffteinsi_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	  fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
	  fprintf ( fp, "EXECUTABLE dynamic.pl\n");
	  fprintf ( fp, "ALN_MODE   multiple\n");
	  fprintf ( fp, "OUT_MODE   aln\n");
	  fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
	  fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
	  fprintf ( fp, "PARAM      -method maffteinsi_msa\n");
	  fprintf ( fp, "SEQ_TYPE   S\n");
	  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
	  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
	  vfclose (fp);}

	
	//MAFFT runs through dynamic.pl
	//Guide Tree mode can be selected
	sprintf (list[n][0], "maffteinsi_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	  fprintf ( fp, "DOC: dynamic [%s]\n", "runs mafft from dynamic.pl");
	  fprintf ( fp, "EXECUTABLE dynamic.pl\n");
	  fprintf ( fp, "ALN_MODE   multiple\n");
	  fprintf ( fp, "OUT_MODE   pair\n");
	  fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
	  fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
	  fprintf ( fp, "PARAM      -method maffteinsi_pair\n");
	  fprintf ( fp, "SEQ_TYPE   S\n");
	  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
	  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
	  vfclose (fp);}
	
   



    sprintf (list[n][0], "dynamic_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Will use a method depending on the provides number of sequences");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}
    
    sprintf (list[n][0], "3dcoffee_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Will use a method depending on the provides number of sequences");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method 3dcoffee_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}
    
    sprintf (list[n][0], "expresso_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Will use a method depending on the provides number of sequences");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method expresso_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}

    sprintf (list[n][0], "accurate_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Will use a method depending on the provides number of sequences");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method accurate_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}
    
    sprintf (list[n][0], "psicoffee_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: dynamic [%s]\n", "Will use a method depending on the provides number of sequences");
    fprintf ( fp, "EXECUTABLE dynamic.pl\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    -seq&bnsp\n");
    fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
    fprintf ( fp, "PARAM      -method psicoffee_msa\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
    fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
    vfclose (fp);}






    sprintf (list[n][0], "clustaloNF_pair");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: clustalo.pl [%s]\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "EXECUTABLE clustalo.pl\n");
    fprintf ( fp, "ALN_MODE   pairwise\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "PARAM1 all \n");
    fprintf ( fp, "IN_FLAG    &bnsp\n");
    fprintf ( fp, "OUT_FLAG   >\n");
    fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", CLUSTALO_4_TCOFFEE);
    vfclose (fp);}

    sprintf (list[n][0], "clustaloNF_msa");
    sprintf (list[n][1], "%s", vtmpnam(NULL));
    n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "DOC: clustalo.pl [%s]\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "EXECUTABLE clustalo.pl\n");
    fprintf ( fp, "ALN_MODE   multiple\n");
    fprintf ( fp, "OUT_MODE   aln\n");
    fprintf ( fp, "IN_FLAG    &bnsp\n");
    fprintf ( fp, "OUT_FLAG   \n");
    //fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
    fprintf ( fp, "PARAM1 all \n");
    fprintf ( fp, "PARAM      &bnsp2");
    fprintf ( fp, "PARAM      \n");
    fprintf ( fp, "SEQ_TYPE   S\n");
    fprintf ( fp, "ADDRESS    %s\n", CLUSTALO_ADDRESS);
    fprintf ( fp, "PROGRAM    %s\n", CLUSTALO_4_TCOFFEE);
    vfclose (fp);}


	sprintf (list[n][0], "clustalw2_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: clustalw [%s]\n", CLUSTALW2_ADDRESS);
	fprintf ( fp, "EXECUTABLE clustalw2\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    %sINFILE=\n",CWF);
	fprintf ( fp, "OUT_FLAG   %sOUTFILE=\n",CWF);
	fprintf ( fp, "PARAM      %sOUTORDER=INPUT %sNEWTREE=SCRATCH_FILE %salign\n",CWF,CWF,CWF);
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", CLUSTALW2_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", CLUSTALW2_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "clustalw2_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC clustalw[%s]\n", CLUSTALW2_ADDRESS);
	fprintf ( fp, "EXECUTABLE clustalw2\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    %sINFILE=\n",CWF);
	fprintf ( fp, "OUT_FLAG   %sOUTFILE=\n", CWF);
	fprintf ( fp, "PARAM      %sOUTORDER=INPUT %sNEWTREE=SCRATCH_FILE %salign\n",CWF,CWF,CWF);
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", CLUSTALW2_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", CLUSTALW2_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "clustalw_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC: clustalw [%s]\n", CLUSTALW_ADDRESS);
	fprintf ( fp, "EXECUTABLE clustalw2\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    %sINFILE=\n", CWF);
	fprintf ( fp, "OUT_FLAG   %sOUTFILE=\n",CWF);
	fprintf ( fp, "PARAM      %sOUTORDER=INPUT %sNEWTREE=SCRATCH_FILE %salign\n",CWF,CWF,CWF);
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", CLUSTALW_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", CLUSTALW_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "clustalw_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC clustalw[%s]\n", CLUSTALW_ADDRESS);
	fprintf ( fp, "EXECUTABLE clustalw2\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    %sINFILE=\n", CWF);
	fprintf ( fp, "OUT_FLAG   %sOUTFILE=\n", CWF);
	fprintf ( fp, "PARAM      %sOUTORDER=INPUT %sNEWTREE=SCRATCH_FILE %salign  >/dev/null&bnsp2>/dev/null\n\n",CWF,CWF,CWF);
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", CLUSTALW_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", CLUSTALW_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "uppNF_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE upp.pl\n");
	fprintf ( fp, "DOC MSA [%s]\n", MSA_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "PARAM1 all \n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", UPP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", UPP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "uppNF_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE upp.pl\n");
	fprintf ( fp, "DOC UPP [%s]\n", UPP_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   \n");
	//fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "PARAM1 all \n");
	fprintf ( fp, "PARAM      &bnsp2");
	fprintf ( fp, "PARAM      \n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", UPP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", UPP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "upp_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE upp.pl\n");
	fprintf ( fp, "DOC MSA [%s]\n", UPP_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "PARAM1 one \n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", UPP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", UPP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "upp_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE upp.pl\n");
	fprintf ( fp, "DOC MSA [%s]\n", UPP_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "PARAM1 one \n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   \n");
	//fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "PARAM      &bnsp2");
	fprintf ( fp, "PARAM      \n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", UPP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", UPP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "msa_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE msa.pl\n");
	fprintf ( fp, "DOC MSA [%s]\n", MSA_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MSA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MSA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "msa_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE msa.pl\n");
	fprintf ( fp, "DOC MSA [%s]\n", MSA_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   \n");
	//fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "PARAM      &bnsp2");
	fprintf ( fp, "PARAM      \n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MSA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MSA_4_TCOFFEE);
	vfclose (fp);}

	
		sprintf (list[n][0], "dca_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE dca.pl\n");
	fprintf ( fp, "DOC DCA [%s]\n", DCA_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", DCA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DCA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "dca_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE dca.pl\n");
	fprintf ( fp, "DOC DCA [%s]\n", DCA_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   \n");
	//fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "PARAM      &bnsp2");
	fprintf ( fp, "PARAM      \n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", DCA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DCA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "dialigntx_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC dialign-tx [%s]\n", DIALIGNTX_ADDRESS);
	fprintf ( fp, "EXECUTABLE dialign-tx\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	if ( isdir (DIALIGNTX_DIR))
		fprintf ( fp, "PARAM1 %s \n", DIALIGNTX_DIR);
	else
		fprintf ( fp, "PARAM1 %s \n", get_mcoffee_4_tcoffee());
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   &bnsp\n");
	fprintf ( fp, "PARAM      >/dev/null&bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", DIALIGNTX_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DIALIGNTX_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "dialigntx_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC dialign-tx [%s]\n", DIALIGNTX_ADDRESS);
	fprintf ( fp, "EXECUTABLE dialign-tx\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	if ( isdir (DIALIGNTX_DIR))
		fprintf ( fp, "PARAM1 %s \n", DIALIGNTX_DIR);
	else
		fprintf ( fp, "PARAM1 %s \n", get_mcoffee_4_tcoffee());
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   &bnsp\n");
	fprintf ( fp, "PARAM      >/dev/null&bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", DIALIGNTX_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DIALIGNTX_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "dialignt_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC dialign-tx [%s]\n", DIALIGNT_ADDRESS);
	fprintf ( fp, "EXECUTABLE dialign-tx\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	if ( isdir (DIALIGNT_DIR))
		fprintf ( fp, "PARAM1 %s \n", DIALIGNT_DIR);
	else
		fprintf ( fp, "PARAM1 %s \n", get_mcoffee_4_tcoffee());
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   &bnsp\n");
	fprintf ( fp, "PARAM      >/dev/null&bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", DIALIGNT_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DIALIGNT_4_TCOFFEE);

	vfclose (fp);}

	sprintf (list[n][0], "dialignt_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC dialign-tx [%s]\n", DIALIGNT_ADDRESS);
	fprintf ( fp, "EXECUTABLE dialign-tx\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	if ( isdir (DIALIGNT_DIR))
		fprintf ( fp, "PARAM1 %s \n", DIALIGNT_DIR);
	else
		fprintf ( fp, "PARAM1 %s \n", get_mcoffee_4_tcoffee());
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   &bnsp\n");
	fprintf ( fp, "PARAM      >/dev/null&bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", DIALIGNT_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", DIALIGNT_4_TCOFFEE);

	vfclose (fp);}

	sprintf (list[n][0], "poa_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Partial Order Graph Alignment [%s]\n", POA_ADDRESS);
	fprintf ( fp, "EXECUTABLE poa\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "PARAM1 -toupper \n");
	fprintf ( fp, "IN_FLAG    -read_fasta&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -clustal&bnsp\n");
	if (file_exists (POA_DIR, POA_FILE1))
		fprintf ( fp, "PARAM      %s/%s&bnsp2>/dev/null\n",POA_DIR,POA_FILE1);
	else
		fprintf ( fp, "PARAM      %s/%s&bnsp2>/dev/null\n", get_mcoffee_4_tcoffee(), POA_FILE1);
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", POA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",POA_4_TCOFFEE);

	vfclose (fp);}

	sprintf (list[n][0], "poa_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Partial Order Graph Alignment [%s]\n", POA_ADDRESS);
	fprintf ( fp, "EXECUTABLE poa\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "PARAM1 -toupper \n");
	fprintf ( fp, "IN_FLAG    -read_fasta&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -clustal&bnsp\n");
	if (file_exists (POA_DIR, POA_FILE1))
		fprintf ( fp, "PARAM      %s/%s&bnsp2>/dev/null\n",POA_DIR,POA_FILE1);
	else
		fprintf ( fp, "PARAM      %s/%s&bnsp2>/dev/null\n", get_mcoffee_4_tcoffee(), POA_FILE1);
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", POA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",POA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "msaprobs_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC probcons [%s]\n", MSAPROBS_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "EXECUTABLE msaprobs\n");
	fprintf ( fp, "ADDRESS    %s\n",MSAPROBS_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",MSAPROBS_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "msaprobs_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC probcons [%s]\n", MSAPROBS_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "EXECUTABLE msaprobs\n");
	fprintf ( fp, "ADDRESS    %s\n",MSAPROBS_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",MSAPROBS_4_TCOFFEE);
	vfclose (fp);}
	
	sprintf (list[n][0], "probcons_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC probcons [%s]\n", PROBCONS_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if ( strm (retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))fprintf ( fp, "EXECUTABLE probconsRNA\n");
	else fprintf ( fp, "EXECUTABLE probcons\n");
	fprintf ( fp, "ADDRESS    %s\n", PROBCONS_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",PROBCONS_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "probcons_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC probcons [%s]\n", PROBCONS_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	if ( strm (retrieve_seq_type(), "DNA") || strm (retrieve_seq_type(), "RNA"))fprintf ( fp, "EXECUTABLE probconsRNA\n");
	else fprintf ( fp, "EXECUTABLE probcons\n");
	fprintf ( fp, "ADDRESS    %s\n", PROBCONS_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",PROBCONS_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "probconsRNA_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC probcons [%s]\n", PROBCONSRNA_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "EXECUTABLE probconsRNA\n");
	fprintf ( fp, "ADDRESS    %s\n", PROBCONSRNA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",PROBCONSRNA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "probconsRNA_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC probcons [%s]\n", PROBCONSRNA_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "EXECUTABLE probconsRNA\n");
	fprintf ( fp, "ADDRESS    %s\n", PROBCONSRNA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n",PROBCONSRNA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "muscle_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Muscle [%s]\n", MUSCLE_ADDRESS);
	fprintf ( fp, "EXECUTABLE muscle\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -in&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -out&bnsp\n");
	fprintf ( fp, "PARAM      -quiet&bnsp-maxmb&bnsp0&bnsp>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MUSCLE_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MUSCLE_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "muscle_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Muscle [%s]\n", MUSCLE_ADDRESS);
	fprintf ( fp, "EXECUTABLE muscle\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -in&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -out&bnsp\n");
	fprintf ( fp, "PARAM      -quiet&bnsp-maxmb&bnsp0&bnsp>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MUSCLE_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MUSCLE_4_TCOFFEE);
	vfclose (fp);}



	sprintf (list[n][0], "mus4_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Muscle [%s]\n", MUS4_ADDRESS);
	fprintf ( fp, "EXECUTABLE mus4\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    --input&bnsp\n");
	fprintf ( fp, "OUT_FLAG   --output&bnsp\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MUS4_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MUS4_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "mus4_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC Mus4 [%s]\n", MUS4_ADDRESS);
	fprintf ( fp, "EXECUTABLE mus4\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    --input&bnsp\n");
	fprintf ( fp, "OUT_FLAG   --output&bnsp\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", MUS4_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", MUS4_4_TCOFFEE);
	vfclose (fp);}



	sprintf (list[n][0], "t_coffee_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE t_coffee\n");
	fprintf ( fp, "DOC T-Coffee [%s]\n", TCOFFEE_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -infile&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
	fprintf ( fp, "PARAM      -n_core=1\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", TCOFFEE_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", TCOFFEE_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "t_coffee_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE t_coffee\n");
	fprintf ( fp, "DOC T-Coffee [%s]\n", TCOFFEE_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -infile&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -outfile&bnsp\n");
	fprintf ( fp, "PARAM      -n_core=1\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", TCOFFEE_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", TCOFFEE_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "pcma_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC PCMA [%s]\n", PCMA_ADDRESS);
	fprintf ( fp, "EXECUTABLE pcma\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", PCMA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", PCMA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "pcma_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC PCMA [%s]\n", PCMA_ADDRESS);
	fprintf ( fp, "EXECUTABLE pcma\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", PCMA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", PCMA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "kalign_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE kalign\n");
	fprintf ( fp, "DOC kalign [%s]\n", KALIGN_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -i&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -o&bnsp\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", KALIGN_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", KALIGN_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "kalign_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE kalign\n");
	fprintf ( fp, "DOC kalign [%s]\n", KALIGN_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -i&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -o&bnsp\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", KALIGN_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", KALIGN_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "amap_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE amap\n");
	fprintf ( fp, "DOC amap [%s]\n", AMAP_ADDRESS);
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", AMAP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", AMAP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "amap_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE amap\n");
	fprintf ( fp, "DOC amap [%s]\n", AMAP_ADDRESS);
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", AMAP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", AMAP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "proda_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC proda [%s]\n", PRODA_ADDRESS);
	fprintf ( fp, "EXECUTABLE proda\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", PRODA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", PRODA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "proda_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC proda [%s]\n", PRODA_ADDRESS);
	fprintf ( fp, "EXECUTABLE proda\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", PRODA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", PRODA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "prank_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC prank [%s]\n", PRANK_ADDRESS);
	fprintf ( fp, "EXECUTABLE tc_generic_method.pl\n");
	fprintf ( fp, "ALN_MODE  pairwise\n");
	fprintf ( fp, "OUT_MODE  aln\n");
	fprintf ( fp, "PARAM -method=%s -mode=seq_msa -tmpdir=%s\n",(getenv("PRANK_4_TCOFFEE"))?getenv("PRANK_4_TCOFFEE"):PRANK_4_TCOFFEE, get_tmp_4_tcoffee());
	fprintf ( fp, "IN_FLAG -infile=\n");
	fprintf ( fp, "OUT_FLAG -outfile=\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", PRANK_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", PRANK_4_TCOFFEE);
	vfclose (fp);}



	sprintf (list[n][0], "fsa_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC fsa [%s]\n", FSA_ADDRESS);
	fprintf ( fp, "EXECUTABLE fsa\n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", FSA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", FSA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "fsa_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC fsa [%s]\n", FSA_ADDRESS);
	fprintf ( fp, "EXECUTABLE fsa\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    &bnsp\n");
	fprintf ( fp, "OUT_FLAG   >\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", FSA_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", FSA_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "tblastx_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC tblastx [%s]\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "EXECUTABLE tc_generic_method.pl\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "PARAM  -mode=tblastx_msa\n");
	fprintf ( fp, "OUT_MODE  L\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", NCBIBLAST_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "tblastpx_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC tblastpx [%s]\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "EXECUTABLE tc_generic_method.pl\n");
	fprintf ( fp, "ALN_MODE   multiple\n");
	fprintf ( fp, "PARAM  -mode=tblastpx_msa\n");
	fprintf ( fp, "OUT_MODE  L\n");
	fprintf ( fp, "IN_FLAG    -infile=\n");
	fprintf ( fp, "OUT_FLAG   -outfile=\n");
	fprintf ( fp, "PARAM      &bnsp2>/dev/null\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", NCBIBLAST_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "plib_msa");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC t_coffee [%s]\n", "t_coffee");
	fprintf ( fp, "EXECUTABLE plib_msa\n");
	fprintf ( fp, "ALN_MODE  multiple\n");
	fprintf ( fp, "OUT_MODE  fL\n");
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", "www.tcoffee.org");
	fprintf ( fp, "PROGRAM    %s\n", "t_coffee");
	vfclose (fp);}


	sprintf (list[n][0], "em");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;
	if (method==NULL || lstrstr (method,"em@"))
	{

		fp=vfopen (list[n-1][1], "w");
		if (method)
		{
			char **l2;
			l2=string2list2 ( method, "@");
			fprintf ( fp, "PARAM -method=%s -mode=seq_msa -tmpdir=%s\n",l2[2], get_tmp_4_tcoffee());
			fprintf ( fp, "ALN_MODE  %s\n", l2[3]);
			free_char (l2, -1);
		}
		fprintf ( fp, "EXECUTABLE tc_generic_method.pl\n");
		fprintf ( fp, "OUT_MODE  aln\n");
		fprintf ( fp, "IN_FLAG -infile=\n");
		fprintf ( fp, "OUT_FLAG -outfile=\n");
		fprintf ( fp, "SEQ_TYPE   S\n");
		fprintf ( fp, "ADDRESS    %s\n", ADDRESS_BUILT_IN);
		fprintf ( fp, "PROGRAM    %s\n", PROGRAM_BUILT_IN);
		vfclose (fp);
	}

	sprintf (list[n][0], "consan_pair");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "DOC consan (sfold) RNA pairwise sequence aligner [%s]\n", CONSAN_ADDRESS);
	fprintf ( fp, "EXECUTABLE fasta_seq2consan_aln.pl \n");
	fprintf ( fp, "ALN_MODE   pairwise\n");
	fprintf ( fp, "OUT_MODE   aln\n");
	fprintf ( fp, "IN_FLAG    -i&bnsp\n");
	fprintf ( fp, "OUT_FLAG   -o&bnsp\n");
	fprintf ( fp, "PARAM      -d&bnsp%s&bnsp2>/dev/null\n",get_mcoffee_4_tcoffee());
	fprintf ( fp, "SEQ_TYPE   S\n");
	fprintf ( fp, "ADDRESS    %s\n", CONSAN_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", CONSAN_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "RNAplfold");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE RNAplfold \n");
	fprintf ( fp, "ALN_MODE   predict\n");
	fprintf ( fp, "SEQ_TYPE   RNA\n");
	fprintf ( fp, "ADDRESS    %s\n", RNAPLFOLD_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", RNAPLFOLD_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "HMMtop");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE HMMtop \n");
	fprintf ( fp, "ALN_MODE   predict\n");
	fprintf ( fp, "SEQ_TYPE   PROTEIN\n");
	fprintf ( fp, "ADDRESS    %s\n", HMMTOP_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", HMMTOP_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "GOR4");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE GORIV \n");
	fprintf ( fp, "ALN_MODE   predict\n");
	fprintf ( fp, "SEQ_TYPE   PROTEIN\n");
	fprintf ( fp, "ADDRESS    %s\n", GOR4_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", GOR4_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "wublast_client");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE wublast.pl \n");
	fprintf ( fp, "ALN_MODE   predict\n");
	fprintf ( fp, "SEQ_TYPE   PROTEIN\n");
	fprintf ( fp, "ADDRESS    %s\n", EBIWUBLASTc_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", EBIWUBLASTc_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "blastpgp_client");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE blastpgp.pl \n");
	fprintf ( fp, "ALN_MODE   predict\n");
	fprintf ( fp, "SEQ_TYPE   PROTEIN\n");

	fprintf ( fp, "ADDRESS    %s\n", EBIBLASTPGPc_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", EBIBLASTPGPc_4_TCOFFEE);
	vfclose (fp);}


	sprintf (list[n][0], "local_ncbiblast");
	sprintf (list[n][1], "%s", vtmpnam(NULL));
	n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
	fprintf ( fp, "EXECUTABLE blastall \n");
	fprintf ( fp, "ALN_MODE   predict\n");
	fprintf ( fp, "SEQ_TYPE   PROTEIN\n");

	fprintf ( fp, "ADDRESS    %s\n", NCBIBLAST_ADDRESS);
	fprintf ( fp, "PROGRAM    %s\n", NCBIBLAST_4_TCOFFEE);
	vfclose (fp);}

	sprintf (list[n][0], "famsa_msa"); 
	sprintf (list[n][1], "%s", vtmpnam(NULL)); 
	n++;
	if (method==NULL || strm (method, list[n-1][0])){
		fp=vfopen (list[n-1][1], "w"); 
		fprintf ( fp, "DOC famsa [%s]\n", "FAMSA"); 
		fprintf ( fp, "EXECUTABLE famsa\n"); 
		fprintf ( fp, "ALN_MODE multiple\n"); 
		fprintf ( fp, "OUT_MODE aln\n"); 
		fprintf ( fp, "IN_FLAG &bnsp\n"); 
		fprintf ( fp, "OUT_FLAG &bnsp\n"); 
		fprintf ( fp, "PARAM -t 1 &bnsp2>/dev/null\n"); 
		fprintf ( fp, "SEQ_TYPE S\n"); 

		fprintf ( fp, "ADDRESS %s\n", "famsa"); 
		fprintf ( fp, "PROGRAM %s\n", "famsa"); 
		vfclose (fp);
	}

 	sprintf (list[n][0], "famsa_pair"); 
	sprintf (list[n][1], "%s", vtmpnam(NULL)); n++;
	if (method==NULL || strm (method, list[n-1][0])){
		fp=vfopen (list[n-1][1], "w"); 
		fprintf ( fp, "DOC famsa [%s]\n", "FAMSA"); 
		fprintf ( fp, "EXECUTABLE famsa\n"); 
		fprintf ( fp, "ALN_MODE pairwise\n"); 
		fprintf ( fp, "OUT_MODE aln\n"); 
		fprintf ( fp, "IN_FLAG &bnsp\n"); 
		fprintf ( fp, "OUT_FLAG &bnsp\n"); 
		fprintf ( fp, "PARAM &bnsp2>/dev/null\n"); 
		fprintf ( fp, "SEQ_TYPE S\n"); 
		
		fprintf ( fp, "ADDRESS %s\n", "famsa"); 
		fprintf ( fp, "PROGRAM %s\n", "famsa"); 
		vfclose (fp);
	}
	
	//programatically add something to the configuraion file
	//Used to cause the creation fo a large number of temp files
	
	for (a=0; a<n; a++)
	  {
	    if (strm (method, list[a][0]))
	      printf_file (list[a][1], "a", "METHOD %s \n", list[a][0]); 
	  }
	
	list[n]=NULL;
	return list;
}



Constraint_list * plib_msa (Constraint_list *CL)
{
  int *list=(int*)vcalloc ((CL->S)->nseq,sizeof(int));
  int s,a,ns;
  int n;

  n=get_int_variable ("N_4_PLIB");
  n=(!n)?10:n;

  for (ns=s=a=0; a<n-1 && s!=-1; a++)
    {
      list[ns++]=s;
      CL=add_seq2cl (s,CL);
      s=cl2worst_seq (CL,list,ns);
      if (s!=-1)fprintf (stderr, "\n\tMaster Sequence: %s", (CL->S)->name[s]);
    }
  return CL;
}

int cl2worst_seq (Constraint_list *CL, int *list, int ns)
{
  Sequence *S=CL->S;
  int **score;
  int *bsp;
  int  bs, a, b, s1, s2,r1;
  int len_normalize=1;

  score=declare_int (S->nseq, 2);
  for (a=0; a<S->nseq; a++)score[a][0]=a;
  for (s1=0; s1<S->nseq; s1++)
    for (r1=1; r1<=S->len[s1]; r1++)
      for (b=1; b<CL->residue_index[s1][r1][0]; b+=ICHUNK)
	{
	  s2=CL->residue_index[s1][r1][b+SEQ2];
	  score[s1][1]+=CL->residue_index[s1][r1][b+WE];
	  score[s2][1]+=CL->residue_index[s1][r1][b+WE];
	}

  if (len_normalize)for (a=0; a<S->nseq; a++)score[a][1]/=strlen (S->seq[a]);

  //get rid of sequences already used
  for (a=0; a<ns; a++)
    {
      score[list[a]][1]=-1;
    }
  for (b=0,a=0; a<S->nseq;a++)
    {
      if (score[a][1]==-1);
      else
	{
	  score[b][0]=score[a][0];
	  score[b][1]=score[a][1];
	  b++;
	}
    }
  if (b!=0)
    {
      bsp=flash_sort_int (score,2,1,0, b-1);
      bs=bsp[0];
    }
  else
    bs=-1;


  free_int (score, -1);
  return bs;

}
Constraint_list *add_seq2cl(int s, Constraint_list *CL)
{
  static char *master;
  static char *seq;
  static char *lib;
  static int dumped;


  if (!dumped)
    {
      int a;
      seq   =vtmpnam(NULL);

      for (a=0; a<(CL->S)->nseq; a++)
	printf_file (seq, "a", ">%s\n%s\n", (CL->S)->name[a], (CL->S)->seq[a]);
      dumped=1;
    }
  master=vtmpnam(NULL);
  lib   =vtmpnam(NULL);
  printf_file (master, "w", ">%s\n", (CL->S)->name[s]);

  printf_system ("cp %s seq",seq);
  printf_system ("cp %s master",master);



  printf_system ("t_coffee -seq %s -master %s -out_lib %s -lib_only -quiet", seq, master,lib);

  return read_constraint_list_file(CL,lib);
}




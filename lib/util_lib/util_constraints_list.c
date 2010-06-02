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
int compare_constraint_list_entry ( const void*vx, const void*vy);
int compare_constraint_list_entry4bsearch ( const void*vx, const void*vy);

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
Constraint_list *make_test_lib (Constraint_list *CL);


Constraint_list *fork_line_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode,Job_TC *job, FILE *local_stderr);	
Constraint_list *fork_cell_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode,Job_TC *job, FILE *local_stderr);
Constraint_list *nfork_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode,Job_TC *job, FILE *local_stderr);
Constraint_list *fork_subset_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, FILE *local_stderr);
int job2first_seq(Job_TC *job);
Constraint_list *produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode)	
{
  Job_TC *job=NULL;
  FILE *local_stderr;
  int njob;
 
  store_seq_type (S->type);
  if ( CL==NULL)CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?vtmpfile():NULL, NULL);
  local_stderr=(CL->local_stderr!=NULL)?CL->local_stderr:stderr;
  
  CL->local_stderr=vfopen("/dev/null", "w");
  job=queue2heap(method2job_list ( method,S,weight, CL->lib_list,CL->DM, CL));
  njob=queue2n(job)+1;
  
  if ( get_nproc()==1 || njob==1)return  nfork_produce_list (CL, S, method, weight, mem_mode,job, local_stderr);
  else if (strstr ( CL->multi_thread, "jobcells"))return fork_cell_produce_list (CL, S, method, weight, mem_mode,job,local_stderr);
  else if (strstr ( CL->multi_thread, "joblines"))return fork_line_produce_list (CL, S, method, weight, mem_mode,job, local_stderr);
  else if (strstr ( CL->multi_thread, "jobs"))return fork_subset_produce_list (CL, S, method, weight, mem_mode,job, local_stderr); //Recommended default
  else return nfork_produce_list (CL, S, method, weight, mem_mode,job, local_stderr);
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

Constraint_list *fork_subset_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, FILE *local_stderr)
{
  //forks lines of the matrix
  int a,b;
  Job_TC *heap,*end,*start, ***jl;
  TC_method *M;
  int n_elements_in, n_new_elements;
  char **pid_tmpfile;
  int   *pid_list;
  int pid,npid, njob;
  int nproc, max_nproc, submited;
  int cseq, seq, nlines;
  int n_aln;
  
  max_nproc=nproc=get_nproc();
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
  
  job=queue2heap (job_list2multi_thread_job_list (job, CL->multi_thread, CL));
  heap=job;
  njob=queue2n (job);
  n_elements_in=CL->ne;
  M=(job->param)->TCM;
  if (M)M->PW_CL=method2pw_cl ( M, CL);
  pid_tmpfile=vcalloc (MAX(njob,get_nproc())+1, sizeof (char*));
  pid_list   =vcalloc (MAX_N_PID, sizeof (int *));
  
  fprintf ( local_stderr, "\n\tMulti Core Mode: %d processors [subset]\n", get_nproc());
  
  jl=split_job_list(job,get_nproc());
  a=npid=0;
  
  a=npid=0;
  while (jl[a])
    {
 
      start=job=jl[a][0];
      end=jl[a][1];
      pid_tmpfile[a]=vtmpnam(NULL);
      pid=vvfork(NULL);
    
      if (pid==0)//child process
	{
	  FILE *fp;
	  int done, todo, t;
	  
	  initiate_vtmpnam(NULL);
	  fp=vfopen (pid_tmpfile[a],"a");
	  todo=0;
	  while (job!=end){todo++;job=job->c;}
	  job=start;
	  
	  done=0;
	  while (job!=end)
	    {
	      if (a==0)output_completion ( local_stderr,done,todo,1, "Submit   Job");
	      job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",duplicate_constraint_list4lib_computation (CL),submit_lib_job, retrieve_lib_job, CL->multi_thread );	
	      
	      job=submit_job (job);
	      retrieve_job (job);
	      t=((job->io)->CL)->ne*((job->io)->CL)->entry_len;
	      for (b=0; b<t; b++)fprintf ( fp, "%d ", ((job->io)->CL)->L[b]);
	      job=job->c;
	      done++;
	    }
	  vfclose (fp);
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
      CL=read_constraint_list_raw_file (CL,pid_tmpfile[pid_list[pid]]);
      remove(pid_tmpfile[pid_list[pid]]);
    }
  
  vfree (pid_list);
  vfree (pid_tmpfile);
  
  job=heap;
  while (job)	job=delete_job (job);
  n_new_elements=CL->ne - n_elements_in;
  compact_list (CL, n_elements_in,n_new_elements, "best");
  compact_list (CL, 0, CL->ne, "default");
  CL->local_stderr=local_stderr;
  free_queue  (heap);
  return CL;
}


Constraint_list *fork_line_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, FILE *local_stderr)
{
  //forks lines of the matrix
  int a,b;
  Job_TC *heap;
  TC_method *M;
  int n_elements_in, n_new_elements;
  char **pid_tmpfile;
  int   *pid_list;
  int pid,npid, njob;
  int nproc, max_nproc, submited;
  int cseq, seq, nlines;
  
  max_nproc=nproc=get_nproc();
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
  
  
  job=queue2heap (job_list2multi_thread_job_list (job, CL->multi_thread, CL));
  
  heap=job;
  M=(job->param)->TCM;
  if (M)M->PW_CL=method2pw_cl ( M, CL);
  n_elements_in=CL->ne;    
  
  /* Cf. parse method for possible out_mode flags*/
  
  njob=queue2n(job)+1;
  pid_tmpfile=vcalloc (njob, sizeof (char*));
 
  pid_list   =vcalloc (MAX_N_PID, sizeof (int *));
  fprintf ( local_stderr, "\n\tMulti Core Mode: %d processors [jobline]\n", get_nproc());
 
  
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
		    job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",duplicate_constraint_list4lib_computation (CL),submit_lib_job, retrieve_lib_job, CL->multi_thread );	
		    job=submit_job (job);
		    retrieve_job (job);
		    constraint_list2raw_file ((job->io)->CL,pid_tmpfile[npid], "a");
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
		if (submited>max_nproc)
		  {
		    //wait for nproc
		    for (a=0; a<nproc; a++)
		      {
			int index;
			local_stderr=output_completion ( local_stderr,npid,nlines,1, "Processed   Job");
			pid=vwait(NULL);
			index=pid_list[pid];
			CL=read_constraint_list_raw_file (CL,pid_tmpfile[index]);
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
      CL=read_constraint_list_raw_file (CL,pid_tmpfile[index]);
      remove (pid_tmpfile[index]);
    }
  
  vfree (pid_list);
  vfree (pid_tmpfile);
  
  n_new_elements=CL->ne - n_elements_in;
  compact_list (CL, n_elements_in,n_new_elements, "best");
  compact_list (CL, 0, CL->ne, "default");
  CL->local_stderr=local_stderr;
  
  free_queue  (heap);
  
  return CL;
}

Constraint_list *fork_cell_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, FILE *local_stderr)
    {
      //forks cells of the matrix
      int a,b, n;
      Job_TC *heap;
      TC_method *M;
      int n_elements_in, n_new_elements;
      int *pid_list;
      char **pid_tmpfile;
      int pid,npid, njob;
      int nproc, max_nproc;
      int submited;
      
      max_nproc=nproc=get_nproc();
            
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

      
    job=queue2heap (job_list2multi_thread_job_list (job, CL->multi_thread, CL));
    
    heap=job;
    M=(job->param)->TCM;
    if (M)M->PW_CL=method2pw_cl ( M, CL);
    n_elements_in=CL->ne;    
    
    /* Cf. parse method for possible out_mode flags*/
    
    njob=queue2n(job)+1;
    pid_tmpfile=vcalloc (njob, sizeof (char*));
    pid_list   =vcalloc (MAX_N_PID, sizeof (int *));
    
    fprintf ( local_stderr, "\n\tMulti Core Mode: %d processors:\n", get_nproc());
    npid=0;
    submited=0;
    while (job)
      {
	job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",duplicate_constraint_list4lib_computation (CL),submit_lib_job, retrieve_lib_job, CL->multi_thread );	
	pid_tmpfile[npid]=vtmpnam(NULL);
	pid=vvfork (NULL);
	if ( pid==0)
	  {
	    initiate_vtmpnam (NULL);
	    job=submit_job (job);
	    retrieve_job (job);
	    constraint_list2raw_file ((job->io)->CL,pid_tmpfile[npid], "w");
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
		  CL=read_constraint_list_raw_file (CL,pid_tmpfile[index]);
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
	CL=read_constraint_list_raw_file (CL,pid_tmpfile[index]);
	remove (pid_tmpfile[index]);
      }
    vfree (pid_list);
    vfree (pid_tmpfile);
    
    n_new_elements=CL->ne - n_elements_in;
    compact_list (CL, n_elements_in,n_new_elements, "best");
    compact_list (CL, 0, CL->ne, "default");
    CL->local_stderr=local_stderr;
    
    free_queue  (heap);
    
    return CL;
    }
Constraint_list *nfork_produce_list   ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode, Job_TC *job, FILE *local_stderr)
    {
    int b;
    int n_aln;
    Job_TC *heap;
    TC_method *M;
    int n_elements_in, n_new_elements;
        
    
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

    
    job=queue2heap (job_list2multi_thread_job_list (job, CL->multi_thread, CL));
    
    heap=job;
    n_aln=queue2n (job);
    M=(job->param)->TCM;

    if (M)M->PW_CL=method2pw_cl ( M, CL);

   
    n_elements_in=CL->ne;    
    
    
    /* Cf. parse method for possible out_mode flags*/
    
   
       
    b=0;

    while (job)
      {
	local_stderr=output_completion ( local_stderr, b, n_aln+1,1, "Submit   Job");
	job=print_lib_job (job, "io->CL=%p control->submitF=%p control->retrieveF=%p control->mode=%s",duplicate_constraint_list4lib_computation (CL),submit_lib_job, retrieve_lib_job, CL->multi_thread );	
	

	job=submit_job (job);

	if (retrieve_job (job))	
	  {
	    
	    CL=merge_constraint_list    ((job->io)->CL, CL, "default");
	    free_constraint_list4lib_computation ( (job->io)->CL);	
	    
	  }
	job=job->c;
	b++;
      }
    job=heap;
    fprintf ( local_stderr, "\n");
    
    while (job)	job=delete_job (job);
  
    n_new_elements=CL->ne - n_elements_in;
    compact_list (CL, n_elements_in,n_new_elements, "best");
    compact_list (CL, 0, CL->ne, "default");
    CL->local_stderr=local_stderr;
    
    free_queue  (heap);
    
    return CL;
    }


Job_TC *job_list2multi_thread_job_list (Job_TC* ojob, char *mt, Constraint_list *CL)
{
  FILE *fp=NULL;
  int mtv, n, nj;
  char *met_file, *seq_file, *lib_file, *lib_list;
  char common[1000], command[1000], T_file[100];
  Job_TC *njob, *heap;
  TC_method *TCM;
 
  
  mtv=(mt==NULL)?0:atoi (mt);
  if ( !mtv || mtv==1)return ojob;
  
  HERE ("***"); myexit (0);
  heap=ojob;
  nj=(queue2n(ojob)/mtv)+1;
  nj=(nj==0)?1:nj;
  
  met_file=vtmpnam (NULL);  
  TC_method2method_file ((ojob->param)->TCM,met_file=vtmpnam (NULL));  
  
  TCM=method_file2TC_method (method_name2method_file ("tc2"));
  sprintf (T_file,"%s",(CL->S)->template_file);
  
  HERE ("2");

  seq_file=vtmpnam (NULL);
  sprintf ( common, "%s -in M%s S%s -lib_only ", TCM->executable, met_file, seq_file);
  if ( T_file[0] && check_file_exists (T_file))
    {
      strcat (common, " -template_file ");strcat (common, T_file); strcat ( common , " ");
    }

  njob=vcalloc ( 1, sizeof (Job_TC));
  njob->jobid=-1;
  n=0;

  HERE ("3");

  while ( ojob && n<nj)
    {
      if (n==0)
	{
	  lib_file=vtmpnam ( NULL);
	  lib_list=vtmpnam ( NULL);
	  fp=vfopen (lib_list, "w");
	  sprintf ( command, "%s -out_lib=%s -lib_list=%s -quiet", common, lib_file, lib_list);
	  njob->c=print_lib_job (NULL, "param->TCM=%p param->method=%s param->aln_c=%s io->in=%s io->out=%s",TCM, "tc2", command,seq_file, lib_file);
	  
	  njob=queue_cat (njob, njob->c);
	}

      fprintf ( fp, "%s\n", (ojob->param)->seq_c);
      ojob=ojob->c;
      if ( (++n)==nj){vfclose (fp);n=0;}
    }
  HERE ("4");
  if (fp && n)vfclose (fp);
  free_queue  (heap);
  return njob;
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
      static char *log_output;
      static int log;

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
  static int   header;
  static int log;
  static char *file, *log_file;
  static int set;
  
  if ( set && log_file==NULL && l==NULL) return 0;
  if (!set ){log_file=get_string_variable ("method_log");if (log_file && strm (log_file, "no"))log_file=NULL; set=1;}
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
  static char *l;
  static int log;
  
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

      seq_list2in_file ( M, (io->CL)->S, p->seq_c, io->in);
      printf_system ("%s ::IGNORE_FAILURE::", p->aln_c);
      add_method_output2method_log (NULL,NULL, NULL, NULL, io->out);
      if (!evaluate_sys_call_io (io->out,p->aln_c, "") || (strm (M->out_mode, "aln") && !(is_aln (io->out) || is_seq(io->out))) ) 
	{
	  job->status=EXIT_FAILURE;
	  //myexit (EXIT_FAILURE);
	  return job;
	}
    }
  else if ( strm2 (M->out_mode, "fA", "fL"))
    {
      io->CL= seq2list(job);
      if (!io->CL)
	{
	  add_warning (stderr, "\nFAILED TO EXECUTE:%s [SERIOUS:%s]", p->aln_c, PROGRAM);     
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
	  fname=vcalloc ( 1000, sizeof (char));
	  bufS=vcalloc ( S->nseq*10, sizeof (char));
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
		  add_warning ( stderr, "\nWARNING: %s is not part of the sequence dataset \n", list[i]);
		  continue;
		}
	      }
	    sprintf ( bufS, "%s", list[1]);	    
	    for ( i=2; i<n+2; i++) {strcat (bufS, " ");strcat ( bufS, list[i]);}
	    	    
	   
	    bufA=make_aln_command (method, in=vtmpnam(NULL),out=vtmpnam(NULL));
	    
	    if (strrchr(bufA, '>')==NULL)strcat (bufA,TO_NULL_DEVICE);
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
	char buf[1000];
	
	sprintf (bufS, "%d",S->nseq);
	for (d=0; d< S->nseq; d++)
	  {
	    sprintf ( buf," %d",d);
	    strcat ( bufS, buf);
	  }
	
	bufA=make_aln_command (method, in=vtmpnam(NULL),out=vtmpnam(NULL));
	
	if (strrchr(bufA, '>')==NULL)strcat (bufA,TO_NULL_DEVICE);
	
	if ( check_seq_type ( method, bufS, S))
	  {
	    
	    job->c=print_lib_job (NULL, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ", method, fname, bufA, bufS, in, out, S->template_file);	      
	    
	    job=queue_cat (job, job->c);
	    
	  }
	vfree (bufA);
	
      }
    else if ( strstr(aln_mode, "pairwise"))
      {
	
	int do_mirror, do_self, x, y, id;
	do_mirror=(strstr(aln_mode, "m_"))?1:0;
	do_self=(strstr(aln_mode, "s_"))?1:0;

	
	for (x=0; x< S->nseq; x++)
	  for ( y=(do_mirror)?0:x; y< S->nseq; y++)
	    {
	      
	      id=DM->similarity_matrix[x][y];
	      
	      if ( x==y && !do_self);
	      else if ( !is_in_range(id,method->minid, method->maxid));  
	      else
		{
		  sprintf (bufS, "2 %d %d",x,y);
		  bufA=make_aln_command (method,in=vtmpnam(NULL),out=vtmpnam (NULL));
		  
		  if (strrchr(bufA, '>')==NULL)strcat (bufA, TO_NULL_DEVICE);
		  if (check_seq_type (method, bufS, S))
		    {
		      job->c=print_lib_job (job->c, "param->TCM=%p param->method=%s param->aln_c=%s param->seq_c=%s io->in=%s io->out=%s ",method,fname,bufA, bufS, in, out);
		      job=queue_cat (job, job->c);
		    }
		  else if ( method->seq_type[0]=='P' && hijack_P_jobs)
		    {
		      //Hijack _P_ jobs without enough templates
		      static TC_method *proba_pairM;
		      
		      fprintf (stderr, "\n\t Information: Method %s cannot be applied to [%s vs %s]. Use proba_pair instead", method->executable, (CL->S)->name[x], (CL->S)->name [y]);
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
  
  if ( !buf)buf=vcalloc ( 1000, sizeof (char));
  
  if ( !list || n==0)return list;
  buf=vcalloc ( 1000, sizeof (char));
    
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
  int n=0, a, ml1=0, ml2=0, ml3=0;
  int status;
  
  
  
  

  
  list=produce_method_file (NULL);
  l2=declare_arrayN(3,sizeof (char), 1000, 10, 100);
  
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
      if (l2[a][4][0])l2[a][5]=pg2path (l2[a][1]);
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
  if (!buf) buf=vcalloc ( 1000, sizeof (char));
  else buf[0]='\0';
  
  if (!line)
    {
      line=vcalloc (LONG_STRING+1, sizeof ( char));
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
      line=vcalloc (LONG_STRING+1, sizeof ( char));
      subcommand=vcalloc ( LONG_STRING, sizeof (char));
    }

  m=vcalloc ( 1, sizeof (TC_method));
  
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
      else if ( (p=strstr (line, "IN_FLAG2"    ))) sscanf (p, "IN_FLAG2 %s"   , m->in_flag2);
      else if ( (p=strstr (line, "IN_FLAG"    ))) sscanf (p, "IN_FLAG %s"   , m->in_flag);
      else if ( (p=strstr (line, "OUT_FLAG"   ))) sscanf (p, "OUT_FLAG %s"  , m->out_flag);
      else if ( (p=strstr (line, "OUT_MODE"   ))) sscanf (p, "OUT_MODE %s"  , m->out_mode);
      else if ( (p=strstr (line, "ALN_MODE"   ))) sscanf (p, "ALN_MODE %s"  , m->aln_mode);
      else if ( (p=strstr (line, "SEQ_TYPE"   ))) sscanf (p, "SEQ_TYPE %s"  , m->seq_type);
      else if ( (p=strstr (line, "WEIGHT"     ))) sscanf (p, "WEIGHT %s"  , m->weight);
      else if ( (p=strstr (line, "MATRIX"     ))){ sscanf (p, "MATRIX %s"  , m->matrix);}
      else if ( (p=strstr (line, "GOP"        ))) sscanf (p, "GOP %d"  , &m->gop);
      else if ( (p=strstr (line, "GEP"        ))) sscanf (p, "GEP %d"  , &m->gep);
      else if ( (p=strstr (line, "MAXID"      ))) sscanf (p, "MAXID %d"  , &m->maxid);
      else if ( (p=strstr (line, "MINID"      ))) sscanf (p, "MINID %d"  , &m->minid);
     
    }
  vfclose ( fp);



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
  
  vfclose ( fp);
  return 1;
}

char *make_aln_command(TC_method *m, char *seq, char *aln)
    {
      char *command;
      char buf[1000];
    
      //      sprintf ( buf, "%s %s %s%s %s%s %s", m->executable, m->param1, m->in_flag, seq,m->param2, m->out_flag,aln, m->param);
      sprintf ( buf, "%s %s %s%s %s %s%s %s", m->executable, m->param1, m->in_flag, seq,m->param2, m->out_flag,aln, m->param);
      command=vcalloc ( strlen (buf)+100, sizeof (char));
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


int dump_constraint_list (Constraint_list *CL)
{
  int *L;
  int a, b;
  if (!CL->L) return 0;
  L=CL->L;
  CL->L=NULL;
  CL->fp=tmpfile ();
  
  for ( a=0; a<CL->ne; a++)
    for (b=0; b<CL->entry_len; b++)
      vwrite_clist (CL, a, b, L[a*CL->entry_len+b]);
  
  vfree (CL->L); CL->L=NULL;
  return 1;
}

int vread_clist ( Constraint_list *CL, int a, int b )
    {
    int x;
    
    if ( a>= CL->ne)
      {
	return UNDEFINED_2;
	
      }
    else if (CL->fp)
       {       
	 
	 fseek (CL->fp, a*CL->el_size*CL->entry_len+b*CL->el_size, SEEK_SET);
	 fread (&x, CL->el_size, 1, CL->fp);
	 
	 return x;
       }
    else if ( CL->L)
       {
	 return CL->L[a*CL->entry_len+b];
       }
    else if (CL->M)
       {
       return (CL->M)[a][b];
       }
    else 
       {
       return UNDEFINED_2; 
       
       }
    return UNDEFINED_2;
    }
int vwrite_clist ( Constraint_list *CL, int a, int b, CLIST_TYPE x)
    {

      CL->seq_indexed=0;
      CL->residue_indexed=0;


      if (CL->fp)
	{
	  fseek (CL->fp, a*CL->el_size*CL->entry_len+b*CL->el_size, SEEK_SET);
	  fwrite(&x, CL->el_size, 1, CL->fp);
	}
      else if (!CL->M)
	{
	  int i,l;
	  
	  
	  i=a*CL->entry_len+b;
	  if (CL->L)
	    {
	      Memcontrol *p;
	      p=(Memcontrol *)CL->L;
	      p-=2;
	      l=(int)p[0].size/p[0].size_element;
	      //read_size_int (CL->L,sizeof (int));
	    }
	  else 
	    {
	      l=CL->chunk;
	      (CL->L)=vcalloc (l,sizeof (int));
	    }
	  
	  if (l<=i)
	    {
	      l+=CL->chunk;
	      (CL->L)=vrealloc ( (CL->L),l*sizeof (int));
	    }
	  (CL->L)[i]=x;
#ifdef gagaga
	  if (((a*CL->entry_len)+b)>=(CL->entry_len*CL->max_L_len))
	    {
	      
	      if (!(CL->L))
		{
		  (CL->L)=vcalloc ((CL->chunk+a)*CL->entry_len+b, sizeof (int));
		}
	      else
		{
		  CL->max_L_len+=CL->chunk;
		  (CL->L)=vrealloc ( (CL->L),CL->entry_len*sizeof (int)*CL->max_L_len+b);
		}
	    }
	  
	    (CL->L)[i]=x;
#endif
       }
    return x;
    }


/*********************************************************************/
/*                                                                   */
/*                        INDEXING FUNCTIONS                         */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

Constraint_list *index_res_constraint_list ( Constraint_list *CL, int field)
        {
	/*
	  function documentation: start
	  Constraint_list *index_res_constraint_list ( Constraint_list *CL, int field)
	  This function reorganises the content of the CL->L matrix, so that a single look up gives
	  the constraints associated with any residue
	  
	  1-if CL->residue_indexed=1 return
		 2-CL->residue_index[Seq X][Res Y of Seq X][0]=Z
		 Z=Number of residues matching X(Y)*3+1
		 CL->residue_index[Seq X][Res Y of Seq X][0]=Z
		 CL->residue_index[Seq X][Res Y of Seq X][c+0]=seq W
		 CL->residue_index[Seq X][Res Y of Seq X][c+1]=res V
		 CL->residue_index[Seq X][Res Y of Seq X][c+2]=weight W(V) Vs X(Y)		 
	  
	  NOTE: Works All right with a sequence against itself
	  NOTE: Any modification of CL->L should result in residue_indexed to be reset	  
	  function documentation: end
	*/
	int a, b, s1, s2, r1, r2, w;





	  

	if ( CL->residue_indexed && CL->residue_field==field);
	else
	   {

	     if ( CL->residue_index)
	          {
		    for ( a=0; a< (CL->S)->nseq; a++)
		      for ( b=0; b<= (CL->S)->len[a]; b++)
		          {
			      vfree(CL->residue_index[a][b]);
			      CL->residue_index[a][b]=vcalloc_nomemset (1, sizeof (int));
			      CL->residue_index[a][b][0]=1;
			  }    
		  }
	       else if ( !CL->residue_index)
	          {
		 

		    CL->residue_index=vcalloc_nomemset ( (CL->S)->nseq, sizeof (int**));
		    for ( a=0; a< (CL->S)->nseq; a++)
		          {
			
			    CL->residue_index[a]=vcalloc_nomemset ( ((CL->S)->len[a]+1), sizeof (int*));
			      for ( b=0; b<= (CL->S)->len[a]; b++)
				  {
				    CL->residue_index[a][b]=vcalloc_nomemset (1,sizeof (int));
				  CL->residue_index[a][b][0]=1;
				  }
			  }
		  }
	       for (a=0;a<CL->ne; a++)
	          {
		    
		  s1=vread_clist (CL, a, SEQ1);
		  s2=vread_clist (CL, a, SEQ2);
		  r1=vread_clist (CL, a, R1);
		  r2=vread_clist (CL, a, R2);
		  w=vread_clist (CL, a, field);
		  
	
		  CL->residue_index[s1][r1][0]+=3;
		  CL->residue_index[s1][r1]=vrealloc ( CL->residue_index[s1][r1], CL->residue_index[s1][r1][0]*sizeof (int));
		  CL->residue_index[s1][r1][CL->residue_index[s1][r1][0]-3]=s2;
		  CL->residue_index[s1][r1][CL->residue_index[s1][r1][0]-2]=r2;
		  CL->residue_index[s1][r1][CL->residue_index[s1][r1][0]-1]=w;

		  CL->residue_index[s2][r2][0]+=3;
		  CL->residue_index[s2][r2]=vrealloc ( CL->residue_index[s2][r2], CL->residue_index[s2][r2][0]*sizeof (int));
		  CL->residue_index[s2][r2][CL->residue_index[s2][r2][0]-3]=s1;
		  CL->residue_index[s2][r2][CL->residue_index[s2][r2][0]-2]=r1;
		  CL->residue_index[s2][r2][CL->residue_index[s2][r2][0]-1]=w;
		  
		  }
	   CL->residue_indexed=1;
	   CL->residue_field=field;
	   }
	return CL;
	}
      
Constraint_list *index_constraint_list ( Constraint_list *CL)
        {
	  /*
	  Function Documentation: start
	  Constraint_list *index_constraint_list ( Constraint_list *CL);
	  Indexes a constraint list
	     1-Checks the flag seq_indexed
	     2-if flag set to 0
	       CL->start_index[seq1][seq2] indicatse the first position for seq1 Vs seq2
	       CL->end_index[seq1][seq3]   indicatse the last  position for seq1 Vs seq2
	  Any modif to CL->L should cause the flag 1 to be set to 0;
	  Function Documentation: end
	*/
	int a, csA, csB, sA, sB;


	if ( CL->seq_indexed);
	else
	    {
	    if ( CL->start_index!=NULL)free_int ( CL->start_index,-1);
	    CL->start_index=declare_int ( (CL->S)->nseq, (CL->S)->nseq);
	    
	    if ( CL->end_index!=NULL)free_int ( CL->end_index,-1);
	    CL->end_index=declare_int ( (CL->S)->nseq, (CL->S)->nseq);
	    
	    csA=vread_clist (CL, 0, SEQ1);
	    csB=vread_clist (CL, 0, SEQ2);

	    CL->start_index[csA][csB]=0;
	    CL->start_index[csB][csA]=0;
	    for ( a=1; a<CL->ne; a++)
	        {
		sA=vread_clist (CL, a, SEQ1);
		sB=vread_clist (CL, a, SEQ2);
		if (sA!=csA || sB!=csB)
		   {
		   CL->end_index[csA][csB]=a;
		   CL->end_index[csB][csA]=a;
		   csA=sA;
		   csB=sB;
		   CL->start_index[csA][csB]=a;
		   CL->start_index[csB][csA]=a;
		   }
		}		
	    CL->end_index[csB][csA]=CL->ne;
	    CL->end_index[csA][csB]=CL->ne;
	    CL->seq_indexed=1;
	    }
	return CL;
	}


char ** reindex_constraint_list (char **profile, int np,char **list, int *inL, Sequence *S)
{
  int a, nl=0, l;
  int ***cache, **index, *entry;
  char **nlist, **cons, *nlib_name;
  Sequence *NS;
  Alignment **P, *A;
  Constraint_list *CL, *NCL;
  

  l=inL[0];
  /*1: Pre-Process the profiles*/
  cache=vcalloc (np, sizeof (int**));
  cons=vcalloc  (np, sizeof (char*));
  P=vcalloc (np, sizeof (Alignment *));

  for ( a=0; a< np; a++)
    {
      int e,b,c, **ca;
      A=P[a]=main_read_aln (profile[a], NULL);
      
      ca=cache[a]=declare_int (A->nseq, A->len_aln+1);
      cons[a]=aln2cons_seq_mat ( A, "blosum62mt");
      for (b=0;b<A->nseq; b++)
	{
	  for (e=0,c=0; c<A->len_aln; c++)
	    {
	      if (is_gap(A->seq_al[b][c]));
	      else ca[b][++e]=c+1;
	    }
	}
    }

  /*2: Index The Sequences*/
  index=declare_int (S->nseq,2);
  for (a=0; a<S->nseq; a++)index[a][0]=index[a][1]=-1;
  for (a=0; a< S->nseq; a++)
    {
      char *name;
      int p, b;
      name=S->name[a];
      for (b=0; b< np; b++)
	{
	  if ((p=name_is_in_list(name, (P[b])->name, (P[b])->nseq, MAXNAMES))!=-1) 
	    {
	      index[a][0]=b;
	      index[a][1]=p;
	    }
	}
    }

  /*3: Read the primary Library*/
  
  nlist=declare_char (read_array_size_new ((void *)list), read_array_size_new ((void *)list[0]));
  for (nl=0,a=0; a< l; a++)
    if ( list[a][0]=='L' || list[a][0]=='A')
      nlist[nl++]=list[a];
  
  CL=declare_constraint_list ( S,NULL, NULL, 0,NULL, NULL);
  CL=read_n_constraint_list (nlist,nl,NULL,NULL,"sim",S->type,stdout, CL, "S");   
  vfree (nlist);
  
  NS=fill_sequence_struc (np, cons, profile);
  NCL=declare_constraint_list ( NS,NULL, NULL, 0,NULL, NULL);
  entry=vcalloc ( CL->entry_len, CL->el_size);
  for (a=0; a<CL->ne; a++)
    {
      int s1, s2, r1, r2, ps1, ps2, pps1, pps2;

      entry=extract_entry (entry, a, CL);
      s1=entry[SEQ1]; s2=entry[SEQ2];
      r1=entry[R1]; r2=entry[R2];
      ps1=index[s1][0];ps2=index[s2][0];
      pps1=index[s1][1];pps2=index[s2][1];

      if (ps1==ps2 || ps1==-1 || ps2==-1)continue;
      entry[SEQ1]=ps1;
      entry[SEQ2]=ps2;
    
      entry[R1]=cache[ps1][pps1][r1];
      entry[R2]=cache[ps2][pps2][r2];
    
      add_entry2list (entry, NCL);
      
    }

  compact_list (NCL, 0, NCL->ne, "default");
  nlib_name=vtmpnam(NULL);
  vfclose(save_constraint_list ( NCL, 0, NCL->ne,nlib_name, NULL, "lib",NCL->S));
  
  
  nlist=declare_char (read_array_size_new ((void *)list), read_array_size_new ((void *)list[0]));
  for (nl=0,a=0; a< l; a++)
    if ( list[a][0]!='L' && list[a][0]!='A')
      sprintf (nlist[nl++], "%s",list[a]);
  sprintf (nlist[nl++], "L%s", nlib_name);
  
  free_arrayN ((void *)cache, 3);
  free_arrayN ((void *)cons, 2);
  free_arrayN ((void *)index, 2);
  vfree (entry);
  for (a=0;a<np; a++)free_aln (P[a]);
  vfree (P);
  free_constraint_list_full (CL);
  free_constraint_list_full (NCL);
  
  inL[0]=nl;
  return nlist;
}
  

Constraint_list * progressive_index_res_constraint_list ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  int **index_seq;
  int **index_res;  
  static Constraint_list *NCL;
  
  int a, b, c, d,s,l, max_len;
  int s1, s2, r1, r2;
  int *L, *NL;
  char *seq;
  Sequence *S, *NS;


  
  if ( NCL==NULL)
    {
      NCL=duplicate_constraint_list (CL);
    }
  
  S=CL->S;
  NS=NCL->S;
  
  L=CL->L;
  NL=NCL->L;
  
  
  
  if (NCL->residue_index)
    {
      free_arrayN((void*)NCL->residue_index, 3);
      NCL->residue_index=NULL;
    }

  for ( a=0; a< S->nseq; a++)
    {
      sprintf ( NS->seq[a], "%s", S->seq[a]);
      NS->len[a]=strlen (S->seq[a]);
    }
  
  
  for ( a=0; a<2; a++)
    {
      s=A->order[ls[a][0]][0];
      vfree (NS->seq[s]);
      NS->seq[s]=sub_aln2cons_seq_mat (A, ns[a], ls[a], "blosum62mt");
      NS->len[s]=strlen (NS->seq[s]);
    }
  for ( max_len=0, a=0; a<S->nseq; a++)max_len=MAX(max_len, S->len[a]);
  index_seq=declare_int ( NS->nseq, max_len+1);
  index_res=declare_int ( NS->nseq, max_len+1);
  

  for ( a=0; a< 2; a++)
    {
      for (b=0; b< ns[a]; b++)
	{
	  s=A->order[ls[a][b]][0];
	  seq=A->seq_al[ls[a][b]];
	  l=strlen (seq);

	  for ( d=0, c=0; c<l; c++)
	    {
	      if ( !is_gap(seq[c]))
		{
		  d++;
		  
		  index_res[s][d]=c+1;
		  index_seq[s][d]=A->order[ls[a][0]][0];
		}
	    }
	}
    }

  for ( a=0; a< CL->ne; a++)
    {
      s1=L[a*CL->entry_len+SEQ1];
      s2=L[a*CL->entry_len+SEQ2];
      r1=L[a*CL->entry_len+R1];
      r2=L[a*CL->entry_len+R2];
      
      if ( index_res[s1][r1])
	{
	  NL[a*CL->entry_len+SEQ1]=index_seq[s1][r1];
	  NL[a*CL->entry_len+R1]=index_res[s1][r1];
	}
      else
	{
	  NL[a*CL->entry_len+SEQ1]=s1;
	  NL[a*CL->entry_len+R1]=r1;
	}
      
      if ( index_res[s2][r2])
	{
	  NL[a*CL->entry_len+SEQ2]=index_seq[s2][r2];
	  NL[a*CL->entry_len+R2]=index_res[s2][r2];
	}
      else
	{
	  NL[a*CL->entry_len+SEQ2]=s2;
	  NL[a*CL->entry_len+R2]=r2;
	}
    }

  NCL->ne=CL->ne;
  NCL->residue_indexed=0;
  
  NCL=compact_list (NCL, 0, NCL->ne, "best");
  
  NCL=index_res_constraint_list (NCL, WE);
  
  free_int ( index_res, -1); free_int (index_seq, -1);
  
  return NCL;
}
/*********************************************************************/
/*                                                                   */
/*                         LIST EXTENTION                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *extend_list_pair (Constraint_list *CL,char *store_mode, int s1, int s2)
        {
	static Sequence *S;
	Constraint_list *CLout;
	/*
	  function documentation: start
	  Constraint_list *extend_list_pair (Constraint_list *CL,char *store_mode, int s1, int s2)
	  This function takes a pair of sequences s1, s2 and perrforms the extention
	  It returns the incoming list CL, with CL->L[s1][s2] now extended
	  See main documentation for store_mode
	  function documentation: end
	*/

	if ( S==NULL)S=declare_sequence ((CL->S)->min_len, (CL->S)->max_len,(CL->S)->nseq); 	
	sprintf ( S->name[0], "%s",(CL->S)->name[s1]);
	sprintf ( S->name[1],"%s",(CL->S)->name[s2]);
	S->nseq=2;
	
	CLout=extend_list (CL, store_mode, CL->extend_clean_mode, CL->extend_compact_mode,CL->do_self, S);
	return CLout;		
	}
Constraint_list *extend_list (Constraint_list *CLin, char *store_mode,char *clean_mode, char *compact_mode,int do_self, Sequence *DO_LIST)
	{
	int a, b, c, d, e, f;
	int wA, wC,w, rA, rC, miscA, miscC, misc;
	static int **posA;
	static int **posC;
	int start_ba, end_ba, start_bc, end_bc, start_ac, end_ac;
	int len;
	int lenA=0;
	int lenC=0;
	int *translation;
	Constraint_list *CLout=NULL;


	/*Do not extend if the List is a Matrix*/
	if ( !CLin->L && CLin->M)
	   {
	   CLin->extend_jit=0;
	   return CLin;
	   }
	
	translation=vcalloc ( (CLin->S)->nseq, sizeof (int));
	for ( a=0; a<(CLin->S)->nseq; a++)
	    {
	    translation[a]=name_is_in_list ((CLin->S)->name[a],DO_LIST->name, DO_LIST->nseq, 100);
	    translation[a]++;/* set translation to -1+1=0 if seq not in list*/
	    }
	
	CLout=declare_constraint_list (CLin->S, NULL,NULL,0, strm("disk", store_mode)?tmpfile():NULL, NULL);
        	
	for ( a=0; a<(CLin->S)->nseq-(1-do_self); a++)
		{	        
		fprintf (CLin->local_stderr, "\nSeq %3d: %5d", a+1,CLout->ne);
		for ( c=a+(1-do_self); c<(CLin->S)->nseq; c++)
			{
			if ( translation[a] && translation[c])
			       {			       
			       get_bounds (CLin,a, c, &start_ac, &end_ac); 						      
			       for ( d=start_ac; d<end_ac; d++)
			           {
				   for ( e=0; e< CLin->entry_len; e++)
				       vwrite_clist(CLout,CLout->ne, e, vread_clist(CLin,d, e));
			           CLout->ne++;
				   }
			       
			       for ( b=0; b<(CLin->S)->nseq; b++)
			           {
				   len=strlen ( (CLin->S)->seq[b]);
				   
				   get_bounds (CLin,b, a, &start_ba, &end_ba);  
				   posA=fill_pos_matrix (CLin,start_ba, end_ba, len, posA, &lenA,(b>a));
				   
				   if ((c!=b && a!=b) ||(do_self==1))
				       {
				   
				       get_bounds (CLin, b, c, &start_bc, &end_bc);
				       posC=fill_pos_matrix (CLin, start_bc, end_bc, len, posC, &lenC, (b>c));
				       
				       for (d=1; d<=len; d++)
				           {
					   if ( posA[d][1]==0 || posC[d][1]==0);
					   else
					       {
					       for (e=2; e<=posA[d][1]+1; e+=(CLin->entry_len-4)) 
						   for ( f=2; f<=posC[d][1]+1; f+=(CLin->entry_len-4))
						       {
						       wA   =posA[d][e+1];
						       miscA=posA[d][e+2];

						       wC   =posC[d][f+1];
						       miscC=posC[d][f+2];

						       rA=posA[d][e];
						       rC=posC[d][f];
						       
						       w   =MIN(wA,wC);
						       
						       misc=MAX(miscA, miscC);
						       
						       vwrite_clist( CLout, CLout->ne, SEQ1, a);
						       vwrite_clist( CLout, CLout->ne, SEQ2, c);
						       vwrite_clist( CLout, CLout->ne, R1  ,rA);
						       vwrite_clist( CLout, CLout->ne, R2  ,rC);
						       vwrite_clist( CLout, CLout->ne, WE  , w);
						       vwrite_clist( CLout, CLout->ne, CONS, 1);
						       vwrite_clist( CLout, CLout->ne, MISC,misc);
        					       CLout->ne++;
						       }
					       }
					   }
				       }
				   }
		
			       CLout=compact_list (CLout,0,CLout->ne,"mirror_sum");
			       CLout=clean ( clean_mode,CLout, 0, CLout->ne);
			       }
			}
		}

	
	vfree (translation);
	return CLout;
	}
void get_bounds (Constraint_list *CL, int s1, int s2, int *start, int *end)
	{

	CL=index_constraint_list (CL);
	
	if ( s1>s2)SWAP(s1, s2);
	
	start[0]=CL->start_index[s1][s2];
	end  [0]=CL->end_index  [s1][s2];
	}

	
int ** fill_pos_matrix (Constraint_list *CL, int beg, int end, int slen, int **pos, int *len, int mirrored)
	{
	int small_chunck;
	int a, r1,r2;
	

	small_chunck=2*CL->entry_len;

	if ( pos==NULL)
		{
		pos=declare_int (slen+1, small_chunck);
		for ( a=0; a<=slen; a++)pos[a][0]=small_chunck;
		len[0]=slen+1;
		}
	else if ( len[0]<=slen)
		{
		free_int ( pos, len[0]);
		pos=declare_int (slen+1, small_chunck);
		for ( a=0; a<=slen; a++)pos[a][0]=small_chunck;
		len[0]=slen+1;
		}
	else
		{
		for ( a=0; a<=slen; a++)pos[a][1]=0;
		}
			
			
	
	
	for ( a=beg; a<end; a++)
		{
		
		if (!mirrored)     {r1=vread_clist (CL, a, R1);r2=vread_clist (CL, a, R2);}
		else if ( mirrored){r1=vread_clist (CL, a, R2);r2=vread_clist (CL, a, R1);}

	       if ( ((pos[r1][1])+(CL->entry_len))>pos[r1][0])
			{
			pos[r1]=vrealloc (pos[r1], (pos[r1][0]+small_chunck)*sizeof (int));
			pos[r1][0]+=small_chunck;
			}
		pos[r1][pos[r1][1]+2]=r2;
		pos[r1][pos[r1][1]+3]=vread_clist(CL,a,WE);
		pos[r1][pos[r1][1]+4]=vread_clist(CL,a,MISC);
		pos[r1][1]+=(CL->entry_len-4);		
		}
	return pos;
	}
Constraint_list * evaluate_constraint_list_reference ( Constraint_list *CL)
        {
	    static CLIST_TYPE *entry;
	    int a, b, c, s1, s2, r1, r2, w;
	    int ***max_res;
	    
	    

	    if ( CL->M)
	      {
		CL->max_value=CL->max_ext_value=20;
		
	      }
	    else 
	      {

		CL->max_value=CL->max_ext_value=0;
		max_res=vcalloc ( (CL->S)->nseq, sizeof (int**));
			
		for ( a=0; a< (CL->S)->nseq; a++)
		  {
		    max_res[a]=vcalloc ( strlen ((CL->S)->seq[a])+1, sizeof (int*));
		    for ( b=0; b<=(CL->S)->len[a]; b++)
		      {
		      max_res[a][b]=vcalloc ( (CL->S)->nseq+1, sizeof (int));
		      }
		  }
		
		for ( a=0; a< CL->ne; a++)
		  {
		    entry=extract_entry ( entry, a, CL);
		    s1=entry[SEQ1];
		    s2=entry[SEQ2];
		    r1=entry[R1];
		    r2=entry[R2];
		    w= entry[WE];
		    if ( w==UNDEFINED || ( (CL->moca) && (CL->moca)->forbiden_residues  && ((CL->moca)->forbiden_residues[s1][r1]==UNDEFINED || (CL->moca)->forbiden_residues[s2][r2]==UNDEFINED)));
		    else
		      {
			
			max_res[s1][r1][s2]+=w;
			max_res[s2][r2][s1]+=w;
			CL->max_value=MAX(w, CL->max_value);
		      }
		  }
		
		for ( a=0; a< (CL->S)->nseq; a++)
		  for ( b=1; b<=(CL->S)->len[a]; b++)
		    {
		      for ( c=0; c< (CL->S)->nseq; c++)
			{
			  max_res[a][b][(CL->S)->nseq]+= max_res[a][b][c];
			}
		      CL->max_ext_value=MAX(max_res[a][b][c],CL->max_ext_value);
		    }           
		
		for ( a=0; a< (CL->S)->nseq; a++)
		  {
		    for ( b=0; b<=(CL->S)->len[a]; b++)
		      vfree ( max_res[a][b]);
		    vfree (max_res[a]);
		  }
		CL->max_ext_value=MAX(1,CL->max_ext_value);
		vfree ( max_res);
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
	entry=vcalloc (CL->entry_len, sizeof (int));

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

Constraint_list *add_entry2list ( CLIST_TYPE *entry, Constraint_list *CL)
	{
	  return insert_entry2list (entry, CL->ne++, CL);
	}
Constraint_list* insert_entry2list(CLIST_TYPE * entry, int pos, Constraint_list *CL)
        {
	int a;
	for ( a=0; a< CL->entry_len; a++)
	    vwrite_clist ( CL,pos, a,entry[a]);
	return CL;
	}
CLIST_TYPE* extract_entry(CLIST_TYPE * entry, int pos, Constraint_list *CL)
        {
	int a;
	
	if ( entry==NULL)entry=vcalloc ( CL->entry_len, CL->el_size);
	
	for (a=0; a< CL->entry_len; a++)entry[a]=vread_clist(CL, pos, a);
	return entry;
	}

	
/*********************************************************************/
/*                                                                   */
/*                         SEARCH IN LIST (ARRAY AND FILE)           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE * compare_list (FILE *OUT, Constraint_list *CL1,Constraint_list *CL2)
	{
	int a;
	float nw_score=0;
	float w_score=0;
	int *l;

	CLIST_TYPE  *entry=NULL;
	int p;
	
	entry_len=CL1->entry_len;
	qsort ( (void *)CL1->L, CL1->ne, CL1->entry_len*sizeof (int), compare_constraint_list_entry);
	qsort ( (void *)CL2->L, CL2->ne, CL2->entry_len*sizeof (int), compare_constraint_list_entry);
	
	entry=vcalloc ( CL1->entry_len, CL1->el_size);
	for ( a=0; a<CL1->ne; a++)
		{
		entry=extract_entry (entry,a,CL1);
		if ((l=main_search_in_list_constraint (entry,&p,4,CL2))!=NULL)
			{
			vwrite_clist ( CL2, p,MISC, 1);
			vwrite_clist ( CL1, a,MISC, 1);
			nw_score++;
			w_score+=l[WE];
			}
		}
	fprintf ( OUT, "%-15s:%d pairs (Evaluated matrix), %d pairs in the other (%s)\n", CL2->list_name, CL2->ne, CL1->ne, CL1->list_name);
        fprintf ( OUT, "%-15s:%d pairs\n", CL1->list_name, CL1->ne);
        fprintf ( OUT, "Acurracy=%.2f%%\n", (nw_score/CL1->ne)*MAXID);
        fprintf ( OUT, "Sensitiv=%.2f%%\n\n", (nw_score/CL2->ne)*MAXID);
	return OUT;
	}


CLIST_TYPE *main_search_in_list_constraint ( int *key,int *p,int k_len,Constraint_list *CL)
	{
	  
	
	  CLIST_TYPE *l=NULL;
	  int start, end;

	  CL=index_constraint_list (CL);
	  
	  start=CL->start_index[key[SEQ1]][key[SEQ2]];
	  end  =CL->end_index  [key[SEQ1]][key[SEQ2]];
	  
	
	  	  
	  entry_len=CL->entry_len;
	  l=bsearch (key, (CL->L)+(start*CL->entry_len), (end-start), sizeof (int)*CL->entry_len,  compare_constraint_list_entry4bsearch);
	  p[0]=CL->L-l;
	  return l;
	}
CLIST_TYPE return_max_constraint_list ( Constraint_list *CL, int field)
        {
	CLIST_TYPE max=0;
	int a;
	for ( a=0; a< CL->ne; a++)max=MAX( vread_clist(CL,a,field), max);
	return max;
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
	CL=sort_constraint_list   (CL, start,len);
	

	CL=invert_constraint_list (CL, start,len);
	if ( start+len==CL->ne)
	    {	   
	    while (vread_clist(CL,CL->ne-1, SEQ1)==-1)CL->ne--;
	    }
	

	return CL;
	}

Constraint_list *invert_constraint_list (Constraint_list *CL, int start,int len)
        {
	int a, b, c;
	CLIST_TYPE tp;
	
	
	for ( a=start, b=start+len-1; a<=b; a++, b--)
	    {
	    for (c=0; c< CL->entry_len; c++)
	        {
		tp=vread_clist(CL, a, c);
		vwrite_clist(CL,a, c, vread_clist(CL, b, c));
		vwrite_clist(CL,b, c, tp);
		}
	    }
	return CL;
	}
	
Constraint_list * sort_constraint_list(Constraint_list *CL, int start, int len)
        {

	CL=sort_constraint_list_on_n_fields (CL, start, len, 0, CL->entry_len);

	return CL;
	}

Constraint_list * sort_constraint_list_on_n_fields (Constraint_list *CL, int start, int len, int first_field, int n_fields)
	{
	  entry_len=CL->entry_len;
	  
	  if (CL->fp)
	    {	   
	      rewind( CL->fp);
	      fseek      ( CL->fp, start*CL->el_size*CL->entry_len , SEEK_SET);
	      hsort_list_file ( CL->fp,        len, CL->el_size, CL->entry_len,first_field,n_fields);
	    }
	  else if ( CL->L)
	    {
	      qsort ( (void *)CL->L, CL->ne, CL->entry_len*sizeof (int), compare_constraint_list_entry);
	      //hsort_list_array ((void**)(CL->L)+start, len, CL->el_size, CL->entry_len,first_field,n_fields);
	    }
	  return CL;
	}
int compare_constraint_list_entry4bsearch ( const void*vx, const void*vy)
{
  int a;
  const int *x=vx, *y=vy;
  for (a=SEQ1; a<=R2; a++)
    {
      if (x[a]<y[a])return -1;
      else if (x[a]>y[a]) return 1;
    }
  return 0;
}
int compare_constraint_list_entry ( const void*vx, const void*vy)
{
  int a;
  const int *x=vx, *y=vy;
  
  for (a=SEQ1; a<=WE; a++)
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
Constraint_list*  fork_read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source);
Constraint_list* nfork_read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source);

Constraint_list* read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source)
{
  
  if ( get_nproc()==1 || n_list<=2)return nfork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode,type,local_stderr, CL, seq_source);
  else if ( strstr (CL->multi_thread, "methods"))
    return fork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode,type,local_stderr, CL, seq_source);
  else
    return nfork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode,type,local_stderr, CL, seq_source);
}
Constraint_list* fork_read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source)
    {
    int a, b;
    Sequence *S;
    char **tmp_list;
    int*proclist;
    int nproc, ns;

    nproc=get_nproc();
    
    proclist=vcalloc (MAX_N_PID, sizeof (int));
    tmp_list=vcalloc (n_list+1, sizeof (char*));
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
	return nfork_read_n_constraint_list(fname,n_list, in_mode,mem_mode,weight_mode, type,local_stderr, CL, seq_source);
        }   
      
    if (!CL)CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?tmpfile():NULL, NULL);    

    if (CL->ne)
      {
	constraint_list2raw_file(CL,tmp_list[n_list], "w");
	CL->ne=0;
      }

    CL->local_stderr=local_stderr;
    fprintf ( local_stderr, "\n\tMulti Core Mode: %d processors:\n", nproc);
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
	    if (CL->ne>in)constraint_list2raw_file (CL,tmp_list[a], "w");
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
		    CL=read_constraint_list_raw_file (CL,tmp_list[b]);
		    compact_list (CL, 0, CL->ne, (strm (weight_mode, "cons")?"cons":"default"));
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
	    CL=read_constraint_list_raw_file (CL,tmp_list[a]);
	    compact_list (CL, 0, CL->ne, (strm (weight_mode, "cons")?"cons":"default"));
	  }
	ns--;
      }
    fprintf ( local_stderr, "\n\n\tAll Methods Retrieved\n");
    
    if (tmp_list[n_list] && check_file_exists (tmp_list[n_list]))
	  {
	    CL=read_constraint_list_raw_file (CL,tmp_list[n_list]);
	    compact_list (CL, 0, CL->ne, (strm (weight_mode, "cons")?"cons":"default"));
	  }
    
    CL->local_stderr=local_stderr;
    CL=evaluate_constraint_list_reference (CL);
    vfree (proclist);
    vfree (tmp_list);
    return CL;
    }
Constraint_list* nfork_read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr, Constraint_list *CL, char *seq_source)
    {
    int a, b;
    Sequence *S;

    
    if (!(CL->S) && (S=read_seq_in_n_list (fname, n_list,type,seq_source))==NULL)
	{
	fprintf ( stderr, "\nNO SEQUENCE WAS SPECIFIED[FATAL]\n");
        myexit(EXIT_FAILURE);
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

        }   
    
    if (!CL)CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?tmpfile():NULL, NULL);    
    CL->local_stderr=local_stderr;
    fprintf ( CL->local_stderr,"\nREAD/MAKE LIBRARIES:[%d]\n",n_list );

    CL=read_constraint_list (CL, fname[0], in_mode, mem_mode,weight_mode);
    compact_list (CL, 0, CL->ne, "default");
    for ( a=1; a< n_list; a++)
        {	
	CL=read_constraint_list (CL, fname[a], in_mode, mem_mode,weight_mode);	
	if (strm (weight_mode, "cons"))compact_list (CL, 0, CL->ne, "cons");
	else
	  compact_list (CL, 0, CL->ne, "default");
	}
    CL->local_stderr=local_stderr;
    
    CL=evaluate_constraint_list_reference (CL);
    
    return CL;
    }
Constraint_list* read_constraint_list(Constraint_list *CL,char *in_fname,char *in_mode, char *mem_mode,char *weight_mode)
        {
	Sequence *SL=NULL, *TS=NULL;
	int a;
	Constraint_list *SUBCL=NULL;
	static char *read_mode;
	char *fname;
	 
	fname=in_fname;	
	if ( !read_mode)read_mode=vcalloc ( STRING, sizeof (char));
	
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
	      CL->L=NULL;
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
		SUBCL=aln_file2constraint_list ( fname,SUBCL,weight_mode);
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

#define is_seq_source(Symbol,Mode,SeqMode)            (Symbol==Mode && (SeqMode==NULL || strm (SeqMode, "ANY") || (SeqMode[0]!='_' && strchr (SeqMode,Symbol)) || (SeqMode[0]=='_' && !strchr (SeqMode,Symbol)))) 
Sequence * read_seq_in_n_list(char **fname, int n, char *type, char *SeqMode)
        {
	int nseq=0;
	int a, b;
	Alignment *A;
	char **sequences=NULL;
	char **seq_name=NULL;
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
	       	
	       buf=name2type_name(fname[a]);mode=buf[0];lname=buf+1;
	      
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
			add_warning ( stderr, "\nWarning: Could not use PDB: %s", lname);
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
		
		    read_seq_in_list (lname,&nseq,&sequences,&seq_name);             
		    S1=fill_sequence_struc ( nseq, sequences, seq_name);
		  
		  for ( b=0; b< S1->nseq; b++)sprintf ( S1->file[b], "%s", lname);
		  nseq=0;free_char (sequences, -1); free_char ( seq_name, -1);
		  sequences=NULL;
		  seq_name=NULL;
		  S1=seq2unique_name_seq (S1);

                  if ((S=merge_seq( S1, S))==NULL){fprintf ( stderr, "\nSequence Error in %s [FATAL:%s]\n",lname,PROGRAM); myexit(EXIT_FAILURE);} 
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
char * expand_constraint_list_file ( char *file)
{
  char *new_file;
  FILE *IN, *OUT;
  int a, b, c, n;
  char **list;
  static char *buf;

  if ( !grep_function ( "'+BLOCK+'", file))return file;

  new_file=vtmpnam (NULL);
  IN=vfopen ( file,"r");
  OUT=vfopen (new_file, "w");
 
  while ( (c=fgetc (IN))!=EOF)
    {
      ungetc (c, IN);
      buf=vfgets (buf, IN);
      if ( !strstr (buf, "+BLOCK+"))
	fprintf (OUT, "%s", buf);
      else
	{
	  list=string2list (buf);
	  n=atoi (list[2]);
	  
	  for (a=0; a< n; a++)
	    {
	      fprintf ( OUT, "%5d %5d ",atoi(list[3])+a, atoi(list[4])+a);
	      for (b=5; b<atoi(list[0]); b++) fprintf ( OUT, "%s ", list[b]);
	      fprintf (OUT, "\n");
	    }
	  free_char (list, -1);
	}
    }
  vfclose (IN); vfclose (OUT);
  return new_file;
}
	  

Constraint_list * make_test_lib(Constraint_list *CL)
{
  int a, b, c, l1;
  Sequence *S;
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
	      vwrite_clist( CL,CL->ne,0,a);
	      vwrite_clist( CL,CL->ne,1,b);
	      vwrite_clist( CL,CL->ne,2,c+1);
	      vwrite_clist( CL,CL->ne,3,c+1);
	      vwrite_clist( CL,CL->ne,4,100);
	      CL->ne++;
	    }
	}
    }
  return CL;
}
Constraint_list * read_constraint_list_file(Constraint_list *CL, char *fname)
        {
	int a, c,e,n,z;
	int seq_len, sn;
	int s1, s2;
	FILE *fp;
	static char *name;
	char *sequence;
	static char *mat;
	static char *dp_mode;
	int max_nseq=0;
	static int *sn_list;
	static int line=2;
	int list_nseq;
	static CLIST_TYPE *entry;
	Sequence *S;
	Sequence *small_S;
	int seq_1_to_n=0;
	Alignment *B=NULL;
	char *buf;
	int lline;
	char *stripped_file1;
	char *stripped_file;



	stripped_file1=strip_file_from_comments ("!", fname);
	stripped_file=expand_constraint_list_file (stripped_file1);
	small_S=read_seq_in_n_list (&fname, 1,NULL, NULL); 
	
	if ( !CL)
	  {
	    CL=declare_constraint_list ( small_S,NULL, NULL, 0,NULL, NULL);  
	    CL->S=small_S;
	  }

	small_S=read_seq_in_n_list (&fname, 1, (CL->S)->type, NULL);
	
	B=seq2aln (small_S, NULL, 1); 
	B=fix_aln_seq  ( B, (CL->S));
	
	if ( CL->S!=small_S)free_sequence (small_S, B->nseq);
		
	lline=measure_longest_line_in_file (fname)+1;
		
	if ( !mat) mat=vcalloc (STRING, sizeof (char));
	if ( !dp_mode) dp_mode=vcalloc (STRING, sizeof (char));
	fp=vfopen (fname, "r");
	while((c=fgetc(fp))!='#' && c!=EOF)if ( c=='\n')max_nseq++;
	vfclose (fp);
	
	
	buf=vcalloc (lline, sizeof (char));	
	sequence=vcalloc (lline, sizeof (char));
	if ( !name)name=vcalloc ( 100, sizeof (char));
	if ( !entry)entry=vcalloc ( CL->entry_len, CL->el_size);
	if ( !sn_list)sn_list=vcalloc (max_nseq, sizeof (int));
	else 
	  {
	    sn_list=vrealloc (sn_list, max_nseq*sizeof (int));
	  }
	S=CL->S;

	seq_1_to_n=((fp=find_token_in_file (fname, NULL, "SEQ_1_TO_N"))!=NULL);
	vfclose (fp);
	if ( sn_list==NULL)sn_list=vcalloc (max_nseq, sizeof (int));
	
	/*Read Constraint list*/
	fp=vfopen(stripped_file,"r");
	fscanf ( fp, "%d\n", &list_nseq);
	for ( a=0; a<list_nseq; a++)
		{
		fscanf ( fp, "%s %d %s\n", name, &seq_len, sequence);
		line++;
		lower_string (sequence);
		
		if ((sn=name_is_in_list (name,S->name, S->nseq, 100))==-1){continue;}
		else
		  {
		    sn_list[a]=sn;
		  }
		}
	c=0;
	
	while (c!=EOF && (c=fgetc(fp))!=EOF)
	  {
	    
	    ungetc(c, fp);
	    
	    if ( c=='#')
			{
			  fscanf ( fp, "#%d %d\n", &s1, &s2);line++;
			/*Check If the sequence numbering is legal*/
			if ( seq_1_to_n){s1--; s2--;}
			
			if (s1<0 || s2 <0)
			  {
			    
			    myexit(fprintf_error (stderr, "ERROR: Wrong Sequence Numbering in %s\n",fname));
			  }
			

			
			s1=sn_list[s1];
			s2=sn_list[s2];

			while (isdigit((c=fgetc(fp))))
				{

				for ( z=0; z<  CL->entry_len; z++)entry[z]=0;
				ungetc(c, fp);				
				n=0;
				entry[n++]=s1;
				entry[n++]=s2;
				while ( (c=fgetc(fp))!='\n')
				  {
				   
					if ( isspace (c));
					else 
						{
						ungetc(c, fp);
						fscanf ( fp, "%d", &entry[n]);
						n++;
						}
					
					if ( n>CL->entry_len)
						{
						add_warning ( stderr, "\nWARNING:PARSING ERROR #1 (Too many Fields) IN %s AT LINE %d: C=%c n=%d\n", fname,line, c,n);
						for ( e=2; e<LIST_N_FIELDS; e++)
							fprintf ( stderr, "%d ", entry[e]);
					
						myexit (EXIT_FAILURE);
						}
					}
				if (c=='\n')line++;
				 
				if ( n<=CONS)entry[CONS]=1;
				
				  
				/*Check The legality of the entry*/
				if ( n>0 && n<3)
					{
					add_warning ( stderr, "\nWARNING:PARSING ERROR #2 IN %s (Not enough Fields) AT LINE %d: C=%c\n", fname,line-1, c);
					for ( e=2; e<LIST_N_FIELDS; e++)
						fprintf ( stderr, "%d ",entry[e]);
					
					myexit (EXIT_FAILURE);
					}
				
				entry[R1]=(B->seq_cache)?B->seq_cache[entry[SEQ1]][entry[R1]]:entry[R1];
				entry[R2]=(B->seq_cache)?B->seq_cache[entry[SEQ2]][entry[R2]]:entry[R2];
				
				if ( entry[R1] && entry[R2])
				  {
				    if ( entry[R1]<=0 || entry[R1]>(CL->S)->len[s1])
				      {
					fprintf ( stderr, "\nERROR: Seq1=%d (len=%d, name=%s), Seq2=%d (len=%d, name=%s), Res1 %d, Res2 %d\n", entry[SEQ1]+1,(CL->S)->len[s1],(CL->S)->name[s1], entry[SEQ2]+1,(CL->S)->len[s2],(CL->S)->name[s2],entry[R1], entry[R2]);
					fprintf ( stderr, "\nERROR: Library %s, line %d, Field 1: Bad residue numbering (%d)[FATAL:%s]\n", fname, line-1,entry[R1], PROGRAM);
					myexit (EXIT_FAILURE);
				      }
				    else if (entry[R2]<=0 || entry[R2]>(CL->S)->len[s2])
				      {
					fprintf ( stderr, "\nERROR: Seq1=%d (len=%d, name=%s), Seq2=%d (len=%d, name=%s), Res1 %d, Res2 %d\n", entry[SEQ1]+1,(CL->S)->len[s1],(CL->S)->name[s1], entry[SEQ2]+1,(CL->S)->len[s2],(CL->S)->name[s2],entry[R1], entry[R2]);
					
					fprintf ( stderr, "\nERROR: Seq1: %d, Seq2 %d, Res1 %d, Res2 %d\n", entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]);
					fprintf ( stderr, "\nERROR: Library %s, line %d, Field 2: Bad residue numbering (%d)[FATAL:%s]\n", fname, line-1, entry[R2],PROGRAM);
					myexit (EXIT_FAILURE);
				      }
				    fscanf ( fp, "\n");
				    if ( (entry[SEQ1]>entry[SEQ2])|| (entry[SEQ1]==entry[SEQ2] && entry[R1]>entry[R2]))
				      {
					SWAP(entry[SEQ1],entry[SEQ2]);
					SWAP(entry[R1], entry[R2]);
				      }
				  
				    for ( z=0; z< CL->entry_len; z++)vwrite_clist( CL,CL->ne, z, entry[z]);
				    
				    CL->ne++;
				  }
				}
			 ungetc ( c, fp);
			 
			 }
	    else if ( c=='!' || c=='C' || c=='\n' || c=='\r'){while ((c=fgetc(fp))!='\n' && c!=EOF && c!='\r');}
	    else
	      {
		fprintf ( stderr, "\n\n PARSING ERROR 3 IN %s AT LINE %d: [%c] \n[read_constraint_list_file]", fname,line,c); 			
		while ((c=fgetc(fp))!='\n' && c!=EOF)fprintf ( stderr, "%c", c);
		fprintf ( stderr, "\n");
		printf_system ( "cp %s faulty_library.tc_lib", fname);
		myexit (EXIT_FAILURE);
	      }
	    if ( c==EOF)ungetc(c, fp);
	  }
	
	free_aln (B);
	vfree(buf);
	vfree(sequence);
	vfclose (fp);	 
	remove(stripped_file);
	
	return CL;
	} 

Constraint_list * read_constraint_list_raw_file(Constraint_list *CL, char *fname)
{
  FILE *fp;
  int n=0, v;
  static int *entry;

  if ( !entry)entry=vcalloc (CL->entry_len, sizeof (int));
  fp=vfopen (fname, "r");
  
  while (fscanf (fp, "%d ", &v)==1)
    {
      if (n==CL->entry_len)
	{
	  add_entry2list (entry, CL);
	  n=0;
	}
      entry[n++]=v;
    }
   if (n==CL->entry_len)
	{
	  add_entry2list (entry, CL);
	  n=0;
	}
   vfclose (fp);
   return CL;
}
	   
Constraint_list * fast_read_constraint_list_file(Constraint_list *CL, char *in_fname)
        {
	  Sequence *NS;
	  int **index;
	  int *list=NULL;
	  int c;
	  FILE *fp;
	  CLIST_TYPE *entry=NULL;
	  char *buf=NULL, *buf2;
	  char *fname;
	  int i;
	  

	  
	  fname=expand_constraint_list_file (in_fname);
	  
	  if (!CL)
	    {
	      return read_constraint_list_file (CL,fname);
	    }
	  
	  entry=vcalloc (sizeof (int), CL->entry_len);
	  NS=read_seq_in_n_list (&fname, 1,NULL, NULL);
	  index=index_seq_name(NS,CL->S);
	  

	  /*Read Constraint list*/
	  fp=vfopen(fname,"r");
	  i=0;
	  while (i<=NS->nseq)
	    {
	      buf=vfgets ( buf, fp);
	      if (buf[0]!='!')i++;
	    }
	  
	  while ( (c=fgetc (fp))!='#' && c!=EOF);
	  if (c==EOF)
	    {
	      vfclose (fp);
	      vfree (list);
	      free_sequence (NS,-1);
	      vfree (entry);
	      vfree (buf);
	      add_warning (stderr, "Warning: incomplete library [%s]",PROGRAM); 
	      return CL;
	    }
	  ungetc (c, fp);
	  
	  

	  while ((buf2=vfgets ( buf, fp))!=NULL)
	    {
	      if (buf2[0]=='!')continue;
	      buf=buf2;
	      
	      list=string2num_list2 (buf, " #\n");
	      
	      if (buf[0]=='#')
		{
		  sscanf ( buf, "#%d %d", &entry[SEQ1], &entry[SEQ2]);
		  entry[SEQ1]=index[entry[SEQ1]-1][0];
		  entry[SEQ2]=index[entry[SEQ2]-1][0];
		  
		}
	      else
		{
		  sscanf (buf, "%d %d %d %d %d", &entry[R1], &entry[R2], &entry[WE], &entry[CONS], &entry[MISC]);
		  
		  
		  
		  
		  
		  
		  if ( (entry[SEQ1]>entry[SEQ2])|| (entry[SEQ1]==entry[SEQ2] && entry[R1]>entry[R2]))
		    {
		      SWAP(entry[SEQ1],entry[SEQ2]);SWAP(entry[R1], entry[R2]);
		    }
		  add_entry2list (entry, CL);
		}
	    }
	  vfclose (fp);
	  
	  
	  
	  vfree (list);
	  free_sequence (NS,-1);
	  vfree (entry);
	  vfree (buf);
	  
	  return CL;
	} 	
  
int read_seq_in_list ( char *fname,  int *nseq, char ***sequences, char ***seq_name)
	{
	  int a;
	int seq_len, sn;

	FILE *fp;
	char name[1000];
	char *sequence;	
	static int max_nseq;
	static int *sn_list;
	int list_nseq;
	int lline;

	fp=vfopen (fname, "r");
	fp=skip_commentary_line_in_file ('!', fp);
	fscanf (fp, "%d\n", &max_nseq);
	for ( lline=0,a=0; a<max_nseq; a++)
	  {
	    int l=0;
	    fscanf (fp, "%*s %d %*s\n", &l);
	    lline=MAX(lline, l);
	  }
	vfclose (fp);
	sequence=vcalloc (lline+1, sizeof (char));
	
	if ( seq_name[0]==NULL)
		{
		seq_name[0]= declare_char (max_nseq,0);
		sequences[0]=declare_char (max_nseq,0);		
		}
	if ( sn_list==NULL)sn_list=vcalloc ( max_nseq, sizeof (int));
	else sn_list=vrealloc (sn_list, max_nseq*sizeof (int));

	fp=vfopen (fname,"r");
	fp=skip_commentary_line_in_file ('!', fp);
	if (fscanf ( fp, "%d\n", &list_nseq)!=1)return 0;
	for ( a=0; a<max_nseq; a++)
	  {
	    fp=skip_commentary_line_in_file ('!', fp);
	    fscanf ( fp, "%s %d %s\n", name, &seq_len, sequence);
	    lower_string (sequence);
	    if ((sn=name_is_in_list (name, seq_name[0], nseq[0], 100))==-1)
	      {
		seq_name[0][nseq[0]]=vcalloc (strlen (name)+1, sizeof (char));
		sprintf (seq_name[0][nseq[0]], "%s", name);
		sequences[0][nseq[0]]=vcalloc (strlen (sequence)+1, sizeof (char));
		sprintf (sequences[0][nseq[0]], "%s", sequence);
		sn_list[a]=nseq[0];
		nseq[0]++;
	      }
	    else
	      {
		sn_list[a]=sn;
	      }
	  }
	vfclose (fp);
	vfree (sequence);
	return 1;
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

FILE * save_constraint_list ( Constraint_list *CL,int start, int len, char *fname, FILE *fp,char *mode, Sequence *S)
        {
	int a, b;
	static int* translation;

	
	
	
	if ( fp==NULL)
	   {
	   if ( translation!=NULL)vfree(translation);
	   translation=vcalloc ( (CL->S)->nseq+1, sizeof (int));
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
	else if (strm(mode, "binary"))
	   {
	   if ( fp==NULL)fp=vfopen ( fname, "wb");
	   fp=save_constraint_list_bin  (fp, CL, 0, CL->ne, translation);
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
int constraint_list2raw_file ( Constraint_list *CL, char *fname, char *mode)
{
  FILE *fp;
  if ( !CL || !CL->ne || !fname){return 0;}
  
  fp=vfopen (fname,mode);
  fp=save_raw_constraint_list (fp, CL, 0, CL->ne, NULL);
  vfclose (fp);
  return CL->ne;
  }
FILE * save_raw_constraint_list   ( FILE *fp,Constraint_list *CL, int start,int len, int *translation)
{
  int a, b;
  for ( b=0; b<CL->entry_len*CL->ne; b++)
    {
      fprintf ( fp, "%d ", CL->L[b]);
    }

  fprintf (fp, "\n");
  return fp;
}
FILE * save_constraint_list_ascii ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation)
	{	
	int a, b, s1, s2;
        CLIST_TYPE x1, x2;

	if (len==start && CL->cpu!=-1)
	    {
	    fprintf (OUT, "! CPU %d\n",get_time());
	    return OUT;
	    }
	else
	    {
	    
	    s1=translation[vread_clist(CL,start,SEQ1)];
	    s2=translation[vread_clist(CL,start,SEQ2)];
	   
	
	    if ( s1!=-1 && s2!=-1)fprintf ( OUT, "#%d %d\n", s1+1, s2+1);
	    for ( a=start; a<(len+start); a++)
		   {
		   x1=translation[vread_clist(CL,a,SEQ1)];
		   x2=translation[vread_clist(CL,a,SEQ2)];
		   if ( x1==-1 || x2==-1);
		   else 
		       {
		       if ( x1!=s1 || x2!=s2)
			  {
			  s1=x1;
			  s2=x2;
			  fprintf ( OUT, "#%d %d\n", s1+1, s2+1);
			  }
		       for ( b=2; b<CL->entry_len; b++) fprintf ( OUT, "%5d ", vread_clist(CL, a, b));
		       fprintf (OUT, "\n");
		       }
		   }
	    }
	return save_list_footer (OUT, CL);
			
	}
FILE * save_constraint_list_bin ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation)
	{	
	int a, b;

	CLIST_TYPE x1, x2;
	

	if (len==start && CL->cpu!=-1)
	    {
	    
	    return OUT;
	    }
	else
	    {
	    for ( a=start; a<(len+start); a++)
		   {
		   x1=translation[vread_clist(CL,a,SEQ1)];
		   x2=translation[vread_clist(CL,a,SEQ2)];
		   if ( x1==-1 || x2==-1);
		   else 
		       {
		       for ( b=2; b<CL->entry_len; b++)
		           {
			   x1=vread_clist(CL,a,b);
			   fwrite (&x1, CL->el_size, 1, OUT);
			   }
		       }
		   }
	    }
	return OUT;			
	}

/*********************************************************************/
/*                                                                   */
/*                         LIST CONVERTION                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/		
Constraint_list * filter_constraint_list (Constraint_list *CL, int field, int T)
{
  int a,b,c;
  if (!CL || !CL->L)return CL;
  for (a=0, b=0; a<CL->ne; a++)
    {
      if (CL->L[a*CL->entry_len+field]<T);
      else
	{
	  for (c=0; c<CL->entry_len; c++)CL->L[b*CL->entry_len+c]=CL->L[a*CL->entry_len+c];
	  b++;
	}
    }
  CL->ne=b;
  CL->residue_indexed=0;
  return CL;
}
int constraint_list_is_connected ( Constraint_list *CL)
{
  int a, b, c, s1, s2;
  int *connexions;
  int rv =1;
  Sequence*S;

  if (!CL->ne) return 1;
  
  S=CL->S;
  
  connexions=vcalloc ((CL->S)->nseq+1, sizeof (int));
  for ( a=0; a< CL->ne; a++)
    {
      s1=CL->L[a*CL->entry_len+SEQ1];
      s2=CL->L[a*CL->entry_len+SEQ2];
      connexions[s1]++;
      connexions[s2]++;
    }
  for (a=0; a<S->nseq; a++)
    {
      if (!connexions[a])
	{
	  add_warning ( stderr, "ERROR: Sequence %s is not connected\n", (CL->S)->name[a]);
	  rv=0;
	}
    }
  return rv;
}
Constraint_list * nfork_relax_constraint_list (Constraint_list *CL);
Constraint_list * fork_relax_constraint_list (Constraint_list *CL);
Constraint_list * relax_constraint_list (Constraint_list *CL)
{
  if ( get_nproc()==1)return  nfork_relax_constraint_list (CL);
  else if (strstr ( CL->multi_thread, "relax"))return fork_relax_constraint_list (CL);
  else return nfork_relax_constraint_list (CL);
}

Constraint_list * fork_relax_constraint_list (Constraint_list *CL)
{
  int a, s1, s2, r1, r2,n;
  int score;
  int thr;
  int chunk, npid, job,pid;
  FILE *fp;
  char **pid_tmpfile;
  int * pid_list;
  int in;

  in=CL->ne;
  if (!CL || !CL->L)return CL;

  fprintf ( CL->local_stderr, "\nLibrary Relaxation: Multi_proc [%d] ", get_nproc());
  
  if ((chunk=CL->ne/get_nproc())==0)chunk=get_nproc();
  
 
  pid_tmpfile=vcalloc ((CL->ne/chunk)+1, sizeof (char*));
  pid_list   =vcalloc (MAX_N_PID, sizeof (int *));
  
  for (npid=0,job=0; job<CL->ne; job+=chunk)
    {
      pid_tmpfile[npid]=vtmpnam(NULL);
      pid=vvfork (NULL);
      if (pid==0)
	{
	  int s,e;
	  
	  initiate_vtmpnam (NULL);
	  s=job;
	  e=MIN((s+chunk),CL->ne);
	  fp=vfopen (pid_tmpfile[npid], "w");
	  for (a=s; a<e; a++)
	    {
	      if (job==0)output_completion (CL->local_stderr,a,chunk,1, "Submit   Job");
	      s1=CL->L[a*CL->entry_len+SEQ1];
	      s2=CL->L[a*CL->entry_len+SEQ2];
	      
	      r1=CL->L[a*CL->entry_len+R1];
	      r2=CL->L[a*CL->entry_len+R2];
	      score=residue_pair_extended_list_pc (CL,s1, r1,s2, r2);
	      CL->L[a*CL->entry_len+WE]=score;
	      fprintf (fp, "%d %d ", a, score);
	    }
	  vfclose (fp);
	  myexit (EXIT_SUCCESS);
	}
      else
	{
	  pid_list[pid]=npid;
	  //set_pid (pid);
	  npid++;
	}
    }
  
  for (a=0; a<npid; a++)
    {
      int i;
      int j=0;
      
      pid=vwait (NULL);
      fp=vfopen (pid_tmpfile[pid_list[pid]], "r");
      while (fscanf (fp, "%d %d ",&i, &score)==2){CL->L[i*CL->entry_len+WE]=score;j++;}
      vfclose (fp);
      remove(pid_tmpfile[pid_list[pid]]);
    }
  
  vfree (pid_list);
  vfree (pid_tmpfile);
  
  thr=10;
    
  for (n=0,a=0; a< CL->ne; a++)
    {
      score=CL->L[a*CL->entry_len+WE];
      
      if (score<=thr);
      else
	{
	  CL->L[n*CL->entry_len+SEQ1]=CL->L[a*CL->entry_len+SEQ1];
	  CL->L[n*CL->entry_len+SEQ2]=CL->L[a*CL->entry_len+SEQ2];
	  CL->L[n*CL->entry_len+R1]=CL->L[a*CL->entry_len+R1];
	  CL->L[n*CL->entry_len+R2]=CL->L[a*CL->entry_len+R2];
	  CL->L[n*CL->entry_len+WE]=score;
	  n++;
	}
      
    }
 
  CL->L=vrealloc (CL->L, n*CL->entry_len*sizeof (int));
  CL->ne=n;
  
  fprintf ( CL->local_stderr, "\nTotal Relaxation: [%d]--->[%d] Entries\n",in, CL->ne); 
  CL->residue_indexed=0;
  return CL;
  
}		
  
Constraint_list * nfork_relax_constraint_list (Constraint_list *CL)
{
  int a, s1, s2, r1, r2;
  int max, score, n;
  int thr;
  
  if (!CL || !CL->L)return CL;

  fprintf ( CL->local_stderr, "\nLibrary Relaxation:[%d] ", CL->ne);
  for (max=0,a=0; a<CL->ne; a++)
    {
      
      s1=CL->L[a*CL->entry_len+SEQ1];
      s2=CL->L[a*CL->entry_len+SEQ2];

      r1=CL->L[a*CL->entry_len+R1];
      r2=CL->L[a*CL->entry_len+R2];
      
      score=residue_pair_extended_list_pc (CL,s1, r1,s2, r2);
      //HERE ("%d %d",  CL->L[a*CL->entry_len+WE],score);
      CL->L[a*CL->entry_len+WE]=score;
       
    }

  thr=10;
    
  for (n=0,a=0; a< CL->ne; a++)
    {
      score=CL->L[a*CL->entry_len+WE];
      
      if (score<=thr);
      else
	{
	  CL->L[n*CL->entry_len+SEQ1]=CL->L[a*CL->entry_len+SEQ1];
	  CL->L[n*CL->entry_len+SEQ2]=CL->L[a*CL->entry_len+SEQ2];
	  CL->L[n*CL->entry_len+R1]=CL->L[a*CL->entry_len+R1];
	  CL->L[n*CL->entry_len+R2]=CL->L[a*CL->entry_len+R2];
	  CL->L[n*CL->entry_len+WE]=score;
	  n++;
	}
      
    }
  CL->L=vrealloc (CL->L, n*CL->entry_len*sizeof (int));
  CL->ne=n;
  
  fprintf ( CL->local_stderr, "--->[%d]\n", CL->ne); 
  CL->residue_indexed=0;
  return CL;
  
}		
// relax constraint list for gene prediction

Constraint_list * expand_constraint_list_4gp (Constraint_list *CL, int T)
{
  int *L;
  int chunck=500;
  int max=500;
  int n=0;
  int s1, s2, r1, r2,a,b,c,d,w;
  Sequence *S;

  CL=index_res_constraint_list (CL,WE);
  S=CL->S;
  L=vcalloc (max*CL->entry_len, sizeof (int));
  for (a=0; a<S->nseq; a++)//loop sequences
    {
      for (b=0; b<S->len[a];b++)
	{
	  for (c=1; c<CL->residue_index[a][b][0]-3;c+=3)
	    {
	      
	      
	      s2=a;
	      r2=b;
	      s1=CL->residue_index[a][b][c];;
	      r1=CL->residue_index[a][b][c+1];
	      w=residue_pair_extended_list_4gp (CL,s1, r1,s2, r2);
	      if (w>T)
		{
		  if ( n>=max){max+=chunck; L=vrealloc (L, max*CL->entry_len*sizeof (int));}
		  L[n*CL->entry_len+SEQ2]=s2;
		  L[n*CL->entry_len+R2  ]=r2;
		  L[n*CL->entry_len+SEQ1]=s1;
		  L[n*CL->entry_len+R1  ]=r1;
		  L[n*CL->entry_len+WE  ]=w;
		  n++;
		}
	      for (d=c+3; d<CL->residue_index[a][b][0]; d+=3)
		{
		  s2=CL->residue_index[a][b][d];
		  r2=CL->residue_index[a][b][d+1];
		  w=residue_pair_extended_list_4gp (CL,s1, r1,s2, r2);
		  if (w>T)
		    {
		      if ( n>=max){max+=chunck; L=vrealloc (L, max*CL->entry_len*sizeof (int));}
		      L[n*CL->entry_len+SEQ1]=s1;
		      L[n*CL->entry_len+SEQ2]=s2;
		      L[n*CL->entry_len+R1]=r1;
		      L[n*CL->entry_len+R2]=r2;
		      L[n*CL->entry_len+WE]=w;
		      n++;
		    }
		}
	    }
	}
    }
  vfree (CL->L);
  CL->L=L;
  CL->residue_indexed=0;
  CL->ne=n;

  CL=compact_list (CL, 0,n, "best");
  CL->residue_indexed=0;
  
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
  int a, s1, s2, r1, r2,n;
  int score;
  int thr;
  int chunk, npid, job,pid;
  FILE *fp;
  char **pid_tmpfile;
  int * pid_list;
  int in;

  in=CL->ne;
  if (!CL || !CL->L)return CL;

  
  if ((chunk=CL->ne/get_nproc())==0)chunk=get_nproc();
  
 
  pid_tmpfile=vcalloc ((CL->ne/chunk)+1, sizeof (char*));
  pid_list   =vcalloc (MAX_N_PID, sizeof (int *));
  
  for (npid=0,job=0; job<CL->ne; job+=chunk)
    {
      pid_tmpfile[npid]=vtmpnam(NULL);
      pid=vvfork (NULL);
      if (pid==0)
	{
	  int s,e;
	  
	  initiate_vtmpnam (NULL);
	  s=job;
	  e=MIN((s+chunk),CL->ne);
	  fp=vfopen (pid_tmpfile[npid], "w");
	  for (a=s; a<e; a++)
	    {
	     
	      s1=CL->L[a*CL->entry_len+SEQ1];
	      s2=CL->L[a*CL->entry_len+SEQ2];
	      
	      r1=CL->L[a*CL->entry_len+R1];
	      r2=CL->L[a*CL->entry_len+R2];
	      score=residue_pair_extended_list_4gp (CL,s1, r1,s2, r2);
	      CL->L[a*CL->entry_len+WE]=score;
	      fprintf (fp, "%d %d ", a, score);
	    }
	  vfclose (fp);
	  myexit (EXIT_SUCCESS);
	}
      else
	{
	  pid_list[pid]=npid;
	  //set_pid (pid);
	  npid++;
	}
    }
  
  for (a=0; a<npid; a++)
    {
      int i;
      int j=0;
      
      pid=vwait (NULL);
      fp=vfopen (pid_tmpfile[pid_list[pid]], "r");
      while (fscanf (fp, "%d %d ",&i, &score)==2){CL->L[i*CL->entry_len+WE]=score;j++;}
      vfclose (fp);
      remove(pid_tmpfile[pid_list[pid]]);
    }
  
  vfree (pid_list);
  vfree (pid_tmpfile);
  
  thr=0;
    
  for (n=0,a=0; a< CL->ne; a++)
    {
      score=CL->L[a*CL->entry_len+WE];
      
      if (score<=thr);
      else
	{
	  CL->L[n*CL->entry_len+SEQ1]=CL->L[a*CL->entry_len+SEQ1];
	  CL->L[n*CL->entry_len+SEQ2]=CL->L[a*CL->entry_len+SEQ2];
	  CL->L[n*CL->entry_len+R1]=CL->L[a*CL->entry_len+R1];
	  CL->L[n*CL->entry_len+R2]=CL->L[a*CL->entry_len+R2];
	  CL->L[n*CL->entry_len+WE]=score;
	  n++;
	}
      
    }
 
  CL->L=vrealloc (CL->L, n*CL->entry_len*sizeof (int));
  CL->ne=n;
  
  CL->residue_indexed=0;
  return CL;
  
}		
  
Constraint_list * nfork_relax_constraint_list_4gp (Constraint_list *CL)
{
  int a, s1, s2, r1, r2;
  int max, score, n;
  int thr;
  int tot;
  

  
 
  if (!CL || !CL->L)return CL;
  
 
   
  for (max=0,a=0; a<CL->ne; a++)
    {
      
      s1=CL->L[a*CL->entry_len+SEQ1];
      s2=CL->L[a*CL->entry_len+SEQ2];

      r1=CL->L[a*CL->entry_len+R1];
      r2=CL->L[a*CL->entry_len+R2];
      
      score=residue_pair_extended_list_4gp (CL,s1, r1,s2, r2);
      
      //HERE ("%d %d",  CL->L[a*CL->entry_len+WE],score);
      CL->L[a*CL->entry_len+WE]=score;
       
    }

  thr=0;
    
  for (tot=0,n=0,a=0; a< CL->ne; a++)
    {
      score=CL->L[a*CL->entry_len+WE];
      if (score<=thr);
      else
	{
	  CL->L[n*CL->entry_len+SEQ1]=CL->L[a*CL->entry_len+SEQ1];
	  CL->L[n*CL->entry_len+SEQ2]=CL->L[a*CL->entry_len+SEQ2];
	  CL->L[n*CL->entry_len+R1]=CL->L[a*CL->entry_len+R1];
	  CL->L[n*CL->entry_len+R2]=CL->L[a*CL->entry_len+R2];
	  CL->L[n*CL->entry_len+WE]=score;
	  n++;
	}
      
    }
  CL->L=vrealloc (CL->L, n*CL->entry_len*sizeof (int));
  CL->ne=n;

  CL->residue_indexed=0;
  return CL;
  
}		
int constraint_list2avg ( Constraint_list *CL)
{
  int a;
  long tot=0;
  for (a=0; a<CL->ne; a++)tot+= CL->L[a*CL->entry_len+WE];
  return (tot)/(CL->ne);
}
int constraint_list2fraction_covered ( Constraint_list *CL)
{
  int a;
  long tot=0, tot2=0;
  int **w;
  w=declare_int ((CL->S)->nseq+2, (CL->S)->max_len+10);
  for (a=0; a<(CL->S)->nseq; a++)tot+=strlen ((CL->S)->seq[a]);
  
  for (a=0; a<CL->ne; a++)
    {
      int s1, s2,r1, r2;
      s1=CL->L[a*CL->entry_len+SEQ1];
      s2=CL->L[a*CL->entry_len+SEQ2];
      r1=CL->L[a*CL->entry_len+R1];
      r2=CL->L[a*CL->entry_len+R2];
      
      if (!w[s1][r1]){tot2++; w[s1][r1]=1;}
      if (!w[s2][r2]){tot2++; w[s2][r2]=1;}
    }
  free_int (w,-1);
  tot=(tot2*100)/tot;
  tot=tot*constraint_list2avg (CL);
  return tot;
}


Constraint_list * shrink_constraint_list (Constraint_list *CL)
{
  int a, b, n, tot;
  Constraint_list *CL2;
  Alignment *A, *B;
  int *ns, **ls;

  ns=vcalloc (2, sizeof (int));
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
  vfree(CL->L);
  CL->L=CL2->L;
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
  
  
  if ( !weight)weight=vcalloc (MAX(2,A->len_aln), sizeof (int));
  
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
  else if ( strncmp ( weight_mode, "len",3)==0)
    {
      weight[1]=A->len_aln;
    }
  else if ( strnm ( weight_mode, "sim", 3) || strm (weight_mode, "default"))
    {
      
      ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-", (strm (weight_mode, "default"))?NULL:(weight_mode+3));
	  if (ref_weight == 0)
		  ref_weight = 1;
      weight[1]=ref_weight;

    }
  else if ( strnm ( weight_mode, "subset", 6))
    {
      ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-",NULL);
      weight[1]=ref_weight;
    }
  
  else if ( strncmp (weight_mode, "winsim", 6)==0)
    {
      weight=get_seq_winsim ( A->seq_al[s1], A->seq_al[s2], "-", weight_mode+6, weight);
    }
  else if (  strncmp ( weight_mode, "cdna", 4)==0)
    {
      ref_weight=get_seq_sim ( A->seq_al[s1], A->seq_al[s2], "-", weight_mode+4);
      col=vcalloc ( A->len_aln+1, sizeof (int));
      if (A->cdna_cache)
	for ( a=0; a<=A->len_aln; a++)col[a]=A->cdna_cache[0][a];
      else
	for ( a=0; a<=A->len_aln; a++)col[a]=1;
      for ( c=0; c< A->len_aln; c++)weight[c]=ref_weight*col[c];
      vfree (col);
    }
  else if ( strm ( weight_mode, "pdb"))
    {
      if ( !(A->CL) || !(A->CL)->T)
	{
	  myexit(fprintf ( stderr, "\nCould not find the PDB structure"));
	}
    }
  else if ( strm (weight_mode, "overaln"))
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




Constraint_list *aln2constraint_list      (Alignment *A, Constraint_list *CL,char *in_weight_mode)
	{
	Constraint_list *CLB=NULL;
	int a, b, c,nres1, nres2;
	int *weight=NULL;
	int s1, s2;
	int fixed_nres1, fixed_nres2;
	int do_pdb=0;
	int pdb_weight=0;
	int set_misc;
	char*alp=NULL;
	char *p, *s;
	char weight_mode [100];
	int *top_seq, *sindex;
	int use_top=0;
	
	
	sprintf ( weight_mode , "%s", (!in_weight_mode || strm (in_weight_mode, "default"))?"sim":in_weight_mode);
	
	if ( !A)
	  return CL;	
	
	if ( !CL)
	  {
	    Sequence *S;
	    S=aln2seq (A);
	    CL=declare_constraint_list (S,NULL, NULL, 0,NULL, NULL);  
	    CL->S=S;
	  }
	CLB=(Constraint_list *)A->CL;
	
	
	do_pdb=(strstr ( weight_mode, "pdb"))?1:0;
	if ( (p=strstr (weight_mode, "_subset_")))
	  {
	    alp=strchr (weight_mode, '_')+1;
	    p[0]='\0';
	  }
	
	A=reorder_aln (A, (CL->S)->name, (CL->S)->nseq);
	sindex=vcalloc (A->nseq, sizeof(int));
	for (a=0; a<A->nseq; a++)sindex[a]=name_is_in_list (A->name[a], (CL->S)->name, (CL->S)->nseq, 100);
	
	
	top_seq=vcalloc (A->len_aln, sizeof (int));
	for (c=0; c<A->len_aln; c++)
	  {
	    for (a=0; a<A->nseq; a++)
	      {
		if (sindex[a]!=-1 && !is_gap(A->seq_al[sindex[a]][c])){top_seq[c]=sindex[a];a=A->nseq;}
	      }
	  }
	
	for ( a=0; a<A->nseq-1; a++)
	  {
	    for (set_misc=0,b=a+1; b< A->nseq; b++)
	      {	
		s1=sindex[a];
		s2=sindex[b];
		if ( s1==-1 || s2==-1)
		  {
		    if ( getenv4debug ("DEBUG_LIBRARY"))
		      fprintf ( stderr, "\n[DEBUG_LIBRARY:aln2constraint_list]Could use a pair of constraints");
		  }
		else if ( s1!=-1 && s2!=-1)
		  {
		    int use_pair;
		    
		    weight=seqpair2weight (a, b, A, CL, weight_mode, weight);
		    
		    for (nres1=A->order[a][1], nres2=A->order[b][1], c=0; c< A->len_aln; c++)
		      {
			int isgop1, isgop2;
			if (use_top && top_seq[c]!=a)continue; //make nr dataset for MSAs
			
			isgop1=is_gop(c, A->seq_al[a]);
			isgop2=is_gop(c, A->seq_al[b]);
			nres1+=!is_gap(A->seq_al[a][c]);
			nres2+=!is_gap(A->seq_al[b][c]);
			
			if ( strm ( weight_mode, "pdb") && CLB)
			  {
			    
			    pdb_weight=MAX(0,(CLB->evaluate_residue_pair)(CLB,0, nres1,1,nres2));
			  }
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
			    
			    fixed_nres1=(!A->seq_cache)?nres1:A->seq_cache[s1][nres1];
			    fixed_nres2=(!A->seq_cache)?nres2:A->seq_cache[s2][nres2];
			    
			    if ( fixed_nres1==-1 || fixed_nres2==-1)
			      {
				fprintf ( stderr, "\nPB: Sequence %s, Residue %d : Cache=%d",A->name[a], nres1,fixed_nres1 );
						fprintf ( stderr, "\nPB: Sequence %s, Residue %d : Cache=%d",A->name[b], nres2,fixed_nres2 );
						myexit(EXIT_FAILURE);
			      }
			    
			    if ( fixed_nres1 && fixed_nres2)
			      {
				vwrite_clist (CL,CL->ne, SEQ1, s1);
				vwrite_clist (CL,CL->ne, SEQ2, s2);
				vwrite_clist (CL,CL->ne, R1,fixed_nres1);
				vwrite_clist (CL,CL->ne, R2,fixed_nres2);
				
				if (do_pdb)
				  {
				    
				    vwrite_clist (CL,CL->ne, WE,(NORM_F/MAXID)*pdb_weight );
				  }
				else
				  {
				    
				    vwrite_clist (CL,CL->ne, WE,(NORM_F/MAXID)*((weight[0]==FORBIDEN)?weight[1]:weight[c]) );
				  }
				vwrite_clist (CL,CL->ne, CONS,1);
				if (!set_misc)
				  {
				    vwrite_clist (CL,CL->ne, MISC,A->len_aln);
				    set_misc=1;
				  }
				else 
				  {
				    vwrite_clist (CL,CL->ne, MISC,0);
				  }
				CL->ne++;
			      }
			  }
		      }
		  }
	      }
	  }
	vfree (top_seq);
	vfree (sindex);
	vfree (weight);
	if (A->A) 
	  {
	    return aln2constraint_list (A->A, CL, weight_mode);
	  }
	else
	  return CL;
	}

double **list2mat (Constraint_list *CLin,int s1,int s2, double *min, double *max)
        {
	double ** mat;
	int a, r1, r2;
	int min_def=0;
	Constraint_list *CL;
	static Sequence *S;

	
	int row, column;
	if ( S==NULL)S=declare_sequence ((CLin->S)->min_len, (CLin->S)->max_len,(CLin->S)->nseq); 	
	sprintf ( S->name[0], "%s",(CLin->S)->name[s1]);
	sprintf ( S->name[1],"%s",(CLin->S)->name[s2]);
        S->nseq=2;
	
        row   =(CLin->S)->len[s1];
	column=(CLin->S)->len[s2];
	
	if ( CLin->extend_jit)	
	    CL=extend_list(CLin,"mem",CLin->extend_clean_mode, CLin->extend_compact_mode, CLin->do_self, S);
	else
	    CL=CLin;


	min[0]=max[0];
	mat=declare_double ( row, column);

	for ( a=0; a<CL->ne; a++)
	    {
	    r1=vread_clist(CL,a,R1)-1;
	    r2=vread_clist(CL,a,R2)-1;
	    if ( vread_clist(CL,a,SEQ1)==s1 &&vread_clist(CL,a,SEQ2)==s2)
		{
		mat[r1][r2]=(double)vread_clist(CL,a,WE);
		if (min_def==0)
		   {
		   min_def=1;
		   min[0]=mat[r1][r2];
		   max[0]=mat[r1][r2];
		   }
		else
		   {
		   min[0]=(min[0]<mat[r1][r2])?min[0]:mat[r1][r2];
		   max[0]=(max[0]>mat[r1][r2])?max[0]:mat[r1][r2];
		   }
		}
	    else if (vread_clist(CL,a,SEQ2)==s1 &&vread_clist(CL,a,SEQ1)==s2) 
		{
		mat[r2][r1]=(double)vread_clist(CL,a,WE);
			if (min_def==0)
		   {
		   min_def=1;
		   min[0]=mat[r2][r1];
		   max[0]=mat[r2][r1];
		   }
		else
		   {
		   min[0]=(min[0]<mat[r2][r1])?min[0]:mat[r2][r1];
		   max[0]=(max[0]>mat[r2][r1])?max[0]:mat[r2][r1];
		   }
		}
	    }
	return mat;
	}

Constraint_list * constraint_list2bin_file(Constraint_list *clist)
        {
	int a,b;
	
	clist->fp=tmpfile();
	for ( a=0; a< clist->ne; a++)
	    for ( b=0; b<clist->entry_len; b++)
	        {
		fwrite (&clist->L[a*clist->entry_len+b],clist->el_size, 1,clist->fp);
		}
	return clist;
	}

FILE * bin_file2constraint_list ( Constraint_list *CL, FILE *fp, char *name)
        {
	int a, b, s1, s2;
	CLIST_TYPE *entry;
	
	if ( fp==NULL)fp=vfopen ( name, "w");
	entry=vcalloc ( CL->entry_len, CL->el_size);
	fprintf ( fp, "%d\n", (CL->S)->nseq);
	for ( a=0; a< (CL->S)->nseq; a++)fprintf (fp, "%s %d %s\n", (CL->S)->name[a], (CL->S)->len[a], (CL->S)->seq[a]);
	

	rewind ( CL->fp);
	fread(entry, CL->el_size, CL->entry_len, CL->fp);
	s1=entry[SEQ1];
	s2=entry[SEQ2];
	fprintf (fp, "#%d %d\n", s1, s2);
	for ( b=2; b< CL->entry_len; b++)fprintf (fp, "%5d ",entry[b]);
	fprintf (fp, "\n");
	for ( a=1; a< (CL->ne); a++)
	    {
	    fread(entry, CL->el_size, CL->entry_len, CL->fp);
	    if ( entry[SEQ1]!=s1 || entry[SEQ2]!=s2)
	       {
	       s1=entry[SEQ1];
	       s2=entry[SEQ2];
	       fprintf (fp, "#%d %d\n", s1, s2);
	       }
	    for ( b=2; b< CL->entry_len; b++)fprintf (fp, "%5d ",entry[b]);
	    fprintf (fp, "\n");
	    }
	fprintf (fp, "! CPU %d\n",get_time());
	
	return fp;
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
	int s1, s2, r1, r2, w, a, b;


	tot_weight=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	for ( a=0; a<CL->ne; a++)
	    {
	    r1=vread_clist(CL,a,R1)-1;
	    r2=vread_clist(CL,a,R2)-1;
	    s1=vread_clist(CL,a,SEQ1);
	    s2=vread_clist(CL,a,SEQ2);
	    w=vread_clist(CL,a,WE);
	    tot_weight[s1][r1]+=w;
	    tot_weight[s2][r2]+=w;
	    }
	for (a=0; a<(CL->S)->nseq; a++)
	  for (b=0; b<(CL->S)->len[a]; b++)
	    tot_weight[a][b]/=(CL->S)->nseq;
	
	return tot_weight;
	}

int **list2residue_total_extended_weight ( Constraint_list *CL)
        {
	  /*Returns 
	    tot_extended_weight[nseq][maxlen]
	    where each residue is associated with the total of its weights in CL
	    ####IMPORTANT
            
	          -the numbering of the residues  goes from 1 to L:
	          -the numbering of the sequences goes from 0 to N-1:
	  */

	static int **tot_extended_weight;
	int s1, s2, r1, r2, w;
	
	if (CL->residue_indexed && tot_extended_weight);
	else
	  {
	    if (tot_extended_weight) free_int (tot_extended_weight, -1);
	    if (CL->residue_indexed==0)index_res_constraint_list (CL,WE);
	    

	    tot_extended_weight=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	    
	    for ( s1=0; s1< (CL->S)->nseq-1; s1++)
	      for ( s2=s1+1; s2< (CL->S)->nseq; s2++)
		for (r1=1; r1<=(CL->S)->len[s1]; r1++)
		  for (r2=1; r2<=(CL->S)->len[s2]; r2++)
		    {
		      w=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
		      tot_extended_weight[s1][r1]+=w;
		      tot_extended_weight[s2][r2]+=w;
		    }
	  }
	return tot_extended_weight;
	}
int **list2residue_partial_extended_weight ( Constraint_list *CL)
        {
	  /*Returns 
	    tot_extended_weight[nseq][maxlen]
	    where each residue is associated with the total of its weights in CL
	    ####IMPORTANT
            
	          -the numbering of the residues  goes from 1 to L:
	          -the numbering of the sequences goes from 0 to N-1:
	  */

	int **tot_extended_weight;
	int s1, s2, r1, r2, w1, w2, a;


	tot_extended_weight=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	for ( a=0; a<CL->ne; a++)
	    {
	    r1=vread_clist(CL,a,R1);
	    r2=vread_clist(CL,a,R2);
	    s1=vread_clist(CL,a,SEQ1);
	    s2=vread_clist(CL,a,SEQ2);
	    w1=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
	    w2=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);
	    if ( w1!=w2)fprintf ( stderr, "*");

	    tot_extended_weight[s1][r1]+=w1;
	    tot_extended_weight[s2][r2]+=w2;
	    }
	return tot_extended_weight;
	}
  

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              clean functions                                            */
/*                                                                                         */
/*                                                                                         */
/*                                                                                         */
/*******************************************************************************************/
Constraint_list *clean ( char *clean_mode,Constraint_list *CL,int start, int len)
	{
	
	if ( strm ( clean_mode, "shadow"))   CL=clean_shadow (CL,start,len);	
	else if ( strm5( clean_mode, "","NO","no","No","default"));
 	else	add_warning ( CL->local_stderr, "\nWARNING: The %s CLEANING MODE DOES NOT EXIST\n", clean_mode);
	
	return CL;
	}


Constraint_list * clean_shadow ( Constraint_list *CL, int start, int len)
	{
	int s1, s2, r1, a, b, end;
	int max, min;
		
	s1=vread_clist (CL, start, SEQ1);
	s2=vread_clist (CL, start, SEQ2);
	r1=vread_clist (CL, start, R1);
	
	
	for ( a=start; a<(start+len);)
		{
		
		max=min=vread_clist (CL, a, WE);
		while ( a<CL->ne && vread_clist (CL, a, SEQ1)==s1 && vread_clist (CL, a, SEQ2)==s2 && vread_clist (CL, a, R1)==r1)
			{
			max=(vread_clist (CL, a, WE)>max)?vread_clist (CL, a, WE):max;
			min=(vread_clist (CL, a, WE)<min)?vread_clist (CL, a, WE):min;
			a++;
			}
		end=a;
		
		if ((end-start)>1)
			{
			for ( b=start; b<end; b++)
				if ( vread_clist (CL, b, WE)<max)vwrite_clist (CL, b, SEQ1,-1);
			}
		start=end;
		if ( start<CL->ne)
			{
			s1=vread_clist (CL, start, SEQ1);
			s2=vread_clist (CL, start, SEQ2);
			r1=vread_clist (CL, start, R1);
			}
		}
	CL=sort_constraint_list_inv (CL, start, (CL->ne-start));
	CL=sort_constraint_list (CL,start,(CL->ne-start)  );	

	return CL;
	}
/*********************************************************************/
/*                                                                   */
/*                         LIST FUNCTIONS                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	

static Constraint_list *fast_merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode);
static Constraint_list *slow_merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode);

Constraint_list *merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode)
{


  if ( !ML) 
    {
      return SL;
    }
  else if ( SL->S == ML->S)
    {
      return fast_merge_constraint_list (SL, ML, mode);
    }
  else
    {
      return slow_merge_constraint_list (SL, ML, mode);
    }
}

Constraint_list *fast_merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode)
{
  int l;

  l=ML->ne+SL->ne;

  ML->max_L_len=l;
  ML->seq_indexed=0;
  ML->residue_indexed=0;
  if (ML->ne==0 || !ML->L)ML->L=vcalloc (l*ML->entry_len, sizeof (int));
  else ML->L=vrealloc (ML->L, sizeof (int)*l*ML->entry_len);
  
  memcpy (ML->L+(ML->ne*ML->entry_len), SL->L, SL->ne*sizeof (int)*SL->entry_len);
 
  ML->ne+=SL->ne;
  return ML;
  
}

Constraint_list *slow_merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode)
{
  int a, s1, s2;
  Sequence *S1, *S2;
  int **name_index;
  int **seq_index;
  CLIST_TYPE *entry=NULL;

  if ( !ML)return SL;
  
  S1=SL->S;S2=ML->S;
  
  name_index=index_seq_name(S1,S2);
  
  seq_index=index_seq_res (S1,S2, name_index);
  
  for ( a=0; a< SL->ne; a++)
    {
      
      entry=extract_entry ( entry, a, SL);
      //HERE ("BEF: %d %d %d %d", entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]); 
      s1=entry[SEQ1]; s2=entry[SEQ2];
 
      if ( S1==S2)
	{
	 add_entry2list(entry, ML);
	}
      else if (name_index[s1][0]==-1 || name_index[s2][0]==-1 || !seq_index[s1] || !seq_index[s2]);
      else 
	{
	  int r1, r2;
	  r1=seq_index[s1][entry[R1]-1];
	  r2=seq_index[s2][entry[R2]-1];
	
	  entry[SEQ1]=name_index[s1][0];
	  entry[SEQ2]=name_index[s2][0];
	  if ( r1!=-1 && r2!=-1 && entry[SEQ1]!=entry[SEQ2])
	    {
	      entry[R1]=r1+1; entry[R2]=r2+1;
	      add_entry2list(entry, ML);
	    }
	  //HERE ("AFT: %d %d %d %d", entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]); 
	}
    }
  free_int (name_index, -1);
  free_int ( seq_index, -1);
  vfree ( entry);
  return ML;
}
  
Constraint_list *modify_weight( Constraint_list *CL,int start, int end,  char *modify_mode)
        {
	int a;
	CLIST_TYPE x;
	
	if ( strm(modify_mode, "default"))return CL;
	for ( a=start; a<end; a++)
	    {
	    x=vread_clist(CL, a, WE);
	    
	    if (strm2 (modify_mode,"sc_eq_cons", "we_eq_cons"))
	        if(x!=UNDEFINED)
		    vwrite_clist(CL, a, WE,  vread_clist(CL, a, CONS));
	    
	    if (strm2(modify_mode,"sc_eq_wePcons","sc_eq_consPwe"))
	        if(x!=UNDEFINED)
		     vwrite_clist(CL, a, WE, vread_clist(CL, a, CONS)*x);		
	    }
	return CL;
	}
	
Constraint_list *compact_list (Constraint_list *CL, int start, int len, char *compact_mode)
	{
	int a;
	int r1, r2, rr1, rr2, s1, rs1, s2, rs2, ra;
	CLIST_TYPE x;
	int debug_compact=0;

	if (debug_compact)fprintf ( stderr, "\n[In: %d %s start=%d len=%d", CL->ne, compact_mode, start, len);

	if ( len==0  || strm3(compact_mode, "no", "No", "NO"))return CL;
	else if ( strm2(compact_mode,"mirror","mirror_sum"));
	else if ( strm4(compact_mode, "default","shrink","shrink_best","shrink_worst"))
	        {
		
		for ( a=start; a<(start+len) ; a++)
		    {
		    
		    if ( vread_clist(CL, a, SEQ1)> vread_clist(CL, a, SEQ2) ||\
		       ( vread_clist(CL, a, SEQ1)==vread_clist(CL, a, SEQ2) &&\
			 vread_clist(CL, a, R1)  > vread_clist(CL, a, R2)     ))
			
			
		        {
			s1=vread_clist(CL, a, SEQ1);
			s2=vread_clist(CL, a, SEQ2);
			r1=vread_clist(CL, a, R1);
			r2=vread_clist(CL, a, R2);
			vwrite_clist(CL, a, SEQ1,s2);
			vwrite_clist(CL, a, SEQ2,s1);
			vwrite_clist(CL, a, R1,r2);
			vwrite_clist(CL, a, R2,r1);
			}
		    }
		}
	
	if (debug_compact)fprintf ( stderr, "\n[2: %d %s start=%d len=%d", CL->ne, compact_mode, start, len);

	sort_constraint_list ( CL, start, len);
	
	rs1=vread_clist(CL, start, SEQ1);	
	rs2=vread_clist(CL, start, SEQ2);
	rr1=vread_clist(CL, start, R1);
	rr2=vread_clist(CL, start, R2);
	ra=start;

	if (debug_compact)fprintf ( stderr, "\n[3: %d %s start=%d len=%d", CL->ne, compact_mode, start, len);


	if ( (rs1==rs2) && (rr1==rr2))vwrite_clist(CL, start, SEQ1,-1);		
	for ( a=start+1; a<(start+len); a++)
		{
		s1=vread_clist(CL, a, SEQ1);
		s2=vread_clist(CL, a, SEQ2);
		r1=vread_clist(CL, a, R1);
		r2=vread_clist(CL, a, R2);
		
		//if ( (s1==s2) && (r1==r2))vwrite_clist(CL, a, SEQ1, -1);
		if ( s1==rs1 && s2==rs2 && r1==rr1 && r2==rr2)
			{
			x=vread_clist(CL, ra, WE);
			if (strm ( compact_mode, "shrink"));
			else if (strm( compact_mode,"default"))//default is set to best
			  vwrite_clist(CL, ra, WE,MAX(vread_clist(CL, a, WE),x));
			else if ( strm ( compact_mode,"mirror_sum"))    
			  vwrite_clist(CL, ra, WE, vread_clist(CL, a, WE)+x);
			else if (strm2 ( compact_mode,"best", "shrink_best"))
			  vwrite_clist(CL, ra, WE,MAX(vread_clist(CL, a, WE),x));
			else if (strm2 ( compact_mode, "worst","shrink_worst"))
			  vwrite_clist(CL, ra, WE,MIN(vread_clist(CL, a, WE), vread_clist(CL, a, WE)));
		
			if (  strm(compact_mode, "shrink"));
			else
			    {
			    vwrite_clist(CL, ra, CONS, vread_clist(CL, ra, CONS)+ vread_clist(CL, a, CONS));
			    vwrite_clist(CL, ra, MISC, vread_clist(CL, ra, MISC)+ vread_clist(CL, a, MISC));
			    }
			vwrite_clist(CL,a, SEQ1, -1);
			
			}
		else
			{
			rs1=s1;
			rs2=s2;
			rr1=r1;
			rr2=r2;
			ra=a;
			}
		}
	
	
	sort_constraint_list_inv(CL,0,CL->ne);
	
	sort_constraint_list    (CL,0,CL->ne);
	
	
	if ( strm3 (compact_mode, "consPwe", "wePcons","cons"))
		{
		for ( a=start; a<(start+len); a++)
			{
			if ( strm2(compact_mode,"consPwe", "wePcons"))
			    vwrite_clist(CL, a, WE, vread_clist(CL,a,WE)* vread_clist(CL,a,CONS));
			else if  (strm (compact_mode, "cons"))
			     vwrite_clist(CL, a, WE, vread_clist(CL,a,CONS)*100);
			}
		}
	if (debug_compact)fprintf ( stderr, "....OUT: %d]\n", CL->ne);	
	return CL;
	}


Constraint_list *rescale_list_simple (Constraint_list *CL,int start, int len,int new_min, int new_max)
	{
	int a, min, max;
	double x;
	/*Rescales between 0 and max2
	  Any value above max1 is set to max1 first and then to max2
	*/



	min=max=vread_clist ( CL,start, WE);
	
	for ( a=start; a<(start+len); a++)
		{
		x=(double)vread_clist ( CL,a, WE);
		if ( x>max)max=(int)x;
		if ( x< min)min=(int)x;
		}
	
	fprintf ( CL->local_stderr, "\n[%d-%d]=>[%d-%d]", min, max, new_min, new_max);

	for ( a=start; a<(start+len); a++)
	    {
	      
	    x=vread_clist(CL,a, WE);
	   
	    if ((max-min)==0)x=100;
	    else x=(((x-min)/(max-min))*new_max)+new_min;
	    
	    vwrite_clist(CL, a, WE,(CLIST_TYPE) x);
	    }
	return CL;
	}
Constraint_list *rescale_list (Constraint_list *CL,int start, int len,int max1, int max2)
	{
	int a, min_val, max_val;
	CLIST_TYPE x;
	
	/*Rescales between 0 and max2
	  Any value above max1 is set to max1 first and then to max2
	*/


	min_val=0;
	max_val=max1;

	
	
	for ( a=start; a<len; a++)
		{
		if (vread_clist ( CL,a, WE)>max1)vwrite_clist(CL, a, WE, max1);
		}
	
	for ( a=start; a<len; a++)
	    {
	    x=vread_clist(CL,a, WE);
	    vwrite_clist(CL, a, WE, (((x-min_val)*max2)/(max_val-min_val)));
	    }
	return CL;
	}


Constraint_list* filter_list (Constraint_list *CL, int start, int len,int T)
	{
	int a;
	int field;
	
	if (T==0)return CL;
	
	field=WE;
	if (T>0)field=WE;
	else if ( T<0)
		{
		field=CONS;
		T=-T;
		}
	
	
	for ( a=start; a<len; a++)
	    if (vread_clist(CL, a, field)<=T)vwrite_clist(CL,a,SEQ1,-1);
	
	CL=sort_constraint_list_inv (CL, 0, CL->ne);
	CL=sort_constraint_list     (CL, 0, CL->ne);
	return CL;
	}





Constraint_list *undefine_list (Constraint_list *CL)
      {
      int a, b;
      int undefined_flag;
      
      for ( a=0;a<CL->ne; a++)
          {
	  for ( b=0, undefined_flag=0; b< LIST_N_FIELDS; b++)
	      {

	      if ( vread_clist (CL, a, b)==UNDEFINED)undefined_flag=1;	  
	      if ( undefined_flag)
	          { 
		  for ( b=0; b< LIST_N_FIELDS; b++)
		      if ( b!=SEQ1 && b!=SEQ2 && b!=R1 && b!=R2)
			  {
			  vwrite_clist(CL, a, b, UNDEFINED);
			  }
		  }
	      }
	  }
      return CL;
      }

int ** seq2defined_residues ( Sequence *S, Constraint_list *CL)
{
  int **seq_count;
  int *entry=NULL;
  int a;
  
  seq_count=declare_int (S->nseq, S->max_len+1);
  for (a=0; a< CL->ne; a++)
      {
	entry=extract_entry(entry, a, CL);
	seq_count[entry[SEQ1]][entry[R1]]++;
	seq_count[entry[SEQ2]][entry[R2]]++;
      }
  vfree (entry);
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
void check_seq_pair_in_list(Constraint_list *CLin,int seq1, int seq2)
      {
      int a, s1, s2, r1, r2;
	   
      for ( a=0; a< CLin->ne; a++)
          {
	  s1=vread_clist(CLin,a,SEQ1);
	  s2=vread_clist(CLin,a,SEQ2);
	  if ( s1==seq1 && s2==seq2)
	     {
	     r1=vread_clist(CLin,a,R1);
	     r2=vread_clist(CLin,a,R2);
	     fprintf ( stderr, "\n[%d][%d %d] [%d %d]",a,s1, r1, s2, r2);
	     }
	  }
      }

void print_CL_mem(Constraint_list *CL, char *function)
     {
      fprintf ( stderr,"%s\n", function);
     if ( CL->fp==NULL && CL->L==NULL) fprintf ( stderr, "\n\tNOTHING");
     if ( CL->fp)fprintf ( stderr, "\n\tFILE SET");
     if ( CL->L)fprintf ( stderr, "\n\tMEM SET\n");
     }

int constraint_list_is_sorted ( Constraint_list *CL)
     {
     int a,b, x1, x2;
     for ( a=0; a< CL->ne-1; a++)
         {
	 for ( b=0; b< CL->entry_len; b++)
	     {
	     x1=vread_clist( CL, a, b);
	     x2=vread_clist( CL, a+1,b);
	     if ( x1<x2)break;
	     else if ( x1>x2)
		 {
		 fprintf ( stderr, "\n[%d][%d]=>%d\n[%d][%d]=>%d\n\n",a, b, x1, a+1, b, x2);
		 return 0;
		 }
	     
	     }
	 }
     return 1;
     }
/*********************************************************************/
/*                                                                   */
/*                        PRUNE CONSTRAINT_LIST                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char * list2prune_list_old ( Sequence *S, int **sm)
{
  int a, b, i;
  char *aln, *seq, *file;
  FILE *fp;
  Sequence *subS;
  Alignment *A;
  
  aln=vtmpnam (NULL);
  seq=vtmpnam (NULL);
  file=vtmpnam (NULL);
  
  output_fasta_seq (seq, A=seq2aln (S,NULL,RM_GAP));
  free_aln (A);
  
  printf_system ( "t_coffee %s -in Xblosum62mt -outfile=%s -msa_mode iterative_tree_aln", seq, aln);
  printf_system  ( "t_coffee -other_pg seq_reformat -in %s -action +trim _aln_n5 -output fasta_seq > %s",aln, seq );
  subS=main_read_seq ( seq);

  fp=vfopen (file, "w");
  for ( a=0; a< subS->nseq; a++)
    {
      i=name_is_in_list (subS->name[a], S->name, S->nseq, 100);
      if ( i==-1) continue;
      for ( b=0; b<S->nseq; b++)
	if (i!=b)fprintf ( fp, "\n2 %d %d",i, b);
    }
  vfclose (fp);
  return file;
}

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
  
  HERE ("NS=%d", ns);
  if (ns==0)ns=n;
  else if (ns<0)ns=-(n*ns)/100;
  else if (ns>=n)ns=n;
  
  
  HERE ("NS=%d", ns);
  

  keep=vcalloc (n, sizeof (int));
  used=vcalloc (n, sizeof (int));
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
  

char * list2prune_list_old3 ( Sequence *S, int **sm)
{
  int **keep;
  int a,b,c,n,s1,s2, tot, nseq;

  char *file;
  FILE *fp;
  nseq=S->nseq;
  
  keep=declare_int (nseq, nseq);

  for (a=0; a< nseq; a++)
    for (b=a+1; b< nseq; b++)
      {
	int bc=0,bsim=0;
	s1=sm[a][b];
	for (c=0; c< nseq; c++)
	  {
	    if ( c==a || c==b) continue;
	    s2=MIN(sm[a][c], sm[c][b]);
	    
	    if (s2>s1 && s2>bsim)
	      {
		bsim=s2;
		bc=c;
	      }
	  }
	if ( bsim) 
	  {
	    keep[a][bc]=1;
	    keep[bc][a]=1;
	    keep[b][bc]=1;
	    keep[bc][b]=1;
	  }
      }

  file=vtmpnam (NULL);
  fp=vfopen (file, "w");
  
  for (n=0,tot=0,a=0; a< nseq; a++)
    for ( b=a+1; b< nseq; b++)
      {
	if (keep[a][b])
	  fprintf (fp, "\n 2 %d %d", a, b);
      }
  vfclose (fp);
  return file;
}
  
/*********************************************************************/
/*                                                                   */
/*                        WEIGHT CONSTRAINT_LIST                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

Constraint_list *weight_constraint_list(Constraint_list * CL, char *seq_weight)

    {
      Weights *W;
     
      if ( CL->ne==0)return CL;
      else if ( strm(seq_weight, "t_coffee")) W=compute_t_coffee_weight(CL);
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
      

      


      if (!CL->L)return NULL;
      
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
      
Constraint_list *re_weight_constraint_list(Constraint_list * CL,Weights *W)
    {
      int a;
      float w;
      float *weight;
      int sA, sB;



      weight=W->SEQ_W;

      if (!CL->L)return CL;
 
      
      
      for ( a=0; a< CL->ne; a++)
	   {
	     sA=CL->L[a*CL->entry_len+SEQ1];
	     sB=CL->L[a*CL->entry_len+SEQ2];
	     
	     w=MIN(weight[sA], weight[sB]);
	     
	     CL->L[a*CL->entry_len+WE]*=w;
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
	  DM=vcalloc ( 1, sizeof (Distance_matrix));
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
	  if (!CL || !CL->L || CL->ne==0)
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
  
      
      ns=vcalloc ( 2, sizeof(int));
      l_s=declare_int ( 2, 1);
      ns[0]=ns[1]=1;
      l_s[0][0]=0;
      l_s[1][0]=1;
      
      if (CL->local_stderr && print>0)fprintf ( (CL->local_stderr), "\nCOMPUTE PAIRWISE SIMILARITY [dp_mode: %s] [distance_matrix_mode: %s][Similarity Measure: %s] \n", NCL->dp_mode,mode, sim_mode);	
      
      for (a=0; a< (CL->S)->nseq; a++)
	{
	if (CL->local_stderr && print>0)fprintf ( (CL->local_stderr), "\n\tSeq: %s", (CL->S)->name[a]);
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
		    if ( CL->L)
		      {
			score=(int)(((float)B->score_aln)/(B->len_aln*SCORE_K));
			score=(int)(CL->L && CL->normalise)?((score*MAXID)/(CL->normalise)):(score);
		      }
		    else if ( CL->M)score=id;
		    
		    
		    if ( score>MAXID)score=(CL->L)?sub_aln2sub_aln_score (B, CL, CL->evaluate_mode, ns, l_s):id;
		     
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
		    //HERE ("%s %d %d ->%d", sim_mode, a, b, (int)id);
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
		if (CL->local_stderr && print>1) fprintf (CL->local_stderr, "\n\t%-*s %-*s identity=%3d%% score=%3d", max_name,(CL->S)->name[a], max_name,(CL->S)->name[b], id_score, (int)score);
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

Constraint_list *read_rna_lib ( Sequence *S, char *fname)
{
  Constraint_list *R;
  char **list;
  int n=0,a;
  

  if (check_file_exists (fname))
    {

      list=read_lib_list ( fname, &n);
    }
  else
    {
      X_template *F;

      list=vcalloc (S->nseq, sizeof (char*));
      for ( a=0; a<S->nseq; a++)
	{
	  if ((F=seq_has_template (S, a, "_F_")))
	    {
	      list[n++]=F->template_file;
	    }
	}
    }

  R=declare_constraint_list ( S,NULL, NULL, 0,NULL, NULL); 
      
  for (a=0; a< n; a++)
    {

      if (list[a])R=fast_read_constraint_list_file (R, list[a]);
    }

  R=index_res_constraint_list (R,WE);
  
  return R;
}

Constraint_list * rna_lib_extension ( Constraint_list *CL, Constraint_list *R)
{
  CLIST_TYPE  *entry=NULL;
  int a,b,c,n1,n2, ne,s1, s2, r1, r2,w;
  int list1[100], list2[100];
  
  
  entry=vcalloc ( CL->entry_len, CL->el_size);
  ne=CL->ne;
  

  for ( a=0; a<ne; a++)
    {

      
      extract_entry (entry, a, CL);
      s1=entry[SEQ1];
      s2=entry[SEQ2];
      r1=entry[R1];
      r2=entry[R2];
      w=entry[WE];
      
      n1=n2=0;
      list1[n1++]=r1;
      for (b=1; b<R->residue_index[s1][r1][0]; b+=3)
	{
	  list1[n1++]=R->residue_index[s1][r1][b+1];
	}
      list2[n2++]=r2;
      for (b=1; b<R->residue_index[s2][r2][0]; b+=3)
	{
	  list2[n2++]=R->residue_index[s2][r2][b+1];
	}
        
      for (b=1; b<n1; b++)
	for (c=1; c<n2; c++)
	{

	  
	  entry[R1]=list1[b];
	  entry[R2]=list2[c];
	  add_entry2list ( entry,CL);
	}
    }
  
  return CL;
}

char *** produce_method_file ( char *method)
{
  static char ***list;
  int n=0;
  FILE *fp;

  if (!list)list=declare_arrayN(3, sizeof (char),1000,2, 100);
  

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
  
  sprintf (list[n][0], "profile_pair");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w"); 
  fprintf ( fp, "EXECUTABLE profile_pair\n");
  fprintf ( fp, "EXECUTABLE2 clustalw\n");
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
	  add_warning (stderr, "\n******************** WARNING: ****************************************\nSAP is not installed\nTMalign will be used instead\ntmalign is FASTER than SAP and *almost* as accurate\n**********************************************************************\n");
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
		fprintf ( fp, "EXECUTABLE rna_pair\n");
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
  
  sprintf (list[n][0], "blast_msa");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w"); 
  fprintf ( fp, "DOC: BLAST multiple Aligner [%s]\n", NCBIBLAST_ADDRESS);
  fprintf ( fp, "EXECUTABLE seq_msa\n");
  fprintf ( fp, "EXECUTABLE2 blastpgp\n" );
  fprintf ( fp, "ALN_MODE   multiple\n");
  fprintf ( fp, "OUT_MODE   fL\n");
  fprintf ( fp, "IN_FLAG    -infile=\n");
  fprintf ( fp, "OUT_FLAG   -outfile=\n");
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", NCBIBLAST_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", NCBIBLAST_4_TCOFFEE);
  vfclose (fp);} 	   
  
  sprintf (list[n][0], "clustalw2_pair");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
  fprintf ( fp, "DOC: clustalw [%s]\n", CLUSTALW2_ADDRESS);
  fprintf ( fp, "EXECUTABLE clustalw\n");
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
  fprintf ( fp, "EXECUTABLE clustalw\n");
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
  fprintf ( fp, "EXECUTABLE clustalw\n");
  fprintf ( fp, "ALN_MODE   multiple\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "IN_FLAG    %sINFILE=\n", CWF);
  fprintf ( fp, "OUT_FLAG   %sOUTFILE=\n", CWF);
  fprintf ( fp, "PARAM      %sOUTORDER=INPUT %sNEWTREE=SCRATCH_FILE %salign\n",CWF,CWF,CWF);
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", CLUSTALW_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", CLUSTALW_4_TCOFFEE);
  vfclose (fp);} 	   
  
  sprintf (list[n][0], "mafftdef_pair");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "EXECUTABLE mafft\n");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "ALN_MODE   pairwise\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "PARAM1  \n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   
  

  sprintf (list[n][0], "mafftdef_msa");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "EXECUTABLE mafft\n");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "ALN_MODE   multiple\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, " \n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   
  
  sprintf (list[n][0], "mafft_pair");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
    fprintf ( fp, "EXECUTABLE mafft\n");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "ALN_MODE   pairwise\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "PARAM1 --localpair --maxiterate 1000 \n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   
  
  sprintf (list[n][0], "mafft_msa");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
   fprintf ( fp, "EXECUTABLE mafft\n");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "ALN_MODE   multiple\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "PARAM1 --localpair --maxiterate 1000 \n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "PARAM      \n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   

  sprintf (list[n][0], "mafftjtt_pair");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "EXECUTABLE mafft \n");
  fprintf ( fp, "ALN_MODE   pairwise\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM1 --jtt 250 --localpair --maxiterate 1000 \n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   

  sprintf (list[n][0], "mafftjtt_msa");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "EXECUTABLE mafft \n");
  
  fprintf ( fp, "ALN_MODE   multiple\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM1 --jtt 250 --localpair --maxiterate 1000 \n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   
	
  sprintf (list[n][0], "mafftgins_pair");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
   fprintf ( fp, "EXECUTABLE mafft\n");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "ALN_MODE   pairwise\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "PARAM1 --globalpair --maxiterate 1000 \n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
  vfclose (fp);} 	   

  sprintf (list[n][0], "mafftgins_msa");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w");
   fprintf ( fp, "EXECUTABLE mafft\n");
  fprintf ( fp, "DOC Mafft [%s]\n", MAFFT_ADDRESS);
  fprintf ( fp, "ALN_MODE   multiple\n");
  fprintf ( fp, "OUT_MODE   aln\n");
  fprintf ( fp, "PARAM1 --globalpair --maxiterate 1000 \n");
  fprintf ( fp, "IN_FLAG    &bnsp\n");
  fprintf ( fp, "OUT_FLAG   >\n");
  fprintf ( fp, "PARAM      &bnsp2>/dev/null\n"); 
  fprintf ( fp, "SEQ_TYPE   S\n");
  fprintf ( fp, "ADDRESS    %s\n", MAFFT_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", MAFFT_4_TCOFFEE);
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
  
  sprintf (list[n][0], "prank_msa");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w"); 
  fprintf ( fp, "DOC prank [%s]\n", PRANK_ADDRESS);
  fprintf ( fp, "EXECUTABLE tc_generic_method.pl\n");
  fprintf ( fp, "ALN_MODE  multiple\n");
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
  
  sprintf (list[n][0], "ncbi_netblast");
  sprintf (list[n][1], "%s", vtmpnam(NULL));
  n++;if (method==NULL || strm (method, list[n-1][0])){fp=vfopen (list[n-1][1], "w"); 
  fprintf ( fp, "EXECUTABLE blastcl3 \n");
  fprintf ( fp, "ALN_MODE   predict\n");
  fprintf ( fp, "SEQ_TYPE   PROTEIN\n");
    
  fprintf ( fp, "ADDRESS    %s\n", NCBIWEBBLAST_ADDRESS);
  fprintf ( fp, "PROGRAM    %s\n", NCBIWEBBLAST_4_TCOFFEE);  
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

  list[n]=NULL;
  return list;
}
  

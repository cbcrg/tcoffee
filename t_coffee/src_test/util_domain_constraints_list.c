#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

/*********************************************************************/
/*                                                                   */
/*                         MISCEANELLOUS                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
int ** get_undefined_list (Constraint_list *CL)
         {
	 int **list;
	 int a;
	 CLIST_TYPE x;
	 int *e;
	 list=declare_int ( (CL->S)->nseq+1, (CL->S)->max_len+1);


	 while (e=extract_entry (CL))
	   {
	     x=e[WE];
	     list[e[SEQ1]][e[R1]]=(x==UNDEFINED);
	     list[e[SEQ2]][e[R2]]=(x==UNDEFINED);
	   }
	 return list;
	 }




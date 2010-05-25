/*--------------------------------------------------------------------------
 *
 *  mshell.c -- Memory management utilities
 *  Description: mshell.c contains routines to protect the programmer
 *  from errors in calling memory allocation/free routines. The standard
 *  library calls supported are: malloc, realloc, strdup, free
 *
 *  The interface header mshell.h redefines these standard library
 *  calls as macros. When the client code is compiled, the macros expand
 *  to calls to this module. This module then calls the system memory
 *  routines, ith additional error checking.
 *
 *  The allocation routines in this module add a data structure to the
 *  top of allocated memory blocks which tag them as legal memory blocks.
 *
 *  When the free routine is called, the memory block to be freed is
 *  checked for legality. If the block is not legal, the memory list
 *  is dumped to stderr and the program is terminated.
 *
 *  Compilation options:
 *  MEM_LIST	Link all allocated memory blocks onto an internal list.
 *	The list can be displayed using Mem_Display().
 *  MEM_WHERE	Save the file/line number of allocated blocks in header.
 *	Requires that compiler supports __FILE__ and __LINE__ preprocessor
 *	directives and __FILE__ string have static or global scope.
 *  Problems: If you have any problems with this module, please contact:
 *	Jim Schimandle, Primary Syncretics, 473 Sapena Court, Suite #6,
 *	Santa Clara, CA 95054, Tel: (408)988-3818, Fax: (408)727-9891
 *  Reference: Encapsulating C memory allocation. Jim Schimandle.
 *	Dr. Dobbs Journal. vol 15(8): 24-35
 *
 */

#define __MSHELL__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mshell.h"

/* Constants */
#define MEMTAG	0xa55a		/* value for mh_tag */

/* Structures */
typedef struct memnod			/* Memory block header info	*/
{					/*------------------------------*/
    unsigned int	mh_tag;		/* Special ident tag		*/
    size_t		mh_size;	/* Size of allocation block	*/
#if defined(MEM_LIST)
    struct memnod	*mh_next;	/* Next memory block		*/
    struct memnod	*mh_prev;	/* Previous memory block	*/
#endif
#if defined(MEM_WHERE)
    char		*mh_file;	/* File allocation was from	*/
    unsigned int	mh_line;	/* Line allocation was from	*/
#endif
} MEMHDR;

/*  Alignment macros -- The macro ALIGN_SIZE defines size of the largest
 *  object that must be aligned on your processor. the RESERVE_SIZE macro
 *  calculates the nearest multiple of ALIGN_SIZE that is larger than the
 *  memory block header. Change ALIGN_SIZE to match alignment requirements
 *  of your processor
 */
#define ALIGN_SIZE	sizeof(double)

#define HDR_SIZE	sizeof(MEMHDR)
#define RESERVE_SIZE	(((HDR_SIZE + (ALIGN_SIZE-1)) / ALIGN_SIZE) \
* ALIGN_SIZE)

/*  Conversion macros -- These macros convert the internal pointer
 *  to a memory block header to/from the pointer used by the client code.
 */
#define CLIENT_2_HDR(a)	((MEMHDR *)  (((char *) (a)) - RESERVE_SIZE))
#define HDR_2_CLIENT(a)	((void *)  (((char *) (a)) + RESERVE_SIZE))

/* Local variables */
static unsigned long	mem_size = 0;	/* Amount of memory used */
#if defined(MEM_LIST)
static MEMHDR	*memlist = NULL;	/* List of memory blocks */
#endif


/* Local functions */
static void	mem_tag_err(void *, char *, int);	/* Tag error */
#if defined(MEM_LIST)
static void	mem_list_add(MEMHDR *);		/* Add block to list */
static void	mem_list_delete(MEMHDR *);	/* Delete block from list */
#define Mem_Tag_Err(a)	mem_tag_err(a, fil, lin)
#else
#define Mem_Tag_Err(a)	mem_tag_err(a, __FILE__, __LINE__)
#endif


/**** Functions accessed only through macros ******************************/

/*--------------------------------------------------------------------------
 *
 *  mem_alloc -- Allocate a memory block
 *  Usage:
 *	void *
 *	mem_alloc(
 *	size_t	size
 *	)
 *  Parameters: size	Size of block in bytes to allocate
 *  Return Value: Pointer to allocated memory block, NULL if not
 *	enough memory.
 *  Description: mem_alloc() makes a protected call to malloc()
 *  Notes: Access this routine using the malloc() macro in mshell.h
 *
 */
void *
mem_alloc(
#if defined(MEM_WHERE)
size_t	size,
char	*fil,
int	lin
#else
size_t	size
#endif
)

{
MEMHDR	*p;

    /* Allocate memory block */
    p = calloc(1, RESERVE_SIZE + size);
    if (p == NULL) {
	return NULL;
    }

    /* Init header */
    p->mh_tag = MEMTAG;
    p->mh_size = size;
    mem_size += size;
#if defined(MEM_WHERE)
    p->mh_file = fil;
    p->mh_line = lin;
#endif

#if defined(MEM_LIST)
    mem_list_add(p);
#endif

    /* Return pointer to client data */
    return HDR_2_CLIENT(p);
}

/*--------------------------------------------------------------------------
 *
 *  mem_realloc -- Reallocate a memory block
 *  Usage:
 *	void *
 *	mem_realloc(
 *	void	*ptr,
 *	size_t	size
 *	)
 *  Parameters: ptr	Pointer to current block
 *		size	Size to adjust block to
 *  Return value: Pointer to new memory block; NULL if memory cannot
 *	be reallocated
 *  Description: mem_realloc() makes a protected call to realloc().
 *  Notes: Access this routine using the realloc() macro in mshell.h
 *
 */

void *
mem_realloc(
#if defined(MEM_WHERE)
void	*ptr,
size_t	size,
char	*fil,
int	lin
#else
void	*ptr,
size_t	size
#endif
)

{
    MEMHDR	*p;

    /* Convert client pointer to header pointer */
    p = CLIENT_2_HDR(ptr);

    /* Check for valid block */
    if (p->mh_tag != MEMTAG)
	{
	Mem_Tag_Err(p);
	return NULL;
	}

    /* Invalidate header */
    p->mh_tag = ~MEMTAG;
    mem_size -= p->mh_size;

#if defined(MEM_LIST)
    mem_list_delete(p);		/* Remove block from list */
#endif

    /* Reallocate memory block */
    p = (MEMHDR *) realloc(p, RESERVE_SIZE + size);
    if (p == NULL)
	{
	return NULL;
	}

    /* Update header */
    p->mh_tag = MEMTAG;
    p->mh_size = size;
    mem_size += size;
#if defined(MEM_WHERE)
    p->mh_file = fil;
    p->mh_line = lin;
#endif

#if defined(MEM_LIST)
    mem_list_add(p);		/* Add block to list */
#endif

    /* Return pointer to client data */
    return HDR_2_CLIENT(p);
}

/*--------------------------------------------------------------------------
 *
 *  mem_strdup -- Save a string in dynamic memory
 *  Usage:
 *	char *
 *	mem_strdup(
 *	char	*str
 *	)
 *  Parameters: str	String to save
 *  Return Value: Pointer to allocated string; NULL if not enough memory
 *  Description: mem_strdup() saves specified string in dynamic memory.
 *  Notes: Access this routine using strdup() macro in mshell.h
 *
 */


char *
mem_strdup(
#if defined(MEM_WHERE)
char	*str,
char	*fil,
int	lin
#else
char	*str
#endif
)

{
    char *s;

#if defined(MEM_WHERE)
    s = mem_alloc(strlen(str)+1, fil, lin);
#else
    s = mem_alloc(strlen(str)+1);
#endif

    if (s != NULL)
	{
	strcpy(s, str);
	}


    return s;

}

/*--------------------------------------------------------------------------
 *
 *  mem_free -- Free a memory block
 *  Usage:
 *	void
 *	mem_free(
 *	void	*ptr
 *	)
 *  Parameters: ptr	Pointer to memory to free
 *  Description: mem_free() frees specified memory block. The block must
 *  be allocated using mem_alloc(), mem_realloc() or mem_strdup().
 *  Notes: Access this routine using the free() macro in mshell.h
 *
 */

void
mem_free(
#if defined(MEM_WHERE)
void	*ptr,
char	*fil,
int	lin
#else
void	*ptr
#endif
)

{
    MEMHDR	*p;

    /* Convert client pointer to header pointer */
    p = CLIENT_2_HDR(ptr);

    /* Check for valid block */
    if (p->mh_tag != MEMTAG)
	{
	Mem_Tag_Err(p);
	return;
	}

    /* Invalidate header */
    p->mh_tag = ~MEMTAG;
    mem_size -= p->mh_size;

#if defined(MEM_LIST)
    mem_list_delete(p);		/* Remove block from list */
#endif

    /* Free memory block */
    free(p);
}


/**** Functions accessed directly *****************************************/

/*--------------------------------------------------------------------------
 *
 *  Mem_Used -- Return amount of memory currently allocated
 *  Usage:
 *	unsigned long
 *	Mem_Used(
 *	)
 *  Parameters: None.
 *  Description: Mem_Used() return number of bytes currently allocated
 *	using the memory management system. Value returned is simply the
 *	sum of the size requests to allocation routines; does not reflect
 *	any overhead required by the memory management system.
 *  Notes: None.
 *
 */

unsigned long
Mem_Used(
void)

{
    return mem_size;
}

/*--------------------------------------------------------------------------
 *
 *  Mem_Display -- Display memory allocation list
 *  Usage:
 *	void
 *	Mem_Display(
 *	FILE	*fp
 *	)
 *  Parameters: fp	File to output data to
 *  Description: Mem_Display() displays contents of memory allocation
 *	list. This function is a no-op if MEM_LIST is not defined.
 *  Notes: none.
 *
 */

void
Mem_Display(
FILE	*fp
)

{
#if defined(MEM_LIST)
    MEMHDR	*p;
    int		idx;

#if defined(MEM_WHERE)
    fprintf(fp, "Index    Size    File(Line) - total size %lu\n", mem_size);
#else
    fprintf(fp, "Index    Size - total size %lu\n", mem_size);
#endif

    idx = 0;
    p = memlist;
    while (p != NULL)
	{
	fprintf(fp, "%-5d %6u", idx++, p->mh_size);
#if defined(MEM_WHERE)
	fprintf(fp, "    %s(%d) ", p->mh_file, p->mh_line);
#endif
	if (p->mh_tag != MEMTAG)
	    {
	    fprintf(fp, "INVALID");
	    }
	fprintf(fp, "\n");
	p = p->mh_next;
    }

#else
    fprintf(fp, "Memory list not compiled (MEM_LIST not defined)\n");

#endif
}

/**** Memory list manipulation functions **********************************/

/*------------------------------------------------------------------------*/
/* mem_list_add() -- Add block to list */

#if defined(MEM_LIST)
static void
mem_list_add(
MEMHDR *p
)

{
    p->mh_next = memlist;
    p->mh_prev = NULL;
    if (memlist != NULL)
	{
	memlist->mh_prev = p;
	}
    memlist = p;

#if defined(DEBUG_LIST)
    printf("mem_list_add()\n");
    Mem_Display(stdout);
#endif
}
#endif

/*------------------------------------------------------------------------*/
/* mem_list_delete() -- Delete block from list */

#if defined(MEM_LIST)
static void
mem_list_delete(
MEMHDR	*p
)

{
    if (p->mh_next != NULL)
	{
	p->mh_next->mh_prev = p->mh_prev;
	}
    if (p->mh_prev != NULL)
	{
	p->mh_prev->mh_next = p->mh_next;
	}
    else
	{
	memlist = p->mh_next;
	}

#if defined(DEBUG_LIST)
    printf("mem_list_delete()\n");
    Mem_Display(stdout);
#endif
}
#endif

/**** Error Display *******************************************************/

/*------------------------------------------------------------------------*/
/* mem_tag_err() -- Display memory tag error */
static void
mem_tag_err(
void	*p,
char	*fil,
int	lin
)

{
    fprintf(stderr, "Memory tag error - %p - %s(%d)\n", p, fil, lin);
#if defined (MEM_LIST)
    Mem_Display(stderr);
#endif
    exit(1);
}

/*------------------------------------------------------------------------*/


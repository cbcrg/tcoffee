/*
 *
 *  mshell .h -- Dynamic memory handler interface
 *
 *  Description: mshell.h provides the interface definitions for the
 *  dynamic memory handler.
 *  See mshell.c for complete documentation.
 *
 */


/* Compilation options */
#define MEM_LIST	/* Build internal list */
#define MEM_WHERE	/* Keep track of memory block source */

/* Interface functions */
unsigned long	Mem_Used(void);
void		Mem_Display(FILE *);

/* Interface functions to access only through macros */
#if defined(MEM_WHERE)
void	*mem_alloc(size_t, char *, int);
void	*mem_realloc(void *, size_t, char *, int);
void	 mem_free(void *, char *, int);
char	*mem_strdup(char *, char *, int);
#else
void	*mem_alloc(size_t);
void	*mem_realloc(void *, size_t);
void	 mem_free(void *);
char	*mem_strdup(char *);
#endif

/* Interface macros */
#if !defined(__MSHELL__)
#if  defined(MEM_WHERE)
#define malloc(a)	mem_alloc( (a), __FILE__, __LINE__ )
#define realloc(a, b)	mem_realloc( (a), (b), __FILE__, __LINE__ )
#define free(a)		mem_free( (a), __FILE__, __LINE__ )
#define strdup(a)	mem_strdup( (a), __FILE__, __LINE__ )
#define calloc(a, b)	malloc( ((a)*(b)) )
#else
#define malloc(a)	mem_alloc( (a) )
#define realloc(a, b)	mem_realloc( (a), (b) )
#define free(a)		mem_free( (a) )
#define strdup(a)	mem_strdup( (a) )
#endif
#endif


#ifndef KM_UTIL_H_
#define KM_UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>


/**
 * \brief Wrapper for fopen.
 *
 * \param name_f The name of the file to open.
 * \param mode The mode in which to open the file
 * \return The opend File. On failure an error message is printed and the program exits.
 */
FILE *
my_fopen(char *name_f, char *mode);

/**
 * \brief Wrapper for the malloc function.
 * \param size The size of the memory block to allocate
 * \return Pointer to the newly allocated block of memory. On failure an error message is printed and the program exits.
 */
void *
my_malloc(size_t size);

/**
* \brief Wrapper for the realloc function.
* \param size The size of the memory block to allocate
* \return Pointer to the newly allocated block of memory. On failure an error message is printed and the program exits.
*/
void *
my_realloc(void *p, size_t size);

/**
* \brief Wrapper for the calloc function.
* \param size The size of the memory block to allocate
* \return Pointer to the newly allocated block of memory. On failure an error message is printed and the program exits.
*/
void *
my_calloc ( size_t num, size_t size );

/**
 * \brief Makes a temporary file.
 *
 * In case of failure, the program exits.
 * \param template The template to used for making the temporary file.
 * \param function The function from where this was called (For error message)
 */
char *
my_make_temp_dir(char *template, char *function);





#endif


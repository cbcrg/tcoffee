/* General purpose header file - rf 12/90 */

#ifndef _H_general
#define _H_general



#define pint int			/* cast ints in printf statements as pint */
typedef int Boolean;			/* Is already defined in THINK_C */

#undef TRUE						
#undef FALSE
#define TRUE 1
#define FALSE 0

#define EOS '\0'				/* End-Of-String */
#define MAXLINE 512			/* Max. line length */


#endif /* ifndef _H_general */

#ifndef __MATPLOT_H__
#define __MATPLOT_H__

/* ==== HEADER matplot.h ==== */

/* GL program for visualising general rectangular MxN matrices. */

/* ANSI C+SGI GL, IRIX 5.3, 26. July 1995. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>

/* ---- DEFINITIONS ---- */

#define DXORIG 100  /* the default lower left corner positions */
#define DYORIG 460

/* ---- PROTOTYPES ---- */

/* display_mat: displays a Row x Col matrix Mat in the window Gid,
 * using a simple colour-coded representation (blue for Lowval, red for
 * Upval etc.) The window may be resized, moved or iconified.
 * The window is closed when <Esc> is pressed.
 */
void display_mat(long Gid, double **Mat, unsigned int Row, unsigned int Col, 
		    double Lowval, double Upval);

/* init_matplot: initialises distmat plotting. The matrix size should
 * be Row x Col (which is the min. size in pixels). Title is put on the top
 * of the window by the window manager. Xorig and Yorig specify
 * the lower left corner position of the window: if either of them
 * is negative, then the window will be positioned interactively.
 * Return value: the 'gid' window ID number.
 */
long init_matplot(unsigned int Row, unsigned int Col, char *Title, 
	long Xorig, long Yorig);

/* plot_mat: creates a simple colour-coded dot representation of Mat,
 * which is assumed to be a Row x Col matrix. The colours vary from
 * blue for Lowval values to red for Upval values. Values lower than
 * Lowval will be black, higher than Upval will be white.
 * Plotting is performed in window no. Gid.
 * Call winclose() to dispose of this window.
 */
void plot_mat(long Gid, double **Mat, unsigned int Row, unsigned int Col, 
		    double Lowval, double Upval);

/* ==== END OF HEADER matplot.h ==== */
#endif

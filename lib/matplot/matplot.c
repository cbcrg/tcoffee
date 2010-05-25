/* ==== FUNCTIONS matplot.c ==== */

/* GL program for visualising general rectangular MxN matrices. */

/* ANSI C+SGI GL, IRIX 5.3, 22. Sept. 1995. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <math.h>
#include <gl/gl.h>
#include <gl/device.h>

/* ---- MODULE HEADER ---- */

#include "matplot.h"

/* ---- DEFINITIONS ---- */

#define RGB_BLACK 0xFF000000L

/* ---- FILE_SCOPE VARS ---- */

static int Dblbuffer;

/* ---- PROTOTYPES ---- */

static unsigned long rainbow_ramp(double X, double Lowval, double Upval);
static int round_id(double X);

/* ==== FUNCTIONS ==== */

/* ---- Initialisation ---- */

/* init_matplot: initialises distmat plotting. The matrix size should
 * be Row x Col (which is the min. size in pixels). Title is put on the top
 * of the window by the window manager. Xorig and Yorig specify
 * the lower left corner position of the window: if either of them
 * is negative, then the window will be positioned interactively.
 * Return value: the 'gid' window ID number.
 */
long init_matplot(unsigned int Row, unsigned int Col, char *Title, 
	long Xorig, long Yorig)
{
    const unsigned int MINSIZE=500;
    
    long Gid, Xsize, Ysize, Xmax, Ymax;
    float Xmagnif, Ymagnif;
    unsigned int Rs, Cs;
    
    if (!Row) Row=MINSIZE;
    if (!Col) Col=MINSIZE;
    
    Xmagnif=(float)Col/MINSIZE;
    Ymagnif=(float)Row/MINSIZE;
    Rs=Row; Cs=Col;
    if (Col<=Row)
    {
	if (Xmagnif<1.0)
	{
	    Cs=MINSIZE;
	    Rs=round_id(Row/Xmagnif);
	}
    }
    else
    {
	if (Ymagnif<1.0)
	{
	    Rs=MINSIZE;
	    Cs=round_id(Col/Ymagnif);
	}
    }
    
    foreground();	/* enable signal catch etc. */
    Xmax=getgdesc(GD_XPMAX); Ymax=getgdesc(GD_YPMAX);
    keepaspect(Col, Row);

    /* Xmax,Ymax: the maximal screen coordinates */
    Xmax=getgdesc(GD_XPMAX); Ymax=getgdesc(GD_YPMAX);
    if (Xorig+Cs>Xmax) Xorig=Xmax-Cs;
    if (Yorig+Rs>Ymax) Yorig=Ymax-Rs;
    if (Xorig>=0 && Yorig>=0)
	prefposition(Xorig, Xorig+Cs, Yorig, Yorig+Rs);
    iconsize(84, 67);
    Gid=winopen(Title);	/* create window */
    
    RGBmode();
    /* check if double buffering is available */
    if (Dblbuffer=getgdesc(GD_BITS_NORM_DBL_BLUE))
    { doublebuffer(); gconfig(); }
    else fputs("init_matplot: single-buffer mode\n", stderr);

    /* enable resize */
    winconstraints();
    keepaspect(Col, Row);
    winconstraints();
    
    getsize(&Xsize, &Ysize);    /* scale drawing into window */
    Xmagnif=(float)Xsize/Col;
    Ymagnif=(float)Ysize/Row;
    pushmatrix();
    scale(Xmagnif, Ymagnif, 1.0);
    cpack(RGB_BLACK); clear();          /* clears window to black */
    if (Dblbuffer) { swapbuffers(); cpack(RGB_BLACK); clear(); }
    
    /* queue input events */
    qdevice(ESCKEY); qdevice(WINFREEZE); qdevice(WINTHAW);
    qdevice(REDRAWICONIC); qdevice(WINQUIT);
    
    return(Gid);
}
/* END of init_matplot */

/* ---- Plotting ---- */

/* display_mat: displays a Row x Col matrix Mat in the window Gid,
 * using a simple colour-coded representation (blue for Lowval, red for
 * Upval etc.) The window may be resized, moved or iconified.
 * The window is closed when <Esc> is pressed.
 */
void display_mat(long Gid, double **Mat, unsigned int Row, unsigned int Col, 
		    double Lowval, double Upval)
{
    Device Dev;
    short Val;
    
    plot_mat(Gid, Mat, Row, Col, Lowval, Upval);
    while(1)
    {
	do  /* read the event queue until empty */
	{
	    Dev=qread(&Val);
	    switch(Dev)
	    {
		case REDRAW:	/* redraw if necessary */
		case REDRAWICONIC:
		case WINFREEZE:
		case WINTHAW:
		    plot_mat(Gid, Mat, Row, Col, Lowval, Upval);
		break;
		
		case ESCKEY:	/* both close up shop */
		case WINQUIT:
		    winclose(Gid); return;
		break;
		
		default: break;
	    }
	} while (qtest());
    }	    /* while(1) */
}
/* END of display_mat */

/* plot_mat: creates a simple colour-coded dot representation of Mat,
 * which is assumed to be a Row x Col matrix. The colours vary from
 * blue for Lowval values to red for Upval values. Values lower than
 * Lowval will be black, higher than Upval will be white.
 * Plotting is performed in window no. Gid.
 * Call winclose() to dispose of this window.
 */
void plot_mat(long Gid, double **Mat, unsigned int Row, unsigned int Col, 
		    double Lowval, double Upval)
{
    register unsigned int i,j;
    
    winset(Gid);    /* make this the current window */
   
    reshapeviewport();
    cpack(RGB_BLACK);
    clear();
    
    for (i=0; i<Row; i++)           /* scans Mat */
	for (j=0; j<Col; j++)
	{
	    /* select code color for Mt values */
	    cpack(rainbow_ramp(Mat[i][j], Lowval, Upval));

	    /* plot the point */
	    rectfi(j, Row-i, j+1, Row-i-1);
 	}

    /* flush graphics */
    if (Dblbuffer) swapbuffers(); else gflush();
}
/* END of plot_mat */

/* ---- Colour ---- */

/* rainbow_ramp: converts a double value X to a colour in RGB
 * representation according to the following approx. scheme:
 * 
 * Colour  <--Black  Blue  Cyan   Green  Yellow  Red  White-->
 *                    |	    |       |       |      |
 * X (relative)       0 ...1/4.....1/2.....3/4.....1
 * X (absolute)      Lowval.....................Upval
 * Return value: a 32-bit integer that contains the RGB info
 * in the SGI order: MSB:alpha, blue, green, red:LSB where
 * the alpha byte is always FF (maximal transparency).
 * If Upval==Lowval, then black is returned.
 */
static unsigned long rainbow_ramp(double X, double Lowval, double Upval)
{
    const unsigned long PURE_WHITE=0xFFFFFFFFL, 
			PURE_BLACK=0xFF000000L, 
			GREEN_SHIFT=0x00000100L, 
			BLUE_SHIFT=0x00010000L;

    unsigned long Col=PURE_BLACK;
    
    if (Upval==Lowval) return(PURE_BLACK);
    
    /* rescale: Lowval=0, Upval=1 */
    X=(X-Lowval)/(Upval-Lowval);
    if (X<0) return(PURE_BLACK);    /* outside range */
    if (X>1) return(PURE_WHITE);
    
    /* find ranges, set colours */
    if (X<=0.25)    /* blue...cyan: full blue+green grows */
	Col+=255*BLUE_SHIFT+(char)round_id(1020*X)*GREEN_SHIFT;
    else if (X>0.25 && X<=0.5)	/* cyan..green: full green+blue shrinks */
	Col+=255*GREEN_SHIFT+(char)round_id(1020*(0.5-X))*BLUE_SHIFT;
    else if (X>0.5 && X<=0.75)	/* green..yellow: full green+red grows */
	Col+=255*GREEN_SHIFT+(char)round_id(1020*(X-0.5));
    else			/* yellow..red: full red+green shrinks */
	Col+=255+(char)round_id(1020*(1.0-X))*GREEN_SHIFT;

    return(Col);
}
/* END of rainbow_ramp */

/* round_id: rounds a double to the nearest integer. Does the
 * same job as (int)rint(X) but rint(),  for some reason,  is not 
 * pure ANSI. This is a patch (it took me 3 hrs to realise what
 * went wrong when called rint() in a pure ANSI environment... 
 */
static int round_id(double X)
{
    double Floor, Ceil;
    
    Floor=floor(X); Ceil=ceil(X);
    if (X-Floor<Ceil-X) return((int)Floor);
    else return((int)Ceil);
}
/* END of round_id */

/* ==== END OF FUNCTIONS matplot.c ==== */

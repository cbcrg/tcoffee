#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <cstdint>  // For uint32_t
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
static float dm[]={
0.178091,
1.18065, 0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 0.131916, 0.012391, 0.022599, 0.020358, 0.030727, 0.015315, 0.048298, 0.053803, 0.020662, 0.023612, 0.216147, 0.147226, 0.065438, 0.003758, 0.009621,
0.056591,
1.35583, 0.021465, 0.0103, 0.011741, 0.010883, 0.385651, 0.016416, 0.076196, 0.035329, 0.013921, 0.093517, 0.022034, 0.028593, 0.013086, 0.023011, 0.018866, 0.029156, 0.018153, 0.0361, 0.07177, 0.419641,
0.0960191,
6.66436 ,0.561459, 0.045448, 0.438366, 0.764167, 0.087364, 0.259114, 0.21494, 0.145928, 0.762204, 0.24732, 0.118662, 0.441564, 0.174822, 0.53084, 0.465529, 0.583402, 0.445586, 0.22705, 0.02951, 0.12109,
0.078123,
2.08141, 0.070143, 0.01114, 0.019479, 0.094657, 0.013162, 0.048038, 0.077, 0.032939, 0.576639, 0.072293, 0.02824, 0.080372, 0.037661, 0.185037, 0.506783, 0.073732, 0.071587, 0.042532, 0.011254, 0.028723,
0.0834977,
2.08101, 0.041103, 0.014794, 0.00561, 0.010216, 0.153602, 0.007797, 0.007175, 0.299635, 0.010849, 0.999446, 0.210189, 0.006127, 0.013021, 0.019798, 0.014509, 0.012049, 0.035799, 0.180085, 0.012744, 0.026466,
0.0904123,
2.56819, 0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 0.017255, 0.004911, 0.796882, 0.017074, 0.285858, 0.075811, 0.014548, 0.015092, 0.011382, 0.012696, 0.027535, 0.088333, 0.94434, 0.004373, 0.016741,
0.114468,
1.76606, 0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 0.105877, 0.049776, 0.014963, 0.094276, 0.027761, 0.01004, 0.187869, 0.050018, 0.110039, 0.038668, 0.119471, 0.065802, 0.02543, 0.003215, 0.018742,
0.0682132,
4.98768, 0.452171, 0.114613, 0.06246, 0.115702, 0.284246, 0.140204, 0.100358, 0.55023, 0.143995, 0.700649, 0.27658, 0.118569, 0.09747, 0.126673, 0.143634, 0.278983, 0.358482, 0.66175, 0.061533, 0.199373,
0.234585,
0.0995, 0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 0.016931, 0.003647, 0.002184, 0.005019, 0.00599, 0.001473, 0.004158, 0.009055, 0.00363, 0.006583, 0.003172, 0.00369, 0.002967, 0.002772,0.002686};

double int_logB (int *i, int n)
	{
	static double *array;
	int a;
	
	if ( array==NULL)array=(double*)vcalloc ( 1000, sizeof (double));
	
	for ( a=0; a< n; a++)
		array[a]=(double)i[a];
	return double_logB(array, n);
	}
double float_logB (float *i, int n)
	{
	static double *array;
	int a;
	 
	if ( array==NULL)array=(double*)vcalloc ( 1000, sizeof (double));
	for ( a=0; a< n; a++)
		array[a]=(double)i[a];
	return double_logB(array, n);
	}
	  
double double_logB (double *x, int n)
	{
	double vx=0;
	double result=0;
	int i;
	
	
	for ( i=0; i<n; i++)vx+=x[i];
	for ( i=0; i<n; i++)result+=lgamma2(x[i]);
	return 	result-lgamma2(vx);
	} 	
double *** make_lup_table ( Mixture *D)
	{
	int a, b, c;
	double ***lup;
	
	lup=(double***)vcalloc ( 9, sizeof (double**));
	for ( a=0; a<9; a++)
		{
		lup[a]=(double**)vcalloc ( 20, sizeof (double*));
		for ( b=0; b< 20; b++)
			lup[a][b]=(double*)vcalloc ( 200, sizeof (double));
		}
	
	for ( a=0; a< 9; a++)
		for ( b=0; b< 20; b++)
			for ( c=0; c< 100; c++)
				lup[a][b][c]=lgamma2(D->ALPHA[a][b]+c);
	
	return lup;
	}
	
double  double_logB2(int j, double *n,Mixture *D)
	{
	double vx=0;
	double result=0;
	int i;
	
	static double ***lup;
	
	
	
	if ( lup==NULL)lup=make_lup_table (D);
	
       

	for ( i=0; i<D->n_aa; i++)vx+=(double)n[i]+D->ALPHA[j][i];
	for ( i=0; i<D->n_aa; i++)
	  {
	    
	    
	    result+=lup[j][i][(int)n[i]];
	  }
	return 	result-lgamma2(vx);
	} 	
			
double compute_exponant ( double *n, int j, Mixture *D)
	{
	
	if ( j>=9)fprintf ( stderr, "\nPB: j=%d", j);
	
	return double_logB2(j, n,D)-D->double_logB_alpha[j];
	}


double *compute_matrix_p ( double *n)
        {
	  
	  /*
	    reads in a frquency list of various amino acids:
	    
	    sum freq(aa)=1 (gaps are ignored)
	    
	    aa[1]=x1
	    aa[2]=x2
	    ....

	    Outputs a similar list with frequencies 'Blurred' using a pam250 mt
	  */

	    
	  
	  static int **matrix;
	  double *R;
	  int a, b;
	  double v,min, tot;
	  
	  
	  if ( !matrix) 
	    {
	      matrix=read_matrice ( "pam250mt");	  
	    }
	  
	  R=(double*)vcalloc ( 26, sizeof (double));


	  for ( a=0; a<26; a++)
	    {
	      if (!is_aa(a+'a'))continue;
	      if ( n[a]==0)continue;
	      
	      for ( b=0; b< 26; b++)
		{
		  if (!is_aa(b+'a'))continue;
		  v=n[a]*(matrix[a][b]);
		  if ( v>0)
		    {
		      R[b]+=v+(10*n[a]);
		    }
		}
	    }
	  
	  min=R[0];
	  for ( min=R[0],a=0; a< 26; a++)min=MIN(min,R[a]);
	  for ( tot=0,   a=0; a< 26; a++)         {R[a]-=min;tot+=R[a];}
	  for ( a=0; a< 26; a++)if ( is_aa(a+'a')){R[a]=R[a]*((float)(100)/(float)tot);}
	  return R;
	}
	      

double *compute_dirichlet_p ( double *n)
	{
	  /*
	    Given a list of frequenceies measured for the residues, this function returns 
	    the p_values associated with each residue in the column
	  */
	 
	int a, b;
	double X_LIST[100];
	double sum, log_sum, max;
	static Mixture *D;
	static double *R;



	if (!D)
		{
		D=read_dirichlet (NULL);
		
		D->n_aa=20;
		R=(double*)vcalloc ( D->n_aa, sizeof (double));
		D->double_logB_alpha=(double*)vcalloc (D->N_COMPONENT , sizeof (double));
		
		D->exponant_list=(double*)vcalloc (D->N_COMPONENT , sizeof (double));
		precompute_log_B ( D->double_logB_alpha,D);
		D->alpha_tot=(double*)vcalloc (D->N_COMPONENT , sizeof (double));
		for ( a=0; a<D->N_COMPONENT; a++)
			for ( b=0; b< D->n_aa; b++)
				D->alpha_tot[a]+=D->ALPHA[a][b];
		}
	
	

	for ( D->tot_n=0,a=0; a< D->n_aa; a++)D->tot_n+=(double)n[a];
	max=D->exponant_list[0]=compute_exponant ( n, 0, D);	
	for ( a=1; a<D->N_COMPONENT; a++)
		{
		D->exponant_list[a]=compute_exponant ( n, a,D);
		max= ( max< D->exponant_list[a])?D->exponant_list[a]:max;
		}
	for ( a=1; a<D->N_COMPONENT; a++)D->exponant_list[a]=D->exponant_list[a]-max;
	
	
	for ( sum=0,log_sum=0,a=0; a< D->n_aa; a++)
		{
		sum+=X_LIST[a]=compute_X (n, a,D);
		}
	log_sum=log(sum);

		
	for (a=0; a<D->n_aa; a++)
		{
		R[a]=(log(X_LIST[a])-log_sum);
		}
	
	
	/*
	printf ( "\n[");
	for ( a=0;a< n_aa; a++)printf ("%d ", n[a]);
	printf ("] score=%f",(float) result );
	
	fprintf ( stderr, "\nRESULT=%f", (float)result);
	exit(0);
	*/
	return R;

	}
	
void precompute_log_B ( double *table,Mixture *D)
	{
	int a;
	for ( a=0; a< D->N_COMPONENT; a++)
		{
		table[a]=double_logB ( D->ALPHA[a], D->n_aa);
		}
	}			
double compute_X (double *n,int i,Mixture *D)
	{
	int  j;
	double term1, term2,result;
	double **alpha;
	double *q;

	
	
	alpha=D->ALPHA;
	q=D->DM_Q;

	for (result=0, j=0; j<D->N_COMPONENT; j++)
		{
		term1=exp (D->exponant_list[j])*q[j];
		term2=(alpha[j][i]+(double)n[i])/(D->alpha_tot[j]+D->tot_n);
		result+=term1*term2;
		}
	return result;
	}
Mixture * read_dirichlet ( char *name)
	{
	FILE *fp;
	int a,b, c;
	float f;	
	Mixture *D;


	D=(Mixture*)vcalloc ( 1, sizeof (Mixture));
	
	
	D->N_COMPONENT=9;
	D->ALPHA= (double**)vcalloc (9, sizeof (double*));
	for ( a=0; a< 9; a++)
		D->ALPHA[a]= (double*)vcalloc (20, sizeof (double));
	D->DM_Q= (double*)vcalloc (9, sizeof (double));
	
	if (name!=NULL)
	  {
	    fp=vfopen ( name, "r");
	    for ( a=0; a< 9; a++)
		{
		fscanf(fp, "%f\n", &f);
		D->DM_Q[a]=(double)f;
		fscanf(fp, "%f", &f);
		
		for ( b=0; b<20; b++)
			{
			fscanf(fp, "%f", &f);
			D->ALPHA[a][b]=(double)f;
			}
		fscanf(fp, "\n");
		}
	    for ( a=0; a< 9; a++)
	      {
		fprintf(stderr, "\n%f\n",(float)D->DM_Q[a] );
		
		for ( b=0; b<20; b++)
		  {
		    fprintf(stderr, "%f ", (float)D->ALPHA[a][b]);
		  }
		fprintf(stderr, "\n");
	      }
	    fprintf ( stderr, "\nN_C=%d",D->N_COMPONENT );	
	    vfclose ( fp);
	  }
	else
	  {
	    for (c=0, a=0; a< 9;a++)
	      {
		D->DM_Q[a]=dm[c++];	
		for (b=0; b<20; b++)
		  D->ALPHA[a][b]=dm[c++];
	      }
	  }
	
	return D;
	}
int *dirichlet_code2aa_lu ()
{
  static int *dm;

  if (dm);
  else
    {
      char aal[21];
      int a;
      
      sprintf (aal, "acdefghiklmnpqrstvwy");
      dm=(int*)vcalloc (265, sizeof (int));
      memset (dm, -1, 20*sizeof(int));
      for (a=0; a<20; a++)
	{
	  dm[dirichlet_code (aal[a])]=aal[a];
	}
    }
  return dm;
}
int *aa2dirichlet_code_lu ()
{
  static int *dm;

  if (dm);
  else
    {
      char aal[21];
      int a;
      
      sprintf (aal, "acdefghiklmnpqrstvwy");
      dm=(int*)vcalloc (265, sizeof (int));
      memset (dm, -1, 265*sizeof(int));
      for (a=0; a<20; a++)
	{
	  dm[aal[a]]=dm[(aal[a]-'a')+'A']=dm[(aal[a]-'a')]=dirichlet_code (aal[a]);
	}
    }
  return dm;
}
int dirichlet_code( char aa)
	{
	
	char x;
	
	x=tolower (aa);
	
	if ( (x<'a') || (x>'z'))
		crash ( "CODE UNDEFINED");
	else if ( x<='a')
	    return x-'a';
	else if ( x<='i')
	    return x-('a'+1);
	else if ( x<= 'n')
	    return x-('a'+2);
	else if ( x<='t')
	    return x-('a'+3);
	else if ( x<='w')
	    return x-('a'+4);
	else if ( x=='y')
	    return x-('a'+5);
	else 
	  {
	    crash ("ERROR in dirichlet_code");
	    return 0;
	  }
	return 0;
	
	}


static const double
two52=  4.50359962737049600000e+15, /* 0x43300000, 0x00000000 */
half=  5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
one =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
pi  =  3.14159265358979311600e+00, /* 0x400921FB, 0x54442D18 */
a0  =  7.72156649015328655494e-02, /* 0x3FB3C467, 0xE37DB0C8 */
a1  =  3.22467033424113591611e-01, /* 0x3FD4A34C, 0xC4A60FAD */
a2  =  6.73523010531292681824e-02, /* 0x3FB13E00, 0x1A5562A7 */
a3  =  2.05808084325167332806e-02, /* 0x3F951322, 0xAC92547B */
a4  =  7.38555086081402883957e-03, /* 0x3F7E404F, 0xB68FEFE8 */
a5  =  2.89051383673415629091e-03, /* 0x3F67ADD8, 0xCCB7926B */
a6  =  1.19270763183362067845e-03, /* 0x3F538A94, 0x116F3F5D */
a7  =  5.10069792153511336608e-04, /* 0x3F40B6C6, 0x89B99C00 */
a8  =  2.20862790713908385557e-04, /* 0x3F2CF2EC, 0xED10E54D */
a9  =  1.08011567247583939954e-04, /* 0x3F1C5088, 0x987DFB07 */
a10 =  2.52144565451257326939e-05, /* 0x3EFA7074, 0x428CFA52 */
a11 =  4.48640949618915160150e-05, /* 0x3F07858E, 0x90A45837 */
tc  =  1.46163214496836224576e+00, /* 0x3FF762D8, 0x6356BE3F */
tf  = -1.21486290535849611461e-01, /* 0xBFBF19B9, 0xBCC38A42 */
/* tt = -(tail of tf) */
tt  = -3.63867699703950536541e-18, /* 0xBC50C7CA, 0xA48A971F */
t0  =  4.83836122723810047042e-01, /* 0x3FDEF72B, 0xC8EE38A2 */
t1  = -1.47587722994593911752e-01, /* 0xBFC2E427, 0x8DC6C509 */
t2  =  6.46249402391333854778e-02, /* 0x3FB08B42, 0x94D5419B */
t3  = -3.27885410759859649565e-02, /* 0xBFA0C9A8, 0xDF35B713 */
t4  =  1.79706750811820387126e-02, /* 0x3F9266E7, 0x970AF9EC */
t5  = -1.03142241298341437450e-02, /* 0xBF851F9F, 0xBA91EC6A */
t6  =  6.10053870246291332635e-03, /* 0x3F78FCE0, 0xE370E344 */
t7  = -3.68452016781138256760e-03, /* 0xBF6E2EFF, 0xB3E914D7 */
t8  =  2.25964780900612472250e-03, /* 0x3F6282D3, 0x2E15C915 */
t9  = -1.40346469989232843813e-03, /* 0xBF56FE8E, 0xBF2D1AF1 */
t10 =  8.81081882437654011382e-04, /* 0x3F4CDF0C, 0xEF61A8E9 */
t11 = -5.38595305356740546715e-04, /* 0xBF41A610, 0x9C73E0EC */
t12 =  3.15632070903625950361e-04, /* 0x3F34AF6D, 0x6C0EBBF7 */
t13 = -3.12754168375120860518e-04, /* 0xBF347F24, 0xECC38C38 */
t14 =  3.35529192635519073543e-04, /* 0x3F35FD3E, 0xE8C2D3F4 */
u0  = -7.72156649015328655494e-02, /* 0xBFB3C467, 0xE37DB0C8 */
u1  =  6.32827064025093366517e-01, /* 0x3FE4401E, 0x8B005DFF */
u2  =  1.45492250137234768737e+00, /* 0x3FF7475C, 0xD119BD6F */
u3  =  9.77717527963372745603e-01, /* 0x3FEF4976, 0x44EA8450 */
u4  =  2.28963728064692451092e-01, /* 0x3FCD4EAE, 0xF6010924 */
u5  =  1.33810918536787660377e-02, /* 0x3F8B678B, 0xBF2BAB09 */
v1  =  2.45597793713041134822e+00, /* 0x4003A5D7, 0xC2BD619C */
v2  =  2.12848976379893395361e+00, /* 0x40010725, 0xA42B18F5 */
v3  =  7.69285150456672783825e-01, /* 0x3FE89DFB, 0xE45050AF */
v4  =  1.04222645593369134254e-01, /* 0x3FBAAE55, 0xD6537C88 */
v5  =  3.21709242282423911810e-03, /* 0x3F6A5ABB, 0x57D0CF61 */
s0  = -7.72156649015328655494e-02, /* 0xBFB3C467, 0xE37DB0C8 */
s1  =  2.14982415960608852501e-01, /* 0x3FCB848B, 0x36E20878 */
s2  =  3.25778796408930981787e-01, /* 0x3FD4D98F, 0x4F139F59 */
s3  =  1.46350472652464452805e-01, /* 0x3FC2BB9C, 0xBEE5F2F7 */
s4  =  2.66422703033638609560e-02, /* 0x3F9B481C, 0x7E939961 */
s5  =  1.84028451407337715652e-03, /* 0x3F5E26B6, 0x7368F239 */
s6  =  3.19475326584100867617e-05, /* 0x3F00BFEC, 0xDD17E945 */
r1  =  1.39200533467621045958e+00, /* 0x3FF645A7, 0x62C4AB74 */
r2  =  7.21935547567138069525e-01, /* 0x3FE71A18, 0x93D3DCDC */
r3  =  1.71933865632803078993e-01, /* 0x3FC601ED, 0xCCFBDF27 */
r4  =  1.86459191715652901344e-02, /* 0x3F9317EA, 0x742ED475 */
r5  =  7.77942496381893596434e-04, /* 0x3F497DDA, 0xCA41A95B */
r6  =  7.32668430744625636189e-06, /* 0x3EDEBAF7, 0xA5B38140 */
w0  =  4.18938533204672725052e-01, /* 0x3FDACFE3, 0x90C97D69 */
w1  =  8.33333333333329678849e-02, /* 0x3FB55555, 0x5555553B */
w2  = -2.77777777728775536470e-03, /* 0xBF66C16C, 0x16B02E5C */
w3  =  7.93650558643019558500e-04, /* 0x3F4A019F, 0x98CF38B6 */
w4  = -5.95187557450339963135e-04, /* 0xBF4380CB, 0x8C0FE741 */
w5  =  8.36339918996282139126e-04, /* 0x3F4B67BA, 0x4CDAD5D1 */
w6  = -1.63092934096575273989e-03; /* 0xBF5AB89D, 0x0B9E43E4 */

static const double zero=  0.00000000000000000000e+00;
static double sin_pi(double x) {
    double y, z;
    int n;
    uint32_t ix;

    // Copy the high 32 bits of x into ix safely
    unsigned long long ull;
    memcpy(&ull, &x, sizeof(x));
    ix = static_cast<uint32_t>(ull >> 32) & 0x7FFFFFFF;

    if (ix < 0x3FD00000) return sin(pi * x);
    y = -x; // x is assumed negative

    // Argument reduction
    z = floor(y);
    if (z != y) { // inexact anyway
        y *= 0.5;
        y = 2.0 * (y - floor(y)); // y = |x| mod 2.0
        n = static_cast<int>(y * 4.0);
    } else {
        if (ix >= 0x43400000) {
            y = zero;
            n = 0; // y must be even
        } else {
            if (ix < 0x43300000) z = y + two52; // exact
            memcpy(&ull, &z, sizeof(z)); // Safe bit manipulation
            n = ull & 1;
            y = n;
            n <<= 2;
        }
    }

    switch (n) {
        case 0: y = sin(pi * y); break;
        case 1:
        case 2: y = cos(pi * (0.5 - y)); break;
        case 3:
        case 4: y = sin(pi * (one - y)); break;
        case 5:
        case 6: y = -cos(pi * (y - 1.5)); break;
        default: y = sin(pi * (y - 2.0)); break;
    }
    return -y;
}
#ifdef OLD_SIN_PI
static double sin_pi(double x)
{
        double y,z;
        int n,ix;
	unsigned long long ull;
	memcpy(&ull, &x, sizeof(x));
	
        ix=(*(long long *)&x)>>32;
        ix &= 0x7fffffff;

        if(ix<0x3fd00000) return sin(pi*x);
        y = -x;                /* x is assume negative */

     /*
      * argument reduction, make sure inexact flag not raised if input
      * is an integer
      */
        z = floor(y);
        if(z!=y) {                                /* inexact anyway */
            y  *= 0.5;
            y   = 2.0*(y - floor(y));                /* y = |x| mod 2.0 */
            n   = (int) (y*4.0);
        } else {
             if(ix>=0x43400000) {
                 y = zero; n = 0;                 /* y must be even */
             } else {
                 if(ix<0x43300000) z = y+two52;        /* exact */
                        n=(*(long long *)&x);
                n &= 1;
                 y  = n;
                 n<<= 2;
             }
         }
        switch (n) {
            case 0:   y =  sin(pi*y); break;
            case 1:
            case 2:   y =  cos(pi*(0.5-y)); break;
            case 3:
            case 4:   y =  sin(pi*(one-y)); break;
            case 5:
            case 6:   y = -cos(pi*(y-1.5)); break;
            default:  y =  sin(pi*(y-2.0)); break;
            }
        return -y;
}
#endif

double lgamma2 ( double x)
{
  int s;
  return lgamma_r ( x, &s);
}
double lgamma_r(double x, int *signgamp)
{
        double t, y, z, nadj = 0, p, p1, p2, p3, q, r, w;
	int i, ix, hx, lx;
	
	unsigned long long ull;
	memcpy(&ull, &x, sizeof(x));
	
	// Extracting higher and lower parts of the double's bit representation
	hx = ull >> 32; // Higher 32 bits
	lx = ull & 0xFFFFFFFF; // Lower 32 bits
	
        *signgamp = 1;
        ix = hx&0x7fffffff;
        if(ix>=0x7ff00000) return x*x;
        if((ix|lx)==0) return one/fabs(x);
        if(ix<0x3b900000) {        /* |x|<2**-70, return -log(|x|) */
            if(hx<0) {
                *signgamp = -1;
                return -log(-x);
            } else return -log(x);
        }
        if(hx<0) {
            if(ix>=0x43300000)         /* |x|>=2**52, must be -integer */
                return x/zero;
            t = sin_pi(x);
            if(t==zero) return one/fabs(t); /* -integer */
            nadj = log(pi/fabs(t*x));
            if(t<zero) *signgamp = -1;
            x = -x;
        }

     /* purge off 1 and 2 */
        if((((ix-0x3ff00000)|lx)==0)||(((ix-0x40000000)|lx)==0)) r = 0;
     /* for x < 2.0 */
        else if(ix<0x40000000) {
            if(ix<=0x3feccccc) {         /* lgamma(x) = lgamma(x+1)-log(x) */
                r = -log(x);
                if(ix>=0x3FE76944) {y = one-x; i= 0;}
                else if(ix>=0x3FCDA661) {y= x-(tc-one); i=1;}
                  else {y = x; i=2;}
            } else {
                  r = zero;
                if(ix>=0x3FFBB4C3) {y=2.0-x;i=0;} /* [1.7316,2] */
                else if(ix>=0x3FF3B4C4) {y=x-tc;i=1;} /* [1.23,1.73] */
                else {y=x-one;i=2;}
            }
            switch(i) {
              case 0:
                z = y*y;
                p1 = a0+z*(a2+z*(a4+z*(a6+z*(a8+z*a10))));
                p2 = z*(a1+z*(a3+z*(a5+z*(a7+z*(a9+z*a11)))));
                p  = y*p1+p2;
                r  += (p-0.5*y); break;
              case 1:
                z = y*y;
                w = z*y;
                p1 = t0+w*(t3+w*(t6+w*(t9 +w*t12)));        /* parallel comp */
                p2 = t1+w*(t4+w*(t7+w*(t10+w*t13)));
                p3 = t2+w*(t5+w*(t8+w*(t11+w*t14)));
                p  = z*p1-(tt-w*(p2+y*p3));
                r += (tf + p); break;
              case 2:
                p1 = y*(u0+y*(u1+y*(u2+y*(u3+y*(u4+y*u5)))));
                p2 = one+y*(v1+y*(v2+y*(v3+y*(v4+y*v5))));
                r += (-0.5*y + p1/p2);
            }
        }
        else if(ix<0x40200000) {                         /* x < 8.0 */
            i = (int)x;
            t = zero;
            y = x-(double)i;
            p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
            q = one+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
            r = half*y+p/q;
            z = one;        /* lgamma(1+s) = log(s) + lgamma(s) */
            switch(i) {
            case 7: z *= (y+6.0);        /* FALLTHRU */
            case 6: z *= (y+5.0);        /* FALLTHRU */
            case 5: z *= (y+4.0);        /* FALLTHRU */
            case 4: z *= (y+3.0);        /* FALLTHRU */
            case 3: z *= (y+2.0);        /* FALLTHRU */
                    r += log(z); break;
            }
     /* 8.0 <= x < 2**58 */
        } else if (ix < 0x43900000) {
            t = log(x);
            z = one/x;
            y = z*z;
            w = w0+z*(w1+y*(w2+y*(w3+y*(w4+y*(w5+y*w6)))));
            r = (x-half)*(t-one)+w;
        } else
     /* 2**58 <= x <= inf */
            r =  x*(log(x)-one);
        if(hx<0) r = nadj - r;
        return r;
}
#ifdef OLD_LGAMMA
double lgamma_r(double x, int *signgamp)
{
        double t,y,z,nadj=0,p,p1,p2,p3,q,r,w;
        int i,hx,lx,ix;

        hx=(*(long long *)&x)>>32;
        lx=(*(long long *)&x);

     /* purge off +-inf, NaN, +-0, and negative arguments */
        *signgamp = 1;
        ix = hx&0x7fffffff;
        if(ix>=0x7ff00000) return x*x;
        if((ix|lx)==0) return one/fabs(x);
        if(ix<0x3b900000) {        /* |x|<2**-70, return -log(|x|) */
            if(hx<0) {
                *signgamp = -1;
                return -log(-x);
            } else return -log(x);
        }
        if(hx<0) {
            if(ix>=0x43300000)         /* |x|>=2**52, must be -integer */
                return x/zero;
            t = sin_pi(x);
            if(t==zero) return one/fabs(t); /* -integer */
            nadj = log(pi/fabs(t*x));
            if(t<zero) *signgamp = -1;
            x = -x;
        }

     /* purge off 1 and 2 */
        if((((ix-0x3ff00000)|lx)==0)||(((ix-0x40000000)|lx)==0)) r = 0;
     /* for x < 2.0 */
        else if(ix<0x40000000) {
            if(ix<=0x3feccccc) {         /* lgamma(x) = lgamma(x+1)-log(x) */
                r = -log(x);
                if(ix>=0x3FE76944) {y = one-x; i= 0;}
                else if(ix>=0x3FCDA661) {y= x-(tc-one); i=1;}
                  else {y = x; i=2;}
            } else {
                  r = zero;
                if(ix>=0x3FFBB4C3) {y=2.0-x;i=0;} /* [1.7316,2] */
                else if(ix>=0x3FF3B4C4) {y=x-tc;i=1;} /* [1.23,1.73] */
                else {y=x-one;i=2;}
            }
            switch(i) {
              case 0:
                z = y*y;
                p1 = a0+z*(a2+z*(a4+z*(a6+z*(a8+z*a10))));
                p2 = z*(a1+z*(a3+z*(a5+z*(a7+z*(a9+z*a11)))));
                p  = y*p1+p2;
                r  += (p-0.5*y); break;
              case 1:
                z = y*y;
                w = z*y;
                p1 = t0+w*(t3+w*(t6+w*(t9 +w*t12)));        /* parallel comp */
                p2 = t1+w*(t4+w*(t7+w*(t10+w*t13)));
                p3 = t2+w*(t5+w*(t8+w*(t11+w*t14)));
                p  = z*p1-(tt-w*(p2+y*p3));
                r += (tf + p); break;
              case 2:
                p1 = y*(u0+y*(u1+y*(u2+y*(u3+y*(u4+y*u5)))));
                p2 = one+y*(v1+y*(v2+y*(v3+y*(v4+y*v5))));
                r += (-0.5*y + p1/p2);
            }
        }
        else if(ix<0x40200000) {                         /* x < 8.0 */
            i = (int)x;
            t = zero;
            y = x-(double)i;
            p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
            q = one+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
            r = half*y+p/q;
            z = one;        /* lgamma(1+s) = log(s) + lgamma(s) */
            switch(i) {
            case 7: z *= (y+6.0);        /* FALLTHRU */
            case 6: z *= (y+5.0);        /* FALLTHRU */
            case 5: z *= (y+4.0);        /* FALLTHRU */
            case 4: z *= (y+3.0);        /* FALLTHRU */
            case 3: z *= (y+2.0);        /* FALLTHRU */
                    r += log(z); break;
            }
     /* 8.0 <= x < 2**58 */
        } else if (ix < 0x43900000) {
            t = log(x);
            z = one/x;
            y = z*z;
            w = w0+z*(w1+y*(w2+y*(w3+y*(w4+y*(w5+y*w6)))));
            r = (x-half)*(t-one)+w;
        } else
     /* 2**58 <= x <= inf */
            r =  x*(log(x)-one);
        if(hx<0) r = nadj - r;
        return r;
}
#endif

double **prf2dmx (double **in, double **prf, int len)
{
  int c,r;

  if (!prf)prf=declare_double (len, 20);
  else if (read_array_size_new(prf)<len){free_double (prf,-1); prf=declare_double (len, 20);}

  for (c=0; c<len; c++)
    {
      memcpy  (prf[c],compute_dirichlet_p(in[c]), 20*sizeof (double));
    }
  return prf;
}
double **aln2prf (Alignment *A, int ns, int *ls, int len, double **prf)
{
  
  int free_ls=0;
  int *lu;
  int c,r,d;
  double tot;
  double prior=0.1;
  
  
  if (!A)return prf;
  if (ns==0 || ls==0 || len==0)
    {
      ns=A->nseq;
      ls=(int*)vcalloc (ns, sizeof (int));
      for (r=0; r<ns; r++)ls[r]=r;
      free_ls=1;
    }
  
  lu=aa2dirichlet_code_lu();
  
  if (!prf)prf=declare_double (len, 20);
  else if (read_array_size_new(prf)<len){free_double (prf,-1); prf=declare_double (len, 20);}

  for (c=0; c<len; c++)
    {
      double *dmr;
      memset (prf[c], 0, sizeof (double)*20);
      for (tot=0,r=0; r<ns; r++)
	{
	  d=lu[A->seq_al[ls[r]][c]];
	  if (d>=0)
	    {
	      prf[c][d]++;
	      tot++;
	    }
	}
      if (tot>0)for (r=0;r<20; r++)prf[c][r]/=tot;
    }
  if (free_ls)vfree(ls);
  return prf;
}

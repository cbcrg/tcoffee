/*UTIL MACROS*/
#define SWAP(x,y) {x=x+y;y=x+y; x=y-x; y=y-2*x;}
#define SWAPP(x,y,tp) {tp=y;y=x;x=tp;}

#define MAX(x, y) (((x) >(y)) ? (x):(y))
#define MAX3(x,y,z) (MAX(MAX(x,y),z))
#define MAX4(a,b,c,d) (MAX(MAX(a,b),MAX(c,d)))
#define MIN(x, y) (((x) <(y)) ? (x):(y))
#define FABS(x) (x<0)?(-x):(x)
#define a_better_than_b(x,y,m) ((m==1)?(((x)>(y))?1:0):(((x)<(y))?1:0))
#define best_of_a_b(x,y,m)   ((m==1)?(MAX((x),(y))):(MIN((x),(y))))
#define strm(x,y)            ((strcmp((x),(y))==0)?1:0)
#define strnm(x,y,n)           ((strncmp((x),(y),(n))==0)?1:0)
#define strm2(a,b,c)         (strm(a,b) || strm(a,c))
#define strm3(a,b,c,d)       (strm2(a,b,c) || strm(a,d))
#define strm4(a,b,c,d,e)     (strm2(a,b,c) || strm2(a,d,e))
#define strm5(a,b,c,d,e,f)   (strm2(a,b,c) || strm3(a,d,e,f))
#define strm6(a,b,c,d,e,f,g) (strm3(a,b,c,d) || strm3(a,e,f,g))
#define declare_name(x) x=vcalloc (MAX(FILENAMELEN,L_tmpnam)+1, sizeof (char)) 
#define is_parameter(x) (x[0]=='-' && !isdigit(x[1])) 

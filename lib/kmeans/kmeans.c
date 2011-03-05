/*****
** kmeans.c
** - a simple k-means clustering routine
** - returns the cluster labels of the data points in an array
** - here's an example
**   extern int *k_means(double**, int, int, int, double, double**);
**   ...
**   int *c = k_means(data_points, num_points, dim, 20, 1e-4, 0);
**   for (i = 0; i < num_points; i++) {
**      printf("data point %d is in cluster %d\n", i, c[i]);
**   }
**   ...
**   free(c);
** Parameters
** - array of data points (double **data)
** - number of data points (int n)
** - dimension (int m)
** - desired number of clusters (int k)
** - error tolerance (double t)
**   - used as the stopping criterion, i.e. when the sum of
**     squared euclidean distance (standard error for k-means)
**     of an iteration is within the tolerable range from that
**     of the previous iteration, the clusters are considered
**     "stable", and the function returns
**   - a suggested value would be 0.0001
** - output address for the final centroids (double **centroids)
**   - user must make sure the memory is properly allocated, or
**     pass the null pointer if not interested in the centroids
** References
** - J. MacQueen, "Some methods for classification and analysis
**   of multivariate observations", Fifth Berkeley Symposium on
**   Math Statistics and Probability, 281-297, 1967.
** - I.S. Dhillon and D.S. Modha, "A data-clustering algorithm
**   on distributed memory multiprocessors",
**   Large-Scale Parallel Data Mining, 245-260, 1999.
** Notes
** - this function is provided as is with no warranty.
** - the author is not responsible for any damage caused
**   either directly or indirectly by using this function.
** - anybody is free to do whatever he/she wants with this
**   function as long as this header section is preserved.
** Created on 2005-04-12 by
** - Roger Zhang ([[Email Removed]])
** Modifications
** -
** Last compiled under Linux with gcc-3
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
int  k_means(double **data, int n, int m, int k, double t, double **centroids);
float mk_means(double **d, int n, int dim, int k,double t, double **centroids, int nrounds);
double ** read_data ( char *file, int n, int dim, int len, char **field_list);
int       file2dim     ( char *file, int *n,int *dim, int *len);
double** shuffle_data (double **d, double **sd, int n, int r);
void display_data (double **d, int n, int dim);
double data2evaluate ( double **d, int n, int dim);

void output_data ( double **data, int n, int dim, int len,  char *infile, char *outfile);

		
int main (int argc, char *argv[])
{
  double **data, **sdata;
  
  int dim, len, n;
  int a, b;
  int k=2;
  char *infile=NULL;
  char *outfile=NULL;
  double t=0.0001;
  int mode=1;
  int scan=0;
  char **field_list;
  int nf=0;
  
  field_list=calloc (1000, sizeof (char*));
  srand(time(NULL));
  for (a=1; a<argc;a++)
    {
      if    ( strcmp (argv[a], "-k")==0){k=atoi(argv[++a]);}    //k
      else if ( strcmp (argv[a], "-i")==0){infile=argv[++a];}   //input
      else if ( strcmp (argv[a], "-o")==0){outfile=argv[++a];}  //output
      else if ( strcmp (argv[a], "-t")==0){t=(double)atof(argv[++a]);}//t
      else if ( strcmp (argv[a], "-m")==0){mode=atoi(argv[++a]);}//mode: 1 simple, n, nrounds
      else if ( strcmp (argv[a], "-f")==0){field_list[nf++]=argv[++a];}
      else
	{
	  fprintf ( stdout, "ERROR: %s: -k <nlust> -i <file> -o <file> -f <field1> -f <filed2>\n", argv[a]);
	
	}
    }
  file2dim (infile,&n, &dim, &len);
  data=read_data(infile, n, dim, len,field_list);
  
  if ( k<0)
    {
      k*=-1;
      mode=(mode<=1)?10:mode;
      
      for (a=2; a<k; a++)
	{
	  fprintf (stdout, "K=%d score=%.3f\n", a,mk_means (data, n, dim,a,t,NULL, mode));
	}
      exit (0);
    }
  else
    mk_means (data, n, dim,k,t,NULL, mode);
  
  output_data (data,n, dim,len,infile,outfile);
}

void display_data (double **d, int n, int dim)
{
  int a, b;
  
  fprintf ( stdout, "\n");
  for (a=0; a<n; a++)
    {
      for (b=0; b<dim; b++)fprintf ( stdout, "%d ", (int)d[a][b]);
      fprintf ( stdout, "\n");
    }
  
}
int   file2dim   ( char *file, int *n,int *dim, int *len)
{
  FILE *fp;
  double **data;
  char *buf;
  int c;
  int mdim=0;
  int cdim=0;
  int mlen=0;
  int clen=0;
  char *s1, *s2;
  
  fp=fopen (file, "r");
  while ((c=fgetc(fp)!=EOF))
    {
      clen=0;
      while ( (c=fgetc(fp))!='\n' && c!=EOF)clen++;
      mlen=(clen>mlen)?clen:mlen;
    }
  len[0]=mlen+10;
  close (fp);
  
  n[0]=0;
  fp=fopen (file, "r");
  buf=calloc(len[0]+1, sizeof (char));
  while ((fgets (buf,len[0], fp)))
    {
      if ( buf[0]='#')
	{
	  n[0]++;
	  cdim=0;
	  strtok(buf, ";");//pass #d;
	  strtok(NULL, ";"); //pass exp;
	  strtok(NULL, ";"); //pass recid;
	  while ((s1=strtok(NULL, ";")))
	    {
	      s2=strtok(NULL, ";");
	      if (strstr (s1, "value::"))
		{
		  cdim++;
		}
	    }
	  mdim=(cdim>mdim)?cdim:mdim;
	}
    }
  free (buf);
  n[0];
  dim[0]=mdim;
  close (fp);
  return n[0];
}
double ** read_data ( char *file, int n, int dim, int mlen, char **fl)
{
  FILE *fp;
  double **data;
  char *buf;
  char *s1;
  char *s2;
  int cdim,cn;
  int a,b,c,p;
  int *fi;
  
  fi=calloc (1000,sizeof (int));
  for (a=0; a<1000; a++)fi[a]=-1;
  buf =calloc (mlen+1,sizeof (char));
  data=calloc (n+1, sizeof (double*));
  for (a=0; a<n; a++)data[a]=calloc(dim+1, sizeof (double));
  
  fp=fopen (file, "r");
  cn=0;
  while ((fgets (buf,mlen, fp)))
    {
      
      if ( buf[0]='#')
	{
	  
	  p=cdim=0;
	  strtok(buf, ";"); p++;//pass #d;
	  strtok(NULL, ";");p++; //pass exp;
	  strtok(NULL, ";");p++; //pass recid;
	  while ((s1=strtok(NULL, ";")))
	    {
	      s2=strtok(NULL, ";");
	      p++;
	      if (fi[p]==-1)
		{
		  b=fi[p]=0;
		  while (fl[b]){if (strcmp(s1,fl[b])==0){fi[p]=1;}b++;}
		}
	      if (fi[p]){data[cn][cdim++]=atof(s2);}
	    }
	  cn++;
	}
    }
  free (buf);
  close (fp);
  return data;
}

void output_data ( double **data,int n, int dim, int mlen, char *infile, char *outfile)
{
  FILE *out;
  FILE *in;
  int cn=0;
  char *buf, *s1,*s2;
  
  if (!outfile)out=stdout;
  else out=fopen (outfile, "w");
  in=fopen (infile, "r");
  
  buf =calloc (mlen+1,sizeof (char));
  cn=0;
  while ((fgets (buf,mlen,in)))
    {
      if (buf[0]=='#')
	{
	  
	  fprintf (out,"%s;",strtok(buf , ";"));//pass #d;
	  fprintf (out,"%s;",strtok(NULL, ";"));//pass exp
	  fprintf (out,"%s;",strtok(NULL, ";"));//pass #rec_id;
	  while ((s1=strtok(NULL, ";")) && s1[0]!='\n')
	    {
	      s2=strtok(NULL, ";");
	      if (strcmp (s1, "bin")!=0)
		{
		  fprintf (out, "%s;%s;", s1, s2);
		}
	    }
	  fprintf (out, "bin;%d;", (int)data[cn++][dim]);
	  fprintf ( out, "\n");
	}
    }
  free (buf);
  close (in);
  if (out!=stdout)close (out);
}

	  
	  

double data2evaluate ( double **d, int n, int dim)
{
  int score=0;
  int a,b,c,s;
  n=10000;
  for (a=0; a<n; a++)
    {
      for (b=a+1; b<n-1; b++)
	{
	  for (s=0,c=0; c<dim; c++)
	    {
	      s+=(d[a][c]==d[b][c])?1:0;
	    }
	  s*=s;
	  score+=s;
	}
    }
  score/=n;
  return sqrt(score);
}
float mk_means (double **data, int n, int dim, int k,double t, double **centroids, int nrounds)
{
  double **result;
  double **sdata;
  float score=0;
  int a, b;
  
  if (nrounds==1)return k_means (data,n,dim,k,t,centroids);
  
  result=calloc (n,sizeof (double*));
  for (a=0; a<n; a++)result[a]=calloc (nrounds+1, sizeof (double));
  
  sdata=calloc ( n, sizeof (double*));
  for (a=0; a<nrounds; a++)
    {
      fprintf ( stderr, "Round %d\n", a+1);
      sdata=shuffle_data (data,sdata, n, 10);
      k_means (sdata, n, dim, k, t, centroids);
      for (b=0; b<n; b++)result[b][a]=data[b][dim];
      fprintf ( stderr, "\n");
    }
  k_means (result, n,nrounds,k,t, centroids);
  score=data2evaluate(result, n, nrounds);
  
  display_data (result,n,nrounds+1);
  
  for(a=0; a<n; a++)
    {
      data[a][dim]=result[a][nrounds];
      free (result[a]);
    }
  
  free(result);
  free(sdata);
  return score;
}
	
double** shuffle_data (double **d, double **sd, int n, int r)
{
  int a,b, sn;
  
  
  if (!sd)sd=calloc( n, sizeof (double*));
  for (a=0; a<r; a++)
    {
      int p=rand()%n;
      sn=0;
      for (b=p; b<n; b++)sd[sn++]=d[b];
      for (b=0; b<p; b++)sd[sn++]=d[b];
    }
  return sd;
}



int k_means(double **data, int n, int m, int k, double t, double **centroids)
{
   /* output cluster label for each data point */
  int a;
   int h, i, j; /* loop counters, of course :) */
   int *counts = (int*)calloc(k, sizeof(int)); /* size of each cluster */
   double old_error, error = DBL_MAX; /* sum of squared euclidean distance */
   double **c = centroids ? centroids : (double**)calloc(k, sizeof(double*));
   double **c1 = (double**)calloc(k, sizeof(double*)); /* temp centroids */

   assert(data && k > 0 && k <= n && m > 0 && t >= 0); /* for debugging */
   for (a=0; a<n; a++)data[a][m]=0;
   /****
   ** initialization */

   for (h = i = 0; i < k; h += n / k, i++) {
      c1[i] = (double*)calloc(m, sizeof(double));
      if (!centroids) {
         c[i] = (double*)calloc(m, sizeof(double));
      }
      /* pick k points as initial centroids */
      for (j = m; j-- > 0; c[i][j] = data[h][j]);
   }

   /****
   ** main loop */

   do {
      /* save error from last step */
      old_error = error, error = 0;

      /* clear old counts and temp centroids */
      for (i = 0; i < k; counts[i++] = 0) {
         for (j = 0; j < m; c1[i][j++] = 0);
      }

      for (h = 0; h < n; h++) {
         /* identify the closest cluster */
         double min_distance = DBL_MAX;
         for (i = 0; i < k; i++) {
            double distance = 0;
            for (j = m; j-- > 0; distance += pow(data[h][j] - c[i][j], 2));
            if (distance < min_distance) {
	      
              data[h][m]=i;
	      min_distance = distance;
            }
         }
         /* update size and temp centroid of the destination cluster */
         for (j = m; j-- > 0; c1[(int)data[h][m]][j] += data[h][j]);
         counts[(int)data[h][m]]++;
         /* update standard error */
         error += min_distance;
      }

      for (i = 0; i < k; i++) { /* update all centroids */
         for (j = 0; j < m; j++) {
            c[i][j] = counts[i] ? c1[i][j] / counts[i] : c1[i][j];
         }
      }

   } while (fabs(error - old_error) > t);

   /****
   ** housekeeping */

   for (i = 0; i < k; i++) {
      if (!centroids) {
         free(c[i]);
      }
      free(c1[i]);
   }

   if (!centroids) {
      free(c);
   }
   free(c1);

   free(counts);
   
   return 1;
}

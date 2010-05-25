#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

main ()
{
  fork_main ();
}
nfork_main()
{
  int a, b, c;
  for (a=0; a<20; a++)
    function ();
}
static int n;
fork_main ()
{
  int a;
  int pid;
  
  for ( a=0; a< 20; a++)
    {
      
      pid=fork();
      
      if (pid==0)
	{
	  n++;
	  function();
	  exit (EXIT_SUCCESS);
	}
      else{;}
    }
  fprintf (stderr,"DONE %d", n);
  while (n>=0)
    {
      wait (NULL);
      n--;
    }
}
int function ()
{
  int a, b;
  for ( a=0; a<1000; a++)
    for (b=0; b<100000;b++);
      
  fprintf (stdout, "HERE %d\n", getpid());
}

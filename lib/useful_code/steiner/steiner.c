#include<stdio.h>
#include<stdlib.h>
#include"steiner.h"
#include<limits.h>
#define N 100

/* 

Runs 2-approximation algorithm for problem of finding a rectilinear network
connecting all up-right pairs of a set of n points.

Algorithm details and analysis to appear in
L. Pachter, F. Lam, M. Alexandersson, Picking Alignments from (Steiner) Trees, 
Proceedings of the Sixth Annual International Conference on Computational 
Molecular Biology (RECOMB 2002)

*/


void rsma ( POINT wpts[], int root, int a[], int p, EDGELIST **wHead, EDGELIST 
**wTail, int pass, int xmax, int ymax);

void slidearb ( POINT vpts[], int array[], EDGELIST **vHead, EDGELIST **vTail, 
int cost[][N], int split[][N], int left, int right, int pass, int xmax, int ymax );

void insert ( EDGELIST **head, EDGELIST **tail, EDGELIST *edge);

void printlist( EDGELIST *head );

static int udSort( const void *a, const void *b ) 
{
  if( ((POINT *)a)->y < ((POINT *)b)->y)
    return 1;
  else if ( ( ((POINT *)a)->y == ((POINT *)b)->y ) && ( ((POINT *)a)->x < 
    ((POINT *)b)->x ) )
    return 1; 
  else
    return -1;
}

static int lrSort( const void *a, const void *b ) 
{
  if( ((POINT *)a)->x > ((POINT *)b)->x)
    return 1;
  else if ( ( ((POINT *)a)->x == ((POINT *)b)->x ) && ( ((POINT *)a)->y > 
    ((POINT *)b)->y ) ) 
    return 1;
  else
    return -1;
}


EDGELIST *steiner( POINT *pts, int ptNum )

/* steiner outputs the edge list of a 2-approximation rectilinear network 
on point set pts */

{
  int i;
  EDGELIST *head, *tail;
  head = tail = NULL;

  // First pass - up/down
  // Sort points in descending-y order.
 
  qsort( pts, ptNum, sizeof(POINT), udSort );  
  
  /* For each point, construct a list for the slide of points 
  directed to it */

  int **ptLists = new int* [ptNum];
  int *ptCounts = new int [ptNum];
  int j = 0;

  for( i = 0; i < ptNum; i++ ) {
    ptLists[i] = new int [ptNum];
    ptCounts[i] = 0;
  }

{
  int i = 0;
  while( i < ptNum ) {
    j = i + 1;
    while( j <= ptNum - 1 && pts[j].x > pts[i].x )
      j++;
    if (j == ptNum)
      i++; 
    else {
      ptLists[j][ptCounts[j]] = i;
      ptCounts[j]++;
      i++;
      }
  } 
} 

  /* For each point u, call function rsma to find the Rectilinear Steiner 
  Minimum Arborescence between u and the slide of points pointing to u */

  EDGELIST *tempHead, *tempTail;
  tempHead = tempTail = NULL;

  for( j = 0; j < ptNum; j++ ){
    rsma ( pts, j, ptLists[j], ptCounts[j], &head, &tail, 1, 0, 0 );    
  } 
  
  // end up/down pass
  // Second pass - left/right

  // Find maximum x and y values
  int max1 = 0, max2 = 0;    

  for( i = 0; i < ptNum; i++ ){
    if ( pts[i].x > max1 )   
        max1 = pts[i].x;
    if ( pts[i].y > max2 )
        max2 = pts[i].y;
  }

  for ( i = 0; i < ptNum; i++ )  {
    pts[i].x = max1 - pts[i].x;
    pts[i].y = max2 - pts[i].y;
  }
    
  qsort( pts, ptNum, sizeof(POINT), lrSort );

  for( i = 0; i < ptNum; i++ )
    ptCounts[i] = 0;

{
  int i = 1;
  while( i < ptNum ) {
    j = i - 1;
    while( j >= 0 && pts[j].y > pts[i].y )
      j--;
    if (j == -1)
      i++;
    else {
      ptLists[j][ptCounts[j]] = i;
      ptCounts[j]++;
      i++;
      }
  }     
}


  for( j = 0; j < ptNum; j++) {  
    rsma ( pts, j, ptLists[j], ptCounts[j], &head, &tail, 2, max1, max2 );
  }

  //end second pass

  // clean up
  for( i = 0; i < ptNum; i++ )
     delete[] ptLists[i];
  delete[] ptLists;
  delete[] ptCounts;

  return head;

} //end steiner
  

void rsma( POINT wpts[], int root, int a[], int p, EDGELIST **wHead, EDGELIST **wTail, int pass, 
int xmax, int ymax )  
  
  /* Constructs Rectilinear Steiner Minimum Arborescence between root and slide 
  of points using dynamic programming 
    a is array of subscripts for slide corresponding to wpts[j]
    p is number of points in slide */

{

  int i, j, k, r, q;

  if ( p == 0)
     return;

  int m[p][N], s[p][N];
  for ( i = 0; i < p; i++ ) 
      for ( j = 0; j <= i; j++)  {
	m[i][j] = 0;
	s[i][j] = 0;
      }
 
  if (p == 2) {
     m[0][1] = wpts[a[0]].y - wpts[a[1]].y + wpts[a[1]].x - wpts[a[0]].x;
     s[0][1] = 0;
  } 
  else {    
    for ( r = 1; r <= p-1; r++) {
      for ( i = 0; i < p-r; i++) {
	  j = i + r;
	  m[i][j] = INT_MAX;
	  for ( k = i; k < j; k++) {
		q = m[i][k] + m[k+1][j] + wpts[a[k+1]].x - 
			wpts[a[i]].x + wpts[a[k]].y - wpts[a[j]].y;
		if ( q < m[i][j] ) {
		   m[i][j] = q;
		   s[i][j] = k;
	  	}
	  }
       }
    }  
  } //end else
  

  /* edge construction
  edge1, edge2 connect root to parent node in slide arborescence */

  EDGELIST *edge1 = new EDGELIST;
  EDGELIST *edge2 = new EDGELIST;

  if (pass == 1) {
  edge1->start = wpts[root];
  edge1->end.x = wpts[a[0]].x;
  edge1->end.y = wpts[root].y;
  edge2->start = edge1->end;
  edge2->end.x = wpts[a[0]].x;
  edge2->end.y = wpts[a[p-1]].y;
  }
  else  {
	edge1->start.x = xmax - edge2->end.x;
	edge1->start.y = ymax - edge2->end.y;
        edge1->end.x = xmax - wpts[root].x;
        edge1->end.y = edge1->start.y;
	edge2->start = edge1->end;
        edge2->end.x = xmax - wpts[root].x;
	edge2->end.y = ymax - wpts[root].y;
  }
  
  insert( wHead, wTail, edge1 );
  insert( wHead, wTail, edge2 );

  // find remaining edges of arboresence
  slidearb ( wpts, a, wHead, wTail, m, s, 0, p-1, pass, xmax, ymax );   

}

void slidearb ( POINT vpts[], int array[], EDGELIST **vHead, EDGELIST 
**vTail, int cost[][N], int split[][N], int left, int right, int vpass, 
int xmax, int ymax)

/* constructs edges in Rectilinear Steiner Minimum Arborescence for points in 
array vpts.  "Cost" and "split" are matrices indicating optimal values for the 
RSMA of subsets of points corresponding to intervals along the slide */

{  

  if (left == right)
    return;

  EDGELIST *edge3 = new EDGELIST;
  EDGELIST *edge4 = new EDGELIST;

  if ( vpass == 1) {
  edge3->start.x = vpts[array[left]].x;
  edge3->start.y = vpts[array[right]].y;
  edge3->end.x = vpts[array[left]].x;
  edge3->end.y = vpts[array[split[left][right]]].y; 
  edge4->start = edge3->start;                                       
  edge4->end.x = vpts[array[split[left][right] + 1]].x;                 
  edge4->end.y = vpts[array[right]].y;
  } 
  else {
        edge3->start.x = xmax - vpts[array[split[left][right] + 1]].x;
        edge3->start.y = ymax - vpts[array[right]].y;
        edge3->end.x = xmax - edge3->start.x;
        edge3->end.y = ymax - edge3->start.y;
 	edge4->start.x = xmax - vpts[array[left]].x;
	edge4->start.y = ymax - vpts[array[split[left][right]]].y;
        edge4->end = edge3->end;   
  }

  insert ( vHead, vTail, edge3 );
  insert ( vHead, vTail, edge4 );    

  slidearb ( vpts, array, vHead, vTail, cost, split, left, 
	split[left][right], vpass, xmax, ymax);
  slidearb ( vpts, array, vHead, vTail, cost, split, 
	split[left][right] + 1, right, vpass, xmax, ymax);
 
}


void insert ( EDGELIST **head, EDGELIST **tail, EDGELIST *edge)
{
    if ( ( edge->start.x == edge->end.x)  && ( edge->start.y == edge->end.y ) )
	return;
  
    if( *tail == NULL ){                                                                   
        *head = edge;
      	*tail = edge;
    }
    else {
      	(*tail)->next = edge;   
        *tail = edge;
    }
    (*tail)->next = NULL;


}


void printlist( EDGELIST *head )
{

  EDGELIST *temp;
  temp = head;
  while (temp != NULL)  {
        printf("(%d,  %d)  ", temp->start.x, temp->start.y);
        printf("(%d,  %d) \n", temp->end.x, temp->end.y);
        temp = temp->next;
  }
  printf("\n");

}


int main()
{

  int i, j, s;
  EDGELIST *head, *temp;
  
  scanf("%d", &s);
  POINT pts[s];
  for ( i = 0; i < s; i++ )  
	scanf("%d %d", &pts[i].x, &pts[i].y);

  // find maximum x and y values

  int m = 0, n = 0;   
  for( i = 0; i < s; i++ ){
    if ( pts[i].x > m )
        m = pts[i].x;
    if ( pts[i].y > n )
        n = pts[i].y;
  }

  // find rectilinear network connecting all up-right pairs
  head = steiner(pts, s);

  /* Construct matrices A and B corresponding to edges in network
  A is 1 for each point in the constructed network
  B indicates the starting and ending positions of the 1's in A */

  int A[n+1], B[n+1], count;
  for ( i = 0; i < m+1; i++ ) {
      count = 0;
      for (j = 0; j < n+1; j++ ) {
          A[j] = 0;
 	  B[j] = 0;
      }
      temp = head;
      while ( temp != NULL ) {
          if ( temp->start.x <= i && i <= temp->end.x )  {
             for ( j = (temp->start.y) ; j < (temp->end.y) + 1; j++ )
                 A[j] = 1;
          } 
          temp = temp->next;
      }


      if ( A[0] == 1 )   {
	 B[count] = 0;
	 count++;
      }
      for (j = 0; j < n; j++) {
          if (A[j] == 0 && A[j+1] == 1) {
		B[count] = j+1;
		count++;
	  }
	  if (A[j] == 1 && A[j+1] == 0) {
	    	B[count] = j;
		count++;
	  }
      }
      if ( A[n] == 1 )   {      
         B[count] = n;
         count++;
      }

      printf("%d  ", i);
      for (j = 0; j < count; j++) {
		printf("%d ", B[j]);
      }
      printf("\n");


  }	

  return 0;

}

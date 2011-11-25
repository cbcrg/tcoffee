
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

HaschT * hcreate ( int n_elements,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )
       {
	 HaschT *T;
	 int a;
	 
	 n_elements=n_elements*2+1;
	 
	 T=vcalloc ( 1, sizeof (HaschT));
	 T->ne=n_elements;	 
	 T->p=vcalloc (n_elements,sizeof ( Hasch_entry*));
	 for ( a=0; a<n_elements; a++)
	   {
	     T->p[a]=allocate_hasch_entry(NULL,DECLARE,declare_data, free_data);
	   }
	 return T;
       }
HaschT *hdestroy (HaschT *T,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )


       {
	 int a;
	 Hasch_entry *p, *pp;
	 
	 if ( T==NULL)return NULL;

	 for (a=0; a< T->ne; a++)
	   {
	     p=T->p[a];
	     while (p)
	       {
		 pp=p;
		 p=p->n;
		 allocate_hasch_entry(pp,FREE, declare_data, free_data);
	       }
	   }
	 vfree (T->p);
	 vfree ( T);
	 return NULL;
       }


Hasch_entry* hsearch (HaschT *T, int k, int action, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )


       {
	 /*action: FIND,ADD, REMOVE*/
	 Hasch_entry *p, *pi;
	 int h;
	 
	 

	 /* find the key: k->h*/
	 
	 h=k%T->ne;
	 

	 if ( action==ADD || action==FIND)
	   {	
	     p=pi=T->p[h];
	     while (p && p->k!=k){p=p->n;}
	     if (action==ADD && !p)
	      {
		p=insert_hasch_entry_in_list (pi, NULL, NULL, declare_data, free_data);
		p->k=k;	
	      }
	    else if (action==FIND && !p)p=NULL;
	    return p;
	   } 
	 else if ( action==REMOVE)
	   {
	     allocate_hasch_entry(hsearch ( T, k, FIND, declare_data, free_data), FREE, declare_data, free_data);
	     return NULL;
	   }
       return NULL;
       }


Hasch_entry * extract_hasch_entry_from_list (Hasch_entry *e, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )


  {
    /*extracts entry e and returns p, or next if is NULL*/
    Hasch_entry *p=NULL, *n=NULL;
    
    if (!e);
    else
      {
	p=e->p;
	n=e->n;

	if (p)p->n=n;
	if (n)n->p=p;
	e->p=e->n=NULL;
      }
    return e;
  }

Hasch_entry * insert_hasch_entry_in_list (Hasch_entry *p, Hasch_entry *e, Hasch_entry *n, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )


{
    /*inserts entry e between entry p and entry n and returns e*/

    if (!e)e=allocate_hasch_entry (NULL,DECLARE, declare_data, free_data);
    
    
      
    if (!p && !n);
    else if ( !p)p=n->p;
    else if ( !n)n=p->n;

    e->p=p;
    if (p)p->n=e;
    
    e->n=n;
    if (n)n->p=e;
    
    return e;
  }

Hasch_entry * allocate_hasch_entry (Hasch_entry *e, int action,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )


{
    static Hasch_entry *s;
    Hasch_entry *ns;

    if ( !s)s=vcalloc ( 1, sizeof (Hasch_entry));

    if ( action==DECLARE)
      {
	ns=s->p;
	e=extract_hasch_entry_from_list (s, declare_data, free_data);
	if ( e->free_data)(e->free_data)(e->data);
	e->declare_data=declare_data;
	e->free_data=free_data;
	e->declare_data (e);
	e->k=UNDEFINED;
	s=ns;
      }
    else if ( action==FREE)
      {
	extract_hasch_entry_from_list (e,declare_data, free_data );
	e->k=UNDEFINED;
	if ( e->free_data)e->data=(e->free_data)(e->data);
	e->free_data=NULL;
	e->declare_data=NULL;
	s=insert_hasch_entry_in_list (s, e, NULL, declare_data, free_data);
	
      }
    else if ( action==FREE_STACK)
      {
	while (s)
	  {
	    e=s->p;
	    allocate_hasch_entry (s, FREE, declare_data,free_data);
	    vfree (s);
	    s=e;
	  }
      }
    else crash ("Unknown MODE for allocate_hasch_entry\n");
    return e;
  }
    
/*********************************************************************/
/*                                                                   */
/*                         Get string key                                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/


int string2key (char *s, Char_node *n)
{
  static Char_node *root;

  if ( !root)root=declare_char_node (DECLARE);
  
  if ( n==NULL && s==NULL)
    {
      declare_char_node (FREE_STACK);
    }
  else if (n==NULL)
    {
      return string2key(s, root);
    }
  else if ( s[0]=='\0')
    {
      return n->key;
    }
  else
    {
      return string2key(s+1, (n->c[(int)s[0]])?(n->c[(int)s[0]]):(n->c[(int)s[0]]=declare_char_node (DECLARE)));
    }
  return 0;
}

Char_node * declare_char_node (int action)
{
static struct Char_node **heap;
static int heap_size, free_heap, a;
static int key;  
 if ( action==DECLARE)
    {
      if ( free_heap==0)
	{
	  free_heap=100;
	  
	  heap=vrealloc (heap,(heap_size+free_heap)*sizeof (struct Char_node *));
	  for ( a=heap_size; a<heap_size+free_heap; a++)
	    {
	      (heap[a])=vcalloc ( 1, sizeof ( struct Char_node));
	      (heap[a])->c=vcalloc ( 256, sizeof (Char_node*));
	      (heap[a])->key=key++;
	    }
	  heap_size+=free_heap;
	}
      return heap[heap_size-(free_heap--)];
    }
  else if ( action==FREE_STACK)
    {
      for (a=0; a< heap_size; a++)
	{
	  heap[a]->key=key++;
	  vfree ( heap[a]->c);
	  (heap[a])->c=vcalloc ( 256, sizeof (Char_node*));
	}
      free_heap=heap_size;
      return NULL;
    }
  return NULL;
}
  

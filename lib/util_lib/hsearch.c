
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
  
/* old declare_char_node (too hungry)
Char_node * declare_char_node (int action)
{
  static int key;
  Char_node *cn;
  static Char_node *root;

  if ( action==DECLARE)
    {
      cn=vcalloc (1, sizeof (Char_node));
      cn->key=++key;
      cn->c=vcalloc (256, sizeof (Char_node *));
      
      }
  return cn;
}
*/			 
 /***************************************************************************
  *cr
  *cr            (C) Copyright 1995-2011 The Board of Trustees of the
  *cr                        University of Illinois
  *cr                         All Rights Reserved
  *cr
  ***************************************************************************/
 
 /***************************************************************************
  * RCS INFORMATION:
  *
  *      $RCSfile: hash.c,v $
  *      $Author: johns $        $Locker:  $             $State: Exp $
  *      $Revision: 1.14 $      $Date: 2010/12/16 04:08:55 $
  *
  ***************************************************************************
  * DESCRIPTION:
  *   A simple hash table implementation for strings, contributed by John Stone,
  *   derived from his ray tracer code.
  ***************************************************************************/
#define HASH_LIMIT 0.5
 
 typedef struct hash_node_t {
   int data;                           /* data in hash node */
   char * key;                   /* key for hash lookup */
   struct hash_node_t *next;           /* next node in hash chain */
 } hash_node_t;
 
 /*
  *  hash() - Hash function returns a hash number for a given key.
  *
  *  tptr: Pointer to a hash table
  *  key: The key to create a hash number for
  */
 static int hash(const hash_t *tptr, const char *key) {
   int i=0;
   int hashvalue;
  
   while (*key != '\0')
     i=(i<<3)+(*key++ - '0');
  
   hashvalue = (((i*1103515249)>>tptr->downshift) & tptr->mask);
   if (hashvalue < 0) {
     hashvalue = 0;
   }    
 
   return hashvalue;
 }
 
 /*
  *  hash_init() - Initialize a new hash table.
  *
  *  tptr: Pointer to the hash table to initialize
  *  buckets: The number of initial buckets to create
  */
 VMDEXTERNSTATIC void hash_init(hash_t *tptr, int buckets) {
 
   /* make sure we allocate something */
   if (buckets==0)
     buckets=16;
 
   /* initialize the table */
   tptr->entries=0;
   tptr->size=2;
   tptr->mask=1;
   tptr->downshift=29;
 
   /* ensure buckets is a power of 2 */
   while (tptr->size<buckets) {
     tptr->size<<=1;
     tptr->mask=(tptr->mask<<1)+1;
     tptr->downshift--;
   } /* while */
 
  
   /* allocate memory for table */
   tptr->bucket=(hash_node_t **) calloc(tptr->size, sizeof(hash_node_t *));
 
   return;
 }
 
 /*
  *  rebuild_table() - Create new hash table when old one fills up.
  *
  *  tptr: Pointer to a hash table
  */
 static void rebuild_table(hash_t *tptr) {
   hash_node_t **old_bucket, *old_hash, *tmp;
   int old_size, h, i;
 
   old_bucket=tptr->bucket;
   old_size=tptr->size;
 
   /* create a new table and rehash old buckets */
   hash_init(tptr, old_size<<1);
   for (i=0; i<old_size; i++) {
     old_hash=old_bucket[i];
     while(old_hash) {
       tmp=old_hash;
       old_hash=old_hash->next;
       h=hash(tptr, tmp->key);
       tmp->next=tptr->bucket[h];
       tptr->bucket[h]=tmp;
       tptr->entries++;
     } /* while */
   } /* for */
 
   /* free memory used by old table */
   free(old_bucket);
 
   return;
 }
 
 /*
  *  hash_lookup() - Lookup an entry in the hash table and return a pointer to
  *    it or HASH_FAIL if it wasn't found.
  *
  *  tptr: Pointer to the hash table
  *  key: The key to lookup
  */
 VMDEXTERNSTATIC int hash_lookup(const hash_t *tptr, const char *key) {
   int h;
   hash_node_t *node;
 
 
   /* find the entry in the hash table */
   h=hash(tptr, key);
   for (node=tptr->bucket[h]; node!=NULL; node=node->next) {
     if (!strcmp(node->key, key))
       break;
   }
 
   /* return the entry if it exists, or HASH_FAIL */
   return(node ? node->data : HASH_FAIL);
 }
 
 /*
  *  hash_insert() - Insert an entry into the hash table.  If the entry already
  *  exists return a pointer to it, otherwise return HASH_FAIL.
  *
  *  tptr: A pointer to the hash table
  *  key: The key to insert into the hash table
  *  data: A pointer to the data to insert into the hash table
  */
 VMDEXTERNSTATIC int hash_insert(hash_t *tptr, const char *key, int data) {
   int tmp;
   hash_node_t *node;
   int h;
 
   /* check to see if the entry exists */
   if ((tmp=hash_lookup(tptr, key)) != HASH_FAIL)
     return(tmp);
 
   /* expand the table if needed */
   while (tptr->entries>=HASH_LIMIT*tptr->size)
     rebuild_table(tptr);
 
   /* insert the new entry */
   h=hash(tptr, key);
   node=(struct hash_node_t *) malloc(sizeof(hash_node_t));
   node->data=data;
   
   //duplicate key
   node->key=calloc ( strlen (key)+1, sizeof (char));
   sprintf (node->key, "%s", key);
   //node->key=key;
   node->next=tptr->bucket[h];
   tptr->bucket[h]=node;
   tptr->entries++;
 
   return HASH_FAIL;
 }
 
 /*
  *  hash_delete() - Remove an entry from a hash table and return a pointer
  *  to its data or HASH_FAIL if it wasn't found.
  *
  *  tptr: A pointer to the hash table
  *  key: The key to remove from the hash table
  */
 VMDEXTERNSTATIC int hash_delete(hash_t *tptr, const char *key) {
   hash_node_t *node, *last;
   int data;
   int h;
 
   /* find the node to remove */
   h=hash(tptr, key);
   for (node=tptr->bucket[h]; node; node=node->next) {
     if (!strcmp(node->key, key))
       break;
   }
 
   /* Didn't find anything, return HASH_FAIL */
   if (node==NULL)
     return HASH_FAIL;
 
   /* if node is at head of bucket, we have it easy */
   if (node==tptr->bucket[h])
     tptr->bucket[h]=node->next;
   else {
     /* find the node before the node we want to remove */
     for (last=tptr->bucket[h]; last && last->next; last=last->next) {
       if (last->next==node)
         break;
     }
     last->next=node->next;
   }
 
   /* free memory and return the data */
   data=node->data;
   free(node);
 
   return(data);
 }
 
 
 
 /*
  * hash_destroy() - Delete the entire table, and all remaining entries.
  * 
  */
 VMDEXTERNSTATIC void hash_destroy(hash_t *tptr) {
   hash_node_t *node, *last;
   int i;
 
   for (i=0; i<tptr->size; i++) {
     node = tptr->bucket[i];
     while (node != NULL) { 
       last = node;   
       node = node->next;
       free (last->key);
       free(last);
     }
   }     
 
   /* free the entire array of buckets */
   if (tptr->bucket != NULL) {
     free(tptr->bucket);
     memset((void*)tptr, 0, sizeof(hash_t));
   }
 }
 
 /*
  *  alos() - Find the average length of search.
  *
  *  tptr: Pointer to a hash table
  */
 static float alos(hash_t *tptr) {
   int i,j;
   float alos=0;
   hash_node_t *node;
 
 
   for (i=0; i<tptr->size; i++) {
     for (node=tptr->bucket[i], j=0; node!=NULL; node=node->next, j++);
     if (j)
       alos+=((j*(j+1))>>1);
   } /* for */
 
   return(tptr->entries ? alos/tptr->entries : 0);
 }
 
 
 /*
  *  hash_stats() - Return a string with stats about a hash table.
  *
  *  tptr: A pointer to the hash table
  */
 VMDEXTERNSTATIC char * hash_stats(hash_t *tptr) {
   static char buf[1024];
 
   sprintf(buf, "%u slots, %u entries, and %1.2f ALOS",
     (int)tptr->size, (int)tptr->entries, alos(tptr));
 
   return(buf);
 }
 

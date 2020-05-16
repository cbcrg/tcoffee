
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <climits>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

HaschT * hcreate ( int n_elements,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) )
       {
	 HaschT *T;
	 int a;
	 
	 n_elements=n_elements*2+1;
	 
	 T=(HaschT*)vcalloc ( 1, sizeof (HaschT));
	 T->ne=n_elements;	 
	
	 T->p=(Hasch_entry**)vcalloc (n_elements,sizeof ( Hasch_entry*));
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

    if ( !s)s=(Hasch_entry*)vcalloc ( 1, sizeof (Hasch_entry));

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
	  
	  heap=(Char_node**)vrealloc (heap,(heap_size+free_heap)*sizeof (struct Char_node *));
	  for ( a=heap_size; a<heap_size+free_heap; a++)
	    {
	      (heap[a])=(Char_node*)vcalloc ( 1, sizeof ( struct Char_node));
	      (heap[a])->c=(Char_node**)vcalloc ( 256, sizeof (Char_node*));
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
	  (heap[a])->c=(Char_node**)vcalloc ( 256, sizeof (Char_node*));
	}
      free_heap=heap_size;
      return NULL;
    }
  return NULL;
}
  



/* Create a new hashtable. */
hashtable_t * ht_destroy( hashtable_t *ht)
{
  entry_t *previous = NULL;
  entry_t *next = NULL;
  int a;
  
  for (a=0; a<ht->size; a++)
    {
      next=ht->table[a];
      while (next)
	{
	  previous=next;
	  next=next->next;
	  free(previous);
	}
    }
  free (ht->table);
  free (ht);
  return NULL;
}

hashtable_t *ht_create( int size ) {

	hashtable_t *hashtable = NULL;
	int i;

	if( size < 1 ) return NULL;

	/* Allocate the table itself. */
	if( ( hashtable = (hashtable_t*) malloc( sizeof( hashtable_t ) ) ) == NULL ) {
		return NULL;
	}

	/* Allocate pointers to the head nodes. */
	if( ( hashtable->table=( entry_s **)malloc( sizeof( entry_t * ) * size ) ) == NULL ) {
		return NULL;
	}
	for( i = 0; i < size; i++ ) {
		hashtable->table[i] = NULL;
	}

	hashtable->size = size;

	return hashtable;	
}


/* Hash a string for a particular hash table. */
int ht_hash (hashtable_t *hashtable, char *key )
{
  size_t i = 0;
  unsigned long hash = 0;
  int length=strlen (key);
  while (i != length) {
    hash += key[i++];
    hash += hash << 10;
    hash ^= hash >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;
  return hash  % hashtable->size;
}
int ht_hash_simple( hashtable_t *hashtable, char *key ) {

	unsigned long int hashval=0;
	int i = 0;
	int h;
	/* Convert our string to an integer */
	while( hashval < ULONG_MAX && i < strlen( key ) ) {
		hashval = hashval << 8;
		hashval += key[ i ];
		i++;
	}
	
	h= hashval % hashtable->size;
	
	return h;
}

/* Create a key-value pair. */
entry_t *ht_newpair( char *key, char *value ) {
	entry_t *newpair;

	if( ( newpair = (entry_t*)malloc( sizeof( entry_t ) ) ) == NULL ) {
		return NULL;
	}

	if( ( newpair->key = strdup( key ) ) == NULL ) {
		return NULL;
	}

	if( ( newpair->value = strdup( value ) ) == NULL ) {
		return NULL;
	}

	newpair->next = NULL;

	return newpair;
}

/* Insert a key-value pair into a hash table. */
void ht_set( hashtable_t *hashtable, char *key, char *value ) {
	int bin = 0;
	entry_t *newpair = NULL;
	entry_t *next = NULL;
	entry_t *last = NULL;

	bin = ht_hash( hashtable, key );

	next = hashtable->table[ bin ];

	while( next != NULL && next->key != NULL && strcmp( key, next->key ) > 0 ) {
		last = next;
		next = next->next;
	}

	/* There's already a pair.  Let's replace that string. */
	if( next != NULL && next->key != NULL && strcmp( key, next->key ) == 0 ) {

		free( next->value );
		next->value = strdup( value );

	/* Nope, could't find it.  Time to grow a pair. */
	} else {
		newpair = ht_newpair( key, value );

		/* We're at the start of the linked list in this bin. */
		if( next == hashtable->table[ bin ] ) {
			newpair->next = next;
			hashtable->table[ bin ] = newpair;
	
		/* We're at the end of the linked list in this bin. */
		} else if ( next == NULL ) {
			last->next = newpair;
	
		/* We're in the middle of the list. */
		} else  {
			newpair->next = next;
			last->next = newpair;
		}
	}
}

/* Retrieve a key-value pair from a hash table. */
char *ht_get( hashtable_t *hashtable, char *key ) {
	int bin = 0;
	entry_t *pair;

	bin = ht_hash( hashtable, key );

	/* Step through the bin, looking for our value. */
	pair = hashtable->table[ bin ];
	while( pair != NULL && pair->key != NULL && strcmp( key, pair->key ) > 0 ) {
		pair = pair->next;
	}

	/* Did we actually find anything? */
	if( pair == NULL || pair->key == NULL || strcmp( key, pair->key ) != 0 ) {
		return NULL;

	} else {
		return pair->value;
	}
	
}

hashtable_t * array2hash (char **list)
{
  if (!list)return NULL;
  else return array2hashN(list,-1);
}
hashtable_t * array2hashN (char **list, int nn)
{
  int hn;
  hashtable_t *ht;
  Memcontrol*M;
  
  
  if (!list)return NULL;
  if (nn==-1)nn=arrlen(list);
  
  M=(Memcontrol*)list;
  M-=2;
  hn=M[0].hn;
  ht=M[0].ht;
  
  

  if (nn!=hn && ht)
    {
      ht_destroy(ht);
      ht=NULL;
    }
  
  
  if (!ht)
    {
      int a;
      char index[100];
      ht=ht_create(nn*10);
      
      for (a=0; a<nn; a++)
	{
	  sprintf (index, "%d", a);
	  ht_set (ht,list[a], index);
	}
        M[0].hn=nn;
	M[0].ht=ht;
    }
  
  return ht;
}
hashtable_t *hupdate (char **list)
{
  int hn;
  hashtable_t *ht;
  Memcontrol*M;
  int nn;
  
  if (!list)return NULL;
  nn=arrlen(list);
  
  M=(Memcontrol*)list;
  M-=2;
  hn=M[0].hn;
  ht=M[0].ht;
  
  if (ht)hfree(list);
  return array2hash(list);
}
hashtable_t * hfree (void *list)
{
  Memcontrol*M;

  if (!list) return NULL;
  M=(Memcontrol*)list;
  M-=2;
  
  if (M[0].ht)
    {
      ht_destroy(M[0].ht);
      M[0].ht=NULL;
    }
  return NULL;
}
int name2hindex (char *name, char **list)
{
  return name_is_in_hlist(name, list,-1);
}
int name_is_in_hlist (char *name, char **list, int n)
{
  hashtable_t *ht;
  char *r;
  

  if (!name){hfree(list);return -1;}
  if (!list)return -1;
  ht=array2hashN(list,n);
  
  if ((r=ht_get(ht, name))!=NULL)
    {
      return atoi (r);
    }
  else
    return -1;
}
int name_is_in_hlist_old (char *name, char **list, int n)
{
  static hashtable_t *ht;
  static char **rlist;
  char *r;
  
  if (list !=rlist || list==NULL)
    {
      rlist=list;
      if (ht)ht=ht_destroy(ht);
    }
  if (list==NULL)return -1;
  
  if (!ht)
    {
      int a;
      char index[100];
      
      ht=ht_create (n*10);
      for (a=0;a<n; a++)
	{
	  sprintf (index, "%d", a);
	  ht_set (ht,list[a], index);
	  //HERE ("SET %s %s", list[a], index);
	}
    }
  if ((r=ht_get(ht, name))!=NULL)
    {
      return atoi (r);
    }
  else
    return -1;
}
char *check_hlist_for_dup (char **name, int n)
{
  hashtable_t *ht;
  char *r;
  char index[100];
  int a;
  char *dup=NULL;
  Memcontrol*M;
  
  if (! name)return NULL;
  if (n<=0)n=arrlen (name);
  ht=ht_create (n);
  
  for (a=0;a<n; a++)
    {
      
      if ((r=ht_get(ht, name[a]))!=NULL)
	{
	  dup=strcatf (dup, " --- Duplicated Sequence: %s\n", name[a]);
	}
      else
	{
	  sprintf (index, "%d", a);
	  ht_set (ht,name[a], index);
	}
    }
  
  M=(Memcontrol*)name;
  M-=2;
  if (M[0].ht)hfree(name);
  M[0].ht=ht;
  
  return dup;
}
  

int old_main( int argc, char **argv ) {

	hashtable_t *hashtable = ht_create( 65536 );

	ht_set( hashtable, "key1", "inky" );
	ht_set( hashtable, "key2", "pinky" );
	ht_set( hashtable, "key3", "blinky" );
	ht_set( hashtable, "key4", "floyd" );

	printf( "%s\n", ht_get( hashtable, "key1" ) );
	printf( "%s\n", ht_get( hashtable, "key2" ) );
	printf( "%s\n", ht_get( hashtable, "key3" ) );
	printf( "%s\n", ht_get( hashtable, "key4" ) );

	return 0;
}

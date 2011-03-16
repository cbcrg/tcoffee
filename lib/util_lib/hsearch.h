/*********************************************************************************************/
/*                                                                                           */
/*         STRUCTURES FOR HSEARCH                                                            */
/*                                                                                           */
/*********************************************************************************************/
#define FIND           0
#define ADD            1
#define REMOVE         2
#define DECLARE        3
#define MARK           4
#define UNMARK         5
#define FREE           6
#define FREE_STACK     7
#define FREE_ALL       8
#define FREE_MARK      9
#define INFO           10
     
struct HaschT
{
  int ne;
  struct Hasch_entry **p;
};
typedef struct HaschT HaschT;

struct Hasch_entry
{
  struct Hasch_entry *n;
  struct Hasch_entry *p;
  int k;
  struct Hasch_data  *data;
  struct Hasch_data * (*free_data)(struct Hasch_data *); 
  struct Hasch_data * (*declare_data)(struct Hasch_entry*);
  int tag;
};
typedef struct Hasch_entry Hasch_entry;
struct Char_node
{
 struct Char_node **c;
 int key;
 
};
typedef struct Char_node Char_node;

HaschT * hcreate ( int n_elements,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
HaschT *hdestroy (HaschT *T,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry* hsearch (HaschT *T, int k, int action, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry * extract_hasch_entry_from_list (Hasch_entry *e, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry * insert_hasch_entry_in_list (Hasch_entry *p, Hasch_entry *e, Hasch_entry *n, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry * allocate_hasch_entry (Hasch_entry *e, int action,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );




 
int string2key (char *s, Char_node *n);
Char_node * declare_char_node (int action);


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
  *      $RCSfile: hash.h,v $
  *      $Author: johns $        $Locker:  $             $State: Exp $
  *      $Revision: 1.12 $      $Date: 2010/12/16 04:08:55 $
  *
  ***************************************************************************
  * DESCRIPTION:
  *   A simple hash table implementation for strings, contributed by John Stone,
  *   derived from his ray tracer code.
  ***************************************************************************/
 #ifndef HASH_H
 #define HASH_H
 
 typedef struct hash_t {
   struct hash_node_t **bucket;        /* array of hash nodes */
   int size;                           /* size of the array */
   int entries;                        /* number of entries in table */
   int downshift;                      /* shift cound, used in hash function */
   int mask;                           /* used to select bits for hashing */
 } hash_t;
 
 #define HASH_FAIL -1
 
 #if defined(VMDPLUGIN_STATIC)
 #define VMDEXTERNSTATIC static
 #include "hash.c"
 #else
 
 #define VMDEXTERNSTATIC 
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 void hash_init(hash_t *, int);
 
 int hash_lookup (const hash_t *, const char *);
 
 int hash_insert (hash_t *, const char *, int);
 
 int hash_delete (hash_t *, const char *);
 
 void hash_destroy(hash_t *);
 
 char *hash_stats (hash_t *);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif
 
 #endif

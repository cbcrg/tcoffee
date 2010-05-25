

typedef struct {
  int x, y;
} POINT;

struct pointlist {
  POINT pt;
  struct pointlist *next;
};
typedef struct pointlist PTLIST;

struct edgelist {
  POINT start, end;
  struct edgelist *next;
};
typedef struct edgelist EDGELIST;

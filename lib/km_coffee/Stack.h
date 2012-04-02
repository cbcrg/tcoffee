#ifndef STACK_H_
#define STACK_H_

// #include "km_util.h"

typedef struct
{
	size_t size;
	size_t reserved;
	void **data;
	void *last;
} Stack;


Stack*
Stack_init();

void
push(Stack *stack, void *data);

void*
top(Stack *stack);

void *
pop(Stack *stack);

void
delStack(Stack *stack);
#endif

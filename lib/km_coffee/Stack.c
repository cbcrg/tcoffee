// #include "Stack.h"
#include "km_coffee_header.h"

Stack *
Stack_init()
{
	Stack *stack = my_malloc(sizeof(Stack));
	stack->reserved = 0;
	stack->size=0;
	stack->data=NULL;
	stack->last=NULL;
	return stack;
}

void
push(Stack *stack, void *data)
{
	if (stack->reserved +1 > stack->size)
	{
		stack->reserved += 10;
		stack->data=realloc(stack->data, stack->reserved);
	}
	stack->data[stack->size] = data;
	stack->last=stack->data[stack->size];
	++stack->size;
}

void *
pop(Stack *stack)
{
	--stack->size;
	if (stack->size == 0)
		stack->last=NULL;
	else
		stack->last=stack->data[stack->size-1];
	return stack->data[stack->size];
}

void
delStack(Stack *stack)
{
	free(stack->data);
	free(stack);
}
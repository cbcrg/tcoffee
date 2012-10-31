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


/**
 * \brief Initializes the stack.
 * \return A pointer to a new stack.
 */
Stack*
Stack_init();

/**
 * \brief Adds a new element to the stack.
 * \param stack The stack to which the new element should be added.
 * \param data A pointer to the object to add.
 */
void
push(Stack *stack, void *data);

/**
 * \brief Return the first element on the stack (The last one added).
 * \param stack The stack.
 * \return A pointer ot the topmost element.
 */
void*
top(Stack *stack);

/**
 * \brief Removes the topmost element from the stack.
 *
 * \param stack The stack to delete.
 * \return Pointer to the element just removed from the stack.
 */
void *
pop(Stack *stack);

/**
 * \brief Deletes the stackit self and frees the memory it requires.
 * \param stack The stack to free.
 * \warning Elements still on the stack are not freed!
 */
void
delStack(Stack *stack);
#endif

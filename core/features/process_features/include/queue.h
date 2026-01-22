#ifndef QUEUE_H
#define QUEUE_H

#include <stdlib.h>
#include <stdint.h>

// Define the initial capacity of the queue
#define INITIAL_CAPACITY 16

// Define the Queue structure
typedef struct _queue {
    uint64_t *data;
    size_t front;
    size_t back;
    size_t queue_size;
    size_t capacity;
} Queue;

// Queue function prototypes
void init_queue(Queue *queue);
void clear_queue(Queue *queue);
int is_empty(Queue *queue);
size_t size(Queue *queue);
void enqueue(Queue *queue, uint64_t value);
uint64_t dequeue(Queue *queue);
uint64_t peek(Queue *queue);
void free_queue(Queue *queue);


#endif // QUEUE_H

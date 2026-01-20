#include "../include/common.h"
#include "../include/queue.h"

void init_queue(Queue *queue) {
    queue->capacity = INITIAL_CAPACITY;
    queue->data = (uint64_t *)malloc(queue->capacity * sizeof(uint64_t));
    if (!queue->data) {
        fprintf(stderr, "Failed to allocate memory for queue\n");
        exit(EXIT_FAILURE);
    }
    queue->front = 0;
    queue->back = 0;
    queue->queue_size = 0;
}

void clear_queue(Queue *queue) {
    queue->front = 0;
    queue->back = 0;
    queue->queue_size = 0;
}

int is_empty(Queue *queue) {
    return queue->queue_size == 0;
}

size_t size(Queue *queue) {
    return queue->queue_size;
}

void enqueue(Queue *queue, uint64_t value) {
    if (queue->queue_size == queue->capacity) {
        queue->capacity *= 2;
        queue->data = (uint64_t *)realloc(queue->data, queue->capacity * sizeof(uint64_t));
        if (queue->front > queue->back) {
            for (int i = 0; i < queue->front; i++) {
                queue->data[queue->capacity / 2 + i] = queue->data[i];
            }
            queue->back = queue->capacity / 2 + queue->front;
        }
    }
    queue->data[queue->back] = value;
    queue->back = (queue->back + 1) % queue->capacity;
    queue->queue_size++;
}

uint64_t dequeue(Queue *queue) {
    if (is_empty(queue)) {
        printf("Queue underflow\n");
        exit(EXIT_FAILURE);
    }
    uint64_t value = queue->data[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->queue_size--;
    return value;
}

uint64_t peek(Queue *queue) {
    if (is_empty(queue)) {
        printf("Queue is empty\n");
        exit(EXIT_FAILURE);
    }
    return queue->data[queue->front];
}

void free_queue(Queue *queue) {
    free(queue->data);
    queue->data = NULL;
    queue->front = 0;
    queue->back = 0;
    queue->queue_size = 0;
    queue->capacity = 0;
}
#ifndef MEMORY_H
#define MEMORY_H

#include <stddef.h> // For size_t
#include "common.h"

// Memory management constants
#define FEATURE_BARCODE_BLOCK_SIZE 30000
#define UMI_STORAGE_BLOCK 10000000
#define FEATURE_SEQUENCE_BLOCK_SIZE 100000
#define BARCODE_STORAGE_BLOCK_SIZE 30000
#define CB_COUNTS_BLOCK_SIZE  10000


// Function prototypes for memory management
memory_pool_collection* initialize_memory_pool_collection();
void free_memory_pool_collection(memory_pool_collection *pools);
void free_memory_pool_storage(memory_pool *pool);
storage_block* allocate_storage_block(size_t block_size);
memory_pool* initialize_storage_pool(size_t block_size, size_t blocks_per_pool);
void* allocate_memory_from_pool(memory_pool *pool);
void expand_memory_pool(memory_pool *pool);

#endif // MEMORY_H

#include "../include/memory.h"
#include "../include/globals.h"
// Function definitions for memory management

memory_pool_collection* initialize_memory_pool_collection() {
    memory_pool_collection *pools = (memory_pool_collection *)malloc(sizeof(memory_pool_collection));
    if (!pools) {
        perror("Failed to allocate memory for memory pool collection");
        exit(EXIT_FAILURE);
    }
    pools->feature_counts_pool = initialize_storage_pool(dynamic_struct_sizes.feature_counts, FEATURE_BARCODE_BLOCK_SIZE);
    pools->feature_umi_counts_pool = initialize_storage_pool(dynamic_struct_sizes.feature_umi_counts, UMI_STORAGE_BLOCK);
    pools->feature_sequences_pool = initialize_storage_pool(dynamic_struct_sizes.feature_sequences, FEATURE_SEQUENCE_BLOCK_SIZE);
    pools->unmatched_barcodes_features_block_pool = initialize_storage_pool(dynamic_struct_sizes.unmatched_barcodes_features_block, BARCODE_STORAGE_BLOCK_SIZE);
    pools->cb_counts_pool = initialize_storage_pool(dynamic_struct_sizes.cb_probe_counts, CB_COUNTS_BLOCK_SIZE);
    return pools;
}

void free_memory_pool_collection(memory_pool_collection *pools) {
    free_memory_pool_storage(pools->feature_counts_pool);
    free_memory_pool_storage(pools->feature_umi_counts_pool);
    free_memory_pool_storage(pools->feature_sequences_pool);
    free_memory_pool_storage(pools->unmatched_barcodes_features_block_pool);
    free_memory_pool_storage(pools->cb_counts_pool);
    free(pools);
}

void free_memory_pool_storage(memory_pool *pool) {
    storage_block *current = pool->first_block;
    while (current != NULL) {
        storage_block *next = current->next;
        free(current->storage);
        free(current);
        current = next;
    }
    free(pool);
}

storage_block* allocate_storage_block(size_t block_size) {
    storage_block *new_block = (storage_block *)malloc(sizeof(storage_block));
    if (!new_block) {
        fprintf(stderr, "Failed to allocate new storage block %f GB\n", (double)block_size/1000000000);
        return NULL;
    }
    new_block->storage = (unsigned char *)malloc(block_size);
    if (!new_block->storage) {
        fprintf(stderr, "Failed to allocate new storage block %f GB\n", (double)block_size/1000000000);
        return NULL;
    }
    new_block->next = NULL;
    memset(new_block->storage, 0, block_size);
    return new_block;
}

memory_pool* initialize_storage_pool(size_t block_size, size_t blocks_per_pool) {
    memory_pool *pool = (memory_pool *)malloc(sizeof(memory_pool));
    if (!pool) {
        perror("Failed to allocate memory for memory pool");
        exit(EXIT_FAILURE);
    }
    pool->block_size = block_size;
    pool->blocks_per_pool = blocks_per_pool;
    pool->free_blocks = blocks_per_pool;
    pool->first_block = allocate_storage_block(block_size * blocks_per_pool);
    pool->current_block = pool->first_block;
    pool->next_free = pool->current_block->storage;
    return pool;
}

void* allocate_memory_from_pool(memory_pool *pool) {
    if (pool->free_blocks == 0) {
        expand_memory_pool(pool);
    }
    pool->free_blocks--;
    void *ret = pool->next_free;
    pool->next_free += pool->block_size;
    return ret;
}

void expand_memory_pool(memory_pool *pool) {
    storage_block *new_block = allocate_storage_block(pool->block_size * pool->blocks_per_pool);
    pool->current_block->next = new_block;
    pool->current_block = new_block;
    pool->next_free = pool->current_block->storage;
    pool->free_blocks += pool->blocks_per_pool;
}
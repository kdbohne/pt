#pragma once

#include "common.h"

#include <cstddef>
#include <list>

#ifndef L1_CACHE_LINE_SIZE
#define L1_CACHE_LINE_SIZE 64
#endif

#define ARENA_ALLOC(arena, Type) new (arena.alloc(sizeof(Type))) Type

class alignas(L1_CACHE_LINE_SIZE) MemoryArena
{
public:
    MemoryArena(size_t block_size = 262144) : block_size(block_size) {}
    ~MemoryArena();

    void *alloc(size_t bytes_count);
    void reset();

private:
    MemoryArena(const MemoryArena &) = delete;
    MemoryArena &operator=(const MemoryArena &) = delete;

    size_t block_size;
    size_t current_block_pos = 0, current_alloc_size = 0;
    uint8_t *current_block = nullptr;
    std::list<std::pair<size_t, uint8_t *>> used_blocks, available_blocks;
};

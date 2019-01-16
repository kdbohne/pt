#include "memory.h"

#include <cstdlib>

static void *alloc_aligned(size_t size)
{
    // TODO: HAVE_ALIGNED_MALLOC/HAVE_POSIX_MEMALIGN compilation flags
//    _aligned_malloc(size, L1_CACHE_LINE_SIZE);

    void *ptr;
    if (posix_memalign(&ptr, L1_CACHE_LINE_SIZE, size) != 0)
        ptr = nullptr;

    return ptr;
}

template<typename T>
static T *alloc_aligned(size_t count)
{
    return (T *)alloc_aligned(count * sizeof(T));
}

static void free_aligned(void *ptr)
{
    if (!ptr)
        return;

    // TODO: HAVE_ALIGNED_MALLOC compilation flag
//    _aligned_free(ptr);

    free(ptr);
}

MemoryArena::~MemoryArena()
{
    free_aligned(current_block);

    for (auto &block : used_blocks)
        free_aligned(block.second);
    for (auto &block : available_blocks)
        free_aligned(block.second);
}

void *MemoryArena::alloc(size_t bytes_count)
{
    int align = alignof(std::max_align_t);
    bytes_count = (bytes_count + align - 1) & ~(align - 1);

    // Find space within an existing block or allocate a new block.
    if (current_block_pos + bytes_count > current_alloc_size)
    {
        if (current_block)
        {
            used_blocks.push_back(std::make_pair(current_alloc_size, current_block));
            current_block = nullptr;
            current_alloc_size = 0;
        }

        // Attempt to find a block that has enough space.
        for (auto it = available_blocks.begin(); it != available_blocks.end(); ++it)
        {
            if (it->first >= bytes_count)
            {
                current_alloc_size = it->first;
                current_block = it->second;
                available_blocks.erase(it);
                break;
            }
        }

        // None of the existing blocks have enough space, so allocate a new one.
        if (!current_block)
        {
            current_alloc_size = std::max(bytes_count, block_size);
            current_block = alloc_aligned<uint8_t>(current_alloc_size);
        }

        current_block_pos = 0;
    }

    void *result = current_block + current_block_pos;
    current_block_pos += bytes_count;

    return result;
}

void MemoryArena::reset()
{
    current_block_pos = 0;
    available_blocks.splice(available_blocks.begin(), used_blocks);
}

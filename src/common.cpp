#include "common.h"
#include <cstdio>

void report_(const char *prefix, bool fatal, const char *format, ...)
{
    std::printf("[%s] ", prefix);

    va_list args;
    va_start(args, format);
    std::vprintf(format, args); // TODO LOG
    va_end(args);

    std::printf("\n");

    if (fatal)
        __builtin_trap();
}

#pragma once

#include <stdint.h>

#define UNUSED(x) ((void)(x))

void report_(const char *prefix, bool fatal, const char *format, ...) __attribute__((__format__(__printf__, 3, 4)));

#define error(format, ...) report_("ERROR", false, format, __VA_ARGS__)
#define fatal(format, ...) report_("FATAL", true, format, __VA_ARGS__)
#define assert(expr) \
    ((expr) ? (void)0 : fatal("(%s:%d) Assertion failed: %s", __FILE__, __LINE__, #expr))

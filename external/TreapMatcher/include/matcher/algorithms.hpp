#pragma once

#include "common.hpp"
#include "hasher.hpp"

#include <algorithm>
#include <assert.h>

#define BROTLI_MAX_STATIC_DICTIONARY_MATCH_LEN 37

#define MAX_ZOPFLI_LEN_QUALITY_11 325

size_t FindMatchLengthWithLimit(const uint8_t* s1, const uint8_t* s2, size_t limit);

#define BROTLI_WINDOW_GAP 16
/* Maximum distance, see section 9.1. of the spec. */
size_t MaxBackwardLimit(int lgwin);
#include <matcher/algorithms.hpp>
#include <matcher/memory.hpp>

static const uint32_t kInvalidMatch = 0xfffffff;

size_t FindMatchLengthWithLimit(const uint8_t* s1, const uint8_t* s2, size_t limit)
{
  size_t matched = 0;
  size_t limit2 = (limit >> 3) + 1;  /* + 1 is for pre-decrement in while */
  while (--limit2) {
    if (BROTLI_UNALIGNED_LOAD64(s2) == BROTLI_UNALIGNED_LOAD64(s1 + matched)) {
      s2 += 8;
      matched += 8;
    } else {
      uint64_t x = BROTLI_UNALIGNED_LOAD64(s2) ^ BROTLI_UNALIGNED_LOAD64(s1 + matched);
      size_t matching_bits = (size_t)__builtin_ctzll(x);
      matched += matching_bits >> 3;
      return matched;
    }
  }
  limit = (limit & 7) + 1;  /* + 1 is for pre-decrement in while */
  while (--limit) {
    if (s1[matched] == *s2) {
      ++s2;
      ++matched;
    } else {
      return matched;
    }
  }
  return matched;
}

size_t MaxBackwardLimit(int lgwin) {
  return (1u << lgwin) - BROTLI_WINDOW_GAP;
}
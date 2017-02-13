#include <matcher/memory.hpp>
#include <cstring>

uint32_t BROTLI_UNALIGNED_LOAD32(const void *p) {
  uint32_t t;
  std::memcpy(&t, p, sizeof t);
  return t;
}

uint64_t BROTLI_UNALIGNED_LOAD64(const void *p) {
  uint64_t t;
  memcpy(&t, p, sizeof t);
  return t;
}

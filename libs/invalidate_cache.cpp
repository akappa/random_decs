#include <invalidate_cache.hpp>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <unistd.h>
#include <vector>

#define TURN_OFF_DEAD_STORE_OPT(x) __asm__ volatile("" : "+g"(x));

namespace {
std::vector<size_t> cache_sizes()
{
  auto levels = {
    _SC_LEVEL1_DCACHE_SIZE,
    _SC_LEVEL2_CACHE_SIZE,
    _SC_LEVEL3_CACHE_SIZE, 
    _SC_LEVEL4_CACHE_SIZE
  };
  std::vector<size_t> sizes;
  for (auto level : levels) {
    auto size = sysconf(level);
    if (size > 0) {
      sizes.push_back(size);
    } else {
      break;
    }
  }
  return sizes;
}
}

void wipe_caches()
{
  auto max_cache = cache_sizes().back();
  std::vector<std::uint8_t> bytes(max_cache, 0);
  std::iota(bytes.begin(), bytes.end(), 0UL);
  std::default_random_engine g(42); 
  std::shuffle(bytes.begin(), bytes.end(), g);
  std::uint8_t sink;
  TURN_OFF_DEAD_STORE_OPT(sink)
  for (auto i : bytes) {
    sink = bytes[i];
  }
}

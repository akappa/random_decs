#pragma once

#include <cstdint>
#include <limits>

// Importing basic data types
using std::uint32_t;
using std::uint8_t;
using std::size_t;

class BackwardMatch {
public:
  uint32_t distance;
  uint32_t length_and_code;

  BackwardMatch()
    : distance(std::numeric_limits<uint32_t>::max()),
      length_and_code(std::numeric_limits<uint32_t>::max())
  { }

  BackwardMatch(size_t dist, size_t len)
    : distance(dist), length_and_code(len << 5)
  { }

  size_t match_length() const {
    return length_and_code >> 5;
  }
};
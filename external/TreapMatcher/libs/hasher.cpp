#include <matcher/algorithms.hpp>
#include <matcher/hasher.hpp>
#include <matcher/memory.hpp>

#include <assert.h>
#include <cstring>
#include <iostream>

constexpr size_t HashToBinaryTree::MAX_TREE_COMP_LENGTH;
constexpr size_t HashToBinaryTree::BUCKET_BITS;
constexpr size_t HashToBinaryTree::BUCKET_SIZE;
constexpr size_t HashToBinaryTree::MAX_TREE_SEARCH_DEPTH;
constexpr size_t HashToBinaryTree::MAX_NUM_MATCHES_HASH;

BackwardMatch* HashToBinaryTree::StoreAndFindMatches(
  const uint8_t* const data, const size_t cur_ix, 
  const size_t ring_buffer_mask, const size_t max_length, const size_t max_backward,
  size_t* const best_len, BackwardMatch* matches
)
{
//  if (cur_ix == 2097662) {
//    printf("Here we are\n");
//  }
  const size_t   cur_ix_masked      = cur_ix & ring_buffer_mask;
  const size_t   max_comp_len       = std::min<size_t>(max_length, MAX_TREE_COMP_LENGTH);
  const bool     should_reroot_tree = max_length >= MAX_TREE_COMP_LENGTH;
  const uint32_t key                = HashBytesHash(&data[cur_ix_masked]);
        size_t   prev_ix            = buckets_[key];
  /* The forest index of the rightmost node of the left subtree of the new
     root, updated as we traverse and re-root the tree of the hash bucket. */
        size_t   node_left          = left_child(cur_ix);
  /* The forest index of the leftmost node of the right subtree of the new
     root, updated as we traverse and re-root the tree of the hash bucket. */
        size_t   node_right         = right_child(cur_ix);
  /* The match length of the rightmost node of the left subtree of the new
     root, updated as we traverse and re-root the tree of the hash bucket. */
        size_t   best_len_left      = 0;
  /* The match length of the leftmost node of the right subtree of the new
     root, updated as we traverse and re-root the tree of the hash bucket. */
        size_t   best_len_right     = 0;

  if (should_reroot_tree) {
    buckets_[key] = cur_ix;
  }
  for (auto depth_remaining = MAX_TREE_SEARCH_DEPTH; ; --depth_remaining) {
    const size_t backward       = cur_ix - prev_ix;
    const size_t prev_ix_masked = prev_ix & ring_buffer_mask;
    if (backward == 0 || backward > max_backward || depth_remaining == 0) {
      if (should_reroot_tree) {
        forest_[node_left] = forest_[node_right] = invalid_pos_;
      }
      break;
    }
    {
      const size_t cur_len = std::min<size_t>(best_len_left, best_len_right);
      assert(cur_len <= MAX_TREE_COMP_LENGTH);
      auto len = cur_len + FindMatchLengthWithLimit(
        &data[cur_ix_masked + cur_len],
        &data[prev_ix_masked + cur_len],
        max_length - cur_len
      );
      assert(0 == std::memcmp(&data[cur_ix_masked], &data[prev_ix_masked], len));
      if (matches && len > *best_len) {
        *best_len = len;
        *matches++ = BackwardMatch(backward, len);
      }
      if (len >= max_comp_len) {
        if (should_reroot_tree) {
          forest_[node_left]  = forest_[left_child(prev_ix)];
          forest_[node_right] = forest_[right_child(prev_ix)];
        }
        break;
      }
      if (data[cur_ix_masked + len] > data[prev_ix_masked + len]) {
        best_len_left = len;
        if (should_reroot_tree) {
          forest_[node_left] = prev_ix;
        }
        node_left = right_child(prev_ix);
        prev_ix = forest_[node_left];
      } else {
        best_len_right = len;
        if (should_reroot_tree) {
          forest_[node_right] = prev_ix;
        }
        node_right = left_child(prev_ix);
        prev_ix = forest_[node_right];
      }
    }
  }
  return matches;
}

static char *copy_limited(const char *input, const size_t length, char *buffer, size_t buf_len)
{
  int min = std::min<int>(length, buf_len);
  std::copy(input, input + min, buffer);
  buffer[min < buf_len ? min + 1 : buf_len - 1] = '\0';
  return buffer;
}

size_t HashToBinaryTree::FindAllMatches(
    const uint8_t* data, const size_t ring_buffer_mask, const size_t cur_ix,
    const size_t max_length, const size_t max_backward,
    BackwardMatch* matches
)
{
  BackwardMatch* const orig_matches = matches;
  const size_t cur_ix_masked = cur_ix & ring_buffer_mask;
  size_t best_len = 1;
  const size_t short_match_max_backward = 64;
  size_t stop = cur_ix - short_match_max_backward;
  size_t i;
  if (cur_ix < short_match_max_backward) { stop = 0; }

//  char buffer[16];
//  std::cout << "data = " << copy_limited(reinterpret_cast<const char*>(data), max_length, buffer, 16) << "\n"
//            << "mask = " << ring_buffer_mask << "\n"
//            << "cur_ix = " << cur_ix << "\n"
//            << "max_length = " << max_length << "\n"
//            << "max_backward = " << max_backward << std::endl;

  for (i = cur_ix - 1; i > stop && best_len <= 2; --i) {
    size_t prev_ix = i;
    const size_t backward = cur_ix - prev_ix;
    if (backward > max_backward) {
      break;
    }
    prev_ix &= ring_buffer_mask;
    if (data[cur_ix_masked] != data[prev_ix] ||
        data[cur_ix_masked + 1] != data[prev_ix + 1]) {
      continue;
    }
    {
      const size_t len =
          FindMatchLengthWithLimit(&data[prev_ix], &data[cur_ix_masked],
                                   max_length);
      if (len > best_len) {
        best_len = len;
        *matches++ = BackwardMatch(backward, len);
      }
    }
  }
  if (best_len < max_length) {
    matches = StoreAndFindMatches(data, cur_ix, ring_buffer_mask,
        max_length, max_backward, &best_len, matches);
  }
//  for (auto m_it = orig_matches; m_it != matches; ++m_it) {
//    std::cout << "M: dist = " << m_it->distance
//              << ", length_and_code = " << m_it->length_and_code << std::endl;
//  }
//  std::cout << "----------------------------" << std::endl;

  return matches - orig_matches;
}

void HashToBinaryTree::StoreHashRange(
  const uint8_t *data, const size_t mask, const size_t ix_start,
  const size_t ix_end
)
{
  size_t i = ix_start + 63 <= ix_end ? ix_end - 63 : ix_start;
  for (; i < ix_end; ++i) {
    StoreHash(data, mask, i);
  }
}

void HashToBinaryTree::StoreHash(
  const uint8_t *data, const size_t mask, const size_t ix
) {
  /* Maximum distance is window size - 16, see section 9.1. of the spec. */
  const size_t max_backward = window_mask_ - BROTLI_WINDOW_GAP + 1;
  StoreAndFindMatches(data, ix, mask, MAX_TREE_COMP_LENGTH, max_backward, NULL, NULL);
}

void HashToBinaryTree::StitchToPreviousBlock(
    size_t num_bytes, size_t position, const uint8_t* ringbuffer,
    size_t ringbuffer_mask
)
{
  if (num_bytes >= HashTypeLength() - 1 and position >= MAX_TREE_COMP_LENGTH) {
    /* Store the last `MAX_TREE_COMP_LENGTH - 1` positions in the hasher.
       These could not be calculated before, since they require knowledge
       of both the previous and the current block. */
    const size_t i_start = position - MAX_TREE_COMP_LENGTH + 1;
    const size_t i_end = std::min<size_t>(position, i_start + num_bytes);
    size_t i;
    for (i = i_start; i < i_end; ++i) {
      /* Maximum distance is window size - 16, see section 9.1. of the spec.
         Furthermore, we have to make sure that we don't look further back
         from the start of the next block than the window size, otherwise we
         could access already overwritten areas of the ring-buffer. */
      const size_t max_backward = window_mask_ - std::max<size_t>(BROTLI_WINDOW_GAP - 1, position - i);
      /* We know that i + MAX_TREE_COMP_LENGTH <= position + num_bytes, i.e. the
         end of the current block and that we have at least
         MAX_TREE_COMP_LENGTH tail in the ring-buffer. */
      StoreAndFindMatches(ringbuffer, i, ringbuffer_mask,
          MAX_TREE_COMP_LENGTH, max_backward, NULL, NULL);
    }
  }
}

uint32_t HashToBinaryTree::HashBytesHash(const uint8_t *data)
{
  uint32_t h = BROTLI_UNALIGNED_LOAD32(data) * kHashMul32;
  /* The higher bits contain more mixture from the multiplication,
     so we take our results from there. */
  return h >> (32 - BUCKET_BITS);
}

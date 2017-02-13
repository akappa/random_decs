#pragma once

#include "common.hpp"

#include <vector>

/* A (forgetful) hash table where each hash bucket contains a binary tree of
   sequences whose first 4 bytes share the same hash code.
   Each sequence is MAX_TREE_COMP_LENGTH long and is identified by its starting
   position in the input data. The binary tree is sorted by the lexicographic
   order of the sequences, and it is also a max-heap with respect to the
   starting positions. */
class HashToBinaryTree {
private:
  static constexpr size_t MAX_TREE_COMP_LENGTH{128};
  static constexpr size_t BUCKET_BITS{17};
  static constexpr size_t BUCKET_SIZE{1 << BUCKET_BITS};
  static constexpr size_t MAX_TREE_SEARCH_DEPTH{64};
  static constexpr size_t MAX_NUM_MATCHES_HASH{64 + MAX_TREE_SEARCH_DEPTH};

  /* The window size minus 1 */
  size_t window_mask_;

  /* Hash table that maps the 4-byte hashes of the sequence to the last
     position where this hash was found, which is the root of the binary
     tree of sequences that share this hash bucket. */
  uint32_t buckets_[BUCKET_SIZE];

  /* The union of the binary trees of each hash bucket. The root of the tree
     corresponding to a hash is a sequence starting at buckets_[hash] and
     the left and right children of a sequence starting at pos are
     forest_[2 * pos] and forest_[2 * pos + 1]. */
  std::vector<uint32_t> forest_;

  /* A position used to mark a non-existent sequence, i.e. a tree is empty if
     its root is at invalid_pos_ and a node is a leaf if both its children
     are at invalid_pos_. */
  uint32_t invalid_pos_;

  bool is_dirty_;

  size_t left_child(size_t pos)
  {
    return 2 * (pos & window_mask_);
  }

  size_t right_child(size_t pos)
  {
    return left_child(pos) + 1;
  }

  /* kHashMul32 multiplier has these properties:
     * The multiplier must be odd. Otherwise we may lose the highest bit.
     * No long streaks of ones or zeros.
     * There is no effort to ensure that it is a prime, the oddity is enough
       for this use.
     * The number has been tuned heuristically against compression benchmarks. */
  static const uint32_t kHashMul32 = 0x1e35a7bd;


public:
  HashToBinaryTree()  { reset(); }
  HashToBinaryTree(size_t lg_win, size_t position, size_t bytes, bool is_last)
    : window_mask_((1u << lg_win) - 1u), invalid_pos_(0 - window_mask_)
  {
    for (auto i = 0; i < BUCKET_SIZE; i++) {
      this->buckets_[i] = invalid_pos_;
    }
    auto num_nodes = (position == 0 && is_last) ? bytes : window_mask_ + 1;
    if (num_nodes > forest_size()) {
      this->forest_.resize(num_nodes * 2);
    }
    this->is_dirty_ = false;
  }

  size_t forest_size() { return this->forest_.size() / 2; }

  static size_t HashTypeLength() { return 4; }

  static size_t StoreLookahead() { return MAX_TREE_COMP_LENGTH; }

  static size_t MaxMatches() { return MAX_NUM_MATCHES_HASH; }

  void reset() { is_dirty_ = true; }

  bool dirty() const { return is_dirty_; }

  /* Stores the hash of the next 4 bytes and in a single tree-traversal, the
     hash bucket's binary tree is searched for matches and is re-rooted at the
     current position.
  
     If less than MAX_TREE_COMP_LENGTH data is available, the hash bucket of the
     current position is searched for matches, but the state of the hash table
     is not changed, since we can not know the final sorting order of the
     current (incomplete) sequence.
  
     This function must be called with increasing cur_ix positions. */
  BackwardMatch* StoreAndFindMatches(
    const uint8_t* const data, const size_t cur_ix, 
    const size_t ring_buffer_mask, const size_t max_length, const size_t max_backward,
    size_t* const best_len, BackwardMatch* matches
  );

  /* Finds all backward matches of &data[cur_ix & ring_buffer_mask] up to the
     length of max_length and stores the position cur_ix in the hash table.
  
     Sets *num_matches to the number of matches found, and stores the found
     matches in matches[0] to matches[*num_matches - 1]. The matches will be
     sorted by strictly increasing length and (non-strictly) increasing
     distance. */
  size_t FindAllMatches(
    const uint8_t* data, const size_t ring_buffer_mask, const size_t cur_ix,
    const size_t max_length, const size_t max_backward,
    BackwardMatch* matches
  );

  void StoreHashRange(
    const uint8_t *data, const size_t mask, const size_t ix_start,
    const size_t ix_end
  );

  /* Stores the hash of the next 4 bytes and re-roots the binary tree at the
     current sequence, without returning any matches.
     REQUIRES: ix + MAX_TREE_COMP_LENGTH <= end-of-current-block */
  void StoreHash(
    const uint8_t *data, const size_t mask, const size_t ix
  );

  void StitchToPreviousBlock(
    size_t num_bytes, size_t position, const uint8_t* ringbuffer,
    size_t ringbuffer_mask
  );

  static uint32_t HashBytesHash(const uint8_t *data);
};
#include <matcher/ring_buffer.hpp>

#include <cstring>
using std::memcpy;

void RingBufferSetup(int lgwin, int lgblock, RingBuffer* rb)
{
  /* Allocate at least lgwin + 1 bits for the ring buffer so that the newly
     added block fits there completely and we still get lgwin bits and at least
     read_block_size_bits + 1 bits because the copy tail length needs to be
     smaller than ring-buffer size. */
  int window_bits              = 1 + std::max<int>(lgwin, lgblock);
  int tail_bits                = lgblock;
  *(uint32_t*)&rb->size_       = 1u << window_bits;
  *(uint32_t*)&rb->mask_       = (1u << window_bits) - 1;
  *(uint32_t*)&rb->tail_size_  = 1u << tail_bits;
  *(uint32_t*)&rb->total_size_ = rb->size_ + rb->tail_size_;
}


void CopyInputToRingBuffer(
  RingBuffer* ringbuffer, const size_t input_size, const uint8_t* input_buffer
)
{
  RingBufferWrite(input_buffer, input_size, ringbuffer);
  /* TL;DR: If needed, initialize 7 more bytes in the ring buffer to make the
     hashing not depend on uninitialized data. This makes compression
     deterministic and it prevents uninitialized memory warnings in Valgrind.
     Even without erasing, the output would be valid (but nondeterministic).

     Background information: The compressor stores short (at most 8 bytes)
     substrings of the input already read in a hash table, and detects
     repetitions by looking up such substrings in the hash table. If it
     can find a substring, it checks whether the substring is really there
     in the ring buffer (or it's just a hash collision). Should the hash
     table become corrupt, this check makes sure that the output is
     still valid, albeit the compression ratio would be bad.

     The compressor populates the hash table from the ring buffer as it's
     reading new bytes from the input. However, at the last few indexes of
     the ring buffer, there are not enough bytes to build full-length
     substrings from. Since the hash table always contains full-length
     substrings, we erase with dummy zeros here to make sure that those
     substrings will contain zeros at the end instead of uninitialized
     data.

     Please note that erasing is not necessary (because the
     memory region is already initialized since he ring buffer
     has a `tail' that holds a copy of the beginning,) so we
     skip erasing if we have already gone around at least once in
     the ring buffer.

     Only clear during the first round of ring-buffer writes. On
     subsequent rounds data in the ring-buffer would be affected. */
  if (ringbuffer->pos_ <= ringbuffer->mask_) {
    /* This is the first time when the ring buffer is being written.
       We clear 7 bytes just after the bytes that have been copied from
       the input buffer.

       The ring-buffer has a "tail" that holds a copy of the beginning,
       but only once the ring buffer has been fully written once, i.e.,
       pos <= mask. For the first time, we need to write values
       in this tail (where index may be larger than mask), so that
       we have exactly defined behavior and don't read uninitialized
       memory. Due to performance reasons, hashing reads data using a
       LOAD64, which can go 7 bytes beyond the bytes written in the
       ring-buffer. */
    std::memset(ringbuffer->buffer_ + ringbuffer->pos_, 0, 7);
  }
}

/* Push bytes into the ring buffer. */
void RingBufferWrite(
    const uint8_t *bytes, size_t n, RingBuffer* rb
) {
  if (rb->pos_ == 0 && n < rb->tail_size_) {
    /* Special case for the first write: to process the first block, we don't
       need to allocate the whole ring-buffer and we don't need the tail
       either. However, we do this memory usage optimization only if the
       first write is less than the tail size, which is also the input block
       size, otherwise it is likely that other blocks will follow and we
       will need to reallocate to the full size anyway. */
    rb->pos_ = (uint32_t)n;
    RingBufferInitBuffer(rb->pos_, rb);
    memcpy(rb->buffer_, bytes, n);
    return;
  }
  if (rb->cur_size_ < rb->total_size_) {
    /* Lazily allocate the full buffer. */
    RingBufferInitBuffer(rb->total_size_, rb);
    /* Initialize the last two bytes to zero, so that we don't have to worry
       later when we copy the last two bytes to the first two positions. */
    rb->buffer_[rb->size_ - 2] = 0;
    rb->buffer_[rb->size_ - 1] = 0;
  }
  {
    const size_t masked_pos = rb->pos_ & rb->mask_;
    /* The length of the writes is limited so that we do not need to worry
       about a write */
    RingBufferWriteTail(bytes, n, rb);
    if (masked_pos + n <= rb->size_) {
      /* A single write fits. */
      memcpy(&rb->buffer_[masked_pos], bytes, n);
    } else {
      /* Split into two writes.
         Copy into the end of the buffer, including the tail buffer. */
      memcpy(&rb->buffer_[masked_pos], bytes, std::min<size_t>(n, rb->total_size_ - masked_pos));
      /* Copy into the beginning of the buffer */
      memcpy(&rb->buffer_[0], bytes + (rb->size_ - masked_pos), n - (rb->size_ - masked_pos));
    }
  }
  rb->buffer_[-2] = rb->buffer_[rb->size_ - 2];
  rb->buffer_[-1] = rb->buffer_[rb->size_ - 1];
  rb->pos_ += (uint32_t)n;
  if (rb->pos_ > (1u << 30)) {
    /* Wrap, but preserve not-a-first-lap feature. */
    rb->pos_ = (rb->pos_ & ((1u << 30) - 1)) | (1u << 30);
  }
}

void RingBufferInitBuffer(const uint32_t buflen, RingBuffer* rb)
{
  static const size_t kSlackForEightByteHashingEverywhere = 7;
  rb->data_.resize(2 + buflen + kSlackForEightByteHashingEverywhere);
  size_t i;
  rb->cur_size_ = buflen;
  rb->buffer_ = rb->data_.data() + 2;
  rb->buffer_[-2] = rb->buffer_[-1] = 0;
  for (i = 0; i < kSlackForEightByteHashingEverywhere; ++i) {
    rb->buffer_[rb->cur_size_ + i] = 0;
  }
}

void RingBufferWriteTail(const uint8_t *bytes, size_t n, RingBuffer* rb)
{
  const size_t masked_pos = rb->pos_ & rb->mask_;
  if (masked_pos < rb->tail_size_) {
    /* Just fill the tail buffer with the beginning data. */
    const size_t p = rb->size_ + masked_pos;
    memcpy(&rb->buffer_[p], bytes, std::min<size_t>(n, rb->tail_size_ - masked_pos));
  }
}

void RingBufferInit(RingBuffer* rb) {
  rb->cur_size_ = 0;
  rb->pos_ = 0;
  rb->buffer_ = 0;
}
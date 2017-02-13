#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

using std::uint32_t;
using std::uint8_t;

/* A RingBuffer(window_bits, tail_bits) contains `1 << window_bits' bytes of
   data in a circular manner: writing a byte writes it to:
     `position() % (1 << window_bits)'.
   For convenience, the RingBuffer array contains another copy of the
   first `1 << tail_bits' bytes:
     buffer_[i] == buffer_[i + (1 << window_bits)], if i < (1 << tail_bits),
   and another copy of the last two bytes:
     buffer_[-1] == buffer_[(1 << window_bits) - 1] and
     buffer_[-2] == buffer_[(1 << window_bits) - 2]. */
typedef struct RingBuffer {
  /* Size of the ring-buffer is (1 << window_bits) + tail_size_. */
  uint32_t size_;
  uint32_t mask_;
  uint32_t tail_size_;
  uint32_t total_size_;

  uint32_t cur_size_;
  /* Position to write in the ring buffer. */
  uint32_t pos_;
  /* The actual ring buffer containing the copy of the last two bytes, the data,
     and the copy of the beginning as a tail. */
  std::vector<uint8_t> data_;
  /* The start of the ring-buffer. */
  uint8_t *buffer_;
} RingBuffer;

void RingBufferInit(RingBuffer* rb);

void RingBufferSetup(int lgwin, int lgblock, RingBuffer* rb);

/*
   Copies the given input data to the internal ring buffer of the compressor.
   No processing of the data occurs at this time and this function can be
   called multiple times before calling WriteBrotliData() to process the
   accumulated input. At most input_block_size() bytes of input data can be
   copied to the ring buffer, otherwise the next WriteBrotliData() will fail.
 */
void CopyInputToRingBuffer(
  RingBuffer* ringbuffer, const size_t input_size, const uint8_t* input_buffer
);

/* Push bytes into the ring buffer. */
void RingBufferWrite(
    const uint8_t *bytes, size_t n, RingBuffer* rb
);

/* Allocates or re-allocates data_ to the given length + plus some slack
   region before and after. Fills the slack regions with zeros. */
void RingBufferInitBuffer(const uint32_t buflen, RingBuffer* rb);

void RingBufferWriteTail(const uint8_t *bytes, size_t n, RingBuffer* rb);
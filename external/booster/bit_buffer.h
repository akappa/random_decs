/* ========================================================
   routines for writing bits in a internal memory buffer
   ======================================================== */

#define _bb_MIN_BUFFER_SIZE 10  // initial buffer size;

typedef struct {
  int32_t pos;            // current position in buffer
  int32_t offset;         // number of free bits in current position
  int32_t size;           // capacity of buffer           
  uint32_t *buffer;       // pointer to the buffer 
} _bb_bit_buffer;


// ============= error functions ===========

static void _bb_fatal_error(char *s)
{
  fprintf(stderr,"%s",s); exit(1);
}

static void _bb_out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f); exit(1);
}



// ****** creation of a minimal bit_biffer
static void _bb_init_bit_buffer(_bb_bit_buffer *b)
{
  b->pos=0;                // beginning of buffer 
  b->offset=32;            // free bits 
  b->size=_bb_MIN_BUFFER_SIZE;
  b->buffer = (uint32_t *) malloc(b->size*sizeof(uint32_t));
  if(b->buffer==NULL)
    _bb_out_of_mem("_bb_init_bit_buffer");
  b->buffer[0]=0;          // clean current pos;
}

// ******* double the size of a bit_buffer
static void _bb_enlarge_bit_buffer(_bb_bit_buffer *b)
{
  assert(b->size>0);
  assert(b->buffer!=NULL);
  b->size*=2;
  b->buffer = (uint32_t *) realloc(b->buffer,b->size*sizeof(uint32_t));
  if(b->buffer==NULL)
    _bb_out_of_mem("_bb_enlarge_bit_buffer");
}

// ***** deallocate bit buffer
static void _bb_free_bit_buffer(_bb_bit_buffer *b)
{

  if(b->buffer==NULL) 
    _bb_fatal_error("Illegal free in _bb_free_bit_buffer()\n");
  
  free(b->buffer);
  b->buffer=NULL;
  b->pos=b->size=b->offset=0;
}


/* ****************************************************************
    write bits [0,n-1] of v (MSB first) in buffer b
    **************************************************************** */
void _bb_write_bits32(_bb_bit_buffer *b, int n, uint32_t v)
{
  int extra;

  assert(n>0 && n<=32);
  // make sure bits n+1 ... 31 are zero
  if(n<32)
    v &= ((1<<n)-1);
     
  // make sure there is enough space
  if(b->pos+1>=b->size) {
    _bb_enlarge_bit_buffer(b);
    assert(b->pos+1 < b->size);
  }
  assert(b->offset!=0);

  // write bits: easy if there are enough free bits 
  if(n<=b->offset) {
    b->buffer[b->pos] |= v << (b->offset-n);
    b->offset -= n;
    if(b->offset==0) {
      b->buffer[++b->pos]=0;  // increment b->pos
      b->offset=32;
    }
  }
  // if not enough space, write in pos and pos+1
  else {
    extra = n-b->offset;
    b->buffer[b->pos] |= v >> extra;      // write b->offset bits
    b->pos++;                             // go to  next pos
    b->buffer[b->pos] =  v <<(32-extra);  // write extra bits to next pos 
    b->offset = 32-extra;
  }
} 

// ******* write up to 64 bits
void _bb_write_bits64(_bb_bit_buffer *b, int n, uint64_t v)
{
  assert(n>0 && n<=64);
  if(n<=32)
    _bb_write_bits32(b,n,(uint32_t) v);
  else {
    _bb_write_bits32(b,n-32,(uint32_t) (v>>32));
    _bb_write_bits32(b,32,(uint32_t) v);
  }
}

// ***** output bit buffer to file f via fbit_write()
static void _bb_bit_buffer2file(_bb_bit_buffer *b, FILE *f,
				void (* fbit_write32)(FILE *, int, uint32_t))
{
  int n;

  for(n=0;n<b->pos;n++) {
    fbit_write32(f,32,b->buffer[n]);
    // fprintf(stderr,"%08x  ",b->buffer[n]);
  }
  if(b->offset<32) 
    fbit_write32(f,32-b->offset,b->buffer[b->pos]>>b->offset);
}




#if 0
// ******* write gammacode for n to bit_buffer b ***********
static inline int _bb_write_gammacode(_bb_bit_buffer *b, uint32_t n)
{
  int m,bits;

  assert(n>0);
  bits=0; m=n;
  while (m) { bits++; m >>= 1; }
  assert(bits>0);
  // b is the number of digits in n
  if(bits>1)
    _bb_write_bits(0,bits-1);
  _bb_write_bits(n,bits);
  return 2*bits-1;  
}

// ******** write bit buffer to file (byte aligned) *********
// ******** return number of bytes written
static int _bb_write_to_file(_bb_bit_buffer *b, FILE *f)
{
  int n,written;

  // compute how many bytes are used
  n=b->pos;
  if(b->offset<32)
    n = n+1;
  // write
  written = fwrite(b->buffer,sizeof(int32_t), n,f);
  if(written!=n) {
    perror("(_bb_write_to_file)");
    exit(1);
  }

  return n;
}
#endif

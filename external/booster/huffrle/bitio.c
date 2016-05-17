/*! \file
    Routines handling the input/output at the bit level
 */

#include <stdlib.h>
#include "huffrle.h"



/* ************************************************************
   "high level" bit input/output procedure
   ************************************************************ */

/*!
   write n using the (N-1)xN scheme. The most significant bit indicates 
   if the representation extend to the successive bytes. 
   Return the total number of bits written */
int fwriteN_1xN(FILE *f, uint32 n, int N)
{
  uint32 t, power;
  int bits=0;

  assert(N<=16 && N>1);
  power=(1<<(N-1));
  do {
    t = n & (power-1);       // takes last N-1 bits
    n = n >> (N-1);
    if(n>0) t |= power;    // set Nth bit
    if(f!=NULL) fbit_write(f,N,t);
    bits += N;
  } while(n>0);
  return bits;
}

/*!
    Read an integer using the N-1xN scheme */
uint32 freadN_1xN(FILE *f, int N)
{
  int i,flag;
  uint32 ris,n,power;

  assert(N<=16 && N>1);
  ris=0;
  flag=1;
  power=1<<(N-1);

  for(i=0;flag!=0;i+=N-1) {
    assert(i+N-1<=32);          // ris is 32 bit
    n= fbit_read(f,N);          // read N bit
    flag= n & power;            // continuation bit on?  
    n &= (power-1);             // set the Nth bit to off  
    ris |= n << i;              // add N-1 bits (from most to less signif.)
  }
  return ris;
}


/*!
   write t (with 0<=t<n) to file f using log-skewed encoding.
   This encoding never uses more than ceil(log n) bits, but it can 
   use less bits when t is small 
   \param f: output file
   \param t: value to be written in f
   \param n: max. value+1 allowed for t
   \return number of bits used to write t;
*/
int fwrite_logskewed(FILE *f, uint32 t, uint32 n)
{
  int k,bits_written=0;
  uint32 u,delta;
  
  assert(n>0);
  assert(t<n);
  k=31; u=1<<k;          // u is 2^k 
  while(1) {
    assert(n>0);
    assert(t<n);
    while(n<u) {
      u >>=1; k--;
    }
    // now u contains the most significant bit of n and 0 elsewere
    if(n==u) {
      if(k>0) { // if k=0, n=2^k=1 t=0 and there is nothing to do;  
        if(f!=NULL) fbit_write(f,k,t);
	bits_written += k;
      } 
      break;
    }
    else {  
      delta = n-u;
      assert(delta>0);
      if(t>=delta) {           // encode t using k+1 bits
        if(f!=NULL) {
	  fbit_write(f,1,1);
	  fbit_write(f,k,t-delta);
	}
	bits_written += k+1;	
        break;
      }
      else {
        if(f!=NULL) 
	  fbit_write(f,1,0);   // write 0
	bits_written++;
        n=delta;               // encode(t,delta);
      }
    }
  }
  return bits_written;
}

/*!
   write t (with a<=t<b) to file f using at most ceil(log(b-a)) bits
   and possibly using less bits when b-t is small 
   \param f: output file
   \param t: value to be written in f
   \param a: min value allowed for t
   \param b: max value+1 allowed for t
   \return number of bits used to write t;
*/
int fwrite_top_range(FILE *f, int32 t, int32 a, int32 b)
{

  assert(t>=a && t < b);
  return fwrite_logskewed(f, b-1-t, b-a);
}

/*!
   write t (with a<=t<b) to file f using at most ceil(log(b-a)) bits
   and possibly using less bits when t-a is small 
   \param f: output file
   \param t: value to be written in f
   \param a: min value allowed for t
   \param b: max value+1 allowed for t
   \return number of bits used to write t;
*/
int fwrite_bottom_range(FILE *f, int32 t, int32 a, int32 b)
{

  assert(t>=a && t < b);
  return fwrite_logskewed(f, t-a, b-a);
}

       

/*!
  read an integer t, (with 0<=t<n) from file f using log-skewed decoding 
*/
uint32 fread_logskewed(FILE *f, uint32 n)
{
  int k,bit;
  uint32 u,delta,t=0;

  assert(n>0);
  k=31; u=1<<k;          // u is 2^k 
  while(1) { 
    assert(n>0);
    while(n<u) {
      u >>=1; k--;
    }
    // now u contains the most significant bit of n and 0 elsewere
    if(n==u) {
      if(k>0)  // if k=0, n=1 and there is nothing to read  
        t += fbit_read(f,k); 
      break;
    }
    else {  
      delta = n-u;
      assert(delta>0);
      bit=fbit_read(f,1);
      if(bit==1) {
        t += delta + fbit_read(f,k);
	break;
      }
      else {                  // bit==0
        n=delta;              // decode(delta);
      }
    }
  }
  return t;
}

/*!
  read from file f an integer t, (with a<=t<b) written
  using fwrite_bottom_range()
*/
int32 fread_bottom_range(FILE *f, int32 a, int32 b)
{
  assert(a<=b);
  return a + fread_logskewed(f, (uint32) b-a);
}

/*!
  read from file f an integer t, (with a<=t<b) written
  using fwrite_top_range()
*/
int32 fread_top_range(FILE *f, int32 a, int32 b)
{
  assert(a<=b);
  return b-1-fread_logskewed(f, (uint32) b-a);
}



/*!
  write a 32 bit uint most significat byte first*/
void uint32_write(FILE *f, uint32 u)
{  
  assert(f!=NULL);
  fbit_write(f, 8, (u>>24) & 0xffL);
  fbit_write(f, 8, (u>>16) & 0xffL);
  fbit_write(f, 8, (u>> 8) & 0xffL);
  fbit_write(f, 8,  u      & 0xffL);
}



/* *****************************************************
   Low-level procedures to read and write bits using a Bit_buffer. 
   Unread/unwritten bits of Bit_buffer are the most significant ones.
   ***************************************************** */
uint32 Bit_buffer;     /* 32 bit buffer */
int  Bit_buffer_size;  /* number of unread/unwritten bits in Bit_buffer */

/*!
   Initializes the Bit_buffer to zero */
void init_bit_buffer(void)
{
  Bit_buffer= (uint32) 0;
  Bit_buffer_size=0;
}

/*!
   Write n (<= 24) bits taken from vv. The content of 
   Bit_buffer is flushed out until it contains <8 bits
   /param n: number of bits to be written 
   /param the least significant bits of v are the bits to be written */
__inline__ static void fbit_write24(FILE *f, int n, uint32 v)  
{
  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);
  assert( v < 1u << (n+1) );

  /* ------- add n bits to Bit_buffer -------- */
  Bit_buffer_size += n;       // add first, to compute the correct shift
  Bit_buffer |= (v << (32 - Bit_buffer_size));  // compact to end of the buffer

  /* ------- flush Bit_buffer as much as possible ----- */
  while (Bit_buffer_size>=8) {            
    if( putc((Bit_buffer>>24),f) == EOF) {
      perror("Unable to write bits (bit_write24)");
      exit(1);
    }
    Bit_buffer <<= 8;                       
    Bit_buffer_size -= 8;                 
  }                                           
} 

/*!
   Write to f n bits taken from v (possibly n > 24) */
void fbit_write(FILE *f, int n, uint32 v)
{  
  assert(f!=NULL);
  assert(n <= 32);
  if (n > 24){
    fbit_write24(f, n-24, (v>>24) & 0xffL);
    fbit_write24(f, 24, v & 0xffffffL);
  } else {
    fbit_write24(f,n,v);
  }
}

/*!
  flush the content of Bit_buffer, padding the unused bits with 0's */
void fbit_flush(FILE *f)
{
  assert(f!=NULL);
  if(Bit_buffer_size!=0) {
    assert(Bit_buffer_size<8);
    fbit_write24(f, 8 - (Bit_buffer_size%8), 0);  // pad with zero !
  }
}



// ###### bit input routines

/*! 
    Read n (<=24) bits from file f (and Bit_buffer) */
static uint32 fbit_read24(FILE *f, int n)
{
  uint32 t,u;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);

  /* --- read groups of 8 bits until size>= n --- */
  while(Bit_buffer_size<n) {
    t = (uint32) getc(f);
    if(t==(uint32) EOF) {
      perror("Unable to read more bits (fbit_read)");
      exit(1);
    }
    Bit_buffer |= (t << (24-Bit_buffer_size));
    Bit_buffer_size += 8;
  }
  /* ---- write n top bits in u ---- */
  u = Bit_buffer >> (32-n);
  /* ---- update buffer ---- */
  Bit_buffer <<= n;
  Bit_buffer_size -= n;
  return u;
}


/*!
   Read n bits from Bit_buffer */

uint32 fbit_read(FILE *f, int n)
{  
  uint32 u = 0;

  assert(n <= 32);
  if (n > 24){
    u =  fbit_read24(f, n-24)<<24;
    u |= fbit_read24(f, 24);
    return u;
  } else 
    return fbit_read24(f,n);
}

/*!
  Read an uint32 from f and the Bit_buffer */
uint32 uint32_read(FILE *f)
{
  uint32 u;

  u =  fbit_read24(f,8)<<24;
  u |= fbit_read24(f,8)<<16;
  u |= fbit_read24(f,8)<<8;
  u |= fbit_read24(f,8);
  return u;
}





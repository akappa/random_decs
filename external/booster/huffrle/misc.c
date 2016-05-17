/*! \file
    Miscellaneaous routines 
 */
#include "huffrle.h"



// ###### error messages

void _huff_fatal_error(char *s)
{
  fprintf(stderr,"%s\n",s);
  exit(1);
}

void _huff_out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f);
  exit(1);
}


// ###### bit manipulation

/*!
  compute # bits to represent u 
  \param u: a 32 bit unsigned integer */
int _huff_bin_digits(uint32 u)    
{
  int i = 1;
  uint32 r = 1;
  
  while((i<=32) && (r<u)){
    r=2*r+1;
    i = i+1;
  }
    
  assert(i<=32);
  return i;
}

/*!
   Convert the uint32 variable n in binary writing the binary
   digits to buffer[] and returning the number of digits written
   \param n: the 32-bit uint to be converted.
   \param buffer[]: an already allocated buffer (size 32) 
   \return number of bits written to buffer[]
*/
int _huff_uint32_to_binstring(uint32 n, uint8 *buffer)
{
  int k,written;
  uint32 mask;

  assert(n>0);
  k=31;
  mask=1<<k;
  // ----- now we search where the digits begins ------
  while((mask&n)==0) {
    k--;             // we maintain the invariant mask=1<<k
    mask >>= 1;      // shift right mask
    if(k<0)
      _huff_fatal_error("Something is really wrong! (uint32_to_binstring)");
  }
  // ----- now we write the digits --------------------
  written=k+1;        // number of digits written in buffer[]
  while (k>=0) {      // write digits
    buffer[k]=n&mask?1:0;
    mask >>=1;
    k--;
  }
  return written;
}

// ----- this function prints any char in a readable form
void _huff_fpretty_putchar(FILE *f, int c)
{

  if(c>=32 && c<127)      // printable char
    fprintf(f,"    %c", c);
  else if(c=='\n')
    fprintf(f,"   \\n");        // \n
  else if(c=='\t')
    fprintf(f,"   \\t");        // \t
  else
    fprintf(f," %04x", c);      // print hex code
}

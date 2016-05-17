/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   if_defs.c
   Ver 1.0 (6 March 2006) Giovanni Manzini (manzini@mfn.unipmn.it)

   Inversion frequecies encoding/decoding for use in the compression booster
   Comment: the experimental results are not encouraging:
     1) poor if used without boosting (in accordance with theory)
     2) worse than wavelet+rle if used with the booster. I have tried 
        also fibonacci and ternary codes to no avail.
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "booster.h"


// ====== external global variables ===============
extern FILE *_boost_Infile;            // input file for decompression
extern FILE *_boost_Outfile;           // output file for compression
extern int _boost_Verbose;             // verbosity level
extern int _boost_Infile_size;         // input file size


// compute the # of bits in the binary representation 
// of n, that is 1+\floor(log_2 n)
static inline int bin_digits(int n) {
  int b = 0;

  assert(n>0);
  while (n) { b++; n >>= 1; }
    return b;
}
// compute the len of the gamma code for n
static int gamma_len(int n)
{
  int len;

  assert(n>0);
  len = 2*bin_digits(n)-1;
  assert(len>0);
  return len;
}
static int ternary_len(int n)
{
  int i=0;

  assert(n>0);
  while(n>0) {
    i++;
    n = n/3;
  }
  return 1+2*i;
}
static int fibonacci_len(int n)
{
  int i=2,f1,f2,f3;

  assert(n>0);
  f1=1;f2=2;
  while(n>=f2) { 
    i++;
    f3=f1+f2;
    f1=f2; f2=f3;
  }
  return i;
}



static int fwrite_bottom_range(FILE *f, int32 t, int32 a, int32 b);


int tot_len, tot_bitmap, tot_occ, tot_gap;

/* *********************************************************************
   compute the full cost of the wavelet tree for t[0,n), where
   each t[i] is in the range [0,asize)
   the cost includes the encoding of n, the bitmap to mark which symbols
   do appear in t[], and finally the wavelet tree itself
   ********************************************************************* */ 
int if_full_cost(symbol *t, int n, int asize)
{
  int i,j,remaining,tot,written,gap,mf,maxocc,cost=0;
  int new_asize, occ[ALPHA_SIZE]={0}, alpha_map[ALPHA_SIZE];

  // count the # of occ of each symbol 
  for(i=0;i<n;i++) 
    occ[t[i]]++;                      

  // virtually remap symbols in [0,new_asize)
  // and compute the most frequent symbol
  mf=asize;                            // illegal value
  new_asize=0;                         // # chars seen so far
  for(maxocc=0,i=0;i<asize;i++) {
    if(occ[i] != 0) {                  // symbol i is used?
      alpha_map[i]=new_asize++;        // remap symbol i
      if(occ[i]>maxocc) {
	maxocc=occ[i];
	mf = i;
      }
    }
  }
  assert(new_asize<=asize);
  assert(mf<asize);

  // output bitmap of symbols present in t[0..n-1]
  cost += asize;
  tot_bitmap += asize;

  // output segment length using gamma coding
  assert(n>0);
  assert(n<_boost_Infile_size+2);
  tot = gamma_len(n);
  cost += tot;
  tot_len += tot;

  // if there is only one symbol we are done
  if(new_asize==1)
    return cost;

  // encode #occ of symbols in segment 
  remaining = n; tot=0;
  for(written=i=0;written<new_asize-1;i++) {
    if(occ[i]==0) continue;
    assert(occ[i]>0);
    assert(occ[i]<=remaining);
    tot += fwrite_bottom_range(NULL,occ[i],1,remaining);
    remaining -= occ[i];
    written++;
  } 
  tot_occ += tot; 
  cost += tot;

  // now encode gaps
  for(tot=i=0;i<asize;i++) {
    if(occ[i]==0 || i==mf) 
       continue;
    written=gap=0;
    for(j=0;j<n;j++) {
      if(t[j]==i) {
	tot += gamma_len(gap+1);
	gap=0;
	written++;
      }
      else if(t[j]>i || t[j]==mf)
	gap++;
    }
    assert(written==occ[i]);
  }
  tot_gap += tot;
  cost += tot;
  return cost;
}




int if_encode(symbol *t, int n, int asize)
{
  return if_full_cost(t, n, asize);
}

int if_estart(void)
{  
  tot_gap=tot_bitmap=tot_len=tot_occ=0;
  return 0;
}

int if_estop(void)
{
  if(_boost_Verbose>0) {
    fprintf(stderr,"Detailed costs for the optimal partition:\n    ");
    fprintf(stderr,"lengths: %d, bitmaps: %d, occs: %d gaps: %d\n",
	    (tot_len+7)/8, (tot_bitmap+7)/8,(tot_occ+7)/8,(tot_gap+7)/8);
  }
  return 0;
}




// initialize the struct charaterizing this algorithm
void init_if_gamma_compressor(base_compressor *b)
{
  b->name = "if coding (gamma)";
  b->bound = NULL;
  b->cost = if_full_cost;
  b->merge = NULL;
  b->encode = if_encode;
  b->decode = NULL;
  b->enc_start = if_estart;
  b->enc_stop = if_estop; 
  b->dec_start = NULL;
  b->dec_stop = NULL; 
}



#if 0
static void fbit_write(FILE *f, int n, uint32 v)
{
  return;
}
#endif



/*!
   write t (with 0<=t<n) to file f using log-skewed encoding.
   This encoding never uses more than ceil(log n) bits, but it can 
   use less bits when t is small 
   \param f: output file
   \param t: value to be written in f
   \param n: max. value+1 allowed for t
   \return number of bits used to write t;
*/
static int fwrite_logskewed(FILE *f, uint32 t, uint32 n)
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
        if(f!=NULL) /*fbit_write(f,k,t)*/;
	bits_written += k;
      } 
      break;
    }
    else {  
      delta = n-u;
      assert(delta>0);
      if(t>=delta) {           // encode t using k+1 bits
        if(f!=NULL) {
	  /*fbit_write(f,1,1);
	    fbit_write(f,k,t-delta)*/;
	}
	bits_written += k+1;	
        break;
      }
      else {
        if(f!=NULL) 
	  /*fbit_write(f,1,0)*/;   // write 0
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
#if 0
static int fwrite_top_range(FILE *f, int32 t, int32 a, int32 b)
{

  assert(t>=a && t < b);
  return fwrite_logskewed(f, b-1-t, b-a);
}
#endif

/*!
   write t (with a<=t<b) to file f using at most ceil(log(b-a)) bits
   and possibly using less bits when t-a is small 
   \param f: output file
   \param t: value to be written in f
   \param a: min value allowed for t
   \param b: max value+1 allowed for t
   \return number of bits used to write t;
*/
static int fwrite_bottom_range(FILE *f, int32 t, int32 a, int32 b)
{

  assert(t>=a && t < b);
  return fwrite_logskewed(f, t-a, b-a);
}


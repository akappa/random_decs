/* ==============================================================
   Cost functions (for the moment) for 2-pass arithmetic coding 
   the idea is that first we encode the # occ of each symbol
   then apply arithmetic coding decrementig the count. In this way
   ac coding should get a compression very close to H_0 (so I did
   not implement ac for the moment). I have also inserted a separate
   accounting for the bits used for: 1) segment lengths, 2) encoding
   of the number of occurrences, 3) the code itslef for ac 
   (only estimated).

   The results are not encouraging. For bible:

   Command line: bst -a4 -vv large/bible
   Decompression not supported!
   File size: 4047392 bytes
   ...
   Alphabet size: 63
   No partitioning: 2197256 bytes
   Boosting: *** 2pass Arith. Coding (bound only) ***
   ...
   Detailed costs for the optimal partition:
      lengths: 14932, occ's: 146473, ac code: 909595
   Output file size: 1071041 bytes
      Number of segments: 7557
   Actual compressed size: 1071000 bytes
   ...    

   The problem is that plain 1-pass ac coding (with def. parameters)
   uses 789890 bytes overall. My impression is that there is always a
   local omogeneity that 1pass ac is expoiting and 2pass is missing. 
   ============================================================== */ 

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
static inline int gamma_len(int n)
{
  int len;

  assert(n>0);
  len = 2*bin_digits(n)-1;
  assert(len>0);
  return len;
}
static int fwrite_bottom_range(FILE *f, int32 t, int32 a, int32 b);



// bound the oputput size using cost of encoding # of occ + entropy 
double arit3_bound(stats *s, int alpha_size)
{
  int i, remaining;
  double tot=0;

  assert(s->count_tot>0);
  // encode length of segment
  tot += fwrite_bottom_range(NULL,s->count_tot,1,_boost_Infile_size+2);
  // encod #occ of symbols in segment 
  remaining = s->count_tot;
  for(i=0;i<alpha_size && remaining>0;i++) {
    tot += fwrite_bottom_range(NULL,s->count[i],0,remaining+1);
    remaining -= s->count[i];
  } 
  assert(remaining==0);

  for(i=0;i<alpha_size;i++) {
    // tot += gamma_len(1+s->count[i])/2;
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i])/log(2.0);
    }
  }
  if(_boost_Verbose>1 && (s->count_tot == _boost_Infile_size+1))
    fprintf(stderr,"Estimate for unpartitioned BWT: %f bytes\n",tot/8);
  return tot;
}


int tot_len, tot_occ, tot_code;


int arit3_cost(symbol *t, int n, int asize)
{
  int tot, remaining, cost=0, i, occ[256]={0};

  // count number of occ of each char
  for(i=0;i<n;i++)
    occ[t[i]]++;
  
  // encoding of the string length
  assert(n>0);
  assert(n<_boost_Infile_size+2);
  tot = fwrite_bottom_range(NULL,n,1,_boost_Infile_size+2);
  tot_len += tot;
  cost += tot;

  // encode #occ of symbols in segment 
  remaining = n; tot=0;
  for(i=0;i<asize && remaining>0;i++) {
    assert(occ[i]>=0);
    assert(occ[i]<=remaining);
    tot += fwrite_bottom_range(NULL,occ[i],0,remaining+1);
    remaining -= occ[i];
  } 
  assert(remaining==0);
  tot_occ += tot;
  cost += tot;

  tot=0;
  for(i=0;i<asize;i++) {

    if(occ[i]!=0) {
      assert(occ[i]>0);
      tot += 
       (int) ceil(occ[i]*log(((double) n)/occ[i])/log(2.0));
    }
  }
  tot_code += tot;
  cost += tot;
  if(_boost_Verbose>1 && (n == _boost_Infile_size+1))
    fprintf(stderr,"Estimate for unpartitioned BWT: %d bytes\n",(cost+7)/8);
  return cost;
}

int arit3_encode(symbol *t, int n, int asize)
{
  return arit3_cost(t, n, asize);
}

int arit3_estart(void)
{  
  tot_code=tot_occ=tot_len=0;
  return 0;
}

int arit3_estop(void)
{
  if(_boost_Verbose>0) {
    fprintf(stderr,"Detailed costs for the optimal partition:\n    ");
    fprintf(stderr,"lengths: %d, occ's: %d, ac code: %d\n",
	    (tot_len+7)/8, (tot_occ+7)/8,(tot_code+7)/8);
  }
  return 0;
}




// initialize the struct charaterizing this algorithm
void init_2pass_ac_compressor(base_compressor *b)
{
  b->name = "2pass Arith. Coding (bound only)";
  b->bound = arit3_bound;
  b->cost = arit3_cost;
  b->merge = NULL;
  b->encode = arit3_encode;
  b->decode = NULL;
  b->enc_start = arit3_estart;
  b->enc_stop = arit3_estop; 
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

       

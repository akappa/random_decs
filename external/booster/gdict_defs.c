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

inline int bits (int n) {
  int b = 0;
  while (n) { b++; n >>= 1; }
    return b;
}


// bound the oputput size using gonzalo's strange cost 
double gonzalo_cost(stats *s, int alpha_size)
{
  int i, numchar=0, bits_n;
  double tot=0, cost;

  assert(s->count_tot>0);
  bits_n = bits(_boost_Infile_size);
  for(i=0;i<alpha_size;i++)
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i])/log(2.0);
      numchar++;
    }
  assert(tot>=0);
  if(numchar == 1) 
    cost = alpha_size * bits_n + 8 + 32;
  else
    cost = alpha_size * bits_n + 32 + 64 + (numchar + 2) * 8 + 
           numchar*bits(numchar) + numchar*6*32 + 1.3 * tot;

  if(_boost_Verbose>1 && (s->count_tot == _boost_Infile_size+1))
    fprintf(stderr,"Estimate for unpartitioned BWT: %f bytes\n",cost/8);
  return cost;
}


// initialize the struct charaterizing this algorithm
void init_gdict_compressor(base_compressor *b)
{
  b->name = "Gonzalo's Dict.";
  b->bound = gonzalo_cost;
  b->cost = NULL;
  b->merge = NULL;
  b->encode = NULL;
  b->decode = NULL;
  b->enc_start = NULL;
  b->enc_stop = NULL; 
  b->dec_start = NULL;
  b->dec_stop = NULL; 
}



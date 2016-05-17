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
extern double _boost_CmdLinePar1;      // optional command line parameter #1
extern int _boost_CmdLinePar1_flag;    // !=0 if parameter #1 was defined
extern double _boost_CmdLinePar2;      // optional command line parameter #2
extern int _boost_CmdLinePar2_flag;    // !=0 if parameter #2 was defined
extern double _boost_CmdLinePar3;      // otional command line parameter #3
extern int _boost_CmdLinePar3_flag;    // !=0 if parameter #3 was defined


static void uint32_write(FILE *f, uint32 u);
static uint32 uint32_read(FILE *f);


// bound the oputput size using the modified 0-order entropy 
double h0star(stats *s, int alpha_size)
{
  int i, numchar=0;
  double tot=0, mu;

  assert(s->count_tot>0);
  mu = _boost_CmdLinePar1_flag ? _boost_CmdLinePar1 : 1.0;
  for(i=0;i<alpha_size;i++)
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i]);
      numchar++;
    }
  assert(tot>=0);
  if(numchar == 1) 
    return mu*log(alpha_size) + (1+log(s->count_tot+1));
  return mu * numchar * log(alpha_size) + tot;
}

// ---- dummy compression algorithm
int dummy_compr_aux(symbol *t, int n, int alpha_size)
{
  int i;

  assert(alpha_size>0);
  assert(n>0);
  uint32_write(_boost_Outfile,n);
  for(i=0;i<n;i++)
    putc(t[i],_boost_Outfile);
  return 8*(n+4);
}



int dummy_decompr_aux(symbol *t, int n, int alpha_size)
{
  int i,m;

  assert(alpha_size>0);
  m = uint32_read(_boost_Infile);
  assert(m<=n);
  for(i=0;i<m;i++)
    t[i] = getc(_boost_Infile);
  return m;
}


// initialize the struct charaterizing this algorithm
void init_Dummy_compressor(base_compressor *b)
{
  b->name = "Dummy";
  b->bound = h0star;
  b->cost = NULL;
  b->merge = NULL;
  b->encode = dummy_compr_aux;
  b->decode = dummy_decompr_aux;
  b->enc_start = NULL;
  b->enc_stop = NULL; 
  b->dec_start = NULL;
  b->dec_stop = NULL; 
}


// === auxiliary functions ==========

/* write an uint32 to f */
static void uint32_write(FILE *f, uint32 u)
{

  assert(f!=NULL);
  putc((u>>24) & 0xffL,f);
  putc((u>>16) & 0xffL,f);
  putc((u>> 8) & 0xffL,f);
  putc(      u & 0xffL,f);

}

/* Read an uint32 from f */
static uint32 uint32_read(FILE *f)
{
  uint32 u;

  u =  ((uint32) getc(f)) <<24;
  u |= ((uint32) getc(f)) <<16;
  u |= ((uint32) getc(f))  <<8;
  u |= ((uint32) getc(f));
  return u;
}

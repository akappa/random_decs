#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "booster.h"
#include "huffrle_aux.h"


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


// bound the oputput size using the modified 0-order entropy 
double huff_h0star(stats *s, int alpha_size)
{
  int i, numchar=0;
  double mu,tot=0;

  assert(s->count_tot>0);
  mu = _boost_CmdLinePar1_flag ? _boost_CmdLinePar1 : 12.0;
  for(i=0;i<alpha_size;i++)
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i]);
      numchar++;
    }
  assert(tot>=0);
  if(numchar == 1) 
    return mu * log(alpha_size) + (1+log(s->count_tot+1));
  return mu * numchar * log(alpha_size) +  tot;
}


int huffrle_compr_aux(symbol *t, int n, int alpha_size)
{
  plain_text p;
  int bits; 

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose; 

  bits = huffrle_compr(&p, RLE_ALL, _boost_Outfile);
  return bits;
}

int huffrle_decompr_aux(symbol *t, int n, int alpha_size)
{
  plain_text p;

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose;

  return huffrle_decompr(&p, RLE_ALL, _boost_Infile);
}


int huffrle_cost(symbol *t, int n, int alpha_size)
{
  plain_text p;

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose;

  return huffrle_compr(&p, RLE_ALL, NULL);
}

static int huffrle_init_buffer(void)
{
  init_bit_buffer();
  return 0;
}

static int huffrle_flush_buffer(void)
{
  fbit_flush(_boost_Outfile);
  return 0; 
}



// initialize the struct charaterizing this algorithm
void init_Huffrle_compressor(base_compressor *b)
{
  b->name = "Huffrle (RLE + Huffman)";
  b->bound = huff_h0star;
  b->cost = huffrle_cost;
  b->merge = NULL;
  b->encode = huffrle_compr_aux;
  b->decode = huffrle_decompr_aux;
  b->enc_start = huffrle_init_buffer;
  b->enc_stop = huffrle_flush_buffer; 
  b->dec_start = huffrle_init_buffer;
  b->dec_stop = NULL; 
}


// ======= rle0 + huffman


int huffrle0_cost(symbol *t, int n, int alpha_size)
{
  plain_text p;

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose;

  return huffrle_compr(&p, 0, NULL);
}


int huffrle0_compr_aux(symbol *t, int n, int alpha_size)
{
  plain_text p;
  int bits; 

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose; 

  bits = huffrle_compr(&p, 0, _boost_Outfile);
  return bits;
}


int huffrle0_decompr_aux(symbol *t, int n, int alpha_size)
{
  plain_text p;

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose;

  return huffrle_decompr(&p, 0, _boost_Infile);
}


// initialize the struct charaterizing this algorithm
void init_Huffrle0_compressor(base_compressor *b)
{
  b->name = "Huffrle0 (RLE0 + Huffman)";
  b->bound = huff_h0star;
  b->cost = huffrle0_cost;
  b->merge = NULL;
  b->encode = huffrle0_compr_aux;
  b->decode = huffrle0_decompr_aux;
  b->enc_start = huffrle_init_buffer;
  b->enc_stop = huffrle_flush_buffer; 
  b->dec_start = huffrle_init_buffer;
  b->dec_stop = NULL; 
}


// ======= static huffman without rle


int huffnorle_cost(symbol *t, int n, int alpha_size)
{
  plain_text p;

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose;

  return huffrle_compr(&p, NO_RLE, NULL);
}


int huffnorle_compr_aux(symbol *t, int n, int alpha_size)
{
  plain_text p;
  int bits; 

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose; 

  bits = huffrle_compr(&p, NO_RLE, _boost_Outfile);
  return bits;
}


int huffnorle_decompr_aux(symbol *t, int n, int alpha_size)
{
  plain_text p;

  p.text = t;
  p.tsize = n;
  p.asize = alpha_size;
  Huffrle_Verbose = _boost_Verbose;

  return huffrle_decompr(&p, NO_RLE, _boost_Infile);
}


// initialize the struct charaterizing this algorithm
void init_Huffnorle_compressor(base_compressor *b)
{
  b->name = "Huffnorle (Huffman without RLE)";
  b->bound = huff_h0star;
  b->cost = huffnorle_cost;
  b->merge = NULL;
  b->encode = huffnorle_compr_aux;
  b->decode = huffnorle_decompr_aux;
  b->enc_start = huffrle_init_buffer;
  b->enc_stop = huffrle_flush_buffer; 
  b->dec_start = huffrle_init_buffer;
  b->dec_stop = NULL; 
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   wv_defs.c
   Ver 1.0 (6 March 2006) Giovanni Manzini (manzini@mfn.unipmn.it)

   Wavelet tree encoding/decoding for use in the compression booster
   The internal nodes of the wavelet tree are encoded using rle + a 
   generic prefix encoding of the integers. At the momento only 
   gamma code is implemented but others can be easily added since the
   prefix encoder is a parameter of the wavewlet encoding function.  
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include "booster.h"
#include "bit_buffer.h"

// ====== external global variables ===============
extern FILE *_boost_Infile;            // input file for decompression
extern FILE *_boost_Outfile;           // output file for compression
extern int _boost_Verbose;             // verbosity level
extern int _boost_Infile_size;         // input file size


// ======= local global variables =================
static int tot_len;      // number of bits used for encoding segments' lengths 
static int tot_bitmap;   // number of bits used for encoding segments' bitmaps 
static int tot_inodes;   // number of bits used for encoding wt internal nodes
static int Verbose;      // "local" verbosity level

// ----- wavelet tree node encoded via rle ------------ 
typedef struct {
  int last;            // bit of the pending run 
  int run;             // length of pending run
  _bb_bit_buffer bb;   // bit buffer containing the encoding of previous runs 
} wv_rle_node;



// ======= prototypes =============================
static void fatal_error(char *s);
static int wv_rle_build(symbol *t, int n, int asize, FILE *f,
			uint64_t (* prefix_enc)(uint32_t), 
                        int (* prefix_len)(uint32_t));
static int wv_rle_inodes_enc(symbol *t,int n,int *map,int asize, FILE *f,
			     uint64_t (* prefix_enc)(uint32_t),
			     int (* prefix_len)(uint32_t));
static int wv_rle_decode(symbol *t, int limit, int asize, FILE *f,
			 uint32_t (* pfxcode_read)(FILE *));
static void wv_rle_inodes_dec(symbol *t, int n, int *alpha_map, int asize,
			      FILE *f, uint32_t (*pfxcode_read)(FILE *));


// ==== bit I/O
static void init_bit_buffer(void);
static void fbit_flush(FILE *f);
static void fbit_write(FILE *f, int n, uint32_t v);
static void fbit_write64(FILE *f, int n, uint64_t v);
static uint32_t fbit_read(FILE *f, int n);


// =================================================================
// gamma encoding functions. If you want to use a different
// prefix encoding (eg delta encoding) you must supply similar 
// a pair of similar functions, with exactly the same data types.
// See wv_rle_build_full() for their usage. The only limitation is 
// that a codeword can be at most 64 bit long
// =================================================================

// compute the len of the gamma code for n
static int gammacode_len(uint32_t n)
{
  int b=0;

  assert(n>0);
  while (n) { b++; n >>= 1; }
  assert(b>0);
  return 2*b-1;
}

// compute the codeword for n using gamma coding
// (leading zeros are not a problem since this function is always  
// used combined with gammacode_len 
static uint64_t gammacode_word(uint32_t n)
{
  // printf("--> %10d  ",n); // !!!!!!!!!!!!
  return (uint64_t) n;   
}

// get a gammacode from a file. To get bits from the file 
// this function can use only the procedure
//   uint32_t fbit_read(FILE *f, n)   
// to read up to 32 bits from f
static uint32_t gammacode_read(FILE *f)
{
  int b; uint32_t v;

  // read a sequence of b zeroes
  for(b=0; ;b++) 
    if(fbit_read(f,1)!=0) break;
  assert(b<32);
  // get the binary digits of v
  if(b==0) v = 1;
  else v = (1<<b) | fbit_read(f,b);
  // fprintf(stderr,"<-- %10d  ",v); // !!!!!!!!!
  return v;
}

  


// ======= coding/decoding and cost estimate functions =====

// evaluate the cost of wv+rle+gamma on t[0,n)
static int wv_rle_gamma_cost(symbol *t, int n, int asize)
{
  return wv_rle_build(t,n,asize,NULL,gammacode_word,gammacode_len);
}

// actual wv+rle+gamma compression on t[0,n)
static int wv_rle_gamma_encode(symbol *t, int n, int asize)
{
  return wv_rle_build(t,n,asize,_boost_Outfile,gammacode_word,gammacode_len);
}


// decode a segment of length at most n
// return the number of symbols decoded
static int wv_rle_gamma_decode(symbol *t, int limit, int asize)
{
  return wv_rle_decode(t,limit,asize,_boost_Infile,gammacode_read);
}

static int wv_rle_gamma_enc_start(void)
{  
  Verbose = _boost_Verbose;
  init_bit_buffer();
  tot_inodes=tot_bitmap=tot_len=0;
  return 0;
}

static int wv_rle_gamma_enc_stop(void)
{
  fbit_flush(_boost_Outfile);
  if(_boost_Verbose>1) {
    fprintf(stderr,"Detailed costs (in bits) for the optimal partition:\n");
    fprintf(stderr,"   lengths: %d, bitmaps: %d, inodes: %d\n",
	    (tot_len), (tot_bitmap),(tot_inodes));
  }
  return 0;
}

// initialize the struct characterizing this algorithm
void init_wv_rle_gamma_compr(base_compressor *b)
{
  b->name = "Wavelet tree + RLE + gamma";
  b->bound = NULL;
  b->cost = wv_rle_gamma_cost;
  b->merge = NULL;
  b->encode = wv_rle_gamma_encode;
  b->decode = wv_rle_gamma_decode;
  b->enc_start = wv_rle_gamma_enc_start;
  b->enc_stop = wv_rle_gamma_enc_stop; 
  b->dec_start = wv_rle_gamma_enc_start;
  b->dec_stop = NULL; 
}


/* *********************************************************************
   decode a wavelet tree read from file f and write the decoded string
   to t[]. The decode string must be of length at least limit and is
   over the alphabet [0,asize).
   ********************************************************************* */
static int wv_rle_decode(symbol *t, int limit, int asize, FILE *f,
			uint32_t (* pfxcode_read)(FILE *))
{
   int new_asize,i,n,alpha_remap[ALPHA_SIZE];

  // ----- read the symbol bitmap  
  // set new_asize as the number of distinct symbols
  // remap from [0..new_asize) to  [0..asize)
  new_asize=0;                         // # chars seen so far
  for(i=0;i<asize;i++) {               // scan bitmap on file 
    if(fbit_read(f,1) != 0) {          // symbol i is used?
      alpha_remap[new_asize++] = i;    // image of symbol i
    }
  }
  assert(new_asize<=asize);

  // ---- get segment length using gamma coding
  n = pfxcode_read(f);
  if(Verbose>1) 
    fprintf(stderr,"Segment size: %d (limit: %d)\n",n,limit);
  assert(n>0);
  if(n>limit)
    fatal_error("Plain text buffer overflow (wv_rle_decode)\n");


  // ---- decode wavelet tree internal nodes

  // first the easy case in which there are no internal nodes
  if(new_asize==1) {       
    for(i=0;i<n;i++)
      t[i] = alpha_remap[0];
    return n;
  }
  
  // now the hard case in which there are internal nodes
  wv_rle_inodes_dec(t,n,alpha_remap,new_asize,f,pfxcode_read);
  return n;    
}


static void wv_rle_inodes_dec(symbol *t, int n, int *alpha_map, int asize,
			      FILE *f, uint32_t (*pfxcode_read)(FILE *))
{
  int node_size[2*ALPHA_SIZE]={0}; // size of string associated to node 
  int i,u,len,tot_len,c,nodes_to_go,cur;

  // -------- decode root node
  u=1;                     // node number (root)
  i=0;                     // symbols written so far
  c=0;                     // first symbol is zero by construction
  len=pfxcode_read(f)-1;   // remove the initial 0
  assert(len<n);
  while(1) {            
    node_size[2*u+c] += len;  // update size of child
    for(;len>0;len--,i++)     // write len c's
      t[i]= (symbol) c;
    if(i==n) break;           // end of root node
    c = 1-c;                  // change symbol
    len = pfxcode_read(f);    // lenght of next run
    if(len+i>n) {
      fprintf(stderr,"n=%d len+i=%d\n",n,len+i);
      fatal_error("Invalid wavelet tree root node (wv_rle_inodes_dec)\n");
    }
  }    
  assert(i==n);
  assert(node_size[2]+node_size[3]==n);


  // -------- decode internal non-root nodes
  cur=0;                        // next symbol to be processed
  nodes_to_go=0;                // # of nodes till the end of the level
  for(u=2;u<asize;u++) {
    // {for(i=0;i<n;i++) fprintf(stderr,"%4d",t[i]);}// !!!!!!!!!!!

    if(nodes_to_go==0) {        // if we reached the end of a tree level
      cur=0; nodes_to_go=u;     // reset the current symbol to 0 
    }
    i=0;
    c=0;                        // first symbol is zero by construction
    tot_len=0;
    len=pfxcode_read(f)-1;      // remove initial 0
    while(1) {
      node_size[2*u+c] += len;  // update size of child
      tot_len += len;           // update # of symbols seen
      for(;len>0;i++) {         // update len occurrences of cur
	if(t[i]==cur) {
	  if(c) t[i] += 1;      // update cur
	  len--;
	}
	else if(t[i]>cur)
	  t[i] += 1;            // update symbols greater than cur
	assert(i<n);
      }
      // check if we are done with this node
      if(tot_len==node_size[u]) {
	for( ;i<n;i++)
	  if(t[i]>cur) t[i]++; // update remaining symbols
	break;
      }
      c = 1-c;                 // change symbol
      len=pfxcode_read(f);     // get next run
      if(len+tot_len>node_size[u]) {
	fprintf(stderr,"Error processing node %d: ",u);
	fatal_error("Invalid wavelet tree node (wv_rle_inodes_dec)\n");
      }
    }  
    assert(node_size[2*u]+node_size[2*u+1]==node_size[u]);	
    assert(node_size[2*u]!=0);
    assert(node_size[2*u+1]!=0);
     
    cur+=2;           // update cur (symbol affected by next internal node )
    nodes_to_go--;    // number of nodes missing till the end of the level
  }
  
  // ----------- remap symbols with the "nodes_to_go" correction
  for(i=0;i<n;i++) {
    assert(t[i]<asize);
    t[i] = (symbol) alpha_map[(t[i]+nodes_to_go)%asize];
  }
  return;
}




/* *********************************************************************
   Build the full rle wavelet tree for t[0,n), where each t[i] is in the 
   range [0,asize). The wavelet tree includes the encoding of n, the bitmap 
   to mark which symbols do appear in t[], and finally the internal nodes.
   If f!=NULL the wavelet tree is written to f, otherwise the function
   only returns the size (in bits) of the wavelet tree.

   Internal nodes are encoded using rle and the prefix encoding of the 
   integers indicated by the functions prefix_enc and prefix_len.
   prefix_enc(n) returns the encoding of n, while prefix_len(n)
   returns the length of the encoding.
   ********************************************************************* */ 
static int wv_rle_build(symbol *t, int n, int asize, FILE *f,
			uint64_t (* pfxcode_word)(uint32_t), 
			int (* pfxcode_len)(uint32_t))
{
  int new_asize, i,tot,cost=0;
  int occ[ALPHA_SIZE]={0}, alpha_map[ALPHA_SIZE];

  // set new_asize as the number of distinct symbols
  // remap the old symbols in [0..new_asize-1]
  // (also count the # of occ of each symbol (not used for the moment) 
  for(i=0;i<n;i++) 
    occ[t[i]]++;                       // count occ
  new_asize=0;                         // # chars seen so far
  for(i=0;i<asize;i++) {               // scan alphabet 
    if(occ[i] != 0) {                  // symbol i is used?
      alpha_map[i]=new_asize++;        // remap symbol i
    }
  }
  assert(new_asize<=asize);

  // output bitmap of symbols present in t[0..n-1]
  if(f!=NULL) {
    for(i=0;i<asize;i++) 
      if(occ[i]!=0) fbit_write(f,1,1);
      else          fbit_write(f,1,0);	
  }
  cost += asize;
  tot_bitmap += asize;

  // output segment length using gamma coding
  assert(n>0);
  assert(n<_boost_Infile_size+2);
  tot = pfxcode_len(n);
  if(f!=NULL) 
    fbit_write64(f,tot,pfxcode_word(n));
  cost += tot;
  tot_len += tot;

  // compute cost of wavelet tree internal nodes
  if(new_asize>1) {
    //for(i=0;i<n;i++) fprintf(stderr,"%4d",alpha_map[t[i]]);
    //fprintf(stderr,"\n"); // !!!!!!!!!!!!!!!!!
    tot = wv_rle_inodes_enc(t,n,alpha_map,new_asize,f,
			    pfxcode_word,pfxcode_len);
    cost += tot;
    tot_inodes += tot;
  }
  return cost;
}


/* ******************************************************************
   build the internal nodes of the wavelet tree for 
     map[t[0]], map[t[1]], ... , map[t[n-1]]
   each map[t[i]] is in the range [0,asize). If f!=NULL the internal 
   nodes are written to f, otherwise we only compute the cost of the 
   encoding. See wv_rle_build for details
   ****************************************************************** */
static int wv_rle_inodes_enc(symbol *t,int n,int *map,int asize, FILE *f,
			uint64_t (* pfxcode_word)(uint32_t), 
			int (* pfxcode_len)(uint32_t))
{
  wv_rle_node wvt[ALPHA_SIZE];   // internal nodes of wavelet tree
  int i,u,p,rc,c,b,len,cost=0;

  assert(asize<=ALPHA_SIZE);
  // init internal nodes [1..asize-1] 
  for(i=1;i<asize;i++) {
    wvt[i].last=0;       // start with a dummy 0
    wvt[i].run=1;
    if(f!=NULL)
      _bb_init_bit_buffer(&(wvt[i].bb));
  }
  for(i=0;i<n; ) {
    // get next run of equal symbols
    for(c=map[t[i]],rc=1; (++i<n) && c==map[t[i]]; rc++);
    assert(c>=0 && c<asize);
    assert(rc>0);
    assert(rc<n);
    u = c+asize;       // u is the leaf corresponding to c
    while(u>1) {       // add run to u's ancestors 
      p = u>>1;        // p is u's parent
      b = (u&1);       // b is the bit correspondif to u
      assert(p>0 && p<asize);
      if(wvt[p].last==b)
	wvt[p].run += rc;         // continue the previous run 
      else {
	assert(wvt[p].run>0);
	len = pfxcode_len(wvt[p].run);  // flush previous run
        if(f!=NULL) {                   // write bits to buffer
	  _bb_write_bits64(&(wvt[p].bb),len,pfxcode_word(wvt[p].run));
	}
	cost += len;                   // update total cost 
	wvt[p].last=b;                 // start new run
	wvt[p].run=rc;
      }
      u = p;          // go up one level
    } // end while
    // here i>=n or map[t[i]]!=c
    // if i<n map[t[i]] is the start of the next run
  } // end for
  for(i=1;i<asize;i++) {
    len = pfxcode_len(wvt[i].run);      // flush last run of each node
    if(f!=NULL) {                       // write last run and write buffer
      _bb_write_bits64(&(wvt[i].bb),len,pfxcode_word(wvt[i].run));
      _bb_bit_buffer2file(&(wvt[i].bb),f,fbit_write);
      _bb_free_bit_buffer(&(wvt[i].bb));
    }
    cost += len; 
  }
  return cost;
}

/* ==================================================================
   Low-level procedures to read and write bits using a Bit_buffer. 
   Unread/unwritten bits of Bit_buffer are the most significant ones.
   ================================================================== */
static uint32_t Bit_buffer;    /* 32 bit buffer */
static int  Bit_buffer_size;   /* # of unread/unwritten bits in Bit_buffer */
static uint64_t Bit_buffer_written; // total number of bits written

/*!
   Initializes the Bit_buffer to zero */
static void init_bit_buffer(void)
{
  Bit_buffer= (uint32_t) 0;
  Bit_buffer_size=0;
  Bit_buffer_written=0;
}


/* ************************************************************
   write to f n (<=32) bits taken from v
   Bit_buffer is flushed out until it contains <8 bits
   ************************************************************ */
static void fbit_write(FILE *f, int n, uint32_t v)
{
  assert(Bit_buffer_size<8);
  assert(n>0 && n<=32);

  // make sure bits n+1 ... 31 of v are zero
  if(n<32)
    v &= ((1<<n)-1);

  if(Bit_buffer_size+n<=32) {
    /* ------- add n bits to Bit_buffer -------- */
    Bit_buffer_size += n;                         // compute the correct shift
    Bit_buffer |= (v << (32 - Bit_buffer_size));  // compact to end of buffer
    /* ------- flush Bit_buffer as much as possible ----- */
    while (Bit_buffer_size>=8) {
      Bit_buffer_written += 8;
      if( putc((Bit_buffer>>24),f) == EOF) {
	perror("Unable to write bits (fbit_write/wv_defs.c)");
	exit(1);
      }
      Bit_buffer <<= 8;                       
      Bit_buffer_size -= 8;                 
    }
  }
  else { // n+Bit_buffer_size>32
    int delta=n+Bit_buffer_size-32;
    Bit_buffer |= (v>>delta); // now bit_buffer has 32 bits
    Bit_buffer_written += 32;
    if( putc((Bit_buffer>>24) & 0xFF,f) == EOF)
	perror("Unable to write bits (fbit_write/wv_defs.c)");
    if( putc((Bit_buffer>>16) & 0xFF,f) == EOF) 
	perror("Unable to write bits (fbit_write/wv_defs.c)");
    if( putc((Bit_buffer>>8) & 0xFF,f) == EOF) 
	perror("Unable to write bits (fbit_write/wv_defs.c)");
    if( putc((Bit_buffer) & 0xFF,f) == EOF) {
	perror("Unable to write bits (fbit_write/wv_defs.c)");
	exit(1);
    }
    // ---- insert the remaining delta bits in the buffer
    Bit_buffer = v <<(32-delta);
    Bit_buffer_size = delta;
  } 
}

static void fbit_write64(FILE *f, int n, uint64_t v)
{
  if(n>32) {
    fbit_write(f,n-32,(uint32_t) (v>>32));
    n = 32;
  }
  fbit_write(f,n,(uint32_t) v);
}


/*!
  flush the content of Bit_buffer, padding the unused bits with 0's */
static void fbit_flush(FILE *f)
{
  assert(f!=NULL);
  if(Bit_buffer_size!=0) {
    assert(Bit_buffer_size<8);
    fbit_write(f, 8 - (Bit_buffer_size%8), 0);  // pad with zero !
  }
  // fprintf(stderr,">>> Bits written: %lld (%lld bytes)\n",
  //	  Bit_buffer_written, Bit_buffer_written/8);
}



// ###### bit input routines

// read n (<= 32) bits from file f
static uint32_t fbit_read(FILE *f, int n)
{
  uint32_t w,delta,v;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=32);

  // easy case: we already have enough bits
  if(n<=Bit_buffer_size) {
    v = (Bit_buffer >> (32-n) );
    Bit_buffer <<= n;      // discard n bits from buffer
    Bit_buffer_size -= n;
    return v;
  }

  // if there are available bits place them
  if(Bit_buffer_size>0)
    v = (Bit_buffer >> (32-n));
  else 
    v=0;
  delta = n-Bit_buffer_size; // number of missing bits
  assert(delta>0);
  Bit_buffer_size=0; // no more available bits      
  Bit_buffer=0;      // now the bit buffer is in a coherent state

  // add blocks of 8 bits
  while(delta>=8) {
    w = (uint32_t) getc(f);
    if(w==(uint32_t) EOF) {
      perror("Unable to read more bits (fbit_read)");
      exit(1);
    }
    assert(w<256);
    delta -= 8;
    v |= (w<<delta);
  }
  assert(delta<8);

  if(delta>0) {
    w = (uint32_t) getc(f);
    if(w==(uint32_t) EOF) {
      perror("Unable to read more bits (fbit_read)");
      exit(1);
    }
    v |= (w>>(8-delta));
    Bit_buffer_size = 8-delta;
    Bit_buffer = w << (32-Bit_buffer_size);
  }
  return v;
}


#if 0
// read n (<=64) bits from f (NOT USED)
static uint64_t fbit_read64(FILE *f, int n)
{
  uint64_t v;

  assert(n>0 && n<=64);
  if(n>32) {
    v = fbit_read(f,n-32);
    v = v << 32;
    v |= fbit_read(f,32);
  }
  else
    v = fbit_read(f,n);

  return v;
}
#endif


static void fatal_error(char *s)
{
  fprintf(stderr,"%s",s); exit(1);
}

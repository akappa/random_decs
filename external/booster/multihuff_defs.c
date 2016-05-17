/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   multihuff_defs.c
   Ver 1.0 (June 15 2005)

   multiple table huffman coding/decoding for use by the 
   compression booster.
   This is essentially the multi-huffman code by Julian Seward's
   bzip2 tool with minimum adaptations. 

   G. Manzini
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "booster.h"
#include "bzlib_booster.h"


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


// ======= local global variables ======
static FILE *Outfile, *Infile;
static int Verbose;
static int MH_bits_total;

// ======= prototypes ==================
static int multihuff_compr_aux(uchar *in, int len, int asize,int rlechar);
static int multihuff_decompr_aux(uchar *dest, int len, int asize,int rlechar);
static int int_log2(int u);
static void fatal_error(char *s);
static void out_of_mem(char *s);
static void init_bit_buffer(void);
static void bit_write(int n, uint32 v);     // write n<=32 bits from v  
static void bit_flush(void);                // flush the content of buffer
static uint32 bit_read(int n);              // read n<=32 bits

// ========= macros ==========
#define single_bit_read_macro(A) (A=bit_read(1))
#define bit_read_macro(A,B) (A=bit_read(B))

/* -------- arrays used by multihuf ------------ */
static uchar huf_len[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];  // coding and decoding
static int huf_code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];   // coding
static int rfreq[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];      // coding
static int mtf_freq[BZ_MAX_ALPHA_SIZE];                // coding
static uchar huf_minLens[BZ_N_GROUPS];                // decoding
static int huf_limit[BZ_N_GROUPS][BZ_MAX_CODE_LEN];   // decoding
static int huf_base[BZ_N_GROUPS][BZ_MAX_CODE_LEN];    // decoding
static int huf_perm[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];  // decoding

static int mh_init_compr()
{
  Outfile = _boost_Outfile;
  Verbose = _boost_Verbose;
  init_bit_buffer();
  return 0;
}

static int mh_stop_compr()
{
  bit_flush();
  return 0;
}

static int mh_init_decompr()
{
  Infile = _boost_Infile;
  Verbose = _boost_Verbose;
  init_bit_buffer();
  return 0;
}

 
// bound the oputput size using the modified 0-order entropy 
static double mh_h0star(stats *s, int alpha_size)
{
  int i, numchar=0;
  double mu,tot=0;

  assert(s->count_tot>0);
  mu = _boost_CmdLinePar1_flag ? _boost_CmdLinePar1 : 30.0;
  for(i=0;i<alpha_size;i++)
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i]);
      numchar++;
    }
  assert(tot>=0);
  if(numchar == 1) 
    return mu*log(alpha_size) + (1+log(s->count_tot+1));
  return mu*numchar*log(alpha_size) + tot;
}


// ------ multi-table huffman compr/decompr with rle on all symbols ------- 
int multihuff_compr(uchar *in, int len, int alpha_size)
{
  return multihuff_compr_aux(in, len, alpha_size, 1);
}

int multihuff_decompr(uchar *dest, int limit, int alpha_size)
{
  return multihuff_decompr_aux(dest,limit, alpha_size, 1);
}

static int multihuff_cost(uchar *t, int n, int asize)
{
  int bits;
  Outfile = NULL;
  Verbose = _boost_Verbose;
  init_bit_buffer();
  bits =  multihuff_compr_aux(t,n,asize,1);
  // fprintf(stderr,"Input: %d bytes, Output %d bits\n",n,bits); // !!!!!!!!!!
  return bits;
}



// ------ multi-table huffman compr/decompr with rle on 0's ------ 
int multihuff0_compr(uchar *in, int len, int alpha_size)
{
  return multihuff_compr_aux(in, len, alpha_size,0);
}

int multihuff0_decompr(uchar *dest, int limit, int alpha_size)
{
  return multihuff_decompr_aux(dest,limit, alpha_size, 0);
}

static int multihuff0_cost(uchar *t, int n, int asize)
{
  int bits;
  Outfile = NULL;
  Verbose = _boost_Verbose;
  init_bit_buffer();
  bits =  multihuff_compr_aux(t,n,asize,0);
  // fprintf(stderr,"Input: %d bytes, Output %d bits\n",n,bits); // !!!!!!!!!!
  return bits;
}

// initialize the struct charaterizing this algorithm
void init_mh_compressor(base_compressor *b)
{
  b->name = "RLE + Multiple-Table Huffman Coding";
  b->bound = mh_h0star;
  b->cost = multihuff_cost;
  b->merge = NULL;
  b->encode = multihuff_compr;
  b->decode =  multihuff_decompr;
  b->enc_start = mh_init_compr;
  b->enc_stop = mh_stop_compr; 
  b->dec_start = mh_init_decompr;
  b->dec_stop = NULL; 
}

void init_mh0_compressor(base_compressor *b)
{
  b->name = "RLE0 + Multiple-Table Huffman Coding";
  b->bound = mh_h0star;
  b->cost = multihuff0_cost;
  b->merge = NULL;
  b->encode = multihuff0_compr;
  b->decode =  multihuff0_decompr;
  b->enc_start = mh_init_compr;
  b->enc_stop = mh_stop_compr; 
  b->dec_start = mh_init_decompr;
  b->dec_stop = NULL; 
}




/* ********************************************************************
   rle+compression of a string using Huffman with multiple tables 
   input
     int   len          size of mtf sequence
     uchar *in          input mtf sequence
     int   alpha_size   size of the alphabet
     int   rle_flag     rle strategy used (see note)
   output
     the number of bits produced by the encoding
     if Outfile!=NULL output stream is written to Outfile
   Note
     If rle_flag==0 rle is used only for the runs of zeros. A run of 
     length t is encoded writing t+1 in binary (using the symbols 0 and 1) 
     omitting the MSB (this is also called 1/2 encoding). Any other 
     symbol c is encoded with c+1 and without rle.
     If rle_flag!=0 rle is used for all runs of symbols. a run of t c's
     is encoded with t written in binary MSB first using c as 1 and 
     RLE_symbol as 0. This mode is designed for use without MTF,
     but used without partitioning achieves a poor compression and 
     It is too slow for a pratical use with the booster with true-cost
     (in addition I think it has a large overhead for each segment so 
     it is not suitable for being used in many segments; this is not
     surprising since multihuff was designed to work with mtf in one shot). 
   ******************************************************************* */
int multihuff_compr_aux(uchar *in, int len, int alpha_size, int rle_flag)
{
  int v, t, i, j, gs, ge, totc, bt, bc, iter;
  int nSelectors, minLen, maxLen, new_len;
  int nGroups, RLE_symbol, EOF_symbol;
  uint16 cost[BZ_N_GROUPS];
  int  fave[BZ_N_GROUPS];
  uint16* mtfv;
  uchar *selector;

  // --- initialize bit count ---
  MH_bits_total=0;  
  RLE_symbol=alpha_size;       // this is used only when rle_flag!=0
  EOF_symbol=alpha_size+1;     // this is always used
  alpha_size +=2;

  // --- allocate temporary buffer
  mtfv = (uint16 *) malloc((len+1)*sizeof(uint16));
  if(mtfv==NULL) out_of_mem("multiuhf_compr");

  // encode runs of equal symbols; write the result to mtfv[]
  new_len=0;
  if(rle_flag!=0) {
    // encode each run of symbols in binary using RLE_symbol 
    int i, next, count,bits,mask,b;
    for(i=0;i<len;) {
      next=in[i];                       // first symbol of a run
      count=1;                          // initial length of the run
      for(i=i+1;i<len;i++)              // compute run length
	if(in[i]==next) count++;
	else break;
      // --- we have "count" occurrences of "next"
      assert(count>0);
      if(count==1) mtfv[new_len++] = next;
      else {
	bits=int_log2(count);            //  number of bit to write count
	mask = ( 1 << (bits-1) );
	 for (b=bits-1; b>=0; b--) {     // write count in binary ("bits" bit)
	   if ((count & mask) == 0 )     // get b-th bits of count
	     mtfv[new_len++]=RLE_symbol; // RLE_symbol stands for 0
	   else
	     mtfv[new_len++] = next;     // next stands for 1
	   count <<= 1;                  // avoid dangerous right shift of mask
	 }
      }
    }
  }
  else {     // run length encode runs of 0's
    int i, count;
    for(i=0;i<len; ) {
      if(in[i]!=0)
	mtfv[new_len++]=in[i++]+1;    // encode c as c+1
      else {                          // a run of 0 starts here
	count=2;                      // for 1/2 encoding we need (len+1)
	while(++i<len && in[i]==0)
	  count++;
        while(count>1) {              // write 1-2 encoding of count
	  mtfv[new_len++]= (count&1) ? 1:0; // write bit 
	  count = count>>1;
	}
      } // end else
    } // end for
  }  // end rle of rlechar

  // end of segment: encode an EOF_symbol
  mtfv[new_len++] = EOF_symbol;
  if(Verbose > 3) 
    fprintf(stderr,"block size after rle encoding: %d\n",new_len);
   
  // init mtf_freq[]
  for(i=0;i<alpha_size;i++) mtf_freq[i]=0;
  for(i=0;i<new_len;i++) mtf_freq[mtfv[i]]++;
  // init huf_len[][]
  for (t = 0; t < BZ_N_GROUPS; t++)
    for (v = 0; v < alpha_size; v++)
      huf_len[t][v] = BZ_GREATER_ICOST;
  // alloc selector[]
  selector = (uchar *) malloc((1+new_len/BZ_G_SIZE)*sizeof(uchar));
  if(selector==NULL) out_of_mem("multihuf_compr");
  
  /*--- Decide how many coding tables to use ---*/
  assert(new_len > 0);
  if (new_len < 200)  nGroups = 2; else
    if (new_len < 600)  nGroups = 3; else
      if (new_len < 1200) nGroups = 4; else
	if (new_len < 2400) nGroups = 5; else
	  nGroups = 6;

  /*--- Generate an initial set of coding tables ---*/
  // each table uses BZ_LESSER_ICOST for a group of consecutive 
  // chars (gs to ge) and BZ_GREATER_ICOST for the others chars
  { 
    int nPart, remF, tFreq, aFreq;
    
    nPart = nGroups;
    remF  = new_len;
    gs = 0;
    while (nPart > 0) {
      tFreq = remF / nPart;
      ge = gs-1;
      aFreq = 0;
      while (aFreq < tFreq && ge < alpha_size-1) {
	ge++;
	aFreq += mtf_freq[ge];
      }
      if (ge > gs 
	  && nPart != nGroups && nPart != 1 
	  && ((nGroups-nPart) % 2 == 1)) {
	aFreq -= mtf_freq[ge];
	ge--;
      }
      if (Verbose > 3)
	fprintf(stderr,"      initial group %d, [%d .. %d], "
		"has %d syms (%4.1f%%)\n",
		nPart, gs, ge, aFreq, 
		(100.0 * (float)aFreq) / (float)(new_len) );
      for (v = 0; v < alpha_size; v++)
	if (v >= gs && v <= ge) 
	  huf_len[nPart-1][v] = BZ_LESSER_ICOST; else
	    huf_len[nPart-1][v] = BZ_GREATER_ICOST;
      nPart--;
      gs = ge+1;
      remF -= aFreq;
    }
  }

  /*--- 
    Iterate up to BZ_N_ITERS times to improve the tables.
    ---*/
  for (iter = 0; iter < BZ_N_ITERS; iter++) {
    for (t = 0; t < nGroups; t++) fave[t] = 0;
    for (t = 0; t < nGroups; t++)
      for (v = 0; v < alpha_size; v++)
	rfreq[t][v] = 0;
    nSelectors = 0;
    totc = 0;
    gs = 0;
    while (True) {
      /* Set group start & end marks. --*/
      if (gs >= new_len) break;
      ge = gs + BZ_G_SIZE - 1;      // size is at most BZ_G_SIZE
      if (ge >= new_len) ge = new_len-1;
      /*-- 
	Calculate the cost of this group as coded
	by each of the coding tables.
	--*/
      for (t = 0; t < nGroups; t++) cost[t] = 0;
      if (nGroups == 6) {
	register uint16 cost0, cost1, cost2, cost3, cost4, cost5;
	cost0 = cost1 = cost2 = cost3 = cost4 = cost5 = 0;
	for (i = gs; i <= ge; i++) { 
	  uint16 icv = mtfv[i];
	  cost0 += huf_len[0][icv];
	  cost1 += huf_len[1][icv];
	  cost2 += huf_len[2][icv];
	  cost3 += huf_len[3][icv];
	  cost4 += huf_len[4][icv];
	  cost5 += huf_len[5][icv];
	}
	cost[0] = cost0; cost[1] = cost1; cost[2] = cost2;
	cost[3] = cost3; cost[4] = cost4; cost[5] = cost5;
      } else {
	for (i = gs; i <= ge; i++) { 
	  uint16 icv = mtfv[i];
	  for (t = 0; t < nGroups; t++) cost[t] += huf_len[t][icv];
	}
      }
      /*-- 
	Find the coding table which is best for this group,
	and record its identity in the selector table.
	--*/
      bc = 999999999; bt = -1;
      for (t = 0; t < nGroups; t++)
	if (cost[t] < bc) { bc = cost[t]; bt = t; };
      totc += bc;
      fave[bt]++;
      selector[nSelectors++] = bt;
      /*-- 
	Increment the symbol frequencies for the selected table.
	--*/
      for (i = gs; i <= ge; i++)
	rfreq[bt][ mtfv[i] ]++;
      gs = ge+1;    // consider next group 
    }
    if (Verbose >3) {
      fprintf(stderr,"    pass %d: size is %d, grp uses are ",iter+1,totc/8 );
      for (t = 0; t < nGroups; t++)
	fprintf(stderr, "%d ", fave[t] );
      fprintf (stderr, "\n" );
    }
    /*--
      Recompute the tables based on the accumulated frequencies.
      --*/
    for (t = 0; t < nGroups; t++) {
      BZ2_hbMakeCodeLengths (&(huf_len[t][0]),&(rfreq[t][0]),alpha_size,20);
    }
  }

  /*--- Assign actual codes for the tables. --*/
  for (t = 0; t < nGroups; t++) {
    minLen = 32;
    maxLen = 0;
    for (i = 0; i < alpha_size; i++) {
      if (huf_len[t][i] > maxLen) maxLen = huf_len[t][i];
      if (huf_len[t][i] < minLen) minLen = huf_len[t][i];
    }
    assert(!(maxLen > 20));
    assert(!(minLen < 1));
    BZ2_hbAssignCodes(&(huf_code[t][0]),&(huf_len[t][0]),
		  minLen,maxLen,alpha_size);
  }

  // write coding tables (i.e codeword length).
  assert(nGroups<8);
  bit_write(3,nGroups);
  for (t = 0; t < nGroups; t++) {
    int curr = huf_len[t][0];
    bit_write(5, curr);
    for (i = 0; i < alpha_size; i++) {
      while (curr < huf_len[t][i]) { bit_write(2,2); curr++; /* 10 */ };
      while (curr > huf_len[t][i]) { bit_write(2,3); curr--; /* 11 */ };
      bit_write(1, 0 );
    }
  }
  
  /*--- write selectors and compressed data ---*/
  {
    int sel=0;
    uchar pos[BZ_N_GROUPS], ll_i, tmp2, tmp;
    for (i = 0; i < nGroups; i++) pos[i] = i;
    gs = 0;
    while (True) {
      if (gs >= new_len) break;
      ge = gs + BZ_G_SIZE - 1; 
      if (ge >= new_len) ge = new_len-1;  // establish group boundaries
      assert( selector[sel] < nGroups);
      {
	ll_i=selector[sel];              // get mtf rank for selector
	j = 0;
	tmp = pos[j];
	while ( ll_i != tmp ) {
	  j++; 
	  tmp2 = tmp; tmp = pos[j]; pos[j] = tmp2;
	};
	pos[0] = tmp;
	bit_write(j+1,1);                // write selector mtf rank in unary 
       }
      for (i = gs; i <= ge; i++) {
	assert(mtfv[i]<alpha_size);
	bit_write(huf_len[selector[sel]][mtfv[i]],
		  huf_code[selector[sel]][mtfv[i]]);
      }
      gs = ge+1;
      sel++;
    }
    assert( sel == nSelectors);
  }
  free(selector);
  free(mtfv);
  return MH_bits_total;
}

// ----------------------------------------------------------
// these external variables are required for the use of 
// the macro bit_read_macro() and single_bit_read()
// -----------------------------------------------------------
static uint32 Bit_buffer;       
static int Bit_buffer_size;  


/* *********************************************************
   decode a unary code read from Infile. the output is the
   # of zeroes we see before we see a 1 
   1 --> 0
   01 --> 1
   001 --> 2  
     etc.
   ********************************************************* */
static __inline__ int decode_unary(void)
{
  int t,i=0;

  do {
    t = bit_read(1);
    if(t!=0) break;
    i++;
  } while(1);
  return i;
}


/* ********************************************************************
   This procedure reads a stream of bits from Infile, decodes it using
   multihuffman+unrle, and writes it to dest[]. The decoding ends when
   EOF is encountered. The numer of decoded symbols is returned: it is an 
   error if more the procedure tries to decode more than limit symbols.
   See multihuff_decompr_aux for an explanation of rle_flag
   ******************************************************************** */
int multihuff_decompr_aux(uchar *dest, int limit, int alpha_size, int rle_flag)
{
  int t, i, j, minLen, maxLen, nGroups, written=0;
  int RLE_symbol, EOF_symbol;

  // --- define RLE and EOF symbols
  RLE_symbol=alpha_size;
  EOF_symbol=alpha_size+1;
  alpha_size +=2;

  // get number of groups
  bit_read_macro(nGroups,3);
  /*--- get the coding tables ---*/
  {
    int curr,uc;

    for (t = 0; t < nGroups; t++) {
      bit_read_macro(curr,5);
      for (i = 0; i < alpha_size; i++) {
	while (True) {
	  if (curr < 1 || curr > 20) 
	    fatal_error("multihuf_decompr");
	  single_bit_read_macro(uc);
	  if (uc == 0) break;
	  single_bit_read_macro(uc);
	  if (uc == 0) curr++; else curr--;
	}
        huf_len[t][i] = curr;
      }
    }

    /*--- Create the Huffman decoding tables ---*/
    for (t = 0; t < nGroups; t++) {
      minLen = 32;
      maxLen = 0;
      for (i = 0; i < alpha_size; i++) {
	if (huf_len[t][i] > maxLen) maxLen = huf_len[t][i];
	if (huf_len[t][i] < minLen) minLen = huf_len[t][i];
      }
      BZ2_hbCreateDecodeTables ( 
            &(huf_limit[t][0]), 
            &(huf_base[t][0]), 
            &(huf_perm[t][0]), 
            &(huf_len[t][0]),
            minLen, maxLen, alpha_size
	    );
      huf_minLens[t] = minLen;
    }
  }

   /*------- uncompress data -------*/
  {
     int next,count,rle_sofar,cur_symbol,rank,gSel,to_be_read;
     int zn,zj,zvec, *gLimit, *gPerm, *gBase;
     uchar pos[BZ_N_GROUPS], gMinlen=0;

     gLimit=gPerm=gBase=NULL;  // to avoid annoying compiler warnings
     for (i = 0; i < nGroups; i++) pos[i] = i;
     count=rle_sofar=0; cur_symbol=-1;
     to_be_read=0;
     while (True) {
       if(to_be_read==0) {
	 to_be_read = BZ_G_SIZE;                       
	 rank=decode_unary();    // get mtf rank of new group
	 assert(rank<nGroups);
	 gSel=pos[rank];
	 for(j=rank;j>0;j--)  pos[j]=pos[j-1];
	 pos[0]=(uchar) gSel;
	 // get tables for this group
	 gMinlen = huf_minLens[gSel];                 
	 gLimit = &(huf_limit[gSel][0]);              
	 gPerm = &(huf_perm[gSel][0]);                
	 gBase = &(huf_base[gSel][0]);                
       }
       to_be_read--;
       // get next huffman encoded char
       zn = gMinlen;                  
       // zvec = bit_read(zn);
       bit_read_macro(zvec,zn);
       while (zvec > gLimit[zn]) {                    
	 zn++;                                    
	 // zj=bit_read(1);
	 single_bit_read_macro(zj);
	 zvec = (zvec << 1) | zj;                    
       };                                             
       next = gPerm[zvec - gBase[zn]];                
       // decode next
       assert(next<alpha_size);

       if(rle_flag!=0) {
	 if(next==cur_symbol)
	   count =  (count<<1)+1;      // bit 1
	 else if(next==RLE_symbol)
	   count = count<<1;           // bit 0
	 else {                        // next starts a new run 
	   // ---- there are count occ of runchar
	   if(written+count>limit) 
	     fatal_error("Plain text buffer full! (multihuf_decompr)\n");
	   for(i=count;i>0;i--) 
	     dest[written++] = (uchar) cur_symbol;
	   if(next==EOF_symbol) 
	     break;  // end of segment
	   else {
	     cur_symbol = next;
	     assert(cur_symbol<ALPHA_SIZE);
	     count=1;
	   }
	 }
       } // end RLE_ALL decoding
       else { // rle_flag==0 
	 if(next==0)  
	   count += (1<<rle_sofar++);
	 else if(next==1)
	   count += (2<<rle_sofar++);
	 else {
	   if(count>0) {
	     // write pending run of 0's
	     if(written+count>limit) 
	       fatal_error("Plain text buffer full! (multihuf_decompr)\n");
	     for(i=count;i>0;i--) 
	       dest[written++] = (uchar) 0;
	     count=rle_sofar=0;
	   }
	   if(next==EOF_symbol)
	     break;
	   if(written>=limit)
	     fatal_error("Plain text buffer full! (multihuf_decompr)\n");
	   dest[written++] = (uchar) (next-1);
	 } 
       } // end rle_flag==0 decoding
     }  // end while       
  }
  return written;
}

static void out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f);
  exit(1);
}

static int int_log2(int u)    
{
  int i = 1;
  int r = 1;
  while((i<=32) && (r<u)){
    r=2*r+1;
    i = i+1;
  }
      assert(i<=32);
  return i;
}

static void fatal_error(char *s)
{
  fprintf(stderr,"Fatal error in function %s!\n",s);
  exit(1);
}


/* ********************************************************************
   Low-level procedures to read and write bits using a Bit_buffer. 
   Unread/unwritten bits of Bit_buffer are the most significant ones.
   ******************************************************************** */
static uint32 Bit_buffer;     /* 32 bit buffer */
static int  Bit_buffer_size;  /* number of unread/unwritten bits in Bit_buffer */

/*!
   Initializes the Bit_buffer to zero */
static void init_bit_buffer(void)
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

/*! Write to f n bits taken from v (possibly n > 24) */
void bit_write(int n, uint32 v)
{  
  assert(n <= 32);

  MH_bits_total += n;
  if(Outfile!=NULL) {
    if (n > 24){
      fbit_write24(Outfile, n-24, (v>>24) & 0xffL);
      fbit_write24(Outfile, 24, v & 0xffffffL);
    } else {
      fbit_write24(Outfile,n,v);
    }
  }
}


/*! flush the content of Bit_buffer, padding the unused bits with 0's */
static void bit_flush(void)
{
  if(Outfile!=NULL) {
    if(Bit_buffer_size!=0) {
      assert(Bit_buffer_size<8);
      fbit_write24(Outfile, 8 - (Bit_buffer_size%8), 0);  // pad with zero !
    }
  }
}


// ###### bit input routines
/*!  Read n (<=24) bits from file f (and Bit_buffer) */
__inline__ static uint32 fbit_read24(FILE *f, int n)
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

/*! Read n bits from Bit_buffer */
static uint32 bit_read(int n)
{  
  uint32 u = 0;

  assert(n <= 32);
  if (n > 24){
    u =  fbit_read24(Infile, n-24)<<24;
    u |= fbit_read24(Infile, 24);
    return u;
  } else 
    return fbit_read24(Infile,n);
}

/*-------------------------------------------------------------*/
/*--- Huffman coding low-level stuff                        ---*/
/*---                                             huffman.c ---*/
/*-------------------------------------------------------------*/

/*--
  This file is a part of bzip2 and/or libbzip2, a program and
  library for lossless, block-sorting data compression.

  Copyright (C) 1996-2000 Julian R Seward.  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

  2. The origin of this software must not be misrepresented; you must 
     not claim that you wrote the original software.  If you use this 
     software in a product, an acknowledgment in the product 
     documentation would be appreciated but is not required.

  3. Altered source versions must be plainly marked as such, and must
     not be misrepresented as being the original software.

  4. The name of the author may not be used to endorse or promote 
     products derived from this software without specific prior written 
     permission.

  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
  OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Julian Seward, Cambridge, UK.
  jseward@acm.org
  bzip2/libbzip2 version 1.0 of 21 March 2000

  This program is based on (at least) the work of:
     Mike Burrows
     David Wheeler
     Peter Fenwick
     Alistair Moffat
     Radford Neal
     Ian H. Witten
     Robert Sedgewick
     Jon L. Bentley

  For more information on these sources, see the manual.
--*/


/*---------------------------------------------------*/
#define WEIGHTOF(zz0)  ((zz0) & 0xffffff00)
#define DEPTHOF(zz1)   ((zz1) & 0x000000ff)
#define MYMAX(zz2,zz3) ((zz2) > (zz3) ? (zz2) : (zz3))

#define ADDWEIGHTS(zw1,zw2)                           \
   (WEIGHTOF(zw1)+WEIGHTOF(zw2)) |                    \
   (1 + MYMAX(DEPTHOF(zw1),DEPTHOF(zw2)))

#define UPHEAP(z)                                     \
{                                                     \
   Int32 zz, tmp;                                     \
   zz = z; tmp = heap[zz];                            \
   while (weight[tmp] < weight[heap[zz >> 1]]) {      \
      heap[zz] = heap[zz >> 1];                       \
      zz >>= 1;                                       \
   }                                                  \
   heap[zz] = tmp;                                    \
}

#define DOWNHEAP(z)                                   \
{                                                     \
   Int32 zz, yy, tmp;                                 \
   zz = z; tmp = heap[zz];                            \
   while (True) {                                     \
      yy = zz << 1;                                   \
      if (yy > nHeap) break;                          \
      if (yy < nHeap &&                               \
          weight[heap[yy+1]] < weight[heap[yy]])      \
         yy++;                                        \
      if (weight[tmp] < weight[heap[yy]]) break;      \
      heap[zz] = heap[yy];                            \
      zz = yy;                                        \
   }                                                  \
   heap[zz] = tmp;                                    \
}


/*---------------------------------------------------*/
static void BZ2_hbMakeCodeLengths ( UChar *len, 
                             Int32 *freq,
                             Int32 alphaSize,
                             Int32 maxLen )
{
   /*--
      Nodes and heap entries run from 1.  Entry 0
      for both the heap and nodes is a sentinel.
   --*/
   Int32 nNodes, nHeap, n1, n2, i, j, k;
   Bool  tooLong;

   Int32 heap   [ BZ_MAX_ALPHA_SIZE + 2 ];
   int64_t weight [ BZ_MAX_ALPHA_SIZE * 2 ];
   Int32 parent [ BZ_MAX_ALPHA_SIZE * 2 ]; 

   for (i = 0; i < alphaSize; i++)
      weight[i+1] = ((int64_t) (freq[i] == 0 ? 1 : freq[i])) << 8;

   while (True) {

      nNodes = alphaSize;
      nHeap = 0;

      heap[0] = 0;
      weight[0] = 0;
      parent[0] = -2;

      for (i = 1; i <= alphaSize; i++) {
         parent[i] = -1;
         nHeap++;
         heap[nHeap] = i;
         UPHEAP(nHeap);
      }

      assert( nHeap < (BZ_MAX_ALPHA_SIZE+2));
   
      while (nHeap > 1) {
         n1 = heap[1]; heap[1] = heap[nHeap]; nHeap--; DOWNHEAP(1);
         n2 = heap[1]; heap[1] = heap[nHeap]; nHeap--; DOWNHEAP(1);
         nNodes++;
         parent[n1] = parent[n2] = nNodes;
         weight[nNodes] = ADDWEIGHTS(weight[n1], weight[n2]);
         parent[nNodes] = -1;
         nHeap++;
         heap[nHeap] = nNodes;
         UPHEAP(nHeap);
      }

      assert( nNodes < (BZ_MAX_ALPHA_SIZE * 2));

      tooLong = False;
      for (i = 1; i <= alphaSize; i++) {
         j = 0;
         k = i;
         while (parent[k] >= 0) { k = parent[k]; j++; }
         len[i-1] = j;
         if (j > maxLen) tooLong = True;
      }
      
      if (! tooLong) break;

      for (i = 1; i < alphaSize; i++) {
         j = weight[i] >> 8;
         j = 1 + (j / 2);
         weight[i] = j << 8;
      }
   }
}


/*---------------------------------------------------*/
static void BZ2_hbAssignCodes ( Int32 *code,
                         UChar *length,
                         Int32 minLen,
                         Int32 maxLen,
                         Int32 alphaSize )
{
   Int32 n, vec, i;

   vec = 0;
   for (n = minLen; n <= maxLen; n++) {
      for (i = 0; i < alphaSize; i++)
         if (length[i] == n) { code[i] = vec; vec++; };
      vec <<= 1;
   }
}


/*---------------------------------------------------*/
static void BZ2_hbCreateDecodeTables ( Int32 *limit,
                                Int32 *base,
                                Int32 *perm,
                                UChar *length,
                                Int32 minLen,
                                Int32 maxLen,
                                Int32 alphaSize )
{
   Int32 pp, i, j, vec;

   pp = 0;
   for (i = minLen; i <= maxLen; i++)
      for (j = 0; j < alphaSize; j++)
         if (length[j] == i) { perm[pp] = j; pp++; };

   for (i = 0; i < BZ_MAX_CODE_LEN; i++) base[i] = 0;
   for (i = 0; i < alphaSize; i++) base[length[i]+1]++;

   for (i = 1; i < BZ_MAX_CODE_LEN; i++) base[i] += base[i-1];

   for (i = 0; i < BZ_MAX_CODE_LEN; i++) limit[i] = 0;
   vec = 0;

   for (i = minLen; i <= maxLen; i++) {
      vec += (base[i+1] - base[i]);
      limit[i] = vec-1;
      vec <<= 1;
   }
   for (i = minLen + 1; i <= maxLen; i++)
      base[i] = ((limit[i-1] + 1) << 1) - base[i];
}


/*-------------------------------------------------------------*/
/*--- end                                         huffman.c ---*/
/*-------------------------------------------------------------*/

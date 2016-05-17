/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   rc_defs.c
   Ver 1.0 (27 mar 2006) Giovanni Manzini (manzini@mfn.unipmn.it)

   Range coding encoding/decoding for use in the compression booster
   These routines encode a segment t[0..n-1] using RLE for each
   run of equal symbols followed by order0 arithmetic coding. We use
   use the ESC symbol to introduce symbols not seen before. These new
   symbols are generate using a flat model.

   The actual compression is done using the Carryless rangecoder 
   (c) 1999 by Dmitry Subbotin (Modified into C and extended 2001 
   by Mikael Lundqvist) see file range.c 

   This program is largely based on the source code from the 1987 CACM
   article by Witten, Neal, and Cleary.  
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "booster.h"
#include "range.h"


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

// ======== prototypes ==============
static void fatal_error(char *s);
static uint32 uint32_read(FILE *f);
static void uint32_write(FILE *f, uint32 u);
static int int_log2(int u);
static void init_update_parameters(void);

// ===== default values for updating the order0 model ============ 
#define default_Max_frequency 16383     // maximum allowed frequency count
#define default_Increment 1             // freq. increment for each occurrence

// ===== local global variables =====================     
// order0 model parameters
static int Max_frequency; 
static int Increment;
// encoder and decoder data structures        
static rc_encoder rc;
static rc_decoder rd;


// ============== decompression functions ==========


/* ****************************************************************
   procedure that reads a stream of bits from the decoder data 
   structure rd and decode them using range coding followed by 1/2 unrle.
   The decoded symbols are written to t[]. The decoding stops when 
   the EOF symbol is encountered. No more than n symbols can be decoded.
   **************************************************************** */
static int rc_rle_decompr_aux(symbol *t, int n, int alpha_size)
{
  int r,written,next,count,i,runchar,RLE_symbol,EOF_symbol; 
  Uint *freq_tbl;
  rc_model rcm;

  // init RLE and EOF symbols
  RLE_symbol = alpha_size;
  EOF_symbol = alpha_size+1;
  // freq table for alphabet symbols and EOF
  freq_tbl = (Uint *) malloc(sizeof(Uint)*(alpha_size+2));
  if(freq_tbl==NULL)
    fatal_error("Out_of_memory (rc_rle_decompr_aux)\n");

  // init model 
  r=ModelInit(&rcm,alpha_size+2,freq_tbl,NULL,Increment,Max_frequency,ADAPT);
  if (r != RC_OK)
    fatal_error("Error: maxFreq or totFreq > 65536 (ModelInit-->)\n");

  // start decoding
  written=0;
  next=DecodeSymbol(&rd, &rcm);
  while(next!=EOF_symbol){
    count=1;
    runchar=next;
    next=DecodeSymbol(&rd, &rcm);
    while(next==runchar || next==RLE_symbol) {
      if(next==RLE_symbol)
	count = (count<<1);        // RLE stands for bit 0 
      else
	count = (count<<1)+1;      // runchar stands for bit 1
      next = DecodeSymbol(&rd, &rcm);
    }
    assert(runchar>=0 && runchar<alpha_size);
    assert(count>=1);
    // ---- there are count occ of runchar
    if(written+count>n) 
       fatal_error("Plain text buffer full! (rc_rle_decompr_aux)\n");
    for(i=count;i>0;i--)
      t[written++] = (symbol) runchar;
  }
  free(freq_tbl);
  return written;
}

// ******* as before, but rle is done only on zeroes ********
static int rc_rle0_decompr_aux(symbol *t, int n, int alpha_size)
{
  int r,written,next,count,i,RLE_symbol,EOF_symbol; 
  Uint *freq_tbl;
  rc_model rcm;

  // init RLE and EOF symbols
  RLE_symbol = alpha_size;
  EOF_symbol = alpha_size+1;
  // freq table for alphabet symbols and EOF
  freq_tbl = (Uint *) malloc(sizeof(Uint)*(alpha_size+2));
  if(freq_tbl==NULL)
    fatal_error("Out_of_memory (rc_rle0_decompr_aux)\n");

  // init model 
  r=ModelInit(&rcm,alpha_size+2,freq_tbl,NULL,Increment,Max_frequency,ADAPT);
  if (r != RC_OK)
    fatal_error("Error: maxFreq or totFreq > 65536 (ModelInit-->)\n");

  // start decoding
  written=0;
  next=DecodeSymbol(&rd, &rcm);
  while(next!=EOF_symbol) {
    if(next!=0 && next!=RLE_symbol) {
      assert(next>0 && next<alpha_size);
      if(written>=n) 
	fatal_error("Plain text buffer full! (rc_rle0_decompr_aux)\n");
      t[written++] = (symbol) next;
      next=DecodeSymbol(&rd, &rcm);;
    }
    else {
      count=1;
      while(next==0 || next==RLE_symbol) {
	if(next==RLE_symbol)
	  count = (count<<1) + 1;   // RLE_symbol stands for bit 1 
	else
	  count = (count<<1);       // 0 stands for bit 0
	next = DecodeSymbol(&rd, &rcm);;
      }
      assert(count>1);
      // ---- there are count-1 occ of 0
      if(written+count-1>n) 
	fatal_error("Plain text buffer full! (rc_rle0_decompr_aux)\n");
      for(i=count-1;i>0;i--)
	t[written++] = (symbol) 0;
    }
  }
  free(freq_tbl);
  return written;
}
/* ************ as before but without RLE ************* */
static int rc_decompr_aux(symbol *t, int n, int alpha_size)
{
  int r,written,next; 
  Uint *freq_tbl;
  rc_model rcm;

  // freq table for alphabet symbols and EOF
  freq_tbl = (Uint *) malloc(sizeof(Uint)*(alpha_size+1));
  if(freq_tbl==NULL)
    fatal_error("Out_of_memory (rc_decompr_aux)\n");

  // init model 
  r=ModelInit(&rcm,alpha_size+1,freq_tbl,NULL,Increment,Max_frequency,ADAPT);
  if (r != RC_OK)
    fatal_error("Error: maxFreq or totFreq > 65536 (ModelInit-->)\n");

  written=0;
  while(1) {
    next=DecodeSymbol(&rd, &rcm);
    if(rd.error!=RC_OK) 
      fatal_error("Decoding Error (DecodeSymbol-->)\n");
    if(next==alpha_size) break;
    assert(next>=0 && next<alpha_size);
    if(written>=n) 
       fatal_error("Plain text buffer full! (rc_decompr_aux)\n");
    t[written++] = (symbol) next;
  }
  free(freq_tbl);
  return written;
}


// ================= compression funcions ===================

/* ****************************************************************
   procedure that encodes the segment t[0..n-1] using 1/2 rle 
   followed by range coding via the rc data structure.
   Return the number of bytes written. This value is not 100% accurate
   since from a call to this function to the next some information
   is left "pending" in the status of the current range. However, the 
   total number of produced bytes (over a several successive calls) 
   should be the correct one.
   **************************************************************** */
static int rc_rle_compr_aux(symbol *t, int n, int alpha_size)
{
  int next,count,i,b,bits,mask,r;
  int RLE_symbol,EOF_symbol,already_passed; 
  Uint *freq_tbl;
  rc_model rcm;

  // init RLE and EOF symbols
  RLE_symbol = alpha_size;
  EOF_symbol = alpha_size+1;
  // freq table for alphabet symbols and EOF
  freq_tbl = (Uint *) malloc(sizeof(Uint)*(alpha_size+2));
  if(freq_tbl==NULL)
    fatal_error("Out_of_memory (rc_rle_compr_aux)\n");

  // init model
  r=ModelInit(&rcm,alpha_size+2,freq_tbl,NULL,Increment,Max_frequency,ADAPT);
  if (r != RC_OK) 
    fatal_error("Error: maxFreq or totFreq > 65536 (ModelInit-->)\n");

  already_passed=rc.passed;  // number of bytes already encoded by rc

  for(i=0;i<n;) {
    next=t[i];                        // first symbol of a run
    count=1;                          // initial length of the run
    for(i=i+1;i<n;i++)                // compute run length
      if(t[i]==next) count++;
      else break;
    // --- we have "count" occurrences of "next"
    assert(count>0);
    if(count==1) {
      EncodeSymbol(&rc,&rcm,next);
      if(rc.error!=RC_OK)
	fatal_error("Encoding Error (EncodeSymbol-->)\n");
    }
    else {
      bits=int_log2(count);   	      //  number of bit to write count
      mask = ( 1 << (bits-1) );
      for (b=bits-1; b>=0; b--) {     // write count in binary ("bits" bit)
	if ((count & mask) == 0 ) {           // get b-th bits of count
	  EncodeSymbol(&rc,&rcm,RLE_symbol);  // RLE_symbol stands for 0
	  if(rc.error!=RC_OK)
	    fatal_error("Encoding Error (EncodeSymbol-->)\n");
	} 
	else {
 	  EncodeSymbol(&rc,&rcm,next);        // next stands for 1
	  if(rc.error!=RC_OK)
	    fatal_error("Encoding Error (EncodeSymbol-->)\n");
	}
    	count <<= 1;                  // avoid dangerous right shift of mask
      }
    }
  }
  // end of segment: encode an EOF_symbol
  EncodeSymbol(&rc,&rcm,EOF_symbol);  
  if(rc.error!=RC_OK)
    fatal_error("Encoding Error (EncodeSymbol-->)\n");

  free(freq_tbl);
  return 8*(rc.passed-already_passed);  // not 100% accurate (see above)
}
/* ************************************************************
   As before, but rle is done only on runs of zeros. 
   ************************************************************ */
static int rc_rle0_compr_aux(symbol *t, int n, int alpha_size)
{
  int next,count,i,b,bits,mask,r; 
  int RLE_symbol,EOF_symbol,already_passed; 
  Uint *freq_tbl;
  rc_model rcm;

  // init RLE and EOF symbols
  RLE_symbol = alpha_size;
  EOF_symbol = alpha_size+1;
  // freq table for alphabet symbols and EOF
  freq_tbl = (Uint *) malloc(sizeof(Uint)*(alpha_size+2));
  if(freq_tbl==NULL)
    fatal_error("Out_of_memory (rc_rle0_compr_aux)\n");

  // init model
  r=ModelInit(&rcm,alpha_size+2,freq_tbl,NULL,Increment,Max_frequency,ADAPT);
  if (r != RC_OK) 
    fatal_error("Error: maxFreq or totFreq > 65536 (ModelInit-->)\n");

  already_passed=rc.passed;  // number of bytes already encoded by rc

  // encode t[0,n) + EOF
  for(i=0;i<n;) {
    next=t[i++];                      // first symbol of a run
    if(next!=0) {                     // nonzero symbol
      EncodeSymbol(&rc,&rcm,next);
      if(rc.error!=RC_OK)
	fatal_error("Encoding Error (EncodeSymbol-->)\n");
      continue;
    }
    count=1;                          // initial length of the run of 0's
    while(i<n && t[i]==0) {
      count++; i++;                   // compute run length
    }
    // --- we have "count" occurrences of 0: encode them with 1/2
    assert(count>0);
    count++;                          // encode count+1 excluding MSB 
    bits=int_log2(count);             // number of bits to write count
    assert(bits>=2);                  
    mask = ( 1 << (bits-1) );         // position of the MSB
    count <<=1;                       // discard MSB
    for (b=bits-2; b>=0; b--) {       // write count using 1/2 encoding
      if ((count & mask) == 0 ) {   // get b-th bits of count
	EncodeSymbol(&rc,&rcm,0);   // 0 stands for 0
	if(rc.error!=RC_OK)
	  fatal_error("Encoding Error (EncodeSymbol-->)\n");
      } 
      else {
	EncodeSymbol(&rc,&rcm,RLE_symbol);// RLE_symbol stands for 1
	if(rc.error!=RC_OK)
	  fatal_error("Encoding Error (EncodeSymbol-->)\n");
      }
      count <<= 1;                    // avoid dangerous right shift of mask
    }
  }
  // end of segment: encode an EOF_symbol
  EncodeSymbol(&rc,&rcm,EOF_symbol);  
  if(rc.error!=RC_OK)
    fatal_error("Encoding Error (EncodeSymbol-->)\n");

  free(freq_tbl);
  return 8*(rc.passed-already_passed);  // not 100% accurate (see above)
}
/* ************************************************************
   As before but without RLE
   ************************************************************ */
static int rc_compr_aux(symbol *t, int n, int alpha_size)
{
  int i,r,already_passed; 
  Uint *freq_tbl;
  rc_model rcm;

  // freq table for alphabet symbols and EOF
  freq_tbl = (Uint *) malloc(sizeof(Uint)*(alpha_size+1));
  if(freq_tbl==NULL)
    fatal_error("Out_of_memory (rc_compr_aux)\n");

  // init model
  r=ModelInit(&rcm,alpha_size+1,freq_tbl,NULL,Increment,Max_frequency,ADAPT);
  if (r != RC_OK) 
    fatal_error("Error: maxFreq or totFreq > 65536 (ModelInit-->)\n");

  already_passed=rc.passed;  // number of bytes already encoded by rc

  // encode t[0,n) + EOF
  for(i=0;i<n;i++) {
    assert(t[i]<alpha_size);
    EncodeSymbol(&rc,&rcm,t[i]);
    if(rc.error!=RC_OK)
      fatal_error("Encoding Error (EncodeSymbol-->)\n");
  }
  EncodeSymbol(&rc,&rcm,alpha_size); // encode EOF symbol
  if(rc.error!=RC_OK)
    fatal_error("Encoding Error (EncodeSymbol-->)\n");

  free(freq_tbl);
  return 8*(rc.passed-already_passed); 
}



// =========== cost estimate functions ========== 


/* ************************************************************************
   Cost encoding t[0,n) using 1/2 RLE+range encoding. The cost is an upper 
   estimate since at the end of the encoding we call FinishEncode() that 
   encodes the current status (range) using 4 bytes 
   ************************************************************************ */
static int rc_rle_cost(symbol *t, int n, int alpha_size)
{

  // init update parameters
  init_update_parameters();
  // init encoder in dummy (i.e. no output) mode
  StartEncode(&rc, NULL, NULL, 0);
  // compress
  rc_rle_compr_aux(t,n,alpha_size);
  FinishEncode (&rc);                   // flush buffer
  return 8*rc.passed;                   // return number of bits produced

}
/* ************************************************************
   As before, but rle is done only on runs of zeros. 
   ************************************************************ */
static int rc_rle0_cost(symbol *t, int n, int alpha_size)
{

  // init update parameters
  init_update_parameters();
  // init encoder in dummy (i.e. no output) mode
  StartEncode(&rc, NULL, NULL, 0);
  // compress
  rc_rle0_compr_aux(t,n,alpha_size);
  FinishEncode (&rc);                   // flush buffer
  return 8*rc.passed;                   // return number of bits produced
}
/* ************************************************************
   As before but without RLE
   ************************************************************ */
static int rc_cost(symbol *t, int n, int alpha_size)
{

  // init update parameters
  init_update_parameters();
  // init encoder in dummy (i.e. no output) mode
  StartEncode(&rc, NULL, NULL, 0);
  // compress
  rc_compr_aux(t,n,alpha_size);
  FinishEncode (&rc);                   // flush buffer
  return 8*rc.passed;                   // return number of bits produced
}




static int rc_init_encoding(void)
{
  // write update parameters to _boost_Outfile
  init_update_parameters();
  uint32_write(_boost_Outfile,((Max_frequency-1)<<16)+Increment);
  // init encoder: prepare for writing to _boost_Outfile
  StartEncode(&rc, _boost_Outfile, NULL, 0);
  return 32;                 // number of bits used for update paramters
}


static int rc_stop_encoding(void)
{
  int already_passed;

  already_passed = rc.passed;
  FinishEncode (&rc);                      // flush buffer
  return 8*(rc.passed-already_passed);     // currently this is always 32 ....
}


static int rc_init_decoding(void)
{
  uint32 t;

  // read and check model parameters
  t = uint32_read(_boost_Infile);
  Increment = t & 0xFFFF;
  Max_frequency = 1+((t>>16) & 0xFFFF);
  if(Max_frequency>65536 || Increment>65535)
       fatal_error("Invalid update parameters! (rc_init_encoding)\n");

  // prepare for reading from _boost_Infile
  StartDecode(&rd, _boost_Infile, NULL, 0);

  return 0;
}



/* ***********************************************************************
   bound the oputput size using the modified 0-order entropy
   return |s|H_0^*(s) + \mu log|\Sigma|, where \mu is the -x command 
   line parameter [if -x was not used default is mu=10.0 which is
   (experimentally) a reasonable value for ac]
   ********************************************************************** */
static double rc_h0star(stats *s, int alpha_size)
{
  int i, numchar=0;
  double mu,tot=0;

  assert(s->count_tot>0);
  mu = _boost_CmdLinePar1_flag ? _boost_CmdLinePar1 : 10.0;
  for(i=0;i<alpha_size;i++)
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i]);
      numchar++;
    }
  assert(tot>=0);
  if(numchar == 1) 
    return mu * log(alpha_size) + (1+log(s->count_tot+1));
  return mu * numchar * log(alpha_size) + tot;
}

// initialize the structs for range coding compression
void init_rc_compressor(base_compressor *b)
{
  b->name = "Range Coding (no RLE)";
  b->bound = rc_h0star;
  b->cost = rc_cost;
  b->merge = NULL;
  b->encode = rc_compr_aux;
  b->decode = rc_decompr_aux;
  b->enc_start = rc_init_encoding;
  b->enc_stop = rc_stop_encoding; 
  b->dec_start = rc_init_decoding;
  b->dec_stop = NULL; 
}

void init_rc_rle_compressor(base_compressor *b)
{
  b->name = "Range Coding with RLE";
  b->bound = rc_h0star;
  b->cost = rc_rle_cost;
  b->merge = NULL;
  b->encode = rc_rle_compr_aux;
  b->decode = rc_rle_decompr_aux;
  b->enc_start = rc_init_encoding;
  b->enc_stop = rc_stop_encoding; 
  b->dec_start = rc_init_decoding;
  b->dec_stop = NULL; 
}

void init_rc_rle0_compressor(base_compressor *b)
{
  b->name = "Range Coding with RLE0";
  b->bound = rc_h0star;
  b->cost = rc_rle0_cost;
  b->merge = NULL;
  b->encode = rc_rle0_compr_aux;
  b->decode = rc_rle0_decompr_aux;
  b->enc_start = rc_init_encoding;
  b->enc_stop = rc_stop_encoding; 
  b->dec_start = rc_init_decoding;
  b->dec_stop = NULL; 
}


static void fatal_error(char *s)
{
  fprintf(stderr,"%s",s);
  exit(1);
}

static void init_update_parameters(void)
{
  // init order0 model parameters
  Max_frequency = _boost_CmdLinePar2_flag ? 
                 (int) _boost_CmdLinePar2 : default_Max_frequency;
  Increment = _boost_CmdLinePar3_flag ? 
                 (int) _boost_CmdLinePar3 : default_Increment;
  // check parameters
  if(Max_frequency>65536 || Increment>65535)
       fatal_error("Invalid update parameters! (init_update_parameters)\n");
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

/* write an uint32 to f */
static void uint32_write(FILE *f, uint32 u)
{

  assert(f!=NULL);
  putc((u>>24) & 0xffL,f);
  putc((u>>16) & 0xffL,f);
  putc((u>> 8) & 0xffL,f);
  putc(      u & 0xffL,f);

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






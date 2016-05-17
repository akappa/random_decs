/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ac_defs.c
   Ver 1.0 (15 may 2005) Giovanni Manzini (manzini@mfn.unipmn.it)

   Arithmetic coding encoding/decoding for use in the compression booster
   These routines encode a segment t[0..n-1] using RLE for each
   run of equal symbols followed by order0 arithmetic coding. We use
   use the ESC symbol to introduce symbols not seen before. These new
   symbols are generate using a flat model.

   This program is largely based on the source code from the 1987 CACM
   article by Witten, Neal, and Cleary.  
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
#include <stdio.h>
#include <stdlib.h>
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


// ===== default values for updating the order0 model ============ 
#define default_Max_frequency 16383     // maximum allowed frequency count
#define default_Increment 1             // freq. increment for each occurrence


// ===== local global variables =====================     
static FILE *Infile, *Outfile;
//static int Use_rle;             // set to !=0 if you want RLE
static int No_of_symbols;    	// # of symbols (including EOF and ESC)
static int No_of_chars;         // # of chars
static int Ac_bits_total;       // # of bits emitted by ac
// special symbols
static int ESC_symbol;		// Index of ESC symbol
static int RLE_symbol;		// Index of RLE symbol
static int EOF_symbol;    	// Index of EOF symbol
// order0 model parameters
static int Max_frequency; 
static int Increment;         

// prototypes for ac encoding 
static void ac_symb_output(int ch);
static void encode_symbol(int symbol,int cum_freq[] );
static void start_outputing_bits( void );
static void bit_plus_follow( int bit );
static void start_encoding( void );
static void done_outputing_bits( void );
static void done_encoding( void );
// prototypes for ac decoding
static void start_decoding( void );
static int ac_symb_input(void);
static void start_inputing_bits( void );
static int  decode_symbol(int cum_freq[] );
// misc ac prototypes
static void init_models(void);
static void update_model(int symbol);
static void update_flat_model( int symbol );
static void alloc_freq_tables(void); 
static void free_freq_tables(void);
static void init_extra_symbols(void);
static int int_log2(int u);
static void fatal_error(char *s);
static uint32 uint32_read(FILE *f);
static void uint32_write(FILE *f, uint32 u);
static void init_update_parameters(void);

/* ****************************************************************
   procedure that reads a strem of bits from _boost_Infile and decode them
   using order0 (with ESCapes for new symbols) followed by 1/2 unrle.
   The decoded symbols are written to t[]. The decoding stops when 
   the EOF symbol is encountered. No more than n symbols
   should be decoded.
   **************************************************************** */
static int ac_rle_decompr_aux(symbol *t, int n, int alpha_size)
{
  int written,next,count,i,runchar;

  Infile = _boost_Infile;
  No_of_chars = alpha_size;
  init_extra_symbols();
  alloc_freq_tables();
  init_models();

  written=0;
  next=ac_symb_input();
  while(next!=EOF_symbol){
    count=1;
    runchar=next;
    next = ac_symb_input();
    while(next==runchar || next==RLE_symbol) {
      if(next==RLE_symbol)
	count = (count<<1);        // RLE stands for bit 0 
      else
	count = (count<<1)+1;      // runchar stands for bit 1
      next = ac_symb_input();
    }
    assert(runchar>=0 && runchar<alpha_size);
    assert(count>=1);
    // ---- there are count occ of runchar
    if(written+count>n) 
       fatal_error("Plain text buffer full! (ac_decompr_aux)\n");
    for(i=count;i>0;i--)
      t[written++] = (symbol) runchar;
  }
  free_freq_tables();
  return written;
}
// ******* as before, but rle is done only on zeroes ********
static int ac_rle0_decompr_aux(symbol *t, int n, int alpha_size)
{
  int written,next,count,i;

  Infile = _boost_Infile;
  No_of_chars = alpha_size;
  init_extra_symbols();
  alloc_freq_tables();
  init_models();

  written=0;
  next=ac_symb_input();
  while(next!=EOF_symbol) {
    if(next!=0 && next!=RLE_symbol) {
      assert(next>0 && next<alpha_size);
      if(written>=n) 
	fatal_error("Plain text buffer full! (ac_rle0_decompr_aux)\n");
      t[written++] = (symbol) next;
      next=ac_symb_input();
    }
    else {
      count=1;
      while(next==0 || next==RLE_symbol) {
	if(next==RLE_symbol)
	  count = (count<<1) + 1;   // RLE_symbol stands for bit 1 
	else
	  count = (count<<1);       // 0 stands for bit 0
	next = ac_symb_input();
      }
      assert(count>1);
      // ---- there are count-1 occ of 0
      if(written+count-1>n) 
	fatal_error("Plain text buffer full! (ac_rle0_decompr_aux)\n");
      for(i=count-1;i>0;i--)
	t[written++] = (symbol) 0;
    }
  }
  free_freq_tables();
  return written;
}
/* ************ as before but without RLE ************* */
static int ac_decompr_aux(symbol *t, int n, int alpha_size)
{
  int written,next;

  Infile = _boost_Infile;
  No_of_chars = alpha_size;
  init_extra_symbols();
  alloc_freq_tables();
  init_models();

  written=0;

  while(1) {
    next=ac_symb_input();
    if(next==EOF_symbol) break;
    assert(next>=0 && next<alpha_size);
    if(written>=n) 
       fatal_error("Plain text buffer full! (ac_decompr_aux)\n");
    t[written++] = (symbol) next;
  }
  free_freq_tables();
  return written;
}



/* ****************************************************************
   procedure that encodes the segment t[0..n-1] using 1/2 rle 
   followed by order0 arithmetic coding (with ESCapes for new symbols)
   Write the output bit stream to _boost_Outfile and return the number 
   of bits in the output stream. This value is not 100% accurate
   since from a call to this function to the next some information
   is left "pending" in the status of the current range [low,high]
   with respect to the full range [0,top_value]. However, the total
   number of produced bits (over a several successive calls) should be 
   the correct one (except for <8 bits to complete a byte)
   **************************************************************** */
static int ac_rle_compr_aux(symbol *t, int n, int alpha_size)
{
  int next,count,i,b,bits,mask; 

  Ac_bits_total=0;
  No_of_chars = alpha_size;
  init_extra_symbols();		      // assign codes to RLE, EOF, and ESC  
  alloc_freq_tables();		      // allocate frequency tables 
  init_models();		      // initialize models

  for(i=0;i<n;) {
    next=t[i];                        // first symbol of a run
    count=1;                          // initial length of the run
    for(i=i+1;i<n;i++)                // compute run length
      if(t[i]==next) count++;
      else break;
    // --- we have "count" occurrences of "next"
    assert(count>0);
    if(count==1) ac_symb_output(next);
    else {
      bits=int_log2(count);   	      //  number of bit to write count
      mask = ( 1 << (bits-1) );
      for (b=bits-1; b>=0; b--) {     // write count in binary ("bits" bit)
	if ((count & mask) == 0 )     // get b-th bits of count
	  ac_symb_output(RLE_symbol); // RLE_symbol stands for 0
	else    	       
	  ac_symb_output(next);       // next stands for 1
    	count <<= 1;                  // avoid dangerous right shift of mask
      }
    }
  }
  // end of segment: encode an EOF_symbol
  ac_symb_output(EOF_symbol);
  free_freq_tables();
  return Ac_bits_total;  // not 100% accurate (see above)
}
/* ************************************************************
   As before, but rle is done only on runs of zeros. 
   ************************************************************ */
static int ac_rle0_compr_aux(symbol *t, int n, int alpha_size)
{
  int next,count,i,b,bits,mask; 

  Ac_bits_total=0;
  No_of_chars = alpha_size;
  init_extra_symbols();		      // assign codes to RLE, EOF, and ESC  
  alloc_freq_tables();		      // allocate frequency tables 
  init_models();		      // initialize models

  for(i=0;i<n;) {
    next=t[i++];                      // first symbol of a run
    if(next!=0) {                     // nonzero symbol
      ac_symb_output(next);
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
      if ((count & mask) == 0 )       // get b-th bits of count
	ac_symb_output(0);            // 0 stands for 0
      else    	       
	ac_symb_output(RLE_symbol);   // RLE_symbol stands for 1
      count <<= 1;                    // avoid dangerous right shift of mask
    }
  }
  // end of segment: encode an EOF_symbol
  ac_symb_output(EOF_symbol);
  free_freq_tables();
  return Ac_bits_total;  // not 100% accurate (see above)
}
/* ************************************************************
   As before, but without rle. NOTE: since rle is not used,
   there is no need for the RLE_symbol. We can remove it
   and save some space. I have not done it since preliminary
   reuslts with this compressor are not encouraging.
   ************************************************************ */
static int ac_compr_aux(symbol *t, int n, int alpha_size)
{
  int i; 

  Ac_bits_total=0;
  No_of_chars = alpha_size;
  init_extra_symbols();		      // assign codes to RLE, EOF, and ESC  
  alloc_freq_tables();		      // allocate frequency tables 
  init_models();		      // initialize models

  for(i=0;i<n;i++) 
    ac_symb_output(t[i]);
  // end of segment: encode an EOF_symbol
  ac_symb_output(EOF_symbol);
  free_freq_tables();
  return Ac_bits_total;  // not 100% accurate (see above)
}

static int ac_init_encoding(void)
{
  Outfile = _boost_Outfile;
  assert(Outfile!=NULL);
  init_update_parameters();
  // write update parameters to Outfile
  uint32_write(Outfile,(Max_frequency<<16)+Increment);
  // init ac encoder
  start_outputing_bits();
  start_encoding();
  return 32;                 // number of output bits
}

static int ac_flush_buffer(void)
{
  Ac_bits_total=0;                    // clear bit count
  done_encoding();
  done_outputing_bits();              // flush the buffer
  return Ac_bits_total;
}

static int ac_init_decoding(void)
{
  uint32 t;

  Infile = _boost_Infile;
  assert(Infile!=NULL);
  // read and check model parameters
  t = uint32_read(Infile);
  Increment = t & 0xFFFF;
  Max_frequency = (t>>16) & 0xFFFF;
  if(Max_frequency>65535 || Increment>65535)
       fatal_error("Invalid update parameters! (ac_init_encoding)\n");
  // init decoder
  start_inputing_bits();
  start_decoding();
  return 0;
}


/* ************************************************************************
   Procedure to estimate the # of bits to encode t[0...n-1] 
   Start from scratch the arithmetic encoding of t[0...n-1], and simulate
   the encoding without producing any output but counting the number
   of produced bits (via AC_bits_total)  
   ************************************************************************ */
static int ac_rle_cost(symbol *t, int n, int alpha_size)
{
  int next,count,i,b,bits,mask; 

  init_update_parameters();
  start_outputing_bits();             // init encoder
  start_encoding();
  Ac_bits_total=0;                    // clear bit count
  Outfile = NULL;                     // no output
  No_of_chars = alpha_size;
  init_extra_symbols();		      // assign codes to RLE, EOF, and ESC  
  alloc_freq_tables();		      // allocate frequency tables 
  init_models();		      // initialize models

  for(i=0;i<n;) {
    next=t[i];                        // first symbol of a run
    count=1;                          // initial length of the run
    for(i=i+1;i<n;i++)                // compute run length
      if(t[i]==next) count++;
      else break;
    // --- we have "count" occurrences of "next"
    assert(count>0);
    if(count==1) ac_symb_output(next);
    else {
      bits=int_log2(count);   	      //  number of bit to write count
      mask = ( 1 << (bits-1) );
      for (b=bits-1; b>=0; b--) {     // write count in binary ("bits" bit)
	if ((count & mask) == 0 )     // get b-th bits of count
	  ac_symb_output(RLE_symbol); // RLE_symbol stands for 0
	else    	       
	  ac_symb_output(next);       // next stands for 1
    	count <<= 1;                  // avoid dangerous right shift of mask
      }
    }
  }
  // end of segment: encode an EOF_symbol
  ac_symb_output(EOF_symbol);
  free_freq_tables();
  done_encoding();
  return Ac_bits_total;
}
/* ********** as before, but without rle ************ */
static int ac_rle0_cost(symbol *t, int n, int alpha_size)
{
  int i,next,count,bits,mask,b; 

  init_update_parameters();
  start_outputing_bits();             // init encoder
  start_encoding();
  Ac_bits_total=0;                    // clear bit count
  Outfile = NULL;                     // no output
  No_of_chars = alpha_size;
  init_extra_symbols();		      // assign codes to RLE, EOF, and ESC  
  alloc_freq_tables();		      // allocate frequency tables 
  init_models();		      // initialize models

  for(i=0;i<n;) {
    next=t[i++];                      // first symbol of a run
    if(next!=0) {                     // nonzero symbol
      ac_symb_output(next);
      continue;
    }
    count=1;                          // initial length of the run of 0's
    while(i<n && t[i]==next) {
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
      if ((count & mask) == 0 )       // get b-th bits of count
	ac_symb_output(0);            // 0 stands for 0
      else    	       
	ac_symb_output(RLE_symbol);   // RLE_symbol stands for 1
      count <<= 1;                    // avoid dangerous right shift of mask
    }
  }
  // end of segment: encode an EOF_symbol
  ac_symb_output(EOF_symbol);
  free_freq_tables();
  done_encoding();
  return Ac_bits_total;
}
/* ********** as before, but without rle ************ */
static int ac_cost(symbol *t, int n, int alpha_size)
{
  int i; 

  init_update_parameters();
  start_outputing_bits();             // init encoder
  start_encoding();
  Ac_bits_total=0;                    // clear bit count
  Outfile = NULL;                     // no output
  No_of_chars = alpha_size;
  init_extra_symbols();		      // assign codes to RLE, EOF, and ESC  
  alloc_freq_tables();		      // allocate frequency tables 
  init_models();		      // initialize models

  for(i=0;i<n;i++) 
    ac_symb_output(t[i]);
  // end of segment: encode an EOF_symbol
  ac_symb_output(EOF_symbol);
  free_freq_tables();
  done_encoding();
  return Ac_bits_total;
}




/* ***********************************************************************
   bound the oputput size using the modified 0-order entropy
   return |s|H_0^*(s) + \mu log|\Sigma|, where \mu is the -x command 
   line parameter [if -x was not used default is mu=10.0 which is
   (experimentally) a reasonable value for ac]
   ********************************************************************** */
static double ac_h0star(stats *s, int alpha_size)
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


// initialize the struct for arithmetic coding
void init_ac_compressor(base_compressor *b)
{
  b->name = "Arithmetic Coding (no RLE)";
  b->bound = ac_h0star;
  b->cost = ac_cost;
  b->merge = NULL;
  b->encode = ac_compr_aux;
  b->decode = ac_decompr_aux;
  b->enc_start = ac_init_encoding;
  b->enc_stop = ac_flush_buffer; 
  b->dec_start = ac_init_decoding;
  b->dec_stop = NULL; 
}

// initialize the struct for rle+ac
void init_ac_rle_compressor(base_compressor *b)
{
  b->name = "Arithmetic Coding with RLE";
  b->bound = ac_h0star;
  b->cost = ac_rle_cost;
  b->merge = NULL;
  b->encode = ac_rle_compr_aux;
  b->decode = ac_rle_decompr_aux;
  b->enc_start = ac_init_encoding;
  b->enc_stop = ac_flush_buffer; 
  b->dec_start = ac_init_decoding;
  b->dec_stop = NULL; 
}

// initialize the struct for rle+ac
void init_ac_rle0_compressor(base_compressor *b)
{
  b->name = "Arithmetic Coding with RLE0";
  b->bound = ac_h0star;
  b->cost = ac_rle0_cost;
  b->merge = NULL;
  b->encode = ac_rle0_compr_aux;
  b->decode = ac_rle0_decompr_aux;
  b->enc_start = ac_init_encoding;
  b->enc_stop = ac_flush_buffer; 
  b->dec_start = ac_init_decoding;
  b->dec_stop = NULL; 
}


// ===========================================================
//  arithmetic coding  
// ===========================================================

/* DECLARATIONS USED FOR ARITHMETIC ENCODING AND DECODING */

/* ADAPTIVE SOURCE MODEL */
static int *freq;      	     //!< Symbol frequencies 
static int *cum_freq;  	     //!< Cumulative symbol frequencies
static int *freq_flat; 	     //!< Symbol frequencies flat model
static int *cum_freq_flat;   //!< Cumulative symbol frequencies flat model

/* SIZE OF ARITHMETIC CODE VALUES. */
#define Code_value_bits 16              //!< Number of bits in a code value
typedef long code_value;                //!< Type of an arithmetic code value
#define Top_value (((long)1<<Code_value_bits)-1)      //!< Largest code value

/* HALF AND QUARTER POINTS IN THE CODE VALUE RANGE. */
#define First_qtr (Top_value/4+1)       //!< Point after first quarter
#define Half	  (2*First_qtr)         //!< Point after first half
#define Third_qtr (3*First_qtr)         //!< Point after third quarter

/* THE BIT BUFFER. */
static int buffer;                     // Bits waiting to be input
static int bits_to_go;                 // Number of bits still in buffer
static int garbage_bits;               // Number of bits past end-of-file



// ====================================================
//    ARITHMETIC ENCODING ALGORITHM 
// ====================================================

/* CURRENT STATE OF THE ENCODING. */
static code_value low, high; //Limits of the current code region     
static long bits_to_follow;  //# of opposite bits to output after the next bit

/* **********************************************************************
   send a char to the arith. coder. 
   if ch is new: output the ESC_symbol, encode ch with the flat model, 
   and remove ch from the flat model. Note that we always update the 
   (non-flat) model for the occurrence of ch, but we do not update 
   it for the ESC_symbol (indeed, each time we generate an ESC_symbol
   it is less likely we will see another one in the future).  
   ********************************************************************** */
static void ac_symb_output(int ch)
{
  if(ch==EOF_symbol) {
    encode_symbol(ESC_symbol,cum_freq);    // Encode ESC symbol
    encode_symbol(ch,cum_freq_flat);       // Encode EOF
    return;                                // do not update
  }
  assert(ch>=0 && ch<No_of_chars);
  ch = ch+1;                               // translate to index (obselete!)

  if (freq[ch]==0) {
    encode_symbol(ESC_symbol,cum_freq);    // Encode ESC symbol
    encode_symbol(ch,cum_freq_flat);       // Encode new symbol
    update_flat_model(ch);		   // Update flat model
  }					      
  else 
    encode_symbol(ch,cum_freq);	           // Encode symbol normally
  update_model(ch);                        // Update the model.        
}


// **********************************************************************
// encode symbol with respect to the frequencies in cumu_freq[]
// **********************************************************************
static void encode_symbol(int symbol,int cumu_freq[])
{
  long range;                               /* Size of current code region */

  range = (long)(high-low)+1;
  high = low +                                	/* Narrow the code region   */
    (range*cumu_freq[symbol-1])/cumu_freq[0]-1; /* to that allotted to this */
  low = low +                                	/* symbol.                  */
    (range*cumu_freq[symbol])/cumu_freq[0];

  for (;;) {
    /* Loop to output bits.     */
    if (high<Half) {
      bit_plus_follow(0);                       /* Output 0 if in low half. */
    }
    else if (low>=Half) {                       /* Output 1 if in high half.*/
      bit_plus_follow(1);
      low -= Half;
      high -= Half;                             /* Subtract offset to top.  */
    }
    else if (low>=First_qtr                     /* Output an opposite bit   */
	     && high<Third_qtr) {               /* later if in middle half. */
      bits_to_follow += 1;
      low -= First_qtr;                         /* Subtract offset to middle*/
      high -= First_qtr;
    }
    else break;                                 /* Otherwise exit loop.     */
    low = 2*low;
    high = 2*high+1;                            /* Scale up code range.     */
  }
}

/*! INITIALIZE FOR BIT OUTPUT. */
static void start_outputing_bits( void )
{   buffer = 0;                                 /* Buffer is empty to start */
    bits_to_go= 8;                              /* with.                    */
}

/*! OUTPUT A BIT. */
inline void output_bit( int bit )
{   
  Ac_bits_total++;                          
  buffer >>= 1; if (bit) buffer |= 0x80;    // Put bit in top of buffer
  bits_to_go -= 1;
  if (bits_to_go==0) {                      // Output buffer if it is full
    if(Outfile)
      putc( (char) buffer, Outfile);        // output byte
    bits_to_go = 8;
  }
}

/*! OUTPUT BITS PLUS FOLLOWING OPPOSITE BITS. */
static void bit_plus_follow( int bit )
{   output_bit(bit);                            /* Output the bit.          */
    while (bits_to_follow>0) {
        output_bit(!bit);                       /* Output bits_to_follow    */
        bits_to_follow -= 1;                    /* opposite bits. Set       */
    }                                           /* bits_to_follow to zero.  */
}

/*! START ENCODING A STREAM OF SYMBOLS. */
static void start_encoding( void )
{   low = 0;                                    /* Full code range.         */
    high = Top_value;
    bits_to_follow = 0;                         /* No bits to follow next.  */
}

/*! FLUSH OUT THE LAST BITS. */
static void done_outputing_bits( void )
{ 
  if(Outfile)
    putc( (char) ( buffer >> bits_to_go ), Outfile);
}

/*! FINISH ENCODING THE STREAM. */
static void done_encoding( void )
{   bits_to_follow += 1;                        /* Output two bits that     */
    if (low<First_qtr) bit_plus_follow(0);      /* select the quarter that  */
    else bit_plus_follow(1);                    /* the current code range   */
}                                               /* contains.                */



// ====================================================
//    ARITHMETIC DECODING ALGORITHM 
// ====================================================

/*! INPUT A BIT. */
inline int input_bit( void )
{   int t;
    if (bits_to_go==0) {                         /* Read the next byte if no */
        buffer = getc(Infile);                   /* bits are left in buffer. */
        if (buffer==EOF) {
            garbage_bits += 1;                      /* Return arbitrary bits*/
            if (garbage_bits>Code_value_bits-2) {   /* after eof, but check */
                fprintf(stderr,"Bad input file\n"); /* for too many such.   */
                exit(-1);
            }
        }
        bits_to_go = 8;
    }
    t = buffer&1;                               /* Return the next bit from */
    buffer >>= 1;                               /* the bottom of the byte.  */
    bits_to_go -= 1;
    return t;
}

/* CURRENT STATE OF THE DECODING. */
static code_value value;        //!<  Currently-seen code value

/*! START DECODING A STREAM OF SYMBOLS. */
static void start_decoding( void )
{   int i;
    value = 0;                                  /* Input bits to fill the   */
    for (i = 1; i<=Code_value_bits; i++) {      /* code value.              */
        value = 2*value+input_bit();
    }
    low = 0;                                    /* Full code range.         */
    high = Top_value;
}

/*! get a symbol from the arith. decoder */
static int ac_symb_input(void)
{
  int ch;

  ch = decode_symbol(cum_freq);         // Decode next symbol 
  if(ch==ESC_symbol) {		        // ESC symbol: next symbol is new
    ch = decode_symbol (cum_freq_flat); // Decode new symbol
    if(ch==EOF_symbol) return ch;       // return  EOF, do not update models
    update_flat_model(ch);		// Update freq_flat and cum_freq_flat
  }
  update_model(ch);                     // update main model
  ch = ch-1;                            // Translate to a character (obsolete)
  assert(ch>=0 && ch<No_of_chars);      // check range      
  return ch;
}


/*! INITIALIZE BIT INPUT. */
static void start_inputing_bits( void )
{   bits_to_go = 0;                             /* Buffer starts out with   */
    garbage_bits = 0;                           /* no bits in it.           */
}


/*! DECODE THE NEXT SYMBOL. */
static int decode_symbol( int cumu_freq[] )
{   long range;                 /* Size of current code region              */
    int cum;                    /* Cumulative frequency calculated          */
    int symbol;                 /* Symbol decoded                           */
    range = (long)(high-low)+1;
    cum = (int)                                 /* Find cum freq for value. */
      ((((long)(value-low)+1)*cumu_freq[0]-1)/range);
    for (symbol = 1; cumu_freq[symbol]>cum; symbol++) ; /* Then find symbol. */
    high = low +                                /* Narrow the code region   */
      (range*cumu_freq[symbol-1])/cumu_freq[0]-1;/* to that allotted to this */
    low = low +                                 /* symbol.                  */
      (range*cumu_freq[symbol])/cumu_freq[0];
    for (;;) {                                  /* Loop to get rid of bits. */
        if (high<Half) {
            /* nothing */                       /* Expand low half.         */
        }
        else if (low>=Half) {                   /* Expand high half.        */
            value -= Half;
            low -= Half;                        /* Subtract offset to top.  */
            high -= Half;
        }
        else if (low>=First_qtr                 /* Expand middle half.      */
              && high<Third_qtr) {
            value -= First_qtr;
            low -= First_qtr;                   /* Subtract offset to middle*/
            high -= First_qtr;
        }
        else break;                             /* Otherwise exit loop.     */
        low = 2*low;
        high = 2*high+1;                        /* Scale up code range.     */
        value = 2*value+input_bit();            /* Move in next input bit.  */
    }
    return symbol;
}



// ===========================================================
//  auxiliary procedures 
// ===========================================================

/*! INITIALIZE THE MODELS. */
static void init_models(void)
{
  int i;

  // ------ initialize order 0 model 
  for (i = 0; i<No_of_symbols-1; i++) {         /* Set up initial frequency */
    freq[i] = 0;                                /* counts to be 0 for all   */
    cum_freq[i] = 1;			        /* symbols...               */
  }
  freq[ESC_symbol]=1;                           // ... except for ESC, which
  cum_freq[ESC_symbol]=0;                       // must be the last symbol

  // ------ inizialize flat model 
  for (i = 0; i<(No_of_symbols-1); i++) {       /* Set up initial frequency */
    freq_flat[i] = 1;                           /* counts to be one for all */
    cum_freq_flat[i] = (No_of_symbols - 1)-i-1; /* symbols except for ESC   */
  }                                             // which is not in this model 
  freq_flat[0] = 0;                             // <-- useless? 
}



// update the order 0 model to account for a new occ of symbol
static void update_model( int symbol )
{   int i;                                      /* New index for symbol     */
    if (cum_freq[0]+Increment>Max_frequency) {  /* See if frequency counts  */
        int cum;                                /* are at their maximum.    */
        cum = 0;
        for (i = (No_of_symbols-1); i>=0; i--) {/* If so, halve all the     */
            freq[i] = (freq[i]+1)/2;            /* counts (keeping them     */
            cum_freq[i] = cum;                  /* non-zero).               */
            cum += freq[i];
        }
    }
    i = symbol;
    freq[i] += Increment;                       /* Increment the frequency  */
    while (i>0) {                               /* count for the symbol and */
        i -= 1;                                 /* update the cumulative    */
        cum_freq[i] += Increment;               /* frequencies.             */
    }
}



// the flat model is used for the symbols not yet seen
// when symbol arrives we  assign to it probability 0.
static void update_flat_model( int symbol )
{
    int i;
    i = symbol;				    // zero the frequency	
    freq_flat[i] = 0;			    // count for the symbol and 
    while (i>0){			    // update the cumulative   
      i--;				    // frequencies
      cum_freq_flat[i]--;
    }
}


/* *******************************************************************
   Assign codes to the extra symbols/chars we use (RLE, EOF, ESC)
   NOTE: chars are the elements of the input file (that is what we
   want to encode and decode). symbols are chars + EOF and ESC.
   If there are A chars these are mapped into [1,A], then 
   A+1 is EOF and A+2 is ESC. No_of_symbols is A+3, but since 0
   is not used the actual number of symbols A+2.
   This stargeness should be removed in a future version
   ******************************************************************* */
static void init_extra_symbols(void)
{
  RLE_symbol = No_of_chars; // assign index to RLE_symbol
  No_of_chars ++;
  No_of_symbols = No_of_chars + 1;
  EOF_symbol = No_of_symbols;
  No_of_symbols++;
  ESC_symbol = No_of_symbols;
  No_of_symbols++;
}


static void alloc_freq_tables(void){  
  freq = malloc((No_of_symbols)*sizeof *freq);
  cum_freq = malloc((No_of_symbols)*sizeof *cum_freq);
  freq_flat = malloc((No_of_symbols - 1)*sizeof *freq_flat);
  cum_freq_flat = malloc((No_of_symbols - 1)*sizeof *cum_freq_flat);
}

static void free_freq_tables(void){  
  free(freq);
  free(cum_freq);
  free(freq_flat);
  free(cum_freq_flat);
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

static void init_update_parameters(void)
{
  // init order0 model parameters
  Max_frequency = _boost_CmdLinePar2_flag ? 
                 (int) _boost_CmdLinePar2 : default_Max_frequency;
  Increment = _boost_CmdLinePar3_flag ? 
                 (int) _boost_CmdLinePar3 : default_Increment;
  // check parameters
  if(Max_frequency>65535 || Increment>65535)
       fatal_error("Invalid update parameters! (init_update_parameters)\n");
}

static void fatal_error(char *s)
{
  fprintf(stderr,"%s",s);
  exit(1);
}

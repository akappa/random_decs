/*! \file Huffrle.c
  Procedures for the compression/decompression of sequences using Huffman
  coding optionally preceeded by RLE encoding.
*/ 
#include "huffrle.h"
#include <math.h>


#define MAX_ALPHA 258
#define NO_RLE 257
#define RLE_ALL 256


// ====== shared global variables =====================
int Huffrle_Verbose;  // controls verbosity level

// ========= "local" global variables ==============================
static int Zero, End_of_string;


// ========= compression procedures and data structures ============
typedef struct {
  int new_asize;              // size of the new alphabet
  int old_asize;              // size of the old alphabet
  int minlen;                 // min length of a codeword
  int maxlen;                 // max length of a codeword
  int code_of_len[MAX_ALPHA]; // # of codewords of each length
  int alpha_old[MAX_ALPHA];   // map new->old alphabet
  int alpha_new[MAX_ALPHA];   // map old->new alphabet
  int freq_old[MAX_ALPHA];    // frequency of each char (wrt to old alph.) 
  int freq_new[MAX_ALPHA];    // frequency of each char (wrt to new alph.)
  int code_val[MAX_ALPHA];    // values of code (wrt to new alphabet)
  int code_len[MAX_ALPHA];    // len of code (wrt to new alphabet)
} huff_data;



// ------------------- prototypes --------------------------------
static  int write_hufftree(huff_data *h, FILE *f);
static void stat_lengths(huff_data *h);
static void sort_by_codevalue(huff_data *h, int *sorted_symbols);
static void remap_alphabet(huff_data *h);
static void compute_text_frequencies(plain_text *t, int *freq, int rle_flag);
static void compute_freq_rle_all(plain_text *t, int *freq);
static void compute_freq_rle_single(plain_text *t, int *freq, int rlechar);
static  int write_encoded_text(plain_text *, huff_data *, int rle, FILE *);
static  int write_encoded_rle_single(plain_text *,huff_data *,int rle,FILE *);
static  int write_encoded_rle_all(plain_text *t, huff_data *h, FILE *f);
static void display_huffman_lengths(huff_data *h);
// ------------------------------------------------------------------


/*!
  Compress using Huffman coding after an optional RLE encoding
  the compressed data is written directly to file *f.
  \param t[]: the input string (size n)
  \param n:  size of input string
  \param asize: alphabet size (all symbols are in [0,asize-1]
  \param rleflag: rle procedure to use
  \param f: output file
  \return: numeber of bits used for the encoding  
*/
int huffrle_compr(plain_text *t, int rleflag, FILE *f) 
{
  huff_data h;
  int bits;   

  assert(t->asize<=256);
  assert(t->text!=NULL);

  // ----- define symbols Zero and End_of_string ---------------------
  if(rleflag!=NO_RLE) {
    Zero=t->asize;
    End_of_string = t->asize+1;
  }
  else
    End_of_string = t->asize;
  h.old_asize = End_of_string+1;   // End_of_string is the last symbol 

  // ----- compute frequencies (set freq_old[End_of_string]=1) -------
  compute_text_frequencies(t,h.freq_old,rleflag);

  // remap alphabet "removing" symbols with freq=0. Also remap frequencies 
  // so that in h->freq_new[0] .. h->freq_new[h->new_asize-1] 
  // there are a set of nonzero frequencies
  remap_alphabet(&h);

  // compute huffman codes based on frequencies. we force the frequency
  // of End_of_string to be 0 (rather than 1) to ensure that it gets
  // the largest codeword.
  h.freq_new[h.alpha_new[End_of_string]]=0;
  compute_huffman_code(h.freq_new,h.new_asize,h.code_len,h.code_val);
  h.freq_new[h.alpha_new[End_of_string]]=1;

  if(Huffrle_Verbose>3) display_huffman_lengths(&h);

  // --- write huffman tree to the output file 
  bits = write_hufftree(&h,f);

  // --- write encoded text to the output file 
  bits += write_encoded_text(t, &h, rleflag, f);
  // ---- note: if necessary fbit_flush must be invoked by the caller
  return bits;
}


/*!
   Procedure to write the huffman tree to file *f.
   The encoding of the huffman tree is based on canonical huffman
   coding with some ad-hoc additional tricks to save some extra bits.
   This makes the code difficult to read, but since in the optimal BWT
   compression we encode many different trees it is important to save
   as much as possible. 
   The encoding is done as follows. For each parameter we indicate its range 
   since the range determines the encoding (via fwrite_xxx_range()).
     maxlen: range [1,h->old_asize)
             this is the maximum length of a codeword
     minlen: range [1,maxlen]
             minimum length of a codeword
     # codewords of each length
             range: [0,#chars left]
             (maxlen-minlen) items overall

    
   \param h: huffman data
   \param f: output file
   \return: number of bits written
   */
int write_hufftree(huff_data *h, FILE *f)
{
  int i,chars_left,sorted_symbols[MAX_ALPHA],next,bits=0;

  stat_lengths(h);
  sort_by_codevalue(h,sorted_symbols);
  // -------- write minlen, maxlen  -------------
  bits += fwrite_bottom_range(f,h->maxlen,1,h->old_asize);
  bits += fwrite_bottom_range(f,h->minlen,1,h->maxlen+1);
  // -------- write # char of each len --------------------
  chars_left=h->old_asize;
  for(i=h->minlen;i<=h->maxlen;i++) {
    bits += fwrite_bottom_range(f,h->code_of_len[i],0,chars_left+1);
    chars_left -= h->code_of_len[i];
  }
  // --- write symbols in order of decreasing code value 
  // --- skip the last symbol which is End_of_string
  assert(h->alpha_old[sorted_symbols[h->new_asize-1]]==End_of_string);

  h->code_of_len[h->maxlen]--;  // EOS is not written
  next = h->new_asize-2;        // next char to be written
  for(i=h->maxlen;i>=h->minlen;i--) {
    int left,c,oldc;
    if((left=h->code_of_len[i])!=0) {
      c = h->alpha_old[sorted_symbols[next--]];
      left--;
      assert(c>=left);
      bits += fwrite_top_range(f,c,left,h->old_asize-1);
      while(left>0) {
        oldc = c;
	c = h->alpha_old[sorted_symbols[next--]];
	left--;
	assert(c<oldc && c>=left);
	bits += fwrite_top_range(f,c,left,oldc);
      }
    }
  }
  h->code_of_len[h->maxlen]++;  // restore correct value
  assert(next==-1);

  if(Huffrle_Verbose>3)
    fprintf(stderr,"Bits for representing the Huffman Tree: %d\n",bits);
  return bits;
}


/*!
    Procedure to compute the minimum and maximum lengths of a codeword
    and the number of codewords of a given length
*/
static void stat_lengths(huff_data *h)
{
  int i, j, n;

  n = h->new_asize;              //  number of encoded symbols
  for(i=0;i<MAX_ALPHA;i++) 
    h->code_of_len[i]=0;            // clear code_of_len[]
  h->minlen=MAX_ALPHA; 
  h->maxlen=0;
  for(i=0;i<n;i++) {
    j=h->code_len[i];
    assert(j>0);
    h->code_of_len[j]++;
    if(j<h->minlen) h->minlen=j;
    else if(j>h->maxlen) h->maxlen=j;
  }
  assert(h->minlen>0 && h->maxlen<MAX_ALPHA);
}

/*!
    Procedure to compute a list of symbols ordered by increasing codeword
*/
static void sort_by_codevalue(huff_data *h, int *sorted_symbols)
{
  int i, j, n, base[MAX_ALPHA];

  n=h->new_asize;
  // ------- init base --------------
  base[h->minlen]=0;
  for(i=h->minlen+1;i<=h->maxlen;i++)
    base[i] = base[i-1]+h->code_of_len[i-1];
  assert(base[h->maxlen]+h->code_of_len[h->maxlen]==n);
  // now base[i] contains the # of symbol of len<i

  // build array containing the symbols sorted in increasing code_value
  for(i=0;i<n;i++) {
    j=h->code_len[i];
    sorted_symbols[base[j]++]=i;
  }
  // ---- this is an extra check -----
  for(i=1;i<n;i++) 
    assert(h->code_val[sorted_symbols[i-1]]<h->code_val[sorted_symbols[i]]);
}


/*!
  Procedure to remap the alphabet so that in the new alphabet all
  symbols have freq!=0
  \param h: huffman data
*/
void remap_alphabet(huff_data *h)
{
  int j, i=0;

  for(j=0;j<h->old_asize;j++) 
    if(h->freq_old[j]!=0) {
      h->alpha_new[j]=i;
      h->alpha_old[i]=j;
      h->freq_new[i++]=h->freq_old[j];
    }
  h->new_asize=i;
}



// ------------------- frequencies computation -----------------

/*!
    Procedure to compute the frequecies of the text t with the
    possible use of rle
    \param t: input text
    \param freq[]: array of frequences
    \param rlechar: character whose runs are encoded using rle
*/ 
void compute_text_frequencies(plain_text *t, int *freq, int rleflag)
{
  int i,next;
 
  for(i=0;i<MAX_ALPHA;i++) freq[i]=0;
  if(rleflag==NO_RLE)
    for(i=0;i<t->tsize;i++) { 
      next = t->text[i];
      assert(next<=t->asize);
      freq[next]++; 
    }
  else if(rleflag<256 && rleflag >=0) 
    compute_freq_rle_single(t,freq,rleflag);
  else if(rleflag==RLE_ALL)
    compute_freq_rle_all(t,freq);
  else
    _huff_fatal_error("Invalid rle flag");
  assert(freq[End_of_string]==0);
  freq[End_of_string]=1;
}

/*!
    Procedure to compute the frequencies of a sequence t[0] ... t[n-1]
    which is encoded using rle encoding for each sequence of repeated symbols.
    A sequence of m>0 occ. of char x  is encoded writing m
    in binary, encoding 1's with x and 0's with Zero.
    Note that the encoding of a single char x is a single char x!
    \param t: input text
    \param freq[]: array of frequences
*/ 
void compute_freq_rle_all(plain_text *t, int *freq)
{
  int i, k, j, n, next;
  uint32 count;
  uint8 buffer[32];

  n=t->tsize;
  assert(n>0);
  for(i=0;i<n; ) {
    next=t->text[i]; 
    assert(next<t->asize);
    count=1;
    for(i=i+1;i<n;i++) 
      if(t->text[i]==next) count++;
      else break;
    // --- we have "count" occurrences of "next"
    assert(count>0);
    // fprintf(stderr,"%d count=%d\n",next,count);
    if(count==1) 
      freq[next]++;   // one occ: no rle
    else {
      k=_huff_uint32_to_binstring(count,buffer); // convert count to binary
      assert(buffer[k-1]==1);
      for(j=k-1;j>=0;j--)
	if(buffer[j]) freq[next]++;
	else          freq[Zero]++;
    }    
  }
} 

/*!
    Procedure to compute the frequencies of a sequence t[0] ... t[n-1]
    which is encoded using 1/2 rle encoding for the occurrences
    of rlechar. A sequence of m>0 rlechar's is encoded writing m+1
    in binary, discarding the most significat bit, and encoding 1's
    with Zero and 0's with rlechar. Note that the encoding of a single
    rlechar is a single rlechar, which makes our life easier!.
    \param t: input text
    \param freq[]: array of frequences
    \param rlechar: character whose runs are encoded using rle
*/ 
static void 
compute_freq_rle_single(plain_text *t, int *freq, int rlechar)
{
  int i, k, j, n, next;
  uint32 count;
  uint8 buffer[32];

  n=t->tsize;
  assert(n>0);
  for(i=0;i<n; ) {
    next=t->text[i++];
    assert(next<t->asize);
    if(next!=rlechar) {
      freq[next]++;
      continue;
    }
    // next is equal to rlechar: count occurrences
    count=1;
    while(i<n && t->text[i]==next) {
      count++; i++;
    }
    // we have count occ of next (==rlechar)
    if(count==1) 
      freq[next]++; // one occ
    else {
      k=_huff_uint32_to_binstring(count+1,buffer); // convert count to binary
      assert(k>1&&(buffer[k-1]==1));
      for(j=k-2;j>=0;j--)
	if(buffer[j]==0) freq[next]++;
	else             freq[Zero]++;
    }
  }
}

// ------------ output of compressed data ----------------- 

/*!
    Procedure to write to the output file the actual huffman code of
    the text t with the possible preprocessing of rle.
    \param t: input sequence
    \param h: huffman codes
    \param rleflag: rle procedure to use
    \param f: output file
*/ 
static int 
write_encoded_text(plain_text *t, huff_data *h, int rleflag, FILE *f)
{
  int i,j,n,bits=0;

  n=t->tsize;
  if(rleflag==NO_RLE) 
    for(i=0;i<n;i++) {       // no rle 
      j = h->alpha_new[t->text[i]];
      if(f!=NULL) fbit_write(f,h->code_len[j],h->code_val[j]);
      bits += h->code_len[j];
    } 
  else if(rleflag<256 && rleflag >=0) 
    bits += write_encoded_rle_single(t,h,rleflag,f);
  else if(rleflag==RLE_ALL)
    bits += write_encoded_rle_all(t,h,f);
  else
    _huff_fatal_error("Invalid rle flag");
  // --- write End_of_string -------------
  j = h->alpha_new[End_of_string];
  if(f!=NULL) fbit_write(f,h->code_len[j],h->code_val[j]);
  bits += h->code_len[j];
  return bits;
}


/*!
    Procedure to actually encode a sequence t[0] ... t[n-1] using rle 
    encoding for each sequence of repeated symbols followed by huffman coding
    \param t: input sequence
    \param h: huffman codes
    \param rleflag: rle procedure to use
    \param f: output file
    \return: # of written bits
*/ 
static int 
write_encoded_rle_all(plain_text *t, huff_data *h, FILE *f)
{
  int i, k, j, n, next, len, val, digit, bits=0;
  uint32 count;
  uint8 buffer[32];

  n = t->tsize;
  assert(n>0);
  for(i=0;i<n; ) {
    next=t->text[i]; 
    count=1;
    for(i=i+1;i<n;i++) 
      if(t->text[i]==next) count++;
      else break;
    // --- we have "count" occurrences of "next"
    assert(count>0);
    if(count==1) {
      len = h->code_len[h->alpha_new[next]];
      val = h->code_val[h->alpha_new[next]];
      if(f!=NULL) fbit_write(f,len,val);         // one occ
      bits += len;
    }
    else {
      k=_huff_uint32_to_binstring(count,buffer);  // convert count to binary
      assert(buffer[k-1]==1);
      for(j=k-1;j>=0;j--) {
	digit = buffer[j] ? next : Zero;
	len = h->code_len[h->alpha_new[digit]];
	val = h->code_val[h->alpha_new[digit]];
	if(f!=NULL) fbit_write(f,len,val);
	bits += len;
      }
    }
  }
  return bits;
} 

/*!
    Procedure to actually encode a sequence t[0] ... t[n-1] using 1/2 rle 
    for the occurrences of rlechar followed by huffman coding.
    \param t: input sequence
    \param h: huffman codes
    \param rleflag: rle procedure to use
    \param f: output file
    \return: # of written bits
*/ 
static int 
write_encoded_rle_single(plain_text *t, huff_data *h, int rlechar, FILE *f)
{
  int i, k, j, n, len, val, next, digit, bits=0;
  uint32 count;
  uint8 buffer[32];

  n = t->tsize;
  assert(n>0);
  for(i=0;i<n; ) {
    next=t->text[i++];
    if(next!=rlechar) {
      len = h->code_len[h->alpha_new[next]];
      val = h->code_val[h->alpha_new[next]];
      if(f!=NULL) fbit_write(f,len,val); 
      bits += len;
      continue;
    }
    // next is equal to rlechar: count occurrences
    count=1;
    while(i<n && t->text[i]==next) {
      count++; i++;
    }
    // we have count occ of next (==rlechar)
    if(count==1) {
      len = h->code_len[h->alpha_new[next]];
      val = h->code_val[h->alpha_new[next]];
      if(f!=NULL) fbit_write(f,len,val); 
      bits += len;
    }
    else {
      k=_huff_uint32_to_binstring(count+1,buffer); // convert count to binary
      assert(k>1&&(buffer[k-1]==1));
      for(j=k-2;j>=0;j--) {
	digit = (buffer[j]==0) ? next : Zero;
	len = h->code_len[h->alpha_new[digit]];
	val = h->code_val[h->alpha_new[digit]];
	if(f!=NULL) fbit_write(f,len,val);
	bits += len;
      }
    }
  }
  return bits;
}

/*!
    Display the huffman code length and the optimal code length
    for each symbol (just for debug purposes)
*/  
static void display_huffman_lengths(huff_data *h)
{
  int i,j,totfreq;
  double aux;

  totfreq=0;
  for(i=0;i<h->new_asize;i++)
    totfreq += h->freq_new[i];
  for(i=0;i<h->new_asize;i++) {
    j=h->alpha_old[i];
    fprintf(stderr,"%3d ",j);
    _huff_fpretty_putchar(stderr,j);
    fprintf(stderr,"  Occ: %8d",h->freq_new[i]);
    aux = (1.0 * h->freq_new[i])/totfreq;
    aux=-log(aux)/log(2.0);
    fprintf(stderr,"  OptLen: %7.4f",aux);
    fprintf(stderr,"  CodeLen: %2d",h->code_len[i]);
    fprintf(stderr,"  CodeVal %x\n",h->code_val[i]);
  }
  fprintf(stderr, "Number of symbols encoded: %d\n",totfreq);
}


// =====================================================================
//      decompression procedures and data structures 
// =====================================================================

typedef struct {
  int minlen;                    // minimum length of an huffman code
  int maxlen;                    // maximum length of an huffman code
  int old_asize;                 // size of the uncompressed alphabet
  int num_codewords;             // total number of codewords
  int limit[MAX_ALPHA];          // limit[i] is largest codeword of length i
  int max_rank[MAX_ALPHA];       // max rank of a codeword of length  i
  int sorted_symbols[MAX_ALPHA]; // symbols sorted in increasing code value 
} huff_tree;

// ==================== prototypes =======================================
void read_hufftree(huff_tree *h, FILE *f);
static int decode_text(plain_text *t, int rleflag, huff_tree *h, FILE *f);
static int decode_text_rle_single(plain_text *,int rleflag,huff_tree *,FILE *);
static int decode_text_rle_all(plain_text *t, huff_tree *h, FILE *f);
static int get_codeword(huff_tree *h, FILE *f);
static void display_hufftree_dec(huff_tree *h);
// =================================================================


int huffrle_decompr(plain_text *t, int rlechar, FILE *f)
{
  int written;
  huff_tree h;

  // ----- define symbols Zero and End_of_string ---------------------
  if(rlechar!=NO_RLE) {
    Zero=t->asize;
    End_of_string = t->asize+1;
  }
  else
    End_of_string = t->asize;
  h.old_asize = End_of_string+1;   // End_of_string is the last symbol 

  read_hufftree(&h,f);
  if(Huffrle_Verbose>3)
    display_hufftree_dec(&h);
  written = decode_text(t,rlechar,&h,f);
  return written;
}

/*!
    Procedure to read from f the data encoding the huffman tree,
    initializes all fields of the data structure h
*/
void read_hufftree(huff_tree *h, FILE *f)
{
  int i,n,chars_left,next;
  int code_of_len[MAX_ALPHA];

  // ---------- read maxlen minlen code_of_len[]---------
  h->maxlen = fread_bottom_range(f,1,h->old_asize);
  h->minlen = fread_bottom_range(f,1,h->maxlen+1);
  chars_left=h->old_asize;
  for(i=h->minlen; i<=h->maxlen; i++) {
    code_of_len[i] = fread_bottom_range(f,0,chars_left+1);
    chars_left -= code_of_len[i];
  }
  // ---------- init max_rank[] and num_codewords 
  for(i=h->minlen, n=0; i<=h->maxlen; i++) { 
    n += code_of_len[i];
    h->max_rank[i]=n-1;        // (# codewords of length <=i) - 1   
  }
  h->num_codewords=n;          // total number of codewords
  // ---- init limit[]: limit[i] is the largest code of a symbol of length i  
  h->limit[h->minlen]=code_of_len[h->minlen]-1;  
  for (i=h->minlen+1; i<=h->maxlen; i++) 
    h->limit[i]=((h->limit[i-1]+1)<<1)+code_of_len[i]-1;

  // ---------- read list of symbols in increasing code len---------
  h->sorted_symbols[n-1]=End_of_string; // EOS is not explicitly written
  code_of_len[h->maxlen]--;      // do not count EOS
  next = h->num_codewords-2;     // next char to be written
  for(i=h->maxlen;i>=h->minlen;i--) {
    int left,c,oldc;
    if((left=code_of_len[i])!=0) {
      left--;
      c = fread_top_range(f,left,h->old_asize-1);
      h->sorted_symbols[next--] = c;
      while(left>0) {
        oldc = c;
	left--;
	c = fread_top_range(f,left,oldc);
	h->sorted_symbols[next--] = c;
      }
    }
  }
  code_of_len[h->maxlen]++;   // restore correct value
  assert(next==-1);
}

static int decode_text(plain_text *t, int rleflag, huff_tree *h, FILE *f)
{
  int i,next;

  if(rleflag==NO_RLE) {
    for(i=0; ;i++) {
      next=get_codeword(h,f);
      if(next==End_of_string)
	return i;
      else if(i>=t->tsize) {
	_huff_fatal_error("Plain text buffer full! (decode_text)");
      } else {
	assert(next>=0 && next < t->asize);
	t->text[i] = (uint8) next;
	// fprintf(stdout,"%c",next);
      } 
    }
  }
  else if(rleflag>=0 && rleflag < 256)
    return decode_text_rle_single(t,rleflag,h,f);
  else if(rleflag==RLE_ALL)
    return decode_text_rle_all(t,h,f);
  else
    _huff_fatal_error("Invalid rleflag (decode_text)");
  return -1;
}

static 
int decode_text_rle_single(plain_text *t, int rlechar, huff_tree *h, FILE *f)
{
  int written,next,count,i;

  assert(rlechar>=0 && rlechar<t->asize);
  written=0;
  next=get_codeword(h,f);
  while(next!=End_of_string){
    if(next!=rlechar && next!=Zero) {    // no rle: decode a single symbol
      if(written>=t->tsize)
	_huff_fatal_error("Plain text buffer full! (decode_text_rle_single)");
      assert(next>=0 && next < t->asize);
      t->text[written++] = (uint8) next;      // write next to buffer
      next=get_codeword(h,f);                 // fetch next symbol
      continue;
    }
    // ---- rle encoding of a sequence of rlechar ----
    count=1;
    while(next==rlechar || next==Zero) {
      if(next==Zero)
	count = (count<<1) + 1;   // Zero stands for bit 1 
      else
	count = (count<<1);       // rlechar stands for bit 0
      next = get_codeword(h,f);
    }
    assert(count>1);
    // ---- there are count-1 occ of rlechar
    if(written+count-1 >t->tsize)
       _huff_fatal_error("Plain text buffer full! (decode_text_rle_single)");
    for(i=count-1;i>0;i--)
      t->text[written++] = (uint8) rlechar;
  }
  return written;
}

static int decode_text_rle_all(plain_text *t, huff_tree *h, FILE *f)
{
  int written,next,count,i,rlechar;

  written=0;
  next=get_codeword(h,f);
  while(next!=End_of_string){
    count=1;
    rlechar=next;
    next = get_codeword(h,f);
    while(next==rlechar || next==Zero) {
      if(next==Zero)
	count = (count<<1);        // Zero stands for bit 0 
      else
	count = (count<<1)+1;       // rlechar stands for bit 1
      next = get_codeword(h,f);
    }
    assert(rlechar>=0 && rlechar<t->asize);
    assert(count>=1);
    // ---- there are count occ of rlechar
    if(written+count >t->tsize)
       _huff_fatal_error("Plain text buffer full! (decode_text_rle_single)");
    for(i=count;i>0;i--)
      t->text[written++] = (uint8) rlechar;
  }
  return written;
}



/*!
    Procedure to read and decode a single codeword from file f
*/
static int get_codeword(huff_tree *h, FILE *f)
{
  int i, offset, codeword;

  i=h->minlen;
  codeword = fbit_read(f,i);          // read minlen bits
  while(codeword>h->limit[i]) {
    codeword = (codeword<<1) + fbit_read(f,1); // read one more bit
    i++;
  }
  assert(i<=h->maxlen);
  offset = h->max_rank[i] - (h->limit[i]-codeword);
  assert(offset>=0 && offset < h->num_codewords);
  return h->sorted_symbols[offset];
}

/*! 
  Display the huffman tree used for decompression
*/
static void display_hufftree_dec(huff_tree *h)
{
  int i,j,k,s,codew_of_len;

  i=0;
  for(j=h->minlen;j<=h->maxlen;j++) {
    // compute # of codewords of length j
    if(j==h->minlen)
      codew_of_len = j-1; 
    else
      codew_of_len = h->max_rank[j]-h->max_rank[j-1];
    for(k=0;k<codew_of_len;k++) {
      s = h->sorted_symbols[i++];
      fprintf(stderr,"Symb: %4d  ",s);
      _huff_fpretty_putchar(stderr,s);
      fprintf(stderr,"  CodeLen %3d",j);
      fprintf(stderr,"  CodeVal %x\n",h->limit[j]-codew_of_len+k+1);
    }
  }
}


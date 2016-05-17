/*   Definitions required for using the huffrle.a library
 */


// ---------- struct representing an uncompressed text ------------
typedef struct {
  unsigned char *text;   // array of char (size tsize)
  int tsize;             // size of text 
  int asize;             // alphabet size: chars are in [0,asize-1]
} plain_text;


// ========== prototypes ========================
int huffrle_compr(plain_text *t, int rleflag, FILE *f);
int huffrle_decompr(plain_text *t, int rleflag, FILE *f);
void init_bit_buffer(void);
void fbit_flush(FILE *);


// external variable for controlling the verbosity level 
extern int Huffrle_Verbose;

// flag for rle usage. a value x < 256 forces rle on char x only.
#define NO_RLE 257
#define RLE_ALL 256


/*! \file Huffrle.h
 */

/* ----------- use assertion if DEBUG!=0 ------------- */
#ifndef DEBUG
#define DEBUG 1   /* set DEBUG to 0 to remove assertions and extra checks */
#endif
#if !DEBUG
#define NDEBUG 1  /* do not compile assertions */
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


// --------- type definitions -------------------
typedef unsigned short uint16;
typedef unsigned char uint8;
typedef unsigned int uint32;
typedef unsigned char uchar;
typedef int int32;

// ---------- struct representing an uncompressed text ------------
typedef struct {
  uint8 *text;           // array of char (size tsize)
  int tsize;             // size of text 
  int asize;             // alphabet size: chars are in [0,asize-1]
} plain_text;


// ========== prototypes ========================

// ------------ huffrle.c ---------------
int huffrle_compr(plain_text *t, int rleflag, FILE *f);
int huffrle_decompr(plain_text *t, int rleflag, FILE *f);

// ------------ huffman_codes.c ---------------
void compute_huffman_code(int *freq, int n, int *code_len, int *code_value);

// ------------- misc.c ------------------
void _huff_fatal_error(char *s);
void _huff_out_of_mem(char *f);
int _huff_uint32_to_binstring(uint32 n, uint8 *buffer);
void _huff_fpretty_putchar(FILE *f, int c);

// ------------- bitio.c -----------------
int fwriteN_1xN(FILE *f, uint32 n, int N);
uint32 freadN_1xN(FILE *f, int N);
int fwrite_top_range(FILE *f, int32 t, int32 a, int32 b);
int fwrite_bottom_range(FILE *f, int32 t, int32 a, int32 b);
int32 fread_top_range(FILE *f, int32 a, int32 b);
int32 fread_bottom_range(FILE *f, int32 a, int32 b);
void uint32_write(FILE *f, uint32 u);
void init_bit_buffer(void);
// static? void fbit_write24(FILE *f, int n, uint32 v);
void fbit_write(FILE *f, int n, uint32 v);
void fbit_flush(FILE *f);
// static? uint32 fbit_read24(FILE *f, int n);
uint32 uint32_read(FILE *f);
uint32 fbit_read(FILE *f, int n);

// ============= macros, variables etc. for debug purposes ===============

// print to stderr an integer followed by newline 
#define showdn(x) fprintf(stderr,"%d\n",x)

// external variable for controlling the verbosity level 
extern int Huffrle_Verbose;



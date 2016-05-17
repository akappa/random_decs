/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   unboost.c 
   Ver 1.0   24-may-05 (extracted from booster.c)
   decompress a file produced by booster.c 

   Copyright (C) 2003 Giovanni Manzini (manzini@mfn.unipmn.it)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

   See COPYRIGHT file for further copyright information	   
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <assert.h>
#include <math.h>
#include "booster.h"
#include "compr_defs.h"

/* ===============================================================
   shared global variables: these are the onlt non-static global 
   variables. to avoid conflict they are all prefixes by _boost_
   =============================================================== */
FILE * _boost_Outfile;            // output file
FILE * _boost_Infile;             // input file
int _boost_Infile_size;           // input file size
int _boost_Verbose;               // verbosity level
double _boost_CmdLinePar1;        // command line parameter #1
int _boost_CmdLinePar1_flag;      // !=0 if parameter #1 was defined
double _boost_CmdLinePar2;        // command line parameter #2
int _boost_CmdLinePar2_flag;      // !=0 if parameter #2 was defined
double _boost_CmdLinePar3;        // command line parameter #3
int _boost_CmdLinePar3_flag;      // !=0 if parameter #3 was defined


// ============ prototypes =====================
static void decompress_file(void);
static void fatal_error(char *);
static void out_of_mem(char *f);
static size_t getTime ( void );
static int expand_alphabet(void);
static void open_files_decompression(char *infile_name, char *outfile_name);
static uint32 uint32_read(FILE *f);
static void mtf_decode(symbol *t, int n);

// ------------- local global variables ------------
static int Verbose;
static FILE *Infile;                         // input file;
static FILE *Outfile;                        // output file;
static int Alpha_old[ALPHA_SIZE];            // map new->old alphabet
static int Alpha_bitmap[ALPHA_BITMAP_SIZE];  // bitmap old->new alphabet


/* ***************************************************************
   Decompression of a booster-produced file 
   *************************************************************** */
int main(int argc, char *argv[])
{  
  double end, start;
  char *infile_name=NULL, *outfile_name=NULL;
  int c;
  extern char *optarg;
  extern int optind, opterr, optopt;

  /* ------------- read command line options ----------- */
  Verbose=0;                       // be quiet by default

  if(argc<2) {
    fprintf(stderr, "Usage:\n\t%s infile [-o outfile] [-v]\n",argv[0]);
    fprintf(stderr,"\t-o file        output file name (def. infile.y)\n");
    fprintf(stderr,"\t-v             verbose output\n\n");
    exit(0);
  }
  opterr = 0;
  while ((c=getopt(argc, argv, "o:v")) != -1) {
    switch (c)
    {
      case 'o':
        outfile_name = optarg; break;
      case 'v':
	Verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
    }
  }
  if(optind<argc)
     infile_name=argv[optind];
  // ---- echo command line ---------------
  if(Verbose>1) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }

  // ----- check input parameters--------------
  if(infile_name==NULL) fatal_error("The input file name is required\n");

  // ---- open files and init shared variables ------------
  open_files_decompression(infile_name,outfile_name);
  _boost_Infile = Infile;
  _boost_Verbose = Verbose;

  // ---- do the work ---------------------------
  start = getTime();
  decompress_file();
  end = getTime();

  if(Verbose>0)
    fprintf(stderr,"Elapsed time: %f seconds.\n", end-start);
  fclose(Infile);
  fclose(Outfile);
  return 0;
}



void decompress_file()
{
  int n,i,k,r0,written, occ[ALPHA_SIZE], algo_id, use_mtf,use_bwt;
  int num_seg, new_alpha_size, *rank_next;
  size_t start,end;
  bwt_data b;
  base_compressor *compr_array, *compressor;

  start = getTime();

  // --- get description of algorithm used for compression ----
  compr_array = init_compr_array();// init available compressors
  rewind(Infile); 
  algo_id = getc(Infile);
  // ---  check if BWT_FLAG is set
  if(algo_id & BWT_FLAG) {         
    use_bwt=1; algo_id -= BWT_FLAG; 
  } else
    use_bwt=0;
  // ---  check if MTF_FLAG is set  
  if(algo_id & MTF_FLAG) {         
    use_mtf=1; algo_id -= MTF_FLAG; 
  } else
    use_mtf=0;
  // --- get algorithm
  if(algo_id>=NUM_BASE_COMPRESSORS)
    fatal_error("Invalid algorithm id!\n");
  compressor = &compr_array[algo_id];
  if(Verbose>1)
    fprintf(stderr,"Algo_id=%d, Mtf=%d Bwt=%d\n",algo_id,use_mtf,use_bwt);

  // ---- read size, eof_pos alpha_bitmap and expand alphabet ------------
  b.size = n = uint32_read(Infile);
  if(use_bwt)
    b.eof_pos = uint32_read(Infile);  
  for(i=0;i<(ALPHA_SIZE+31)/32;i++)
    Alpha_bitmap[i] = uint32_read(Infile);
  new_alpha_size = expand_alphabet();
  if(Verbose>1)
    fprintf(stderr,"Alphabet size: %d\n",new_alpha_size);

  // ----- allocate space for BWT
  b.bwt = (uint8 *) malloc(n+1);
  if(b.bwt==NULL)
    out_of_mem("decompress_file");

  // ----- start decompression
  if(compressor->dec_start) 
    (compressor->dec_start)(); // start decoding

  // if bwt then do decompression segment by segment
  if(use_bwt) {
    for(num_seg=written=0;written<n+1;num_seg++) { 
      written += compressor->decode(b.bwt+written,n+1-written,new_alpha_size);
    }
    assert(written==n+1);
    if(Verbose>1)
      fprintf(stderr,"Number of segments: %d\n",num_seg);
  }
  else {
    // if plain order zero compression, then one call suffices
     written = compressor->decode(b.bwt,n,new_alpha_size);
     assert(written==n);
     b.bwt[n]=0;   // dummy symbol to make length n+1
  }

  // ---- end decompression  
  if(compressor->dec_stop) 
    (compressor->dec_stop)(); // stop decoding
  end = getTime();
  
  fprintf(stderr,"Time %lu μs\n", end - start);

  // --- undo mtf if necessary
  if(use_mtf) mtf_decode(b.bwt,n+1);

  // ------ remap bwt and compute freq
  for(i=0;i<ALPHA_SIZE;i++) occ[i]=0;
  for(i=0;i<=n;i++) {
    assert(b.bwt[i]<new_alpha_size);
    b.bwt[i] = Alpha_old[b.bwt[i]];
    if(i!=b.eof_pos)
    occ[b.bwt[i]]++;
  }

  // ----- retrieve text and write it to Outfile -----
  if(use_bwt) {
    rank_next = (int32 *) malloc((n+1)*sizeof(int32));
    if(rank_next==NULL)
      out_of_mem("decompress_file");
    r0 = _bw_bwt2ranknext(&b, occ, rank_next);
    // retrieve t using rank_next + bwt (see _bw_bwt2ranknext)
    k = r0; i=0;
    do {
      k = rank_next[k];
      putc(b.bwt[k],Outfile);
      i++;
    } while(k!=0);
    assert(i==b.size);
    free(rank_next);
  }
  else // no bwt was used: just write b.bwt[] to outilfe
    for(i=0;i<n;i++)
      putc(b.bwt[i],Outfile);

  free(b.bwt);
  free(compr_array);
}


/* *********************************************************************
   transform Alpha_bitmap[] into the mapping new_alpha -> old_alpha 
   this mapping is returned in the global array Alpha_old[]
   return the size of new_alpha
   ********************************************************************* */
static int expand_alphabet()
{
  int i,j,bit,asize=0;

  for(i=j=0;j<ALPHA_SIZE;j++) {
    bit = Alpha_bitmap[j/32] & (1u << (j%32)); 
    if(bit!=0) { 
      Alpha_old[i++] = j;
      asize++;
    }
  }
  return asize;
}


static void mtf_decode(symbol *t, int n)
{
  symbol next, mtf_list[ALPHA_SIZE];
  int i,j;
  
  for(j=0;j<ALPHA_SIZE;j++)
    mtf_list[j]=j;

  for(i=0;i<n;i++) {
    next = mtf_list[t[i]];
    for(j=t[i];j>0;j--)
      mtf_list[j] = mtf_list[j-1];
    mtf_list[0] = next;
    t[i]=next;
  }
}



/* ***************************************************************
   open input and output files; initialize Outfile and Infile
   *************************************************************** */
void open_files_decompression(char *infile_name, char *outfile_name)
{
  FILE *fopen(const char *path, const char *mode);

  /* ------ open input and output files ------ */
  if(infile_name==NULL) 
    fatal_error("Please provide the input file name (open_files)\n");
  else {
    Infile=fopen(infile_name, "rb"); // b is for binary: required by DOS
    if(Infile==NULL) {
      perror(infile_name); exit(1);
    }
  }
  // --- open output file --------------------------
  if(outfile_name==NULL) {
    /* add ".bst" for compression ".y" for decompression */
    outfile_name = (char *) malloc(strlen(infile_name)+3);
    outfile_name = strcpy(outfile_name, infile_name);
    outfile_name = strcat(outfile_name,".y");
  }
  Outfile = fopen(outfile_name,"wb"); // b is for binary: required by DOS
  if(Outfile==NULL) {
    perror(outfile_name);
    exit(1);
  }
}

static void fatal_error(char *s)
{
  fprintf(stderr,"%s",s);
  exit(1);
}

static void out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f);
  exit(1);
}

static size_t getTime ( void )
{
   size_t usertime,systime;
   struct rusage usage;

   getrusage ( RUSAGE_SELF, &usage );

   usertime = (double)usage.ru_utime.tv_sec * 1000000 + usage.ru_utime.tv_usec;

   systime = (double)usage.ru_stime.tv_sec * 1000000 + usage.ru_stime.tv_usec;

   return(usertime+systime);
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

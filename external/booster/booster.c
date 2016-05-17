/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   booster.c 
   Ver 1.0   17-may-04 (bwtopt.c)
   Ver 2.0   28-feb-05 (renamed booster.c)
   bwt optimal partitioning and compression using suffix tree visit

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
static void fatal_error(char *);
static void out_of_mem(char *f);
static double getTime ( void );
static int compact_alphabet(symbol *t, int n);
static void open_files(char *infile_name, char *outfile_name, int use_stdout);
static int *compute_bwt_lcp(bwt_data *b, int need_lcp);
static void uint32_write(FILE *f, uint32 u);
static void mtf_encode(symbol *t, int n);
static void compress_bwt(base_compressor *compressor, 
                          int algo_id, int use_mtf, int flat);
static void compress_plain(base_compressor *compressor, 
                          int algo_id, int use_mtf);
double bwt_partition(bwt_data *b, int *lcp, int asize, 
                            base_compressor *a, float *segment_cost);

// ------------- local global variables ------------
static int Verbose;
static FILE *Infile;                   // input file;
static FILE *Outfile;                  // output file;
static int Infile_size;                // size input file;
static int New_alpha_size;             // alphabet size of input file
static int Freq_new[ALPHA_SIZE];       // freq. of new alphabet symbols
static int Alpha_bitmap[ALPHA_BITMAP_SIZE];  // bitmap old->new alphabet


/* ***************************************************************
   File compression via the compression booster 
   *************************************************************** */
int main(int argc, char *argv[])
{  
  void check_compressor_usage(base_compressor *a,int);
  double end, start;
  char *infile_name=NULL, *outfile_name=NULL;
  int c, algo_id, use_mtf, flat, use_bwt, outfile_is_stdout;
  base_compressor *compr_array;
  compressor_usage usage;
  extern char *optarg;
  extern int optind, opterr, optopt;

  /* ------------- read command line options ----------- */
  _boost_CmdLinePar1_flag=0;       // flag for par #1 (to be set if defined) 
  _boost_CmdLinePar2_flag=0;       // flag for par #2 (to be set if defined) 
  _boost_CmdLinePar3_flag=0;       // flag for par #3 (to be set if defined) 
  Verbose=0;                       // be quiet by default
  algo_id = 7;                     // default range coding
  use_bwt = 1;                     // default use bwt
  use_mtf = 0;                     // default do not use mtf
  outfile_is_stdout=0;             // default do not send output to stdout
  flat = 0;                        // default partition the bwt 
  usage = True_cost;               // default use true cost
  compr_array = init_compr_array();// init available compressors

  if(argc<2) {
    fprintf(stderr, "Usage:\n\t%s infile ",argv[0]);
    fprintf(stderr,"[-fmnv][-a algo_id][-c cost]");
    fprintf(stderr,"[-x val][-y val][-z val][-o file] \n");
    fprintf(stderr,"\t-a algo_id     range: [0,%d] ",NUM_BASE_COMPRESSORS-1);
    fprintf(stderr,"(def. %d=%s)\n", algo_id, compr_array[algo_id].name);
    fprintf(stderr,"\t-c cost        0=Estimate, 1=True (def. True)\n");
    fprintf(stderr,"\t-f             compress flat BWT (do not partition)\n");
    fprintf(stderr,"\t-m             use move-to-front\n");
    fprintf(stderr,"\t-n             do not use BWT at all\n");
    fprintf(stderr,"\t-o file        output file name\n");
    fprintf(stderr,"\t-O             output to standard out\n");
    fprintf(stderr,"\t-v             verbose output\n");
    fprintf(stderr,"\t-x val         extra parameter #1\n");
    fprintf(stderr,"\t-y val         extra parameter #2\n");
    fprintf(stderr,"\t-z val         extra parameter #3\n\n");
    print_compr_array(stderr,compr_array);
    fprintf(stderr,"\n");
    exit(0);
  }
  opterr = 0;
  while ((c=getopt(argc, argv, "a:c:fmno:x:y:z:vO")) != -1) {
    switch (c)
    {
      case 'a':
	algo_id = atoi(optarg); break;
      case 'c':
	usage = (atoi(optarg)==0) ? Bound_only : True_cost; break;
      case 'f':
	flat = 1; break;
      case 'm':
	use_mtf = 1; break;
      case 'n':
	use_bwt = 0; break;
      case 'O':
        outfile_is_stdout = 1; break;
      case 'o':
        outfile_name = optarg; break;
      case 'x':
        _boost_CmdLinePar1 = atof(optarg);
	_boost_CmdLinePar1_flag =1; break;
      case 'y':
        _boost_CmdLinePar2 = atof(optarg);
	_boost_CmdLinePar2_flag =1; break;
      case 'z':
        _boost_CmdLinePar3 = atof(optarg);
	_boost_CmdLinePar3_flag =1; break;
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
  if(algo_id<0 || algo_id>=NUM_BASE_COMPRESSORS || algo_id>=VALID_ALGO_IDS)
    fatal_error("Invalid algorithm id!\n");
  if(usage>1)
    fatal_error("Invalid booster usage!\n");


  // ---- init and check base compressor
  compr_array[algo_id].usage = usage;
  check_compressor_usage(&compr_array[algo_id],flat);

  // ---- open files and init mandatory shared variables ------------
  open_files(infile_name,outfile_name,outfile_is_stdout);
  _boost_Outfile = Outfile;
  _boost_Infile = Infile;
  _boost_Infile_size = Infile_size;
  _boost_Verbose = Verbose;

  // ---- do the work ---------------------------
  start = getTime();
  if(use_bwt)
    compress_bwt(&compr_array[algo_id], algo_id, use_mtf, flat);
  else 
    compress_plain(&compr_array[algo_id], algo_id, use_mtf);
  end = getTime();

  if(Verbose>0)
    fprintf(stderr,"Elapsed time: %f seconds.\n", end-start);
  free(compr_array);
  fclose(Infile);
  fclose(Outfile);
  return 0;
}


/* ***************************************************************
      1. Read Infile remap it and and store it to text[]. 
      2. Optionally compute mtf of text
      3. Write to the output file:
           algo_id | MTF_FLAG | BWT_FLAG      (1 byte)
           input file size                    (4 bytes)
           bitmap for new alphabet            (4*ALPHA_BITMAP_SIZE bytes) 
         Followed by base_compressor() applied to text
   *************************************************************** */
static void compress_plain(base_compressor *compressor, 
                          int algo_id, int use_mtf)
{
  uint8 *text;
  int i,n,extra_bits,main_body=0; 

  // ----- allocate suffix array bwt text -----
  n = Infile_size;                               // length of input text
  text = (uint8 *) malloc(n*sizeof(*text));      // text
  if (! text) 
    out_of_mem("compress_plain");

  // ----- read text and remap it ------
  rewind(Infile); 
  i=fread(text, (size_t) 1, (size_t) n, Infile);
  if(i!=n) fatal_error("Error reading the input file!");
  if(Verbose>0)
    fprintf(stderr,"File size: %d bytes\n",n);
  New_alpha_size=compact_alphabet(text, n);

  // -- use mtf if required
  if(use_mtf) mtf_encode(text,n);

  // ------ write header information -----------------------
  assert(use_mtf==0 || use_mtf==1);
  assert(algo_id>=0 && algo_id<VALID_ALGO_IDS);
  fputc(MTF_FLAG*use_mtf+algo_id,Outfile);
  uint32_write(Outfile,n);
  for(i=0;i<ALPHA_BITMAP_SIZE;i++)
    uint32_write(Outfile,Alpha_bitmap[i]);
  extra_bits = 8*(1+4*(1+ALPHA_BITMAP_SIZE));

  // ----------- do the actual compression! ----------------
  if(compressor->enc_start) 
    extra_bits += compressor->enc_start();  // start encoding
  if(compressor->encode)
    main_body = compressor->encode(text,n, New_alpha_size);
  if(compressor->enc_stop) 
    extra_bits += compressor->enc_stop();   // stop encoding 

  // ---- final comments -----------------
  if(Verbose>0)
    fprintf(stderr, "Output file size: %d bytes\n",
	    (extra_bits+main_body+7)/8);  
  
  // ---- deallocate -------
  free(text);
}


/* ***************************************************************
      1. Read Infile remap it and and store it to text[]. 
      2. Compute the bwt and the lcp arrays (discard text)
      3. Optionally compute mtf of the bwt
      4. Optionally compute an optimal paritioning
      5. Write to the output file:
           algo_id|MTF_FLAG|BWT_FLAG       (1 byte)
           input file size                 (4 bytes)
           bwt eof_pos                     (4 bytes)
           bitmap for new alphabet         (4*ALPHA_BITMAP_SIZE bytes)
         followed by base_compressor() applied to the optimal partition
         or to the whole bwt 
   *************************************************************** */
static void compress_bwt(base_compressor *compressor, 
                          int algo_id, int use_mtf, int flat)
{
  int n, i, k, actual, tot_actual, extra_bits, *lcp=NULL;
  bwt_data b;
  double start, end, tot_est, opt;
  float  *segment_cost=NULL;

  // -- read text, remap alphabet computing New_alpha_size/Alpha_bitmap[] 
  // -- Finally, compute bwt and lcp array
  if(flat) {
    int *sa;
    sa = compute_bwt_lcp(&b,0);  // compute bwt and sa  
    free(sa);                    // we do not need the suffix array    
  } 
  else
    lcp = compute_bwt_lcp(&b,1); // compute bwt and lcp
  n = b.size;
  // to avoid dealing with the eof char we set bwt[eof_pos]=bwt[eof_pos-1]
  assert(b.eof_pos>0);   // the eof symbol cannot be in bwt[0]!
  b.bwt[b.eof_pos] = b.bwt[b.eof_pos-1]; 

  // -- use mtf if required
  if(use_mtf) mtf_encode(b.bwt,n+1);


  // --- for debugging purposes we may want to know the estimated cost of
  // --- each segment. If not comment the following line
#if 0
  if(!flat) {
    segment_cost = (float *) malloc((n+1)*sizeof(float));
    if(segment_cost==NULL)
      fprintf(stderr,"Warning! No enough memory for segment_cost[]");
  }
#endif

  // for debugging purposes: compute the cost without partitioning
  if(Verbose>1) {
    double c=0;

    fprintf(stderr,"Alphabet size: %d\n", New_alpha_size);
    if(compressor->cost && !flat) {
      c =  (compressor->cost)(b.bwt,n+1,New_alpha_size);
      fprintf(stderr,"No partitioning: %d bytes\n",(int) ceil(c/8));
    }
  }

  // ---- compute the optimal partition ---------
  start = getTime();
  if(flat) {
    opt = 0;
    // lcp[0]=n+1;   // one single segment
  } 
  else {
    if(Verbose>1) fprintf(stderr,"Boosting: *** %s ***\n",compressor->name);
    opt = bwt_partition(&b,lcp,New_alpha_size,compressor,segment_cost);
  }
  end=getTime();
  if(Verbose>0)
    fprintf(stderr,"Optimal partition computation: %.2f seconds\n",end-start);

  // ------ write header information -----------------------
  assert(use_mtf==0 || use_mtf==1);
  assert(algo_id>=0 && algo_id<VALID_ALGO_IDS);
  fputc(BWT_FLAG+(MTF_FLAG*use_mtf)+algo_id,Outfile);
  uint32_write(Outfile,n);
  uint32_write(Outfile,b.eof_pos);
  for(i=0;i<ALPHA_BITMAP_SIZE;i++)
    uint32_write(Outfile,Alpha_bitmap[i]);
  extra_bits = 8*(1+4*(2+ALPHA_BITMAP_SIZE));

  // ----------- do the actual compression! ----------------
  actual = tot_actual = 0; tot_est = 0;
  if(compressor->enc_start) 
    extra_bits += compressor->enc_start(); // start encoding

  if(flat) {                 // only one segment 
    k=1;                     
    if(compressor->encode)
      tot_actual = compressor->encode(b.bwt,n+1, New_alpha_size);
  }
  else {                     // use the partition encoded in lcp[]
    for(k=i=0;i<=n; ) { 
      assert(lcp[i]>i);   // bwt[i] -> bwt[lcp[i]-1] is a segment
      if(compressor->encode)
	actual = compressor->encode(b.bwt+i, lcp[i]-i, New_alpha_size);
      tot_actual += actual;
      
      if(Verbose>1 && segment_cost!=NULL) {
	tot_est += segment_cost[i];
	if(Verbose>2) 
	  fprintf(stderr,"Segment %d) [%d, %d]: real=%d est=%f bits\n",
		  k,i,lcp[i]-1,actual,segment_cost[i]);
      }
      i = lcp[i];      // starting point of next segment
      k++;             // increase # of segment 
    }
    assert(i==n+1);
  }
  if(compressor->enc_stop) 
    extra_bits += compressor->enc_stop(); // stop encoding 
  // ------------------------------------------------------


  // ---- final comments -----------------
  if(Verbose>0) {
    fprintf(stderr, "Output file size: %d bytes\n",
	    (extra_bits+tot_actual+7)/8);  
    fprintf(stderr, "  Number of segments: %d\n",k);
    fprintf(stderr, "  Actual compressed size: %d bytes\n",(tot_actual+7)/8);
    fprintf(stderr, "  Booster optimal estimate: %f bytes\n", opt/8);
    if(Verbose>1 && segment_cost!=NULL) 
      fprintf(stderr, "  Sum of segment costs %f bytes\n",tot_est/8);
  }
  // ---- deallocate -------
  if(segment_cost!=NULL) free(segment_cost);
  free(b.bwt);
  if(!flat)
    free(lcp);
}




/* ===================================================================
   functions for mapping and remapping of alphabets
   =================================================================== */

/* ************************************************************************
   remap the text chars so that they are integers in the range
   [0,new_alpha_size-1]. Returns the value new_alpha_size and Freq_new[] 
   Init Alpha_bitmap[] which translates from/to the new alphabet
   ************************************************************************ */
static int compact_alphabet(symbol *t, int n)
{
  int alpha_new[ALPHA_SIZE];
  int freq_old[ALPHA_SIZE];
  int i,j,new_alpha_size=0;
  
  // init freq old
  for(i=0;i<ALPHA_SIZE;i++)
    freq_old[i]=Freq_new[i]=0;
  // init alphabet bitmap
  for(i=0;i<ALPHA_BITMAP_SIZE;i++)
    Alpha_bitmap[i]=0;
  // compute freq old 
  for(i=0;i<n;i++) {
    // assert(t[i]>=0 && t[i]<ALPHA_SIZE);
    if(freq_old[t[i]]++==0) new_alpha_size++;
  }
  // remap alphabet
  for(i=j=0;j<ALPHA_SIZE;j++) 
    if(freq_old[j]!=0) {
      Freq_new[i]=freq_old[j];
      alpha_new[j]=i++;      
      Alpha_bitmap[j/32] |= (1u << (j%32));
    }
  assert(i==new_alpha_size);
  // remap text
  for(i=0;i<n;i++) {
    t[i]=alpha_new[t[i]];
    assert(t[i]<new_alpha_size);
  }
  assert(new_alpha_size<=ALPHA_SIZE);
  return new_alpha_size;  
}




/* ===================================================================
   auxiliary functions
   =================================================================== */


/* ******************************************************************
   Function for computing the bwt and lcp array of the text
   stored in Infile. The lcp array is computed only if need_lcp!=0;
   The bwt is stored in the bwt_data *b data structure, while the
   address of the lcp array is returned explicitly (if need_lcp==0 the 
   procedure returns the address of the sa).  
   The function also remap the text so that all its symbols are in the 
   range [0 .. New_alpha_size-1] and initialized the local global variables
   New_alpha_size, Freq_old[], Freq_new[],
   ****************************************************************** */
#define Plcp_skip 32
#define Large_file_size 160000000
static int *compute_bwt_lcp(bwt_data *b, int need_lcp)
{
  uint8 *text;
  int i,n,overshoot,extra_int32s, *sa;
  double start, end;

  // ----- init ds suffix sort routine -----
  overshoot=init_ds_ssort(500,2000);
  if(overshoot==0)
    fatal_error("ds initialization failed! (compress_file)\n");

  // ----- allocate suffix array bwt text -----
  n = Infile_size;                               // length of input text
  sa=malloc((n+1)*sizeof *sa);                   // suffix array
  b->bwt = (uint8 *) malloc(n+1);                // bwt
  text= (uint8 *) malloc(n+overshoot);           // text
  if (! sa || ! text || ! b->bwt) 
    out_of_mem("compute_bwt_lcp");

  // ----- read text and remap it ------
  rewind(Infile); 
  i=fread(text, (size_t) 1, (size_t) n, Infile);
  if(i!=n) fatal_error("Error reading the input file!");
  if(Verbose>0)
    fprintf(stderr,"File size: %d bytes\n",n);
  New_alpha_size=compact_alphabet(text, n);

  // ----- build suffix array ----------------
  start = getTime();
  ds_ssort(text,sa+1,n);                         // sort suffixes
  _bw_sa2bwt(text, n, sa, b);                    // compute bwt
  end=getTime();
  if(Verbose>0)
    fprintf(stderr,"Suffix array+bwt construction: %.2f seconds\n",end-start);

  // ---- compute bwt + lcp using 6n algorithm ---------
  if(need_lcp) {
    start = getTime();
    extra_int32s = _lcp_sa2lcp_6n(text,b,sa,Freq_new);  
    #if 0
    if(n<Large_file_size)    // for "small" files use lcp6 
      extra_int32s = _lcp_sa2lcp_6n(text,b,sa,Freq_new);  
    else                     // for "large" files use lcp5q
      extra_int32s =_lcp_sa2lcp_5n_q(text,n,sa,Plcp_skip);
    #endif
    end=getTime();

    if(Verbose>0) {
      fprintf(stderr,"lcp construction: %.2f seconds\n",end-start);
      if(Verbose>1) {
	if(n<Large_file_size)
	  fprintf(stderr,"Total memory for lcp6: %.2fn bytes\n",
		  6+(4.0*extra_int32s)/n);
	else fprintf(stderr,"Total memory for lcp5q: %.2fn bytes\n",
		  5+(4.0*extra_int32s)/n);
      }
    }
  }
  free(text); // text no longer needed
  return sa;  // the lcp array has been overwritten to sa[]
}

/* ***********************************************************
   encode (in place) the array t[0,n-1] using mtf
   *********************************************************** */
static void mtf_encode(symbol *t, int n)
{
  symbol old, tmp, mtf_list[ALPHA_SIZE];
  int i,j;
  
  for(j=0;j<ALPHA_SIZE;j++)         // init mtf list
    mtf_list[j]=j;

  for(i=0;i<n;i++) {
    old=mtf_list[0];                  // required for the case j=0
    for(j=0;j<ALPHA_SIZE;j++) {
      if(mtf_list[j]==t[i]) {
	t[i] = j;
	mtf_list[0] = mtf_list[j];  // move to front
        mtf_list[j] = old;
	break;
      }
      else {
        tmp=old; old=mtf_list[j]; mtf_list[j]=tmp;
      }
    }
    if(j>=ALPHA_SIZE)
      fatal_error("Symbol not in mtf list!\n");
  }
}


// this is a very simple procedure to check that what you asked is available
void check_compressor_usage(base_compressor *a, int flat)
{
  if(flat==0) {
    if(a->usage==Bound_only && a->bound==NULL) {
      fprintf(stderr,"Bound_only mode not supported!\n");
      exit(1);
    }
    if(a->usage==True_cost && a->cost==NULL) {
      fprintf(stderr,"True_cost mode not supported!\n");
      exit(1);
    }
    if(a->usage==Advanced && a->merge==NULL) {
      fprintf(stderr,"Advanced mode not supported!\n");
      exit(1);
    }
  }  
  if(a->encode==NULL)
    fprintf(stderr,"Warning actual compression not supported!\n");
  else if(a->decode==NULL)
    fprintf(stderr,"Decompression not supported!\n");
}



/* ***************************************************************
   open input and output files; initializes Outfile, Infile, Infile_size
   *************************************************************** */
static void open_files(char *infile_name, char *outfile_name, int use_stdout)
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
  // ---- store input file length to Input_size (only for compression)
  fseek(Infile,0,SEEK_END);
  Infile_size=ftell(Infile);
  if (Infile_size==-1) {
    perror(infile_name); exit(1);
  } 
  if (Infile_size==0) fatal_error("Input file empty (open_files)\n");

  // --- open output file --------------------------
  if(use_stdout)
    Outfile = stdout;
  else {
    if(outfile_name==NULL) {
      /* ceradte output file name adding the ".bst" suffix */
      outfile_name = (char *) malloc(strlen(infile_name)+5);
      outfile_name = strcpy(outfile_name, infile_name);
      outfile_name = strcat(outfile_name,".bst");
    }
    Outfile = fopen(outfile_name,"wb"); // b is for binary: required by DOS
    if(Outfile==NULL) {
      perror(outfile_name);
      exit(1);
    }
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

static double getTime ( void )
{
   double usertime,systime;
   struct rusage usage;

   getrusage ( RUSAGE_SELF, &usage );

   usertime = (double)usage.ru_utime.tv_sec +
     (double)usage.ru_utime.tv_usec / 1000000.0;

   systime = (double)usage.ru_stime.tv_sec +
     (double)usage.ru_stime.tv_usec / 1000000.0;

   return(usertime+systime);
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


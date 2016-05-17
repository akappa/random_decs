
#ifdef __USE_HUFFMAN__
#define NUM_BASE_COMPRESSORS 13
void init_Huffrle_compressor(base_compressor *);   // huffrle compressor 
void init_Huffrle0_compressor(base_compressor *);  // huffrle0 compressor 
void init_Huffnorle_compressor(base_compressor *b);// huff without rle
#else
#define NUM_BASE_COMPRESSORS 10
#endif

void init_Dummy_compressor(base_compressor *);    // dummy compressor
void init_ac_compressor(base_compressor *b);      // arith coding
void init_ac_rle_compressor(base_compressor *b);  // arith coding
void init_ac_rle0_compressor(base_compressor *b); // arith coding
void init_mh_compressor(base_compressor *);       // multiple table huffman
void init_mh0_compressor(base_compressor *);      // multiple table huffman0
void init_2pass_ac_compressor(base_compressor *b);// 2pass arith coding
void init_wv_rle_gamma_compr(base_compressor *b);       // wv+rle+gamma
void init_if_gamma_compressor(base_compressor *b);// if+gamma
void init_gdict_compressor(base_compressor *);    // gdict
void init_rc_compressor(base_compressor *);       // range coder plain
void init_rc_rle_compressor(base_compressor *);   // rle+range coder
void init_rc_rle0_compressor(base_compressor *);  // rle0+range coder

// +++++ add here init prototypes of other base compressors
// +++++ do not forget to update NUM_BASE_COMPRESSORS!!!!!


 
/* ***************************************************************
   Initialization of the known compressors 
   Returns an array of size NUM_BASE_COMPRESSORS containing the 
   definition of all known compressors
   *************************************************************** */ 
static void out_of_mem(char *s);
base_compressor *init_compr_array()
{
  base_compressor *ca;
  int i=0;

  // alloc compressor array
  ca = (base_compressor *) 
    malloc(NUM_BASE_COMPRESSORS*sizeof(base_compressor));
  if(ca==NULL)
    out_of_mem("init_compr_array");

  // init compressor array
  init_Dummy_compressor(&ca[i++]);
  init_ac_rle_compressor(&ca[i++]);
  init_ac_rle0_compressor(&ca[i++]);
  init_ac_compressor(&ca[i++]);
  init_wv_rle_gamma_compr(&ca[i++]);
  init_mh_compressor(&ca[i++]);
  init_mh0_compressor(&ca[i++]);
  init_rc_rle_compressor(&ca[i++]);
  init_rc_rle0_compressor(&ca[i++]);
  init_rc_compressor(&ca[i++]);
#ifdef __USE_HUFFMAN__ 
  init_Huffrle_compressor(&ca[i++]);
  init_Huffrle0_compressor(&ca[i++]);
  init_Huffnorle_compressor(&ca[i++]);
#endif

  // init_if_gamma_compressor(&ca[i++]);
  // init_2pass_ac_compressor(&ca[i++]);
  // init_gdict_compressor(&ca[i++]);
  // ++++ add here further initializations +++++++	  


  assert(i==NUM_BASE_COMPRESSORS);
  return ca;			  
}

void print_compr_array(FILE *f, base_compressor *ca)
{
  int i;

  fprintf(f,"Available compressors:\n");
  for(i=0;i<NUM_BASE_COMPRESSORS;i++) 
    fprintf(f,"  %2d.  %s\n",i,ca[i].name);
}

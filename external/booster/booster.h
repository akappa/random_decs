/* --- read proptotypes and typedef for ds_ssort --- */
#include "ds_ssort.h"
#include "bwt_aux.h"
#include "lcp_aux.h"

// constats related to alphabet size
#define ALPHA_SIZE _BW_ALPHA_SIZE
#define ALPHA_BITMAP_SIZE ((ALPHA_SIZE+31)/32)

// constants used to build the first byte of the output file
// which identifiesw the compression procedure used
#define BWT_FLAG 128       // 7th (last) bit set if bwt is used
#define MTF_FLAG 64        // 6th bit set if mtf is used
#define VALID_ALGO_IDS 64  // 6 bits available for algo_id

// integer types definitions
typedef unsigned char symbol;
typedef unsigned char UChar;
typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned int UInt32;
typedef int Int32;

// structure representing the occurrence statistics
// for a text fragment  
typedef struct {
  int *count;           // # occ of each char
  int count_tot;        // total # of chars
} stats;


// data structure representing a suffix tree node 
// and therefore a substring (fragment) of the BWT
typedef struct {
  int depth;            // depth in the suffix tree: ignored by compressors
  int start;            // fragment starting position (relative to *Text)
  int end;              // fragment ending position +1 (lenght is end-start)
  stats local_stats;    // statistic for the fragment
  double cost;          // cost (or estimate) of encoding the segment
  // int first_run;        // length of the first run of equal symbols
  // int last_run;         // length of the last run of equal symbols 
} fragment_node;


// enum representing the possible usages of a compressor 
// during the optimal partitioning process
// typedef enum {Bound_only, True_cost, Advanced} booster_usage;
typedef enum {Bound_only, True_cost, Advanced} compressor_usage;


/* ====================================================================
   Note on base_compressor. the functions bound(), cost(), and merge() 
   are called by the optimal partitioning procedure (when booster
   usage is respectively Bound_only, True_cost, Advanced).
   The compression consists in calling enc_start() followed by 
   encode() on each segment, followed by enc_stop(). The encoder 
   must write its output to _boost_Outfile.
   The decompression consists in calling dec_start() followed by 
   decode() on each segment, followed by dec_stop(). Note that each
   segment read by decode() should be "self-delimited" and that decode()
   is called on the next segment untill the right number of symbols
   have been read. The decoder reads its input from _boost_Infile.
   ==================================================================== */

// data structure representing a generic base compressor
typedef struct {
  char *name;                            
     // compressor name
  compressor_usage usage;   
     // how the compressor is used in the optimal partitioning phase  
  double (*bound)(stats *, int asize);        
     // estimate of fragment cost based on occurrence statistics
  int (*cost)(symbol *t, int size, int asize); 
     // cost of compressing t[0,size-1] whith 0<=0 t[i]<asize
  int (*merge)(fragment_node *f, int t, int alpha_size);
     // cost of compressing together fragments f[0], f[1], ... , f[t-1]
  int (*encode)(symbol *t, int size, int alpha_size);  
     // function for actually encoding t[0,size-1] whith 0<=0 t[i]<asize
     // should return # of bits used for encoding 
  int (*decode)(symbol *buffer, int size, int alpha_size);
     // function decoding a fragment. MUST return # of bytes decoded
  int (*enc_start)(void);        // start encoding. return # of bits produced
  int (*enc_stop)(void);         // stop encoding. return # of bits produced
  int (*dec_start)(void);        // start decoding. return value is ignored
  int (*dec_stop)(void);         // stop decoding. return value is ignored
} base_compressor;

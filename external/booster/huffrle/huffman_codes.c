/*! \file Huffman_codes.c
    Procedures for computing the Huffman codes given the (nonzero)
    frequencies of a set of symbols
*/
#include <stdlib.h>
#include <stdio.h>
#include "huffrle.h"


static void build_heap(int *heap, int *freq, int n);
static void heapify(int i, int size, int *heap);
static void build_huffman_tree(int *heap, int n);
static void compute_code_values(int *code_val, int *code_len, int n);
static void compute_code_len(int *len, int *heap, int n);



/*!
    compute the huffman code
    given a set of nonzero frequencies this procedure computes
    the code of each symbol. 
    /param freq[] frequency of each symbol (size n)
    /param n: number of symbols
    /param code_len[] code length of each symbol (size n, to be computed) 
    /param code_value[] numerical value of each symbol (size n, to be computed)
*/
void compute_huffman_code(int *freq, int n, int *code_len, int *code_value)
{
  int i, *heap;

  // start by checking that all symbols have freq!=0
  for(i=0;i<n;i++)
    assert(freq[i]>=0);

  if(n>2) {
    heap = (int *) malloc((2*n+1)*sizeof(int));
    if(heap==NULL) _huff_out_of_mem("compute_huffman_code");
    build_heap(heap,freq, n);             // build huffman heap
    build_huffman_tree(heap, n);          // build tree using the huffman heap
    compute_code_len(code_len, heap, n);  // compute code length of each symbol
    compute_code_values(code_value,code_len,n);
    free(heap);     
  }
  else if(n==2) {
    code_len[0]=code_len[1]=1;
    code_value[0]=0; code_value[1]=1;
  }
  else if(n==1)
    code_len[0]=code_value[0]=0;
  else
    _huff_fatal_error("Illegal number of symbols! (compute_huffman_code)");
}

/*!  
   Procedure for building an huffman heap
   See D. Salomon, Data Compression 2nd edition Sec. 2.5.8
   The huffman heap h[] has length 2n+1 where n is the number of symbols
   Initially, for i=1, ... n h[i] contains the value n+i and
   h[n+i] contains the frequency of symbol i-1  (h[0] is not used)
   Then an heap is build on h[1] ... h[n] using as the value of h[i] the 
   frequency h[h[i].
   \param freq[]: frequency of each symbol (size n)
   \param n: number of symbols with frequency>0 (length of symb)
   \return heap: huffman heap (size 2n+1) 
*/
static
void build_heap(int *heap, int *freq, int n)
{
  int i;

  for (i=1; i<=n; i++) {        
       heap[i]=(n+i);            // heap[i] is a pointer to the frequency 
                                 // of symbol i-1
       heap[heap[i]]=freq[i-1];  // actual frequency of symbol i-1
  }
  for (i=n/2; i>0; i--) {        // build heap 
    heapify(i, n, heap);   
  }
}

/*!  
   heapify procedure for an huffman heap
   \param i: node to be sifted 
   \param size: heap size
   \param heap: array housing the heap
*/
static 
void heapify(int i, int size, int *heap)
{
  int left, right, smallest, tmp;

  if (size>1) {   
    left=2*i;     // left child 
    right=left+1;        // right child

    smallest=i;
    if (left<=size)
      smallest = (heap[heap[left]]<heap[heap[i]]) ? left : smallest;
    if (right<=size)
      smallest = (heap[heap[right]]<heap[heap[smallest]]) ? right : smallest;
    if (smallest!=i) {          // if a smallest child was found move i down 
      tmp=heap[i];
      heap[i]=heap[smallest];
      heap[smallest]=tmp;
      heapify(smallest, size, heap);     // this could be done iteratively
    }
  }
}

/*!  
   Build huffman tree using an huffman heap.
   See D. Salomon, Data Compression 2nd edition Sec. 2.5.8
   \param heap: huffman heap (has size 2n+1)
   \param n: number of symbols
*/
static
void build_huffman_tree(int *heap, int n)
{
  int min, i;

  for (i=n; i>=2; i--) {          // stop when there is only one node left
    min=heap[heap[1]];            // minimum frequency 
    heap[heap[1]]=i;              // replace frequency with pointer 
                                  // (heap[i] will be initialized later)
    heap[1]=heap[i];              // the last node becomes the root
    heapify(1,i-1, heap);         // restore heap heap[1,.., i-1] 
    heap[i]=(min+heap[heap[1]]);  // store in heap[i] the sum of the frequency 
                                  // of the two nodes which are removed
    heap[heap[1]]=i;              // replace frequency with pointer to the 
                                  // new frequency
    heap[1]=i;                    // insert in the root a new node with 
                                  // frequency equal to the sum of the 
                                  // freq. of the removed nodes
    heapify(1,i-1,heap);          // restore heap[i,..,i-1] 
  }
}


/*!
    This procedure compute the numerical value of each code starting
    from the code lengths.
    Symbols with the same code length have consecutive numerical values
    with smaller values given to the smaller symbols.
    \param code_val[]: numerical value of code (size n, to be computed)
    \param code_len[]: length of codeword for symbol i (size n)
    \param n: number of symbols
*/ 
static
void compute_code_values(int *code_val, int *code_len, int n)
{
  int i,j,len,minlen,maxlen;
  int *code_of_len, *limit;  

  assert(n>0);
  // these arrays must be of size > of the max code_length
  code_of_len = calloc(n+1,sizeof(int));
  limit = calloc(n+1,sizeof(int));
  if(!limit || !code_of_len)
    _huff_out_of_mem("compute_code_values");
  // ---------- compute minlen, maxlen and code_of_len ----------
  minlen=maxlen=code_len[0];
  for(i=0;i<n;i++) {
    j=code_len[i];                   // length of symbol i
    code_of_len[j]++;                // one more code of length j 
    if(j>maxlen) maxlen=j;           // update minlen and maxlen
    else if(j<minlen) minlen=j;
  }
  // ---- init limit: limit[i] is the largest code of a symbol of length i  
  limit[minlen]=code_of_len[minlen]-1;  
  for (i=minlen+1; i<=maxlen; i++) 
    limit[i]=((limit[i-1]+1)<<1)+code_of_len[i]-1;
  // ------------- compute the actual codes ------------  
  for (i=0; i<n; i++) {          // rescan al symbols 
    len=code_len[i];             // len of symbol i
    code_val[i] = limit[len]-(--code_of_len[len]);
    // printf("Symb %4d, CodeLen %4d CodeVal %x\n",i,len,code_val[i]); // !!!
  }
  for(i=0;i<=maxlen;i++)
    assert(code_of_len[i]==0);
  free(limit);
  free(code_of_len);
}


/*! 
    compute the code length of each symbol using the huffman heap
    \param len[]: array where length are written (size=n)
    \param heap[]: huffman heap (now containing only one node) (size 2n+1)
    \param n: number of symbols     (size n)
*/
static void compute_code_len(int *len, int *heap, int n)
{
  int i, j, k;

  for (i=2*n; i>n; i--) {
    k=1;                               // length of symbol i-n-1
    for (j=heap[i]; j!=2; j=heap[j])   // follow pointers
      k++;
    len[i-n-1]=k;                      // write len
  }
}


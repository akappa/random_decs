/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   st_visit.c 
   Ver 1.0   17-may-04 (bwtopt1.c)
   Ver 2.0   22-feb-05 (renamed st_visit.c)
   bwt optimal partitioning using suffix tree visit

   Copyright (C) 2003-05 Giovanni Manzini (manzini@mfn.unipmn.it)

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


// set EXTRA_CHEKS to 1 for additional (time consuming) checks
#ifndef NDEBUG
#define EXTRA_CHECKS 0  
#endif

extern int _boost_Verbose;


// ------------- prototypes ---------------
void reduce_stack(int d);
void merge_top_elements(int s, int e);
double compression_cost(fragment_node *s);
// stack related
void push_node(int sa_pos, int lcp, int c);
void enlarge_stack(void);
void init_empty_stack(void);
int free_stack(void);
// statistics 
void init_stats(stats *s, int c);
void merge_stats(stats *s1, stats *s2);
int free_stats(stats *t);
// debugging
void validate_stats(symbol *t, int size,stats *w);
void check_subtree(fragment_node *q);
double pseudo_compression(symbol *t, int size);
void check_opt_partition_cost(bwt_data *b, double opt_estimate);
static void out_of_mem(char *f);
// double comulative_cost(stats *s);



// ------------ local global variables ---------------
static fragment_node *Stack;     // stack containing suffix tree nodes 
static int Top;                  // pointer to the top of the stack
static int Stack_size;           // current size of the stack
static int *End_point;           // endpoints of segments
static float *Segment_cost;      // estimated cost of segments
static symbol *Bwt;              // bwt array
static int Alpha_size;           // # distinct symbols in the BWT
static base_compressor *Algo;    // base compression algorithm
static int Verbose;              // verbosity level


/* ********************************************************************
   This function computes the optiml partition of the BWT with respect
   to the cost function  a->bound() or a->cost() (depending on
   the value of a->usage).
   The partition is computed by means of a post-order suffix tree visit 
   emulated via a single pass over the BWT+LCP array. 
   The optimal partition of the BWT is written to the LCP array and is 
   given by:
        [0, t0-1] [t0, t1-1] [t1, t2-1] .... [tk,n]
   where t0 = lcp[0], t{i+1} = lcp[ti] and lcp[tk]=n+1
   If *segment_cost is !=NULL then in segment_cost[ti] is stored the 
   estimated cost of encoding the segment  [ti,t{i+1}-1]
   ****************************************************************** */
double bwt_partition(bwt_data *b, int *lcp, int asize, base_compressor *a, 
                     float *segment_cost)
{
  void init_global_vars(bwt_data *, int *, int, base_compressor *, float *);
  int i,depth,last_depth,stack_mem_usage;
  double optimal_cost;

  // initi global variables and stack used for suffix tree visit
  init_global_vars(b, lcp, asize, a, segment_cost);
  init_empty_stack();                    
  // insert leaves corresponding to # symbol (eos) and first suffix
  push_node(0, -1, b->bwt[0]);
  push_node(1, 0, b->bwt[1]);
  last_depth = 0;
  // insert all suffixes 
  for(i=2;i<=b->size;i++) {
    depth = lcp[i];
    if(depth<last_depth) {
      reduce_stack(depth);
    }
    push_node(i, depth, b->bwt[i]);
    last_depth=depth;
  }
  reduce_stack(-1);     // reduce stack to suffx tree root 
  assert(Top==0);
  optimal_cost =  Stack[0].cost;
  stack_mem_usage = free_stack();
  if(Verbose>1) 
    fprintf(stderr,"Memory for suffix-tree stack: %d bytes\n",stack_mem_usage);
#if EXTRA_CHECKS
  check_opt_partition_cost(b,optimal_cost);
#endif 
  return optimal_cost;
}



// remove from stack all elements with depth > d
void reduce_stack(int d)
{
  int i,curr_depth;

  // this function should be called only if there is something to do
  assert(Top>=0);  
  assert(Stack[Top].depth>d);

  while(Stack[Top].depth>d) {
    curr_depth=Stack[Top].depth;
    // find first item with depth<d 
    for(i=Top-1;i>0;i--) {
      assert(Stack[i].depth<=curr_depth);
      if(Stack[i].depth<curr_depth)
	break;
    }
    assert(i>0 || d==-1);
    // find optimal encoding for the subtrees rooted at 
    //   Stack[i], Stack[i+1], ..., Stack[Top]
    // and create in Stack[i] a node for the parent-subtree
    merge_top_elements(i,Top);
    Top=i;
    assert(Top>=0);
  }
  assert(Top>0 || d==-1);
}


/* **********************************************************************
   this function check if the leaves of the subtrees rooted at
      Stack[s], Stack[s+1], ... , Stack[e]
   should be compressed together, or not.  
   ********************************************************************** */
void merge_top_elements(int s, int e)
{
  void merge_fragments(fragment_node *s,fragment_node *t);
  int i;
  double merged_cost;

#ifndef NDEBUG
  // --- check consistency of subtrees
  for(i=s;i<=e;i++) 
    check_subtree(&Stack[i]);
#endif

  // create in Stack[s] the union of the segments
  for(i=s+1;i<=e;i++) 
    merge_fragments(&Stack[s], &Stack[i]);

  // Stack[s].cost contains the "splitted" cost now compute "merged" cost 
  merged_cost = compression_cost(&Stack[s]);
  //!!! merged_cost = comulative_cost(&Stack[s].local_stats); // !!!!!!!!!
  // fprintf(stdout,"** [%d %d] merged: %f splitted %f\n", 
  //	 Stack[s].start, Stack[s].end, merged_cost, Stack[s].cost);
  
  // choose whether to merge segments 
  if(merged_cost<=Stack[s].cost) {
    // do the merging
    Stack[s].cost = merged_cost;
    if(Segment_cost!=NULL)
      Segment_cost[Stack[s].start]=merged_cost;
    End_point[Stack[s].start]=Stack[s].end; 
  }
}


#if 0
// ****************************************************************
// questo frammento sara' eventualmente incorporato nella gestione
// della modalita' advanced per RHC
// ****************************************************************
  // --------- update s->first_run  s->last_run
  if(s->last_run==0) {
    // s is a sequence of equal symbols
    assert(s->first_run==s->end-s->start);
    if(Bwt[s->end-1]==Bwt[t->start]) {
      s->first_run += t->first_run;
      s->last_run = t->last_run;
    }
    else 
      s->last_run = (t->last_run==0) ? t->first_run : t->last_run;
  }
  else if(t->last_run==0) {
    // t is a sequence of equal symbols (but s is not)
    assert(t->first_run==t->end-t->start);
    if(Bwt[s->end-1]==Bwt[t->start]) 
      s->last_run += t->first_run;
    else
      s->last_run = t->first_run;
  }
  else {
    // both s and t contain mor than one symbol
    s->last_run = t->last_run;
  }
#endif



// merge two consecutive fragments
void merge_fragments(fragment_node *s,fragment_node *t)
{
  assert(s->end==t->start);

  // --- code for Usage==Advanced to be inserted here

  // ------- update stats, endpoint and total cost ------------
  merge_stats(&(s->local_stats),&(t->local_stats));
  s->end = t->end;
  s->cost += t->cost;
}


/* *****************************************************************
   compute the true or approximated cost (depending on Algo->usage) 
   for encoding the fragment *s with algorithm Algo
   ***************************************************************** */ 
double compression_cost(fragment_node *s)
{

  if(Algo->usage==Bound_only) {
    assert(Algo->bound);
    return (Algo->bound)(&(s->local_stats),Alpha_size);
  }

  assert(Algo->cost!=NULL);
  return (Algo->cost)(Bwt+s->start,s->end-s->start,Alpha_size);

  // the code for usage==Advanced should be included here!
}



/* ================================================================
   Debugging functions
   ================================================================ */

/* **********************************************************************
   verify that the stats in *w corresponds to the actual frequencies 
   measured in t[0..size-1]
   ********************************************************************** */
void validate_stats(symbol *t, int size,stats *w)
{
  int occ[ALPHA_SIZE];
  int i;

  // check that total # of occ equals text size
  if(size!=w->count_tot)
    fprintf(stderr,
            "Wrong size: %d vs %d (validate_stats)\n",size,w->count_tot);   
  // fill occ[] with # of occ in t[0] ... t[size-1] 
  for(i=0;i<Alpha_size;i++) occ[i]=0;
  for(i=0;i<size;i++) {
    assert(t[i]<Alpha_size);
    occ[t[i]]++;
  }
  // check that # occ of each char matches
  for(i=0;i<Alpha_size;i++)
    if(occ[i]!=w->count[i])
      fprintf(stderr,
	      "Mismatch: %d vs %d (validate_stats)\n",occ[i],w->count[i]); 
  return;
}

/* *****************************************************************
   Check that *q correctly represents a leaf cover of a subtree 
   of the suffix tree. *q represents the subtree starting 
   whose first leaf is q->start and whose last leaf is q->end-1.
   The global array End_point[] represents an optimal leaf cover 
   of the subtree w.r.t to the given cost function. We simply check that
   End_point[] defines a partition of the subtree leaves into contiguous
   intervals.
   ****************************************************************** */
void check_subtree(fragment_node *q)
{
  int ep,sp = q->start;    // first leaf in the subtree

  ep=sp;
  while(ep<q->end) {
    // sum += pseudo_compr(Bwt+ep,End_point[ep]-ep);
    assert(End_point[ep]>ep);
    ep = End_point[ep];        // first leaf in the next segment
  }
  assert(ep==q->end);
  assert(ep-sp == (q->local_stats).count_tot);
  // uncomment to check that the stats in q->lstats are correct
  #if EXTRA_CHECKS
  validate_stats(Bwt+sp,ep-sp,&q->local_stats);
  #endif
}



/* ******************************************************************
   reports the true/estimated size of the output of Algo
   with input t[0...size-1].
   this function is used only for debugging purposes and should be
   kept in sync with compression_cost().
   ****************************************************************** */
double pseudo_compression(symbol *t, int size)
{
  double cost;
  stats s; int i;

  if(Algo->usage==Bound_only && Algo->bound!=NULL) {
    s.count = (int *) calloc(Alpha_size,sizeof(int));
    if(s.count==NULL) 
      out_of_mem("pseudo_compression");
    // fill s.count with # of occ in t[0] ... t[size-1] 
    for(i=0;i<size;i++) {
      assert(t[i]<Alpha_size);
      s.count[t[i]]++;
    }
    s.count_tot = size;
    cost = (Algo->bound)(&s,Alpha_size);
    free(s.count);
  }
  else if(Algo->cost!=NULL)
    cost =  (Algo->cost)(t,size,Alpha_size);
  else
    cost = -1;

  return cost;
}


/* ****************************************************************
   check that the cost of each segment in the optimal partition 
   coincides with the cost obtained applying pseudo_compr()
   and that the optimal cost is equal to the sum of the segment costs
   **************************************************************** */
void check_opt_partition_cost(bwt_data *b, double opt_estimate)
{
  double apost_est, apost_est_tot=0;
  int k,n,sp,ep;

  if(Segment_cost==NULL) return;
  n= b->size;
  // check each segment in the optimal solution
  for(k=sp=0;sp<=n;k++) {
    ep=End_point[sp]; 
    assert(ep>sp);     // bwt[sp] -> bwt[sp] is a segment
    apost_est_tot += apost_est = pseudo_compression(b->bwt+sp,ep-sp);
    if(apost_est>0 && fabs(apost_est-Segment_cost[sp])/apost_est>.01)
      fprintf(stderr,"Error in segment %d [%d,%d]: est=%f real=%f bits\n",
	      k,sp,ep-1,Segment_cost[sp],apost_est);
    sp = ep;               // starting point of next segment
  } 
  // check that total cost is equal to sum of segment costs
  if(fabs(apost_est_tot-opt_estimate)>0.01)
    fprintf(stderr,"Error in optimal cost: est=%f real=%f bits\n",
	    opt_estimate,apost_est_tot);
  return;
}


/* *****************************************************************
   compute the cost of encoding a segment whose statistics are
   in *s. In this function the cost is defined as 
     #(dist. chars)log(Alpha_size) + Lambda H_0
   For greater speed, logarithms are to the base e since 
   a multiplicative constant does not affect the optimal partition
   ***************************************************************** */ 
double comulative_cost(stats *s)
{
  int i, numchar=0;
  double tot=0, log_alpha_size,Lambda;

  Lambda = 1.0;  
  log_alpha_size = log(Alpha_size);
  assert(s->count_tot>0);
  for(i=0;i<Alpha_size;i++)
    if(s->count[i]!=0) {
      assert(s->count[i]>0);
      tot += s->count[i]*log(((double) s->count_tot)/s->count[i]);
      numchar++;
    }
  assert(tot>=0);
  if(numchar == 1) 
    return log_alpha_size + Lambda * (1+log(s->count_tot+1));
  return numchar * log_alpha_size + Lambda * tot;
}


/* ============================================================
   functions for handling the stack
   ============================================================ */

// add a node to the stack. the node contains the single character c=bwt[pos]
void push_node(int pos, int lcp, int c)
{
  if(++Top==Stack_size)
    enlarge_stack();
  assert(Top<Stack_size);

  End_point[pos] = pos+1;
  Stack[Top].depth = lcp;
  Stack[Top].start = pos;
  Stack[Top].end = pos+1;
  init_stats(&Stack[Top].local_stats, c);

  // --- code for Usage==Advanced to be inserted here
  // Stack[Top].first_run = 1;
  // Stack[Top].last_run = 0;   // last run coincides with first run

  Stack[Top].cost = compression_cost(&Stack[Top]);

  if(Segment_cost!=NULL)
    Segment_cost[pos]=Stack[Top].cost;

}


/* *****************************************************************
   enlarge the global stack by 1000 positions.
   all fields in the newly allocated stack positions are not 
   initialized with the exception of (Stack[i].local_stats).count 
   which is set to NULL. (A !NULL field is used to denote a count 
   field containing a pointer to Alpha_size integers).
   ***************************************************************** */
void enlarge_stack()
{
  int i,old_size;

  old_size = Stack_size;
  Stack_size += 1000;

  Stack = (fragment_node *) realloc(Stack,Stack_size*sizeof(fragment_node));
  if(Stack==NULL) 
    out_of_mem("enlarge_stack");
  for(i=old_size;i<Stack_size;i++) 
    (Stack[i].local_stats).count = NULL;
}

// count and deallocate the memory used for the stack
int free_stack(void)
{
  int free_stats(stats *t);
  int i,mem=0;

  for(i=Stack_size-1;i>=0;i--) {
    if((Stack[i].local_stats).count!=NULL)
      mem += free_stats(&(Stack[i].local_stats));
  }
  mem += Stack_size*sizeof(fragment_node);
  free(Stack);
  init_empty_stack();
  return mem; 
}


void init_empty_stack(void)
{
  Stack=NULL;
  Stack_size = 0;
  Top=-1;
}


/* ================================================================
   functions for handling statistics
   ================================================================ */

// init *s with a single occurrence of symbol c
void init_stats(stats *s, int c)
{
  int i;

  assert(c>=0 && c<Alpha_size);
  // alloc memory if necessary
  if(s->count==NULL) {
    s->count = (int *) malloc(Alpha_size * sizeof(int));
    if(s->count==NULL)
      out_of_mem("init_stats");
  }
  // init statistics
  for(i=0;i<Alpha_size;i++)
    s->count[i]=0;
  s->count[c] = s->count_tot = 1;
}

/* ******************************************************************
   add the statistics in *s2 to those in *s1
   ****************************************************************** */
void merge_stats(stats *s1, stats *s2)
{
  int i;

  assert(s1->count!=NULL && s2->count!=NULL);
  for(i=0;i<Alpha_size;i++)
    s1->count[i] += s2->count[i];
  s1->count_tot += s2->count_tot;
}


int free_stats(stats *t)
{
  assert(t->count!=NULL);
  free(t->count);
  return Alpha_size*sizeof(int); 
}



/* ========================================================================
   initialization routines
   ======================================================================== */

/* ***********************************************************************
   init static global variables used within this file. 
   *********************************************************************** */ 
void init_global_vars(bwt_data *b, int *lcp, int asize, base_compressor *a, 
		      float *segment_cost)
{
  Bwt = b->bwt;                // array containing the BWT 
  assert(asize<=ALPHA_SIZE);   // check alphabet size is legal
  Alpha_size = asize;          // alphabet size
  End_point = lcp;             // use lcp array for storing fragment endpoints
  Segment_cost = segment_cost; // make segment_cost global
  Algo = a;                    // boosted compressor
  Verbose = _boost_Verbose;    // Verbosity level
}


static void out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f);
  exit(1);
}

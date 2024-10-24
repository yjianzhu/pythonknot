#ifndef BOUND
#define BOUND

#include "poly.h"

#define MAXSTRING 12
#define BIGWEAVE (2*MAXSTRING+2)

/*
---------------------------------------------------------------------------
BOUND.H
  Declarations of public routines found in bound.c
  Macros dealing with weave boundaries
  Public Domain
---------------------------------------------------------------------------
*/
/*
---------------------------------------------------------------------------
  This algorithm models the solved region of a link by a set of simple
  weaves and their associated polynomials, which are called tags.  Due to a
  nice result, any simple weave is uniquely defined by how its inputs are
  matched to its outputs.  A list of the outputs the inputs are matched to,
  in the order of the inputs, is enough to define a simple weave.
  Unfortunatly, there is no similar result about the associated tags.  Since
  I've never seen a workstation capable of storing 9! polynomials, it
  is reasonable to assume no weave of more than 8 inputs ever needs to be
  dealt with.  Weaves of 8 inputs have 16 boundary crossings, so any reduced
  weave with 8 inputs or less can be stored in (8 copies of 4 bits =32 bits).
  The size of an integer nowadays is, you guessed it, 32 bits.  Something of
  type WEAVE is just a simple weave and its tag.  BODY is a compressed
  representation of the simple weave itself.  TAG is the polynomial
  associated with that simple weave.
---------------------------------------------------------------------------
*/

class weave
{
  public:
  int    boundary[2];       /* Representation of the weave, heavily encoded. */
  Poly   tag;        /* Polynomial associated with link represented by weave */
};


/*
---------------------------------------------------------------------------
  Global Variables -- located in bound.c 
---------------------------------------------------------------------------
*/
extern int       list[];                 /* description of first new weave */
extern int       list2[];              /* description of second, if needed */
extern int       old_going_in[];   /* Was *i* an input? *old_going_in[i]*. */
extern int       going_in[];       /* Will *i* be an input? *going_in[i]*. */
extern int       map[];      /* i of old weave becomes map[i] of new weave */
extern int       first;                       /* first boundary crossing to remove */
extern int       second;                     /* second boundary crossing to remove */
extern int       right;                /* Is the crossing being added righthanded? */
extern int       oldcross;        /* number of boundary crossings in the old weave */
extern int       newcross;       /* number of boundary crossings in each new weave */
extern int       oldin;                       /* number of inputs to the old weave */
extern int       newin;                       /* number of inputs to the new weave */

/*
---------------------------------------------------------------------------
  Procedures defined in bound.c
---------------------------------------------------------------------------
*/

/* Manipulate those variables whose values are universal to all weaves at a
   given step.  (These are the global variables declared above).  */
void   b_manip(weave *oldweaves);

/* Add a crossing to a single weave without removing any pair of boundary
   crossings.  Compute list, and perhaps list2, describing the new simple
   weaves.  */
void   b_no_pairs(int *list, int *list2, int *one, int *two);

/* Add a crossing to a simple weave and remove one pair of boundary crossings.
   If the result is one or two simple weaves, well and good.  If the result
   is more complicated than that, b_one_pair is essentially a no-op. */
void   b_one_pair(int *list, int *list2, int *one, int *two);


/*
---------------------------------------------------------------------------
There are three representations of weave boundaries.
1) *list* is an array of ints of size oldcross saying which crossing
   is attached to which other crossing.  If list[i]==j, list[j]==i.
2) The medium-size representation (two ub4's) stores the values of
   *list* for only the inputs (the output can be deduced).  Further, each
   value is stored in 5 bits, with 6 values per ub4.
3) The tiny representation recognizes that the medium-size representation
   is a permutation of the outputs, and permutations can be enumerated.
   The tiny representation is the number for that permutation.
---------------------------------------------------------------------------
*/

/*
---------------------------------------------------------------------------
  The macro b_cross determines if two strings cross in a simple weave in
  standard form.  This formula has been proven to be correct
  (via a large truth table).
---------------------------------------------------------------------------
*/
#define b_cross( inp1, inp2, outp1, outp2) \
 (((inp1)<(inp2))==((inp2)<(outp1))==((outp1)<(outp2))==((inp1)<(outp2)))


/*
---------------------------------------------------------------------------
  Switch the positions of *first* and *second* in the boundary *list*.
---------------------------------------------------------------------------
*/
#define b_switch( list, first, second, temp ) \
if (1) \
{ \
  list[list[first]] = second; \
  list[list[second]] = first; \
  temp = list[first]; \
  list[first] = list[second]; \
  list[second] = temp; \
} else


/*
---------------------------------------------------------------------------
  Should the crossing of *first* and *second* in *list* be righthanded? *i*.
  Boundary crossings are numbered counterclockwise.  Inside a simple weave,
  the string with the lowest input should be an overpass.
---------------------------------------------------------------------------
*/
#define b_right( list, going_in, first, second, i ) \
if (1) \
{ \
  int x = first < list[first] ? first : list[first]; \
  int y = second < list[second] ? second : list[second]; \
  if (x > y) {i = x; x = y; y = i;} \
  i = !(going_in[x] && !going_in[y]); \
} else

#endif  /* BOUND */

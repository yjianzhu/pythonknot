/*
------------------------------------------------------------------------------
  ORDER.H - Bob Jenkins, 1989, Public Domain
  Structures and routines used for choosing an order in which to add crossings
  to the solved region, plus code to make the instructions for how weaves have
  to change as those crossings are added.
------------------------------------------------------------------------------
*/

#ifndef ORDER
#define ORDER

#include "bound.h"
#include "knot.h"

/*  Something of type INSTRUCT can describes manipulations of a weave boundary.
    There are two types of manipulations:  adding a crossing connected by one
    arc, which uses WHICH, PREV, and OVER; and removing a pair of adjacent
    boundary crossings, which requires REDUCTIONS and R.
    Something of type struct INSTRUCT can store one added crossing
        and 9 boundary removals. */

class Instruct
{
  public:
  int  crossing;                                 /* which crossing is added */
  int  which;     /* which arc to add the crossing at (-1 for just removal) */
  int  prev;                 /* is new previous boundary crossing an input? */
  int  over;                             /* is the old string the overpass? */
  int  oldn;          /* number of boundary crossings in the original weave */
  int  newn;             /* number of boundary crossings in the final weave */
  int  reductions;          /* number of boundary crossing pairs to connect */
  int  right;                               /* is this crossing righthanded */
  int  going_in[BIGWEAVE];                    /* which strings are going in */
  int  r0[MAXSTRING+1];               /* first boundary crossing to connect */
  int  r1[MAXSTRING+1];              /* second boundary crossing to connect */
};

/* Procedures defined in order.c */

void   o_make(Link *link, Instruct **list);
                                                        /* make instructions */
/*void   o_show(Instruct *l, int num_crossings);*/
                                                 /* display all instructions */


#endif  /* ifndef ORDER */

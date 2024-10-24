
#ifndef HOMFLY_H
#define HOMFLY_H
#include "poly.h"
#include "knot.h"

// #ifndef KNOT
// class crossing
// {
//   public:
//     dllink *o;                                              /* overpass */
//     dllink *u;                                             /* underpass */
//     int      hand; /* 1 if right handed, -1 if left, 0 if no longer a crossing */
// };


// class Link
// {
//   public:
//     crossing *data;
//     int num_crossings;
// };

// #endif

// #ifndef POLY
// class Term
// {
//   public:
//   int    coef;
//   int    m;
//   int    l;
// };


// class Poly
// {
//   public:
//   Term  *term;
//   int     len;
// };

// #endif



char *homfly_str(char *argv);

Poly *homfly(char *argv);

/* c_homfly: Compute the homfly polynomial for the link */
Poly *c_homfly(Link *link);

#endif // HOMFLY_H

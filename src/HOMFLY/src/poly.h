#ifndef POLY
#define POLY


/*
------------------------------------------------------------------------------
 Something of type TERM is a single term in a polynomial.
 coef is the coefficient of the term.
 L is the power of L in the term.
 M is the power of M in the term.
------------------------------------------------------------------------------
*/
class Term
{
  public:
  long int    coef;
  int    m;
  
  int    l;
};


/*
------------------------------------------------------------------------------
  Something of type POLY is a polynomial of M and L.
  TERM is the array of terms.
  len is the number of terms.
------------------------------------------------------------------------------
*/
class Poly
{
  public:
  Term  *term;
  int     len;
};


/*
------------------------------------------------------------------------------
  Procedures defined in poly.c
------------------------------------------------------------------------------
*/

#define p_init( p) ((p)->len = 0)

#define p_kill( p) \
if (1) \
{ \
  (p)->len = 0; \
} else
    //if ((p)->len) free((char *)(p)->term); \


bool  p_check(Poly *p);                 /* check if a poly is a power of -2 */
void    p_copy(Poly *inp, Poly *outp);                       /* copies a Poly */
char   *p_show(Poly *l);                                 /* displays the Poly */
void    p_add(Poly *inp1, Poly *inp2, Poly *outp);          /* adds two Polys */
void    p_mult(Poly *inp1, Poly *inp2, Poly *outp);     /* multiply two Polys */

/* returns a poly with term added */
void   p_term(int coef, int m, int l, Poly *inp, Poly *outp);

#endif /* POLY */ 

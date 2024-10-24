/*
------------------------------------------------------------------------------
  DLLINK.H
  Declares structures and procedures for the abstract data type DLLINK
  which is implemented in dllink.c .
  Bob Jenkins, 1990.  Public Domain.
------------------------------------------------------------------------------
*/

#ifndef  DLLINK
#define  DLLINK


/*  Something of type DLLINK is a node in a doubly-linked cyclic list.
    C is the data stored in this node.
    A is a pointer to the previous node.
    Z is a pointer to the next node. */

class dllink
{
  public:
  int            c;  /* the crossing number */
  dllink  *a;  /* the previous crossing */
  dllink  *z;  /* the next crossing */
};



/* Procedures defined in dllink.c */

# define l_init( l) (l = ((dllink *)(0)))
void   l_show(dllink *l);                       /* Display all nodes */
void   l_add(dllink *inp, int c, dllink **outp);      /* Add a node */
void   l_del(dllink **l);                           /* Delete a node */

#endif  /* ifndef DLLINK */

/*
---------------------------------------------------------------------------
  MODEL.H
  By Bob Jenkins, August 1990, in association with my Masters Thesis
  Structures and procedures used for modelling complicated weaves
  Public Domain
---------------------------------------------------------------------------
*/

#ifndef MODEL
#define MODEL

# include "bound.h"

/*  When weaves are complicated enough to need modeling, a model of them is
    made using singly linked lists of type node.  Something of type node
    is a single node in a string in a weave. */

class node
{
  public:
  int          self;                     /* number of string m belongs to */
  int          right;      /* the crossing this node is in is righthanded */
  int          over;                          /* this node is an overpass */
  int          correct;             /* is this node where it ought to be? */
  node  *m;                          /* other node in this crossing */
  node  *z;                             /* next node in this string */
  int          o1;                         /* number of original string 1 */
  int          o2;                         /* number of original string 2 */
};


#define BIGMODEL (((MAXSTRING)*(MAXSTRING-1)+2)*sizeof(node))

/* Public procedures defined in model.c */

/* Handle weaves that need to be modeled */
void m_model_weave(int *list, weave *oldweave, weave *newweaves);

#endif /* MODEL */


#ifndef __GAUSS_SOLVER__H__
#define __GAUSS_SOLVER__H__

#include "globaldefinitions.h"

dVector gausSolve(dMatrix const &A, dVector const &b);
dMatrix gausInv(dMatrix const &A);

#endif
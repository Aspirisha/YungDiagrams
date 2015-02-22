#ifndef __LINEAR_PROBLEM_SOLVER__
#define __LINEAR_PROBLEM_SOLVER__

#include "globaldefinitions.h"

bool solveTransportationPotential(dMatrix const &C, dVector const &a, dVector const &b, dMatrix &optX);
// everything less than this we will assume to be 0. 


#endif
#ifndef __GLOBAL_DEFINITIONS__H__
#define __GLOBAL_DEFINITIONS__H__

// uncomment line below to see debug dumps
//#define __DEBUG__MODE__ON__ 

#include <vector>//
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::matrix<double> dMatrix; // matrix of doubles
typedef boost::numeric::ublas::matrix<char> cMatrix;
typedef boost::numeric::ublas::vector<double> dVector;
typedef std::vector<double> line;
typedef enum {EXTRA_SUPPLIER = 1, EXTRA_PUCRHASER, NO_EXTRAS} extraType;

const int max_dimension = 21;
const double machine_zero = 1e-13;
#endif
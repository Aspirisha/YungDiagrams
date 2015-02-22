#include "ioInterface.h"
#include <fstream>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

void __DUMP__VARIABLE__(const dMatrix &A, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << A << endl;
  #endif
}

void __DUMP__VARIABLE__(const boost::numeric::ublas::matrix<char> &A, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << A << endl;
  #endif
}

void __DUMP__VARIABLE__(const dVector &b, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << b << endl;
  #endif
}

void __DUMP__VARIABLE__(int n, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << n << endl;
  #endif
}
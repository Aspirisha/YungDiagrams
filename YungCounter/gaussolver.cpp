#include "gaussolver.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

dVector gausSolve(dMatrix const &A, dVector const &b)
{
  int n = b.size();
  dMatrix B = A;
  dVector r = b;
  dVector x(n);

  std::vector<int> order(n);
  for (int i = 0; i < n; i++)
    order[i] = i;

  for (int i = 0; i < n; i++)
  {
    double notZero = B(i, i);
    for (int j = i + 1; (j < n && !notZero); j++)
    {
      if (B(i, j) != 0) // swap j line with i
      {
        int buf = order[i];
        order[i] = order[j];
        order[j] = buf;
        dVector temp = column(B, i);
        column(B, i) = column(B, j);
        column(B, j) = temp;
        notZero = B(i, i);
      }
    }
    if (!notZero)
      return dVector(0); // no solution cause i-line is zero line

    row(B, i) = row(B, i) / notZero;
    r(i) /= notZero;
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        double coef = B(j, i);
        for (int k = 0; k < n; k++)
          B(j, k) -= coef * B(i, k);
        r(j) -= coef * r(i);
      }
    }
  }
  for (int i = 0; i < n; i++)
  {
    int index = order[i];
    x(index)= r(i);
  }
  return x;
}

dMatrix gausInv(dMatrix const &A)
{
  int n = A.size1();
  dMatrix B = A;
  dMatrix invA = identity_matrix<double>(n);

  for (int i = 0; i < n; i++)
  {
    double notZero = B(i, i);
    for (int j = i + 1; (j < n && !notZero); j++)
    {
      if (B(j, i) != 0) // swap j line with i
      {
        dVector temp = row(B, i);
        row(B, i) = row(B, j);
        row(B, j) = temp;
        temp = row(invA, i);
        row(invA, i) = row(invA, j);
        row(invA, j) = temp;
        notZero = B(i, i);
      }
    }
    if (!notZero)
      return dMatrix(0, 0); // no solution cause i-line is zero line

    row(B, i) = row(B, i) / notZero;
    row(invA, i) = row(invA, i) / notZero;
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        double coef = B(j, i);
        for (int k = 0; k < n; k++)
        {
          B(j, k) -= coef * B(i, k);
          invA(j, k) -= coef * invA(i, k);
        }
      }
    }
  }
  return invA;
}

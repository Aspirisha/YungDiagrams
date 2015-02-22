#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <deque>
#include "linearproblemsolver.h"
#include "ioInterface.h"

using namespace boost::numeric::ublas;

typedef enum {VERT = -1, HOR = 1, ANY = 0} direction;

struct point_t
{
  point_t() : x(0), y (0) {}
  point_t(int _x, int _y) : x(_x), y(_y) {}
  bool operator==(const point_t &p)
  {
    return (p.x == x && p.y == y);
  }
  int x; 
  int y;
};

static void _getBasisCombinations(int start, int m, int n, std::vector<int> &chosen, std::vector<std::vector<int> > &variants);
static void NorthWestCornerMethod(dMatrix const &C, dVector const &a, dVector const &b, dMatrix &optX, double &allProduct, std::vector<point_t> &basisPoints);
static void countPotentials(dMatrix const &C, std::vector<point_t> &basisPoints, dVector &u, dVector &v);
static void countFreePointsMarks(dMatrix const &C, dVector const &u, dVector const &v, std::vector<point_t> const &basisPoints, dMatrix &delta);
static void buildCycle(std::vector<point_t> const &basisPoints, point_t start, int m, int n,  std::deque<point_t> &q);
bool _buildCycle(matrix<char> A, point_t curPoint, int m, int n, direction dir, std::deque<point_t> &q);

bool solveTransportationPotential(dMatrix const &C, dVector const &a, dVector const &b, dMatrix &optX)
{
  // find initial plan via NW corner method
  int m = C.size1();
  int n = C.size2();
  double allProduct = 0;
  std::vector<point_t> basisPoints(m + n - 1);

  NorthWestCornerMethod(C, a, b, optX, allProduct, basisPoints);
  __DUMP__VARIABLE__(optX, "optX");
  double currentTarget = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      currentTarget += optX(i, j) * C(i, j);

  double restProduct = allProduct;
  int x = 0; // current point
  int y = 0;

  while (true)
  {
    dVector u(m);
    dVector v(n);
    countPotentials(C, basisPoints, u, v);

    __DUMP__VARIABLE__(u, "u");
    __DUMP__VARIABLE__(v, "v");
    __DUMP__VARIABLE__(optX, "optX");
    dMatrix delta;
    countFreePointsMarks(C, u, v, basisPoints, delta);
    __DUMP__VARIABLE__(delta, "delta");
    bool allNegative = true;
    point_t maxFreePoint;
    double maxPosDelta = -1e307;

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (delta(i, j) > machine_zero)
        {
          allNegative = false;
          if (maxPosDelta < delta(i, j))
          {
            maxPosDelta = delta(i, j);
            maxFreePoint = point_t(j, i);
          }
        }
      }
    }

    if (allNegative) // optimal solution found
      break;

    std::deque<point_t> q;
    buildCycle(basisPoints, maxFreePoint, m, n, q);
    double lambda = 1e307;
    point_t pointToExclude;
    for (int i = 1; i < q.size(); i += 2)
    {
      if (lambda > optX(q[i].y, q[i].x))
      {
        lambda = optX(q[i].y, q[i].x);
        pointToExclude = q[i];
      }
    }
    for (int i = 0; i < basisPoints.size(); i++)
    {
      if (pointToExclude == basisPoints[i])
        basisPoints[i] = maxFreePoint;
    }
    for (int i = 0; i < q.size(); i++)
      optX(q[i].y, q[i].x) += lambda * (1 - 2 * (i % 2));
    
    currentTarget -= lambda * maxPosDelta;
  }
  return true;
}

void NorthWestCornerMethod(dMatrix const &C, dVector const &a, dVector const &b, dMatrix &optX, double &allProduct, std::vector<point_t> &basisPoints)
{
  int m = C.size1();
  int n = C.size2();
  dVector rest_a = a;
  dVector rest_b = b;

  optX.resize(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      optX(i, j) = 0;


  for (int i = 0; i < m; i++)
    allProduct += a(i);

  double restProduct = allProduct;
  int x = 0; // current point
  int y = 0;
  int ind = 0;
  while (restProduct > machine_zero) // manage machine 0
  {
    if (rest_a(y) < machine_zero && rest_b(x) >= machine_zero)
      y++;
    else if (rest_b(x) < machine_zero && rest_a(y) >= machine_zero)
      x++;
    else if (rest_b(x) < machine_zero && rest_a(y) < machine_zero) // degenerate point
    {
      x++;
      
      if (C(x, y) == 1e307)

        continue;

      basisPoints[ind].x = x;
      basisPoints[ind].y = y;
      y++;
      ind++;
    }
    else
    {
      double d = std::min(rest_a(y), rest_b(x));
      rest_a(y) -= d;
      rest_b(x) -= d;
      restProduct -= d;
      optX(y, x) = d;
      basisPoints[ind].x = x;
      basisPoints[ind].y = y;
      ind++;
    }
  }
}

void countPotentials(dMatrix const &C, std::vector<point_t> &basisPoints, dVector &u, dVector &v)
{
  u(0) = 0;
  std::deque<int> qu;
  std::deque<int> qv;
  int m = C.size1();
  int n = C.size2();
  std::vector<char> def_u(m); // array that shows if u(i) is defined or not
  std::vector<char> def_v(n);
  
  qu.push_back(0);
  def_u[0] = 1;
  for (int i = 1; i < m; i++)
    def_u[i] = 0;
  for (int i = 0; i < n; i++)
    def_v[i] = 0;

  while (!qu.empty() || !qv.empty())
  {
    if (!qu.empty())
    {
      int y = qu[0];
      qu.pop_front();
      for (int ind = 0; ind < basisPoints.size(); ind++)
      {
        if (basisPoints[ind].y == y && def_v[basisPoints[ind].x] == 0)
        {
          v(basisPoints[ind].x) = C(basisPoints[ind].y, basisPoints[ind].x) + u(y); // here may be -
          def_v[basisPoints[ind].x] = 1;
          qv.push_back(basisPoints[ind].x);
        }
      }
    }
    else if (!qv.empty())
    {
      int x = qv[0];
      qv.pop_front();
      for (int ind = 0; ind < basisPoints.size(); ind++)
      {
        if (basisPoints[ind].x == x && def_u[basisPoints[ind].y] == 0)
        {
          u(basisPoints[ind].y) = v(x) - C(basisPoints[ind].y, basisPoints[ind].x); // here may be -
          def_u[basisPoints[ind].y] = 1;
          qu.push_back(basisPoints[ind].y);
        }
      }
    }
  }
}

void countFreePointsMarks(dMatrix const &C, dVector const &u, dVector const &v, std::vector<point_t> const &basisPoints, dMatrix &delta)
{
  int m = C.size1();
  int n = C.size2();
  delta.resize(m, n);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      delta(i, j) = 1;

  for (int i = 0; i < basisPoints.size(); i++)
    delta(basisPoints[i].y, basisPoints[i].x) = 0;

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (delta(i, j))
        delta(i, j) = v(j) - u(i) - C(i, j);
}

void buildCycle(std::vector<point_t> const &basisPoints, point_t start, int m, int n,  std::deque<point_t> &q)
{
  matrix<char> A(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A(i, j) = '0';
  matrix<char> cycleMatrix;
  for (int i = 0; i < basisPoints.size(); i++)
    A(basisPoints[i].y, basisPoints[i].x) = '?';
  A(start.y, start.x) = 'S';
  __DUMP__VARIABLE__(A, "A");

  _buildCycle(A, start, m, n, ANY, q);
}

bool _buildCycle(matrix<char> A, point_t curPoint, int m, int n, direction dir, std::deque<point_t> &q)
{
  __DUMP__VARIABLE__(A, "A");
  if (dir == VERT || dir == ANY)
  {
    for (int i = 0; i < m; i++)
    {
      if (A(i, curPoint.x) == '?') // try
      {
        matrix<char> B = A;
        A(curPoint.y, curPoint.x) == '-' ? B(i, curPoint.x) = '+' :  B(i, curPoint.x) = '-';
        point_t newPoint(curPoint.x, i);
        bool flag = _buildCycle(B, newPoint, m, n, HOR, q);
        if (flag)
        {
          q.push_front(curPoint);
          return true;
        }
      }
      else if (dir != ANY && A(i, curPoint.x) == 'S') // found cycle!
      {
        q.push_front(curPoint);
        return true;
      }
    }
  }

  if (dir == HOR || dir == ANY)
  {
    for (int i = 0; i < n; i++)
    {
      if (A(curPoint.y, i) == '?') // try
      {
        matrix<char> B = A;
        A(curPoint.y, curPoint.x) == '-' ? B(curPoint.y, i) = '+' :  B(curPoint.y, i) = '-';
        point_t newPoint(i, curPoint.y);
        bool flag = _buildCycle(B, newPoint, m, n, VERT, q);
        if (flag)
        {
          q.push_front(curPoint);
          return true;
        }
      }
      else if (dir != ANY && A(curPoint.y, i) == 'S') // found cycle!
      {
        q.push_front(curPoint);
        return true;
      }
    }
  }
  return false;
}
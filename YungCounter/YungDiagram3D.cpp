#include <fstream>
#include <iostream>
#include <algorithm>
#include "YungDiagram3D.h"

using namespace std;

const int YungDiagram3DHandler::maxCellsNumber = 20;
std::vector<size_t> YungDiagram3DHandler::primes;
std::map<ind_pair, size_t> YungDiagram3DHandler::primesToCells;
bool YungDiagram3DHandler::primesAreCounted = false;

YungDiagram3D::YungDiagram3D()
{
  m_cols.insert(col_element(ind_pair(0, 0), 1));
  m_cellsNumber = 1;
  m_number = 1;
  m_rowsX = 1;
  m_rowsY = 1;
  m_numberIsCounted = true;
}

YungDiagram3D::YungDiagram3D(const boost::xint::integer &number)
{
  if (YungDiagram3DHandler::arePrimesCounted() == false)
    YungDiagram3DHandler::countPrimes();

  m_number = number;
  for (size_t i = 0, ind = 0; i < YungDiagram3DHandler::maxCellsNumber; i++)
  {
    for (size_t j = 0; j < i; j++, ind++)
    {
      size_t prime = YungDiagram3DHandler::primesToCells[ind_pair(i, j)];
      size_t h = 0;
      while (number % prime == 0)
      {
        h++;
        number /= prime;
      }
    }
  }
}

void YungDiagram3DHandler::countPrimes()
{
  size_t neededSize = maxCellsNumber * (maxCellsNumber + 1) / 2;
  size_t testsNumber = neededSize * log(neededSize);
  do
  {
    testsNumber *= 10;

    primes.clear();
    vector<bool> isPrime(testsNumber); // if isPrime[i] is 1, then i + 2 is prime number
    fill(isPrime.begin(), isPrime.end(), true);
    for (size_t i = 0; i < testsNumber; i++)
    {
      if (isPrime[i])
      {
        size_t currentPrime = i + 2;
        primes.push_back(currentPrime);
        for (size_t j = 2 * currentPrime; j < testsNumber; j++)
          isPrime[j] = false;
      }
    }
  } while (primes.size() < neededSize);

  for (size_t i = 0, ind = 0; i < maxCellsNumber; i++)
  {
    for (size_t j = 0; j < i; j++, ind++)
    {
      primesToCells.insert(col_element(ind_pair(j, i - j), primes[ind]));
    }
  }
}
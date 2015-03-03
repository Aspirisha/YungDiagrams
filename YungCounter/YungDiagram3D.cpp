#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include "YungDiagram3D.h"

using namespace std;

const int YungDiagram3DHandler::maxCellsNumber = 100;
std::vector<size_t> YungDiagram3DHandler::primes;
std::map<ind_pair, size_t> YungDiagram3DHandler::primesToCells;
bool YungDiagram3DHandler::primesAreCounted = false;

YungDiagram3D::~YungDiagram3D()
{
}

YungDiagram3D::YungDiagram3D()
{
  m_cols.insert(col_element(ind_pair(0, 0), 1));
  m_cellsNumber = 1;
  m_number = 1;
  m_rowsY.push_back(1);
  m_numberIsCounted = true;
}

YungDiagram3D::YungDiagram3D(const boost::xint::integer &number)
{
  if (YungDiagram3DHandler::arePrimesCounted() == false)
    YungDiagram3DHandler::countPrimes();

  m_number = number;
  int y = YungDiagram3DHandler::maxCellsNumber - 1;
  size_t prime = YungDiagram3DHandler::primesToCells[ind_pair(0, y)];
  ++m_number;
  m_cellsNumber = 0;

  m_rowsY.resize(YungDiagram3DHandler::maxCellsNumber + 1);
  for (int i = YungDiagram3DHandler::maxCellsNumber - 1; i >= 0; i--)
  {
    for (int j = 0; j <= i; j++)
    {
      size_t a = 0, b = 0;
      size_t y = i - j;


      if (m_rowsY[j] >= y + 2)
        b = m_cols[ind_pair(j, y + 1)];
      if (m_rowsY[j + 1] >= y + 1)
        b = m_cols[ind_pair(j + 1, y)];
      size_t h = max<>(a, b);

      prime = YungDiagram3DHandler::primesToCells[ind_pair(j, y)];
      while (m_number % prime == 0)
      {
        m_number /= prime;
        h++;
      } 
      m_cols[ind_pair(j, y)] = h;
      if (h && m_rowsY[j] == 0)
        m_rowsY[j] = y + 1;

      m_cellsNumber += h;
    }
  }

  if (m_number > 1)
  {
    cerr << "Number contains too big prime factor. Unable to create such diagram.\n";
    m_number = 1;
    m_cellsNumber = 1;
    m_cols.clear();
    m_cols[ind_pair(0, 0)] = 1;
    m_rowsY.resize(1);
    m_rowsY[0] = 1;
    m_numberIsCounted = true;
    return;
  }

  // downsize vector
  for (int i = 0; i < m_rowsY.size(); i++)
  {
    if (m_rowsY[i] == 0)
    {
      m_rowsY.resize(i);
      break;
    }
  }

  m_numberIsCounted = true;
  m_number = number;
}

boost::xint::integer YungDiagram3D::GetDiagramNumber()
{
  if (m_numberIsCounted)
    return m_number;

  if (YungDiagram3DHandler::arePrimesCounted() == false)
    YungDiagram3DHandler::countPrimes();

  m_number = boost::xint::integer(1);
  for (size_t i = 0; i < m_rowsY.size(); i++)
  {
    for (size_t j = 0; j < m_rowsY[i]; j++)
    {
      size_t delta_h = m_cols[ind_pair(i, j)] - max<>(m_cols[ind_pair(i + 1, j)], m_cols[ind_pair(i, j + 1)]);

      boost::xint::integer power = boost::xint::pow(YungDiagram3DHandler::primesToCells[ind_pair(i, j)], 
        boost::xint::integer(delta_h));
      m_number *= power;
    }
  }

  --m_number;
  m_numberIsCounted = true;
  return m_number;
}

void YungDiagram3DHandler::countPrimes()
{
  size_t neededSize = maxCellsNumber * (maxCellsNumber + 1) / 2;
  size_t testsNumber = neededSize * log(neededSize);
  do
  {
    testsNumber *= 2;

    primes.clear();
    vector<bool> isPrime(testsNumber); // if isPrime[i] is 1, then i + 2 is prime number
    fill(isPrime.begin(), isPrime.end(), true);
    for (size_t i = 0; i < testsNumber; i++)
    {
      if (isPrime[i])
      {
        size_t currentPrime = i + 2;
        primes.push_back(currentPrime);
        for (size_t j = 2 * currentPrime; j < testsNumber; j += currentPrime)
          isPrime[j - 2] = false;
      }
    }
  } while (primes.size() < neededSize);

  cout << primes[maxCellsNumber * (maxCellsNumber - 1) / 2] << endl;
  for (size_t i = 0, ind = 0; i < maxCellsNumber; i++)
  {
    for (size_t j = 0; j <= i; j++, ind++)
    {
      primesToCells[ind_pair(j, i - j)] = primes[ind];
    }
  }

  primesAreCounted = true;
}

void YungDiagram3D::saveToFile(const char *fileName) const
{
  ofstream out(fileName);

  out << m_rowsY.size() << endl;
  for (size_t i = 0; i < m_rowsY.size(); i++)
    out << m_rowsY[i] << " ";
  out << endl;

  for (size_t i = 0; i < m_rowsY.size(); i++)
  {
    for (size_t j = 0; j < m_rowsY[i]; j++)
    {
      out << m_cols[ind_pair(i, j)] << " ";
    }
    out << endl;
  }
}

void YungDiagram3D::readFromFile(const char *fileName)
{
  ifstream in(fileName);

  size_t rowsNumber = 0;
  m_cellsNumber = 0;

  in >> rowsNumber;
  m_rowsY.resize(rowsNumber);
  for (size_t i = 0; i < rowsNumber; i++)
    in >> m_rowsY[i];

  for (size_t i = 0; i < m_rowsY.size(); i++)
  {
    for (size_t j = 0; j < m_rowsY[i]; j++)
    {
      size_t col = 0;
      in >> col;
      m_cols[ind_pair(i, j)] = col;
      m_cellsNumber += col;
    }
  }

  m_numberIsCounted = false;
}

void YungDiagram3D::printToConsole() const
{
  for (size_t i = 0; i < m_rowsY.size(); i++)
  {
    for (size_t j = 0; j < m_rowsY[i]; j++)
    {
      cout << m_cols[ind_pair(i, j)] << " ";
    }
    cout << endl;
  }
}

YungDiagram3D *YungDiagram3DHandler::getRandomWalkDiagram(ProcessType3D type, size_t cellsNumber)
{
  YungDiagram3D *d = new YungDiagram3D;
  static unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(seed);

  /*d->m_cols[ind_pair(0, 0)] = 2;
  d->m_cols[ind_pair(1, 0)] = 2;
  d->m_cols[ind_pair(0, 1)] = 1;
  d->m_cols[ind_pair(0, 2)] = 1;
  d->m_cellsNumber = 6;
  d->m_rowsY.resize(2);
  d->m_rowsY[0] = 3;
  d->m_rowsY[1] = 1;*/


  switch (type)
  {
  case HOOKS:
    for (size_t i = 1; i < cellsNumber; i++)
    {
      vector<double> probs;
      vector<ind_pair> newPoints;
      for (size_t j = 0; j <= d->m_rowsY.size(); j++) // less or equal!
      {
        size_t max_k = 0;
        if (j < d->m_rowsY.size())
          max_k = d->m_rowsY[j];
        for (size_t k = 0; k <= max_k; k++)
        {
          // check if can add a cell to (j, k)
          bool canAddCell = true;
          if (max_k != 0) // else we don't need to check, we can add to an empty cell
          {
            if (j > 0 && d->m_cols[ind_pair(j - 1, k)] == d->m_cols[ind_pair(j, k)])
              canAddCell = false;
            if (k > 0 && d->m_cols[ind_pair(j, k - 1)] == d->m_cols[ind_pair(j, k)])
              canAddCell = false;
          }
          
          if (canAddCell)
          {
            double newProb = 1.0;
            const size_t z = d->m_cols[ind_pair(j, k)]; // cell (j, k, z) is counted
            const size_t y = k + 1; // (x, y, z) are now the cell coordinates over which we can add a cell
            const size_t x = j + 1;

            size_t dx = 0; // dx = max dx : h(x + dx, y) >= 1, i.e. max dx : rows_Y[x - 1 + dx] >= y; doesn't count (x,y,z)
            size_t dy = max_k - k - 1; // dy = max dy : h(x, y + dy) >= 1, i.e.  dy = rows_Y[x - 1] - k - 1, cause we don't need to count(x,y,z) cell here
            size_t dz = 0;

            for (size_t l = j + 1; l < d->m_rowsY.size(); l++)
            {
              if (d->m_rowsY[l] >= y)
                dx++;
            }

            dz = 0; // just goes over 1 .. z. Nobody counts (x,y,z)
            for (size_t l = z; l > 0; l--) // l counts (x, y, z) cell and centre of hook
            {
              newProb *= (l + dx + dy);
              newProb /= (l + dx + dy + 1);
              
              dz++;
              while (d->m_cols[ind_pair(j + dx, k)] <= dz && dx)
                dx--;
              while (d->m_cols[ind_pair(j, k + dy)] <= dz && dy)
                dy--;
            }

            // no more (x,y,z)!
            dx = 0;
            for (size_t l = x - 1; l > 0; l--) // l doesn't count (x,y,z) cell, counts centre of hook
            {
              dz = (d->m_cols[ind_pair(dx, k)] - z - 1); // cell (x, y, z) is not counted
              dy = d->m_rowsY[dx] - y; //  cell (x,y,z) is not counted
              while (d->m_cols[ind_pair(dx, k + dy)] < z && dy) // here wass j
                dy--;
              newProb *= (l + dy + dz);
              newProb /= (l + dy + dz + 1);
              dx++;
            }

            dx = d->m_rowsY.size() - x;
            dy = 0;
            while (d->m_cols[ind_pair(j + dx, dy)] < z && dx)
              dx--;
            for (size_t l = y - 1; l > 0; l--) // l not counts (x,y,z), counts centre of hook
            {
              dz = d->m_cols[ind_pair(j, dy)] - z - 1;
              
              newProb *= (l + dx + dz);
              newProb /= (l + dx + dz + 1);
              dy++;
              while (d->m_rowsY[dx + j] <= dy && dx)
                dx--;
              while (d->m_cols[ind_pair(j + dx, dy)] < z && dx)
                dx--;
            }

            probs.push_back(newProb);
            newPoints.push_back(ind_pair(j, k));
          }
        }
      }
      size_t j = 0;
      std::discrete_distribution<size_t> distrib(probs.size(), 0, probs.size() - 1,
          [&probs, &j](double){
          auto w = probs[j];
          ++j;
          return w;
      });
      size_t randomAncestor = distrib(generator);
      d->m_cols[newPoints[randomAncestor]]++;
      d->m_cellsNumber++;
      if (d->m_rowsY.size() < newPoints[randomAncestor].first + 1)
        d->m_rowsY.push_back(1);
      else if (d->m_rowsY[newPoints[randomAncestor].first] < newPoints[randomAncestor].second + 1)
        d->m_rowsY[newPoints[randomAncestor].first]++;
    }
    break;
  }

  return d;
}

double YungDiagram3DHandler::countHook(YungDiagram3D &d, size_t x, size_t y, size_t z)
{
  size_t dz = d.m_cols[ind_pair(x - 1, y - 1)] - z + 1;
  size_t dx = 0, dy = 0;
  while (d.m_cols[ind_pair(x + dx, y - 1)] >= z)
    dx++;
  while (d.m_cols[ind_pair(x - 1, y + dy)] >= z)
    dy++;

  return (dx + dy + dz);
}

YungDiagram3D *YungDiagram3DHandler::getRandomWalkDiagramFast(ProcessType3D type, size_t cellsNumber)
{
  YungDiagram3D *d = new YungDiagram3D;
  static unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(seed);
  static map<ind_pair, vector<double> > hooks;

  hooks[ind_pair(0, 0)].push_back(1);

  switch (type)
  {
  case HOOKS:
    for (size_t i = 1; i < cellsNumber; i++)
    {
      vector<double> probs;
      vector<ind_pair> newPoints;
      for (size_t j = 0; j <= d->m_rowsY.size(); j++) // less or equal!
      {
        size_t max_k = 0;
        if (j < d->m_rowsY.size())
          max_k = d->m_rowsY[j];
        for (size_t k = 0; k <= max_k; k++)
        {
          // check if can add a cell to (j, k)
          bool canAddCell = true;
          if (max_k != 0) // else we don't need to check, we can add to an empty cell
          {
            if (j > 0 && d->m_cols[ind_pair(j - 1, k)] == d->m_cols[ind_pair(j, k)])
              canAddCell = false;
            if (k > 0 && d->m_cols[ind_pair(j, k - 1)] == d->m_cols[ind_pair(j, k)])
              canAddCell = false;
          }
          
          if (canAddCell)
          {
            double newProb = 1.0;
            const size_t z = d->m_cols[ind_pair(j, k)]; // cell (j, k, z) is counted
            const size_t y = k + 1; // (x, y, z) are now the cell coordinates over which we can add a cell
            const size_t x = j + 1;

            size_t dx = 0; // dx = max dx : h(x + dx, y) >= 1, i.e. max dx : rows_Y[x - 1 + dx] >= y; doesn't count (x,y,z)
            size_t dy = max_k - k - 1; // dy = max dy : h(x, y + dy) >= 1, i.e.  dy = rows_Y[x - 1] - k - 1, cause we don't need to count(x,y,z) cell here
            size_t dz = 0;

            for (size_t l = j + 1; l < d->m_rowsY.size(); l++)
            {
              if (d->m_rowsY[l] >= y)
                dx++;
            }

            dz = 0; // just goes over 1 .. z. Nobody counts (x,y,z)
            for (size_t l = z; l > 0; l--) // l counts (x, y, z) cell and centre of hook
            {
              newProb *= (l + dx + dy);
              newProb /= (l + dx + dy + 1);
              
              dz++;
              while (d->m_cols[ind_pair(j + dx, k)] <= dz && dx)
                dx--;
              while (d->m_cols[ind_pair(j, k + dy)] <= dz && dy)
                dy--;
            }

            // no more (x,y,z)!
            dx = 0;
            for (size_t l = x - 1; l > 0; l--) // l doesn't count (x,y,z) cell, counts centre of hook
            {
              dz = (d->m_cols[ind_pair(dx, k)] - z - 1); // cell (x, y, z) is not counted
              dy = d->m_rowsY[dx] - y; //  cell (x,y,z) is not counted
              while (d->m_cols[ind_pair(dx, k + dy)] < z && dy) // here wass j
                dy--;
              newProb *= (l + dy + dz);
              newProb /= (l + dy + dz + 1);
              dx++;
            }

            dx = d->m_rowsY.size() - x;
            dy = 0;
            while (d->m_cols[ind_pair(j + dx, dy)] < z && dx)
              dx--;
            for (size_t l = y - 1; l > 0; l--) // l not counts (x,y,z), counts centre of hook
            {
              dz = d->m_cols[ind_pair(j, dy)] - z - 1;
              
              newProb *= (l + dx + dz);
              newProb /= (l + dx + dz + 1);
              dy++;
              while (d->m_rowsY[dx + j] <= dy && dx)
                dx--;
              while (d->m_cols[ind_pair(j + dx, dy)] < z && dx)
                dx--;
            }

            probs.push_back(newProb);
            newPoints.push_back(ind_pair(j, k));
          }
        }
      }
      size_t j = 0;
      std::discrete_distribution<size_t> distrib(probs.size(), 0, probs.size() - 1,
          [&probs, &j](double){
          auto w = probs[j];
          ++j;
          return w;
      });
      size_t randomAncestor = distrib(generator);
      d->m_cols[newPoints[randomAncestor]]++;
      d->m_cellsNumber++;
      if (d->m_rowsY.size() < newPoints[randomAncestor].first + 1)
        d->m_rowsY.push_back(1);
      else if (d->m_rowsY[newPoints[randomAncestor].first] < newPoints[randomAncestor].second + 1)
        d->m_rowsY[newPoints[randomAncestor].first]++;
    }
    break;
  }

  return d;
}
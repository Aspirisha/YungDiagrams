#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include "YungDiagram3D.h"
#include "globaldefinitions.h"
#include "LinearProblemSolver.h"

using namespace std;

const int YungDiagram3DHandler::maxCellsNumber = 100;
std::vector<size_t> YungDiagram3DHandler::primes;
std::map<ind_pair, size_t> YungDiagram3DHandler::primesToCells;
bool YungDiagram3DHandler::primesAreCounted = false;

typedef std::pair<boost::xint::integer, boost::xint::integer> num_pair;

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


      a = m_cols[ind_pair(j, y + 1)];
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

boost::xint::integer YungDiagram3D::GetDiagramNumber(bool forceRecount)
{
  if (m_numberIsCounted && !forceRecount)
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

  for (size_t i = 0, ind = 0; i < maxCellsNumber; i++)
  {
    for (size_t j = 0; j <= i; j++, ind++)
    {
      primesToCells[ind_pair(j, i - j)] = primes[ind];
    }
  }

  primesAreCounted = true;
}

// NB we save in format that is easy to read; 
// to display diagram in matlab, DELETE first two lines!
void YungDiagram3D::saveToFile(const char *fileName, bool forMatlab) const
{
  ofstream out(fileName);

  if (!forMatlab)
  {
    out << m_rowsY.size() << endl;
    for (size_t i = 0; i < m_rowsY.size(); i++)
      out << m_rowsY[i] << " ";
    out << endl;
  }

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

double YungDiagram3D::countTransitiveProb(size_t x, size_t y) const
{
  size_t max_k = 0;
  int j = x - 1;
  int k = y - 1;
  if (x - 1 < m_rowsY.size())
    max_k = m_rowsY[x - 1];

  double newProb = 1.0;
  const size_t z = m_cols[ind_pair(x - 1, y - 1)]; // cell (j, k, z) is counted

  size_t dx = 0; // dx = max dx : h(x + dx, y) >= 1, i.e. max dx : rows_Y[x - 1 + dx] >= y; doesn't count (x,y,z)
  size_t dy = max_k - y; // dy = max dy : h(x, y + dy) >= 1, i.e.  dy = rows_Y[x - 1] - k - 1, cause we don't need to count(x,y,z) cell here
  size_t dz = 0;

  for (size_t l = x; l < m_rowsY.size(); l++)
  {
    if (m_rowsY[l] >= y)
      dx++;
  }

  dz = 0; // just goes over 1 .. z. Nobody counts (x,y,z)
  for (size_t l = z; l > 0; l--) // l counts (x, y, z) cell and centre of hook
  {
    newProb *= (l + dx + dy);
    newProb /= (l + dx + dy + 1);
              
    dz++;
    while (m_cols[ind_pair(j + dx, k)] <= dz && dx)
      dx--;
    while (m_cols[ind_pair(j, k + dy)] <= dz && dy)
      dy--;
  }

  // no more (x,y,z)!
  dx = 0;
  for (size_t l = x - 1; l > 0; l--) // l doesn't count (x,y,z) cell, counts centre of hook
  {
    dz = (m_cols[ind_pair(dx, y - 1)] - z - 1); // cell (x, y, z) is not counted
    dy = m_rowsY[dx] - y; //  cell (x,y,z) is not counted
    while (m_cols[ind_pair(dx, k + dy)] < z && dy) // here wass j
      dy--;
    newProb *= (l + dy + dz);
    newProb /= (l + dy + dz + 1);
    dx++;
  }

  dx = m_rowsY.size() - x;
  dy = 0;
  while (m_cols[ind_pair(j + dx, dy)] < z && dx)
    dx--;
  for (size_t l = y - 1; l > 0; l--) // l not counts (x,y,z), counts centre of hook
  {
    dz = m_cols[ind_pair(j, dy)] - z - 1;
              
    newProb *= (l + dx + dz);
    newProb /= (l + dx + dz + 1);
    dy++;
    while (m_rowsY[dx + j] <= dy && dx)
      dx--;
    while (m_cols[ind_pair(j + dx, dy)] < z && dx)
      dx--;
  }

  return newProb;
}

YungDiagram3D *YungDiagram3DHandler::getRandomWalkDiagram(YungDiagram3D::ProcessType3D type, size_t cellsNumber)
{
  YungDiagram3D *d = new YungDiagram3D;
  static unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(seed);

  switch (type)
  {
  case YungDiagram3D::HOOKS:
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
            probs.push_back(d->countTransitiveProb(j + 1, k + 1));
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

bool YungDiagram3D::isCornerCell(size_t x, size_t y) const
{
  size_t max_k = 0;
  if (x < m_rowsY.size())
     max_k = m_rowsY[x];
  bool canAddCell = true;
  if (max_k != 0) // else we don't need to check, we can add to an empty cell
  {
    if (x > 0 && m_cols[ind_pair(x - 1, y)] == m_cols[ind_pair(x, y)])
      canAddCell = false;
    if (y > 0 && m_cols[ind_pair(x, y - 1)] == m_cols[ind_pair(x, y)])
      canAddCell = false;
  }

  return canAddCell;
}

bool YungDiagram3D::canRemoveCell(size_t x, size_t y) const
{
  bool canRemoveCell = true;
  if (m_cols[ind_pair(x + 1, y)] == m_cols[ind_pair(x, y)])
    canRemoveCell = false;
  if (m_cols[ind_pair(x, y + 1)] == m_cols[ind_pair(x, y)])
    canRemoveCell = false;

  return canRemoveCell;
}

void YungDiagram3D::getPredeccessors(std::map<boost::xint::integer, vector<boost::xint::integer> > &predeccessors,
    std::map<boost::xint::integer, vector<ind_pair> > &dif_cells)
{
  if (!m_numberIsCounted)
    GetDiagramNumber();

  predeccessors[m_number] = vector<boost::xint::integer>();
  dif_cells[m_number] = vector<ind_pair>();

  boost::xint::integer cached_number = m_number;

  size_t colsNumber = m_cols.size();
  m_cellsNumber--;
  for (int k = 0; k < m_rowsY.size(); k++)
  {
    for (int l = 0; l < m_rowsY[k]; l++)
    {
      if (!canRemoveCell(k, l))
        continue;

      int resizeType = 0; // means no resize is needed
      if (--m_cols[ind_pair(k, l)] == 0)
      {
        if (l > 0)
        {
          m_rowsY[k]--;
          resizeType = 1;
        }
        else
        {
          m_rowsY.resize(k);
          resizeType = 2;
        }
      }
      boost::xint::integer num = GetDiagramNumber(true);
      predeccessors[cached_number].push_back(num);
      dif_cells[cached_number].push_back(ind_pair(k, l));

      switch (resizeType)
      {
      case 1:
          m_rowsY[k]++;
          break;
      case 2:
          m_rowsY.resize(k + 1);
          break;
      default:
        break;
      }
      m_cols[ind_pair(k, l)]++;

    }
  }
  m_cellsNumber++; // for proper deletion
}

YungDiagram3D *YungDiagram3DHandler::getRandomWalkDiagramFast(YungDiagram3D::ProcessType3D type, size_t cellsNumber, double power)
{
  YungDiagram3D *d = new YungDiagram3D;
  static unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(seed);
  static map<ind_pair, vector<double> > hooks;

  hooks[ind_pair(0, 0)].push_back(1);

  switch (type)
  {
  case YungDiagram3D::RICHARDSON:
    for (size_t i = 1; i < cellsNumber; i++)
    {
      if ((i - 1) % 1000 == 0)
        cout << i << "-cells diagram is generated.\n";
      vector<ind_pair> newPoints;
      for (size_t j = 0; j <= d->m_rowsY.size(); j++) // less or equal!
      {
        size_t max_k = 0;
        if (j < d->m_rowsY.size())
          max_k = d->m_rowsY[j];
        for (size_t k = 0; k <= max_k; k++)
        {
          // check if can add a cell to (j, k)
          if (d->isCornerCell(j, k))
            newPoints.push_back(ind_pair(j, k));
        }
      }
      size_t j = 0;
      std::uniform_int_distribution<size_t> distrib(0, newPoints.size() - 1);

      size_t randomAncestor = distrib(generator);
      d->m_cols[newPoints[randomAncestor]]++;
      d->m_cellsNumber++;
      if (d->m_rowsY.size() < newPoints[randomAncestor].first + 1)
        d->m_rowsY.push_back(1);
      else if (d->m_rowsY[newPoints[randomAncestor].first] < newPoints[randomAncestor].second + 1)
        d->m_rowsY[newPoints[randomAncestor].first]++;
    }
    break;

  case YungDiagram3D::HOOKS:
  case YungDiagram3D::HOOKS_POWERED:

    for (size_t i = 1; i < cellsNumber; i++)
    {
    /*  if ((i - 1) % 1000 == 0)
        cout << i << "-cells diagram is generated.\n";*/
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
          if (d->isCornerCell(j, k))
          {
            double newProb = 1.0;
            const int z = d->m_cols[ind_pair(j, k)]; // cell (j, k, z) is counted

            for (int l = 0; l < z; l++) // l counts (x, y, z) cell and centre of hook
              newProb -= (newProb / (1 + hooks[ind_pair(j, k)][l]));

            for (int l = 0; l < j; l++) // l doesn't count (x,y,z) cell, counts centre of hook
              newProb -= (newProb / (1 + hooks[ind_pair(l, k)][z]));

            for (int l = 0; l < k; l++) // l not counts (x,y,z), counts centre of hook
              newProb -= (newProb / (1 + hooks[ind_pair(j, l)][z]));

            if (type == YungDiagram3D::HOOKS_POWERED)
              newProb = pow(newProb, power);

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


      hooks[newPoints[randomAncestor]].push_back(1);
      // recount hooks that were touched by addin a cell
      const size_t z = d->m_cols[newPoints[randomAncestor]] - 1; // cell (j, k, z) is counted
      const size_t y = newPoints[randomAncestor].second; // (x, y, z) are now the cell coordinates over which we can add a cell
      const size_t x = newPoints[randomAncestor].first;

      for (size_t l = 0; l < z; l++) // l counts (x, y, z) cell and centre of hook
        hooks[ind_pair(x, y)][l]++;

      for (size_t l = 0; l < x; l++) // l doesn't count (x,y,z) cell, counts centre of hook
        hooks[ind_pair(l, y)][z]++;

      for (size_t l = 0; l < y; l++) // l not counts (x,y,z), counts centre of hook
        hooks[ind_pair(x, l)][z]++;



    }
    break;
  }

  return d;
}

void YungDiagram3DHandler::countCoefficientsForKantorovichMetric(boost::xint::integer diagramNum, 
  map<boost::xint::integer, double> &probs, dVector &c, 
  map<boost::xint::integer, vector<boost::xint::integer> > &preds,  map<boost::xint::integer, vector<ind_pair> > &dif_cells,
  map<boost::xint::integer, double> &newProbs)
{
  size_t predsNum = preds[diagramNum].size();
  c.resize(predsNum);
  double accum = 0;

  for (int i = 0; i < predsNum; i++)
  {
    boost::xint::integer predNum = preds[diagramNum][i];
    YungDiagram3D d(predNum);
    double weight = d.countTransitiveProb(dif_cells[diagramNum][i].first + 1, dif_cells[diagramNum][i].second + 1);
    double sumWeights = 0;
    for (int j = 0; j <= d.m_rowsY.size(); j++)
    {
      int max_k = 0;
      if (j < d.m_rowsY.size())
        max_k = d.m_rowsY[j];

      for (int k = 0; k <= max_k; k++)
      {
        if (d.isCornerCell(j, k))
          sumWeights += d.countTransitiveProb(j + 1, k + 1);
      }
    }
    c(i) = probs[predNum] * weight / sumWeights;
    accum += c(i);
  }

  for (int i = 0; i < predsNum; i++)
    c(i) /= accum;
  newProbs[diagramNum] = accum;
}

double YungDiagram3DHandler::countDistance(YungDiagram3D &d1, YungDiagram3D &d2)
{
  vector<boost::xint::integer> neededDiagrams;
  vector<size_t> numbersOfNeededDiagramsOnLevel;
  vector<char> ancestorFlags; 
  map<boost::xint::integer, vector<boost::xint::integer> >predeccessors;
  map<boost::xint::integer, vector<ind_pair> > dif_cells;

  // probably we don't need to count distance between diagrams that are predecessors of the same diagram from
  // te top level. If it's predecessor of first one, it has first bit to be 1, if of the second - then second bit to be 1.
  // BUT, IF it turns out that this diagram is predecessor for both of top level diagrams,
  // we should count it's flag to be 3 - bits for both

  neededDiagrams.push_back(d1.GetDiagramNumber());
  neededDiagrams.push_back(d2.GetDiagramNumber());
  ancestorFlags.push_back(1);
  ancestorFlags.push_back(2);

  numbersOfNeededDiagramsOnLevel.push_back(2);

  int n = d1.m_cellsNumber;
  int startIndexInVector = 0; // indicate start and end of current level diagrams in neededDiagrams
  int endIndexInVector = 1;

  for (int i = n; i > 2; i--) // we assume that we have metric on the second floor
  {
    size_t neededDiagramsNumberOnLevel = 0;
    for (int j = startIndexInVector; j <= endIndexInVector; j++)
    {
      YungDiagram3D d(neededDiagrams[j]);
      d.getPredeccessors(predeccessors, dif_cells);

      for (boost::xint::integer num : predeccessors[neededDiagrams[j]])
      {
        vector<boost::xint::integer>::iterator iter;
        if ((iter = std::find(neededDiagrams.begin(), neededDiagrams.end(), num)) == neededDiagrams.end())
        {
          if (num > 4)
          {
            neededDiagrams.push_back(num);
            neededDiagramsNumberOnLevel++;
            ancestorFlags.push_back(ancestorFlags[j]);
          }
        } 
        else // maybe modify flags
        {
          size_t ind = iter - neededDiagrams.begin();
          ancestorFlags[ind] |= ancestorFlags[j];
        }
      }
    }

    if (i > 3)
      numbersOfNeededDiagramsOnLevel.push_back(neededDiagramsNumberOnLevel);
    startIndexInVector = endIndexInVector + 1;
    endIndexInVector = neededDiagrams.size() - 1;
  }

  map<num_pair, double> distances; // key is pair of diagram numbers, first number is less than second.
  map<boost::xint::integer, double> probs;

  distances.insert(pair<num_pair, double>(pair<size_t, size_t>(2, 4), 1.0));
  distances.insert(pair<num_pair, double>(pair<size_t, size_t>(2, 3), 1.0));
  distances.insert(pair<num_pair, double>(pair<size_t, size_t>(3, 4), 1.0));
  probs[2] = probs[3] = probs[4] = 1.0 / 3;

  endIndexInVector = neededDiagrams.size() - 1;
  size_t solvedTranspProblems = 0;
  

  for (int i = 3, levelIndex = numbersOfNeededDiagramsOnLevel.size() - 1; i <= n; i++, levelIndex--) 
  {
    startIndexInVector = endIndexInVector - numbersOfNeededDiagramsOnLevel[levelIndex] + 1;

    // last is not needed for it's distances will be counted during previous steps
    map<boost::xint::integer, double> newProbs;
    for (int j = startIndexInVector; j < endIndexInVector; j++)
    {
      dVector c1; 
      countCoefficientsForKantorovichMetric(neededDiagrams[j], probs, c1, predeccessors, dif_cells, newProbs);

      for (int k = j + 1; k <= endIndexInVector; k++)
      {
        if (ancestorFlags[k] == ancestorFlags[j] && (ancestorFlags[k] ^ 3)) // boost!? yes, about 2 times
          continue;

        dVector c2;
        vector<size_t> nums2; // numbers of predecessors of d2
        countCoefficientsForKantorovichMetric(neededDiagrams[k], probs, c2, predeccessors, dif_cells, newProbs);
        
        // now we need to solve transportation problem
        size_t s1 = c1.size();
        size_t s2 = c2.size();
        dMatrix costs(s1, s2);
        dMatrix x(s1, s2);
        for (size_t l = 0; l < s1; l++)
        {
          boost::xint::integer d1 = predeccessors[neededDiagrams[j]][l];
          for (size_t p = 0; p < s2; p++)
          {
            boost::xint::integer d2 = predeccessors[neededDiagrams[k]][p];

            if (d1 == d2)
              costs(l, p) = 0;
            else if (d1 < d2)
              costs(l, p) = distances[num_pair(d1, d2)];
            else
              costs(l, p) = distances[num_pair(d2, d1)];
          }
        }

        if (!solveTransportationPotential(costs, c1, c2, x))
          cout << "Potential method error\n";
        else
        {
          solvedTranspProblems++;
          double dist = 0;
          for (size_t l = 0; l < s1; l++)
          {
            for (size_t p = 0; p < s2; p++)
            {
              dist += costs(l, p) * x(l, p);
            }
          }
          if (neededDiagrams[j] < neededDiagrams[k])
            distances.insert(pair<num_pair, double>(num_pair(neededDiagrams[j], neededDiagrams[k]), dist));
          else
            distances.insert(pair<num_pair, double>(num_pair(neededDiagrams[k], neededDiagrams[j]), dist));
        }
      }
    }
    probs = newProbs;
    endIndexInVector = startIndexInVector - 1;
  }

  cout << "During distance finding were solved " << solvedTranspProblems << " transportation problems.\n";

  if (neededDiagrams[0] < neededDiagrams[1])
    return distances[num_pair(neededDiagrams[0], neededDiagrams[1])];
  return distances[num_pair(neededDiagrams[1], neededDiagrams[0])];
}
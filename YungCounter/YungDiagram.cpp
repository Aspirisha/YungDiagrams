#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <array>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include "gaussolver.h"
#include "YungDiagram.h"
#include "LinearProblemSolver.h"

using namespace std;

boost::xint::integer *YungDiagramHandler::partitionsAmount = 0;
long long *YungDiagramHandler::offsets = 0;
double *YungDiagramHandler::probabilities = 0;
size_t *YungDiagramHandler::numbers = 0;
size_t YungDiagramHandler::levelSize = 0;
const int YungDiagramHandler::maxCellsNumber = 600;
double YungDiagramHandler::s_alpha = 0.16;
double YungDiagramHandler::s_gamma = -0.5;

YungDiagram::YungDiagram() : m_cellsNumber(1), m_numberIsCounted(true), 
  m_number(1), m_probability(1), m_ancestorsNumber(0), m_ancestors(0), m_ancestorsAreCounted(false), m_ancestorsColDifferent(0)
{
  m_cols.push_back(1);
}

YungDiagram::~YungDiagram()
{
  delete[] m_ancestors;
}

YungDiagram::YungDiagram(const string &fileName) : m_cellsNumber(0), m_numberIsCounted(false), 
  m_number(0), m_probability(0), m_ancestorsNumber(0), m_ancestors(0), m_ancestorsColDifferent(0)
{
  std::ifstream in(fileName);
  size_t colsNumber = 0;
  in >> m_cellsNumber;
  in >> colsNumber;

  m_cols.resize(colsNumber);

  in >> m_cols[0];
  for (size_t i = 1; i < colsNumber; ++i)
  {
    in >> m_cols[i];
    if (m_cols[i] > m_cols[i - 1])
    {
      std::cout << "Input file contains incorrect diagram.\n";
    }
  }
}

YungDiagram::YungDiagram(const boost::xint::integer &number) :  m_cellsNumber(0), m_numberIsCounted(true), m_number(0),
  m_probability(0), m_ancestorsNumber(0), m_ancestors(0), m_ancestorsColDifferent(0)
{
  if (!YungDiagramHandler::isPartitionsAmountCounted())
    YungDiagramHandler::countPartitionsAmount();

  const boost::xint::integer *partitionsAmount = YungDiagramHandler::getPartitionsAmount();
  const long long *offsets = YungDiagramHandler::getOffsets();
  size_t maxCellsNumber = YungDiagramHandler::getMaxCellsNumber();

  for (size_t i = 0; i < maxCellsNumber; ++i)
  {
    m_number += partitionsAmount[offsets[i] + i];
    if (m_number >= number)
    {
      m_number -= partitionsAmount[offsets[i] + i];
      m_cellsNumber = i + 1;
      break;
    }
  }
  m_number += 1;

  int n = m_cellsNumber;
  int colSize = n;

  while (n)
  {
    if (colSize == 1)
      break;
    m_number += partitionsAmount[offsets[n - 1] + colSize - 2]; // with max k'th col less than current
    if (m_number > number)
    {
      m_number -= partitionsAmount[offsets[n - 1] + colSize - 2]; // O(n)
      colSize--;
    }
    else
    {
      m_cols.push_back(colSize);
      n -= colSize;
      colSize = min(n, colSize);
    }
  }

  for (int i = 0; i < n; i++)
    m_cols.push_back(1);
}

void YungDiagramHandler::countPartitionsAmount()
{
  int size = (maxCellsNumber + 1) * maxCellsNumber >> 1;
  partitionsAmount = new boost::xint::integer[size];

  offsets = new long long[maxCellsNumber];
  for (long long i = 0; i < maxCellsNumber; i++)
    offsets[i] = i * (i + 1) >> 1;

  partitionsAmount[0] = 1;
  //std::cout << partitionsAmount[0] << " ";
  // P(n, k) is number of partitions of n bocks with biggest not more than k
  for (int n = 1; n < maxCellsNumber; n++)
  {
    long long currentOffset = offsets[n];
    partitionsAmount[currentOffset] = 1; // P(n, 1) 
    //std::cout << partitionsAmount[currentOffset] << " ";
    for (int k = 1; k < n; k++)
    {
      partitionsAmount[currentOffset + k] = partitionsAmount[currentOffset + k - 1];
      if (n > k << 1)
        partitionsAmount[currentOffset + k] += partitionsAmount[offsets[n - k - 1] + k];
      else
        partitionsAmount[currentOffset + k] += partitionsAmount[offsets[n - k - 1] + n - k - 1];
      //std::cout << partitionsAmount[currentOffset + k] << " ";
    }
    partitionsAmount[currentOffset + n] = partitionsAmount[currentOffset + n - 1] + 1;
   // std::cout << partitionsAmount[currentOffset + n] << std::endl;
  }
}

void YungDiagramHandler::PrintPartitionsAmount(const string &fileName)
{
  ofstream out(fileName);
  if (!YungDiagramHandler::isPartitionsAmountCounted())
    YungDiagramHandler::countPartitionsAmount();
  
  for (int i = 0, ctr = 0; i < maxCellsNumber; i++)
  {
    for (int j = 0; j <= i; j++, ctr++)
    {
      out << partitionsAmount[ctr] << " ";
    }
    out << endl;
  }
  out.close();
}

boost::xint::integer YungDiagram::GetDiagramNumber() const
{
  if (m_numberIsCounted)
    return m_number;

  m_number = 0;
  if (!YungDiagramHandler::isPartitionsAmountCounted())
    YungDiagramHandler::countPartitionsAmount();

  const boost::xint::integer *partitionsAmount = YungDiagramHandler::getPartitionsAmount();
  const long long *offsets = YungDiagramHandler::getOffsets();

  for (size_t i = 0; i < m_cellsNumber - 1; ++i)
    m_number += partitionsAmount[offsets[i] + i];
  m_number += 1;

  int n = m_cellsNumber;
  for (size_t col : m_cols)
  {
    if (col == 1)
      break;
    m_number += partitionsAmount[offsets[n - 1] + col - 2]; // with max k'th col less than current
    n -= col;
  }

  m_numberIsCounted = true;
  return m_number;
}

double YungDiagram::GetMyProbabilityRichardson()
{
  long double *prevLayerProbs = new long double[1];
  long double *curLayerProbs = new long double[2];
  prevLayerProbs[0] = 1;
  size_t countedLayer = 1;
  while (countedLayer < m_cellsNumber) 
  {

    countedLayer++;
  }
  return 0;
}

void YungDiagram::SaveToFile(const char *fileName) const
{
  std::ofstream out(fileName);

  out << m_cellsNumber << std::endl;
  out << m_cols.size() << std::endl;
  for (size_t col : m_cols)
    out << col << " ";
}

bool YungDiagram::addCell(size_t col)
{
  if (col > m_cols.size())
    return false;

  if (col < m_cols.size())
    m_cols[col]++;
  else
    m_cols.push_back(1);

  m_ancestorsAreCounted = false;
  m_ancestorsNumber = 0;
  m_cellsNumber++;
  m_ancestorsColDifferent.clear();
  m_numberIsCounted = false;

  return true;
}

void YungDiagram::countAncestors()
{
  getAncestorsNumber();
  size_t m_colsNumber = m_cols.size();

  m_ancestors = new size_t[m_ancestorsNumber];

  size_t beforeLast = m_ancestorsColDifferent.size() - 1;
  for (size_t index = 0; index < beforeLast; index++)
  {
    size_t difCols = m_ancestorsColDifferent[index];
    m_cols[difCols]++;
    m_ancestors[index] = YungDiagramHandler::GetSmallDiagramNumber(m_cellsNumber + 1, m_cols);
    m_cols[difCols]--;
  }

  m_cols.push_back(1);
  m_ancestors[m_ancestorsNumber - 1] = YungDiagramHandler::GetSmallDiagramNumber(m_cellsNumber + 1, m_cols);
  m_cols.pop_back();

  m_ancestorsAreCounted = true;
}

size_t YungDiagramHandler::GetSmallDiagramNumber(size_t cellsNumber, const vector<size_t> &cols, bool withExtraCell)
{
  size_t number = 0;
  size_t colsNumber = cols.size();

  if (!partitionsAmount)
  {
    countPartitionsAmount();
  }

  for (size_t i = 0; i < cellsNumber - 1; ++i)
    number += partitionsAmount[offsets[i] + i]._get_digit(0);
  number += 1;

  int n = cellsNumber;
  for (size_t k = 0; k < colsNumber; k++)
  {
    int colSize = cols[k];
    if (colSize == 1)
      break;
    number += partitionsAmount[offsets[n - 1] + colSize - 2]._get_digit(0); // with max k'th col less than current
    n -= colSize;
  }

  return number;
}

/// returns number of first Digram with (n + 1) cell
boost::xint::integer YungDiagramHandler::GetFirstNumberWithNPlusOneCells(size_t n)
{
  if (n >= maxCellsNumber)
  {
    cout << "Too big number of cells te get diagram number.\n";
    return 0;
  }

  YungDiagram d;
  d.m_cols.resize(n + 1);
  d.m_cellsNumber = n + 1;
  for (size_t i = 0; i < n + 1; i++)
    d.m_cols[i] = 1;
  d.m_numberIsCounted = false;
  return d.GetDiagramNumber();
}

void YungDiagram::setColsNumber(size_t colsNumber)
{
  if (m_ancestorsAreCounted)
  {
    delete[] m_ancestors;
    m_probability = 0;
    m_ancestorsAreCounted = false;
  }
  m_cols.resize(colsNumber);
  for (size_t &col : m_cols) // Test how it works !!!!!!!!!!!!!!!!!!!!!!!!!
    col = 0;
}

void YungDiagram::setCellsInCol(size_t colIndex, size_t cellsNumber)
{
  if (m_ancestorsAreCounted)
  {
    delete[] m_ancestors;
    m_probability = 0;
    m_ancestorsAreCounted = false;
  }

  if (colIndex >= m_cols.size())
  {
    cout << "Invalid col index in SetCellsInCol\n";
    return;
  }

  if (colIndex > 0)
  {
    if (m_cols[colIndex - 1] < cellsNumber)
    {
      cout << "Invalid cells number in setCellsInCol: column on the left of changed column can't have more cells.\n";
      return;
    }
  }

  if (colIndex < m_cols.size() - 1)
  {
    if (m_cols[colIndex + 1] > cellsNumber)
    {
      cout << "Invalid cells number in setCellsInCol: column on the right of changed column can't have less cells.\n";
      return;
    }
  }
  size_t prevCellsNumber = m_cols[colIndex];
  m_cols[colIndex] = cellsNumber;
  m_cellsNumber += (cellsNumber - prevCellsNumber);
}

void YungDiagram::resetAncestors()
{
  delete[] m_ancestors;

  m_ancestors = 0;
  m_ancestorsNumber = 0;
  m_ancestorsColDifferent.clear();
}

void YungDiagramHandler::CountProbabilities(ProcessType processType, size_t cellsNumber)
{
  YungDiagram *currentLevelDiagrams = 0;
  YungDiagram *nextLevelDiagrams = new YungDiagram[1];
  nextLevelDiagrams[0].SetMyProbability(1);

  size_t index = 0;
  levelSize = 0;
  size_t nextLevelSize = 1;
  size_t firstDiagramOnNextLevel = 0;
  size_t firstDiagramOnNextNextLevel  = 2;
  double cur_prob = 1;
  for (size_t level = 1; level < cellsNumber; level++)
  {
    currentLevelDiagrams = nextLevelDiagrams;
    levelSize = nextLevelSize;
    firstDiagramOnNextLevel = firstDiagramOnNextNextLevel;

    firstDiagramOnNextNextLevel = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level + 1)._get_digit(0);
    nextLevelSize = firstDiagramOnNextNextLevel - firstDiagramOnNextLevel;
    nextLevelDiagrams = new YungDiagram[nextLevelSize];

    for (size_t i = 0; i < levelSize; ++i) 
    {
      currentLevelDiagrams[i].countAncestors();
      size_t ancestorsNumber = currentLevelDiagrams[i].getAncestorsNumber();
      size_t *ancestors = currentLevelDiagrams[i].getAncestors();

      double divisor = 0;
      double *numerator = new double[ancestorsNumber];
      for (size_t j = 0; j < ancestorsNumber; j++) 
      {
        size_t idx = ancestors[j] - firstDiagramOnNextLevel;
        if (nextLevelDiagrams[idx].m_cellsNumber == 1) // initialize it if it isn't yet
        {
          nextLevelDiagrams[idx].m_cellsNumber = currentLevelDiagrams[i].m_cellsNumber + 1;
          nextLevelDiagrams[idx].m_cols = currentLevelDiagrams[i].m_cols;
          if (j == ancestorsNumber - 1)
            nextLevelDiagrams[idx].m_cols.push_back(1);
          else
            nextLevelDiagrams[idx].m_cols[currentLevelDiagrams[i].m_ancestorsColDifferent[j]]++;
        }

        size_t x = currentLevelDiagrams[i].m_ancestorsColDifferent[j];
        size_t y = nextLevelDiagrams[idx].m_cols[x];
        x++;
        switch (processType) 
        {
        case RICHARDSON:
          numerator[j] = 1;
          break;
        case ALPHA:
          numerator[j] = pow((x * x + y * y), s_alpha);
          break;
        case BETA:
          numerator[j] = 1.0 / pow(nextLevelDiagrams[idx].getAncestorsNumber(), 1);
          break;
        }
        divisor += numerator[j];
      }

      for (size_t j = 0; j < ancestorsNumber; j++) 
      {
        size_t idx = ancestors[j] - firstDiagramOnNextLevel;
        nextLevelDiagrams[idx].m_probability += currentLevelDiagrams[i].m_probability * numerator[j] / divisor;
      }
      delete[] numerator;
      currentLevelDiagrams[i].resetAncestors();
    }
    delete[] currentLevelDiagrams;
  }

  probabilities = new double[nextLevelSize];
  numbers = new size_t[nextLevelSize];
  for (size_t i = 0; i < nextLevelSize; i++)
  {
    probabilities[i] = nextLevelDiagrams[i].m_probability;
    numbers[i] = i + firstDiagramOnNextLevel;
  }

  delete[] nextLevelDiagrams;

  levelSize = nextLevelSize;
}

void YungDiagramHandler::SortByProbability()
{
  sort(0, levelSize - 1);
}

void YungDiagramHandler::sort(int l, int r)
{
  if (r <= l)
    return;

  double pivot = probabilities[r];
  int i = l - 1;
  for (int j = l; j < r; j++)
  {
    if (probabilities[j] < pivot)
    {
      i++;
      double temp = probabilities[i];
      probabilities[i] = probabilities[j];
      probabilities[j] = temp;

      size_t tempInd = numbers[i];
      numbers[i] = numbers[j];
      numbers[j] = tempInd;
    }
  }
  double temp = probabilities[i + 1];
  probabilities[i + 1] = probabilities[r];
  probabilities[r] = temp;
  size_t tempInd = numbers[i + 1];
  numbers[i + 1] = numbers[r];
  numbers[r] = tempInd;

  sort(l, i);
  sort(i + 2, r);
}

void YungDiagramHandler::saveProbabilities(const char *fileName)
{
  FILE *f;
  fopen_s(&f, fileName, "w");
  for (size_t i = 0; i < levelSize; i++)
    fprintf(f, "%u %.16lf\n", numbers[i], probabilities[i]);

  fclose(f);
}

size_t YungDiagram::getAncestorsNumber()
{
  if (m_ancestorsNumber != 0)
    return m_ancestorsNumber;
  
  m_ancestorsNumber = 2; // left column always can be increased
  m_ancestorsColDifferent.clear();
  size_t m_colsNumber = m_cols.size();
  m_ancestorsColDifferent.push_back(0);
  for (size_t i = 1; i < m_colsNumber; i++) 
  {
    if (m_cols[i] < m_cols[i - 1])
    {
      m_ancestorsNumber++;
      m_ancestorsColDifferent.push_back(i);
    }
  }
  m_ancestorsColDifferent.push_back(m_colsNumber);

  return m_ancestorsNumber;
}

YungDiagram *YungDiagramHandler::getRandomDiagram(ProcessType procType, size_t n)
{
  YungDiagram *currentDiagram = new YungDiagram();
  static unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(seed);

  double curProb = 1.0;

  switch (procType)
  {
  case RICHARDSON:
    for (size_t i = 1; i < n; i++)
    {
      size_t ancestorsNumber = currentDiagram->getAncestorsNumber();
      std::uniform_int_distribution<int> distribution(0, ancestorsNumber - 1);
      size_t randomAncestor = distribution(generator);
      currentDiagram->addCell(currentDiagram->m_ancestorsColDifferent[randomAncestor]);
    }
    break;
  case BETA:
    for (size_t i = 1; i < n; i++)
    {
      size_t ancestorsNumber = currentDiagram->getAncestorsNumber();
      vector<double> probs;
       size_t colsNumber = currentDiagram->m_cols.size();
      for (size_t j = 0; j < ancestorsNumber - 1; j++)
      {
        size_t curCol = currentDiagram->m_ancestorsColDifferent[j];
        size_t newAncestorsNumber = ancestorsNumber;
        if (curCol == 0)
        {
          if (colsNumber > 1)
          {
            if (currentDiagram->m_cols[1] == currentDiagram->m_cols[0])
              newAncestorsNumber++;
          }
        }
        else
        {
          if (colsNumber > j + 1)
          {
            if (currentDiagram->m_cols[j] == currentDiagram->m_cols[j + 1])
              newAncestorsNumber++;
          }
          if (currentDiagram->m_cols[j] + 1 == currentDiagram->m_cols[j - 1])
            newAncestorsNumber--;
        }
        probs.push_back(1.0 / pow(newAncestorsNumber, 2));
      }

      size_t newAncestorsNumber = ancestorsNumber;
      if (currentDiagram->m_cols[colsNumber - 1] > 1)
        newAncestorsNumber++;

      probs.push_back(1.0 / pow(newAncestorsNumber, 2));
      size_t j = 0;
      std::discrete_distribution<size_t> distrib(probs.size(), 0, ancestorsNumber - 1,
          [&probs, &j](float){
          auto w = probs[j];
          ++j;
          return w;
      });
      size_t randomAncestor = distrib(generator);
      currentDiagram->addCell(currentDiagram->m_ancestorsColDifferent[randomAncestor]);
    }
    break;
  case ALPHA:
    for (size_t i = 1; i < n; i++)
    {
      size_t ancestorsNumber = currentDiagram->getAncestorsNumber();
      vector<double> probs;
      for (size_t j = 0; j < ancestorsNumber - 1; j++)
      {
        size_t x = currentDiagram->m_ancestorsColDifferent[j];
        size_t y = currentDiagram->m_cols[x] + 1;
        x++; // indexin from 0, but we need from one here
        probs.push_back(pow((x * x + y * y), s_alpha));
      }
      size_t x = currentDiagram->m_cols.size() + 1;
    
      probs.push_back(pow((x * x + 1), s_alpha));
      size_t j = 0;
      std::discrete_distribution<size_t> distrib(probs.size(), 0, ancestorsNumber - 1,
          [&probs, &j](double){
          auto w = probs[j];
          ++j;
          return w;
      });
      size_t randomAncestor = distrib(generator);
      currentDiagram->addCell(currentDiagram->m_ancestorsColDifferent[randomAncestor]);
    }
    break;

  case PLANSHEREL:
    for (size_t i = 1; i < n; i++)
    {
      size_t ancestorsNumber = currentDiagram->getAncestorsNumber();
      vector<double> probs;
      for (size_t j = 0; j < ancestorsNumber; j++)
      {
        size_t x = currentDiagram->m_ancestorsColDifferent[j] + 1;
        size_t y = 1;

        if (j < ancestorsNumber - 1)
          y = currentDiagram->m_cols[x - 1] + 1;

        double newProb = 1.0;
        size_t height = 0; // here probably 1
        size_t width = 0;
        size_t rightColWithHeightCells = currentDiagram->m_cols.size() - 1;
        
        for (size_t k = y - 1; k > 0; k--)
        {
          int dx = rightColWithHeightCells - x + 1;
          newProb *= (k + dx);
          newProb /= (k + dx + 1);
          height++;
          while (currentDiagram->m_cols[rightColWithHeightCells] < height)
            rightColWithHeightCells--;
        }
        for (size_t k = x - 1; k > 0; k--)
        {
          int dy = currentDiagram->m_cols[width] - y;
          newProb *= (k + dy);
          newProb /= (k + dy + 1);
          width++;
        }
        
        newProb *= (i + 1); // here in fact we don't need *(i + 1) for it's the same for all ancestors
        probs.push_back(newProb);
      }
      size_t j = 0;
      std::discrete_distribution<size_t> distrib(probs.size(), 0, ancestorsNumber - 1,
          [&probs, &j](double){
          auto w = probs[j];
          ++j;
          return w;
      });
      size_t randomAncestor = distrib(generator);
      currentDiagram->addCell(currentDiagram->m_ancestorsColDifferent[randomAncestor]);
    }
    break;

  case PLANSHEREL_POWERED_GAMMA:
    for (size_t i = 1; i < n; i++)
    {
      size_t ancestorsNumber = currentDiagram->getAncestorsNumber();
      vector<double> probs;
      for (size_t j = 0; j < ancestorsNumber; j++)
      {
        size_t x = currentDiagram->m_ancestorsColDifferent[j] + 1;
        size_t y = 1;

        if (j < ancestorsNumber - 1)
          y = currentDiagram->m_cols[x - 1] + 1;

        double newProb = 1.0;
        size_t height = 0;
        size_t width = 0;
        size_t rightColWithHeightCells = currentDiagram->m_cols.size() - 1;
        
        for (size_t k = y - 1; k > 0; k--) // we go up by y, so k = height(x) in the beginning and then it goes down till 1; 
        {
          int dx = rightColWithHeightCells - x + 1;
          newProb *= (k + dx);
          newProb /= (k + dx + 1);
          height++;
          while (currentDiagram->m_cols[rightColWithHeightCells] < height)
            rightColWithHeightCells--;
        }
        for (size_t k = x - 1; k > 0; k--)
        {
          int dy = currentDiagram->m_cols[width] - y;
          newProb *= (k + dy);
          newProb /= (k + dy + 1);
          width++;
        }
        
        newProb *= (newProb * (i + 1));
        probs.push_back(pow(newProb, s_gamma));
      }
      size_t j = 0;
      std::discrete_distribution<size_t> distrib(probs.size(), 0, ancestorsNumber - 1,
          [&probs, &j](float){
          auto w = probs[j];
          ++j;
          return w;
      });
      size_t randomAncestor = distrib(generator);
      currentDiagram->addCell(currentDiagram->m_ancestorsColDifferent[randomAncestor]);
    }
    break;
  default:
    cerr << "ERROR: Unknown process type sent to handler random diagram generator.\n";
    break;
  }

  return currentDiagram;
}

vector<size_t> YungDiagramHandler::getRandomWalkFrequencies(ProcessType processType, size_t cellsNumber, size_t bucketsNumber, size_t testsNumber, bool needPermutation)
{
  boost::xint::integer *leftBounds = new boost::xint::integer[bucketsNumber];
  boost::xint::integer maxNumberPlusOne = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(cellsNumber);
  boost::xint::integer minNumber = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(cellsNumber - 1);
  
  leftBounds[0] = minNumber;
  boost::xint::integer intervalSize = (maxNumberPlusOne - minNumber) / bucketsNumber; // may be problems with division: last bucket may become very samll (down to 1 diagram)
  cout << "Min number is " << minNumber << endl;
  cout << "Max number is " << maxNumberPlusOne << endl;
  cout << "Delta is " << (maxNumberPlusOne - minNumber) << endl;
  cout << "bucket interval is " << intervalSize << endl;

  for (size_t i = 1; i < bucketsNumber; i++)
    leftBounds[i] = leftBounds[i - 1] + intervalSize;
  vector<size_t> buckets(bucketsNumber, 0);
  vector<boost::xint::integer> generatedDiagrams;

  for (size_t i = 0; i < testsNumber; i++)
  {
    if (i % 1000 == 0)
      cout << i << " processed.\n";
    YungDiagram *d = YungDiagramHandler::getRandomDiagram(processType, cellsNumber);
    boost::xint::integer number = d->GetDiagramNumber();
    generatedDiagrams.push_back(number);
    delete d;
  }

  if (needPermutation)
  {
    boost::xint::default_random_generator gen;
    size_t minLen = boost::xint::highestbit(minNumber);
    size_t maxLen = boost::xint::highestbit(maxNumberPlusOne);
    size_t deltaLength = maxLen - minLen;
    unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distribution(0, deltaLength);
    set<boost::xint::integer> newNumbers;

    for (size_t i = 0; i < testsNumber; i++)
    {
      size_t bitLen = distribution(generator) + minLen;
      boost::xint::integer num = 0;
      do 
      {
        num = boost::xint::integer::random_by_size(gen, bitLen, true);
      } while (num < minNumber || num >= maxNumberPlusOne || !newNumbers.insert(num).second);
    }
    map<boost::xint::integer, boost::xint::integer> permutation;
    set<boost::xint::integer>::iterator newNumIter = newNumbers.begin();

    size_t reps = 0;
    for (boost::xint::integer num : generatedDiagrams)
    {
      if (permutation.find(num) == permutation.end())
      {
        permutation.insert(std::pair<boost::xint::integer,boost::xint::integer>(num, *newNumIter));
        newNumIter++;
      }
      else
        reps++;
    }
    cout << "There were " << reps << " repeatings.\n";

    for (size_t i = 0; i < testsNumber; i++)
    {
      boost::xint::integer number = permutation.at(generatedDiagrams[i]);

      int p = 0;
      int q = bucketsNumber - 1;
      int m = (p + q + 1) >> 1;
      while (p < q)
      {
        number < leftBounds[m] ? q = (m - 1) : p = m;
        m = (p + q + 1) >> 1;
      }
      buckets[p]++;
    }
  }
  else // no need to make permutation
  {
    map<boost::xint::integer, size_t> freqsOfDiagrams; // lets find out how many same diagrams we have got
    size_t maxSameDiagramsNumber = 1;
    for (size_t i = 0; i < testsNumber; i++)
    {
      boost::xint::integer number = generatedDiagrams[i];

      map<boost::xint::integer, size_t>::iterator it = freqsOfDiagrams.find(number);
      if (it == freqsOfDiagrams.end())
        freqsOfDiagrams.insert(std::pair<boost::xint::integer, size_t>(number, 1));
      else
      {
        it->second++;
        maxSameDiagramsNumber = it->second;
      }

      int p = 0;
      int q = bucketsNumber - 1;
      int m = (p + q + 1) >> 1;
      while (p < q)
      {
        number < leftBounds[m] ? q = (m - 1) : p = m;
        m = (p + q + 1) >> 1;
      }
      buckets[p]++;
    }

    cout << "Maximum amount of same generated diagrams is " << maxSameDiagramsNumber << endl;
    cout << "There were generated " << freqsOfDiagrams.size() << " different diagrams.\n";
  }
  delete[] leftBounds;

  cout << "Finished random walking.\n";
  return buckets;
}

void YungDiagramHandler::saveColumnsToFile(const char* fileName, YungDiagram *d)
{
  ofstream out(fileName);

  for (size_t col : d->m_cols)
    out << col << " ";
  out.close();
}

void YungDiagramHandler::countCoefficientsForKantorovichMetric(size_t diagramNumber, dVector &coefficients, vector<size_t> &predecessorsNums)
{
  YungDiagram d(diagramNumber);
  size_t colsNumber = d.m_cols.size();

  vector<double> a; // c[i] = a[i] / sum(a)
  double sum_coefs = 0;
  coefficients.resize(colsNumber);

  size_t coefsIndex = 0;
  for (int s = 0; s < colsNumber; s++)
  {
    bool cellCanBeRemoved = true;

    if (s != colsNumber - 1)
    {
      if (d.m_cols[s] == d.m_cols[s + 1])
        cellCanBeRemoved = false;
    }

    if (cellCanBeRemoved)
    {
      size_t x = s + 1;
      size_t y = d.m_cols[s];

      double c = 1.0;
      size_t height = 1;
      size_t width = 0;
      size_t rightColWithHeightCells = colsNumber - 1;
          
      for (size_t l = y; l > 0; l--) // we go up by y, so k = height(x) in the beginning and then it goes down till 1; 
      {
        int dx = rightColWithHeightCells - x + 1;
        if (l + dx - 1)
        {
          c *= (l + dx);
          c /= (l + dx - 1);
          height++;
          while (d.m_cols[rightColWithHeightCells] < height)
            rightColWithHeightCells--;
        }
      }
      for (size_t l = x; l > 1; l--)
      {
        int dy = d.m_cols[width] - y;
        if (l + dy - 1)
        {
          c *= (l + dy);
          c /= (l + dy - 1);
          width++;
        }
      }
        
      //c *= c;
      sum_coefs += c;
      coefficients(coefsIndex++) = c;

      d.m_cols[s]--;
      bool needResize = false;
      if (d.m_cols[s] == 0)
      {
        d.m_cols.resize(d.m_cols.size() - 1);
        needResize = true;
      }
      predecessorsNums.push_back(YungDiagramHandler::GetSmallDiagramNumber(d.m_cellsNumber - 1, d.m_cols));

      if (needResize)
        d.m_cols.push_back(0);
      d.m_cols[s]++;
    }
  }

  coefficients.resize(predecessorsNums.size());
  for (int i = 0; i < coefficients.size(); i++)
    coefficients[i] /= sum_coefs;
}

// for small diagrams, who's number in size_t
double YungDiagramHandler::countKantorovichDistance(YungDiagram &d1, YungDiagram &d2)
{
  vector<size_t> neededDiagrams;
  vector<size_t> numbersOfNeededDiagramsOnLevel;
  vector<char> ancestorFlags; 
  // probably we don't need to count distance between diagrams that are predecessors of the same diagram from
  // te top level. If it's predecessor of first one, it has first bit to be 1, if of the second - then second bit to be 1.
  // BUT, IF it turns out that this diagram is predecessor for both of top level diagrams,
  // we should count it's flag to be 3 - bits for both

  neededDiagrams.push_back(d1.GetDiagramNumber()._get_digit(0));
  neededDiagrams.push_back(d2.GetDiagramNumber()._get_digit(0));
  ancestorFlags.push_back(1);
  ancestorFlags.push_back(2);

  numbersOfNeededDiagramsOnLevel.push_back(2);

  int n = d1.m_cellsNumber;
  int startIndexInVector = 0; // indicate start and end of current level diagrams in neededDiagrams
  int endIndexInVector = 1;

  for (int i = n - 1; i > 2; i--) // we assume that we have metric on the second floor
  {
    size_t neededDiagramsNumberOnLevel = 0;
    for (int j = startIndexInVector; j <= endIndexInVector; j++)
    {
      YungDiagram d(neededDiagrams[j]);
      size_t colsNumber = d.m_cols.size();
      d.m_cellsNumber--;
      for (int k = 0; k < colsNumber; k++)
      {
        bool cellCanBeRemoved = true;

        if (k != colsNumber - 1)
        {
          if (d.m_cols[k] == d.m_cols[k + 1])
            cellCanBeRemoved = false;
        }

        if (cellCanBeRemoved)
        {
          d.m_cols[k]--;
          bool needResize = false;
          if (d.m_cols[k] == 0)
          {
            d.m_cols.resize(colsNumber - 1);
            needResize = true;
          }
          size_t num = GetSmallDiagramNumber(d.m_cellsNumber, d.m_cols);
          // check if we haven't added thid diagram as needed already
          vector<size_t>::iterator iter;
          if ((iter = std::find(neededDiagrams.begin(), neededDiagrams.end(), num)) == neededDiagrams.end())
          {
            if (num > 3)
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

          if (needResize)
            d.m_cols.push_back(0);
          d.m_cols[k]++;
        }
      }
      d.m_cellsNumber++; // for proper deletion
    }

    numbersOfNeededDiagramsOnLevel.push_back(neededDiagramsNumberOnLevel);
    startIndexInVector = endIndexInVector + 1;
    endIndexInVector = neededDiagrams.size() - 1;
  }

  map<num_pair, double> distances; // key is pair of diagram numbers, first number is less than second.

  distances.insert(pair<num_pair, double>(pair<size_t, size_t>(2, 3), 1.0));
  endIndexInVector = neededDiagrams.size() - 1;
  size_t solvedTranspProblems = 0;
  for (int i = 3, levelIndex = numbersOfNeededDiagramsOnLevel.size() - 1; i <= n; i++, levelIndex--) 
  {
    startIndexInVector = endIndexInVector - numbersOfNeededDiagramsOnLevel[levelIndex] + 1;

    // last is not needed for it's distances will be counted during previous steps
    for (int j = startIndexInVector; j < endIndexInVector; j++)
    {
      YungDiagram d1(neededDiagrams[j]);

      size_t colsNumber1 = d1.m_cols.size();
      vector<size_t> nums1; // numbers of predecessors of d1
      dVector c1; 
      countCoefficientsForKantorovichMetric(neededDiagrams[j], c1, nums1);

      for (int k = j + 1; k <= endIndexInVector; k++)
      {
        if (ancestorFlags[k] == ancestorFlags[j] && (ancestorFlags[k] ^ 3)) // boost!? yes, about 2 times
          continue;

        dVector c2;
        vector<size_t> nums2; // numbers of predecessors of d2
        countCoefficientsForKantorovichMetric(neededDiagrams[k], c2, nums2);
        
        // now we need to solve transportation problem
        size_t s1 = c1.size();
        size_t s2 = c2.size();
        dMatrix costs(s1, s2);
        dMatrix x(s1, s2);
        for (size_t l = 0; l < s1; l++)
        {
          for (size_t p = 0; p < s2; p++)
          {
            size_t d1 = nums1[l];
            size_t d2 = nums2[p];

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
    endIndexInVector = startIndexInVector - 1;
  }

  cout << "During distance finding were solved " << solvedTranspProblems << " transportation problems.\n";

  if (neededDiagrams[0] < neededDiagrams[1])
    return distances[num_pair(neededDiagrams[0], neededDiagrams[1])];
  return distances[num_pair(neededDiagrams[1], neededDiagrams[0])];
}

double YungDiagramHandler::countKantorovichDistance(size_t n1, size_t n2)
{
  return (countKantorovichDistance(YungDiagram(n1), YungDiagram(n2)));
}

size_t YungDiagramHandler::GetSmallDiagramNumber(size_t colsNumber, ...)
{
  va_list vl;
  va_start(vl, colsNumber);
  YungDiagram d;
  d.m_cols.resize(colsNumber);
  d.m_cellsNumber = 0;
  d.m_numberIsCounted = false;
  for (size_t i = 0; i < colsNumber; i++)
  {
    size_t col = va_arg(vl, size_t);
    d.m_cols[i] = col;
    d.m_cellsNumber += col;
  }

  va_end(vl);
  std::sort(d.m_cols.begin(), d.m_cols.end(), std::greater<size_t>());
  return d.GetDiagramNumber()._get_digit(0);
}

boost::xint::integer YungDiagramHandler::GetFirstNumberWithNCells(size_t n)
{
  return GetFirstNumberWithNPlusOneCells(n - 1);
}

boost::xint::integer YungDiagramHandler::GetLastNumberWithNCells(size_t n)
{
  return GetFirstNumberWithNPlusOneCells(n) - 1;
}

void YungDiagramHandler::printSmallDiagramsPair(size_t n1, size_t n2)
{
  YungDiagram d1(n1);
  YungDiagram d2(n2);

  size_t h = std::max<>(d1.m_cols[0], d2.m_cols[0]);
  size_t w1 = d1.m_cols.size();
  size_t w2 = d2.m_cols.size();

  for (int i = h; i > 0; i--)
  {
    for (int j = 0; j < w1; j++)
    {
      if (d1.m_cols[j] >= i)
        cout << "*";
      else
        cout << " ";
    }
    cout << "  ";
    for (int j = 0; j < w2; j++)
    {
      if (d2.m_cols[j] >= i)
        cout << "*";
      else
        cout << " ";
    }

    cout << endl;
  }
}

void YungDiagramHandler::printSmallDiagram(size_t n)
{
  YungDiagram d(n);

  size_t h = d.m_cols[0];
  size_t w = d.m_cols.size();

  cout << endl;
  for (int i = h; i > 0; i--)
  {
    for (int j = 0; j < w; j++)
    {
      if (d.m_cols[j] >= i)
        cout << "*";
      else
        cout << " ";
    }
    cout << endl;
  }
  
}

void YungDiagramHandler::getBall(size_t num, double r, vector<size_t> &ballDiagrams, vector<double> &ballDistances)
{
  YungDiagram d(num);
  size_t cellsNum = d.getCellsNumber();
  size_t n1 = GetFirstNumberWithNCells(cellsNum)._get_digit(0);
  size_t n2 = GetLastNumberWithNCells(cellsNum)._get_digit(0);

  for (size_t i = n1; i <= n2; i++)
  {
    double dist = countKantorovichDistance(num, i);
    if (dist <= r)
    {
      ballDistances.push_back(dist);
      ballDiagrams.push_back(i);
    }
  }
}

void YungDiagramHandler::getDistancesToSomeDiagrams(size_t num, size_t diagramsNumber, 
                                                    std::vector<size_t> &diagrams, std::vector<double> &distances)
{
  YungDiagram d(num);
  size_t cellsNum = d.getCellsNumber();
  size_t n1 = GetFirstNumberWithNCells(cellsNum)._get_digit(0);
  size_t n2 = GetLastNumberWithNCells(cellsNum)._get_digit(0);

  double delta = 1;
  size_t maxOffset = n2 - n1;
  if (diagramsNumber < maxOffset)
    delta = double(maxOffset) / diagramsNumber;

  for (size_t offset = 0, i = 0; offset <= maxOffset; )
  {
    size_t number = n1 + offset;
    double dist = countKantorovichDistance(num, number);
    cout << "dist = " << dist << endl;
    distances.push_back(dist);
    diagrams.push_back(number);
    ++i;
    offset = size_t(i * delta);
  }
}

// for small diagrams!
void YungDiagramHandler::getLinearCoefficientsEstimationKantorovich(size_t diagramNum, dVector &c, dMatrix &delta, size_t checkedNeighbours)
{
  vector<size_t> diagrams;
  vector<double> dists;

  getLinearCoefficientsEstimationKantorovich(diagramNum, c, delta, checkedNeighbours, diagrams, dists);
}

void YungDiagramHandler::getLinearCoefficientsEstimationKantorovich(size_t diagramNum, dVector &c, dMatrix &delta, 
    size_t checkedNeighbours, std::vector<size_t> &diagrams, std::vector<double> &dists)
{
    // just get all distances between given diagram and all other on the level
  //getBall(diagramNum, 1.1, diagrams, dists);
  getDistancesToSomeDiagrams(diagramNum, checkedNeighbours, diagrams, dists);
  YungDiagram d(diagramNum);
  size_t cellsNum = d.m_cellsNumber;
  dVector ro(diagrams.size());

  cout << "diagrams.size() = " << diagrams.size() << endl;
  for (size_t i = 0; i < diagrams.size(); i++)
    ro(i) = dists[i];
  
  vector<int> rows(cellsNum);
  size_t h = 1;
  for (size_t i = d.m_cols.size(); i > 0; i--)
  {
    while (d.m_cols[i - 1] >= h)
    {
      rows[h - 1] = i;
      h++;
    }
  }

  delta.resize(diagrams.size(), cellsNum);
  for (size_t ind = 0; ind < diagrams.size(); ind++)
  {
    YungDiagram d1(diagrams[ind]);
    h = 1;
    for (int j = d1.m_cols.size(); j > 0; j--) // here int is a must
    {
      while (d1.m_cols[j - 1] >= h)
      {
        delta(ind, h - 1) = fabs((double)rows[h - 1] - j);
        //delta(ind, h - 1) = rows[h - 1] - j;
        h++;
      }
    }
    for (size_t j = h - 1; j < rows.size(); j++)
      delta(ind, j) = rows[j];
  }


  dMatrix dd = boost::numeric::ublas::prod(boost::numeric::ublas::trans(delta), delta);
  dVector b =  boost::numeric::ublas::prod(boost::numeric::ublas::trans(ro), delta);

  cout << "dd.size = " << dd.size1() << " " << dd.size2() << endl;
  cout << "b.size = " << b.size() << endl;
  c = gausSolve(dd, b);
}
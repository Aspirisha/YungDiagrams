#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <array>
#include <string>
#include <map>
#include <set>
#include "YungDiagram.h"

using namespace std;
boost::xint::integer *YungDiagramHandler::partitionsAmount = 0;
long long *YungDiagramHandler::offsets = 0;
double *YungDiagramHandler::probabilities = 0;
size_t *YungDiagramHandler::numbers = 0;
size_t YungDiagramHandler::levelSize = 0;
const int YungDiagramHandler::maxCellsNumber = 600;
double YungDiagramHandler::s_alpha = 0.3;

YungDiagram::YungDiagram() : m_cellsNumber(0), m_numberIsCounted(false), 
  m_number(0), m_probability(0), m_ancestorsNumber(0), m_ancestors(0), m_ancestorsAreCounted(false), m_ancestorsColDifferent(0)
{

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
    if (m_number > number)
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
  for (int n = 1; n < maxCellsNumber; n++)
  {
    long long currentOffset = offsets[n];
    partitionsAmount[currentOffset] = 1;
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
  m_cols[0]++;
  m_ancestors[0] = YungDiagramHandler::GetSmallDiagramNumber(m_cellsNumber + 1, m_cols) - 1;
  m_cols[0]--;

  int index = 1;
  for (size_t difCols : m_ancestorsColDifferent)
  {
    m_cols[difCols]++;
    m_ancestors[index++] = YungDiagramHandler::GetSmallDiagramNumber(m_cellsNumber + 1, m_cols) - 1;
    m_cols[difCols]--;
  }

  m_ancestors[m_ancestorsNumber - 1] = YungDiagramHandler::GetSmallDiagramNumber(m_cellsNumber + 1, m_cols, true) - 1;

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
    if (withExtraCell && k == colsNumber - 1)
      break;

    int colSize = cols[k];
    if (colSize == 1)
      break;
    number += partitionsAmount[offsets[n - 1] + colSize - 2]._get_digit(0); // with max k'th col less than current
    n -= colSize;
  }

  return number;
}

/// returns number of first Digram with (n + 1) cell
boost::xint::integer YungDiagramHandler::GetMaxNumberWithNCells(size_t n)
{
  if (n >= maxCellsNumber)
  {
    cout << "Too big number of cells te get diagram number.\n";
    return 0;
  }

  YungDiagram d;
  d.setColsNumber(n + 1);
  for (size_t i = 0; i < n + 1; i++)
    d.setCellsInCol(i, 1);
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

void YungDiagramHandler::CountRichardsonProbabilities(size_t cellsNumber)
{
  YungDiagram *currentLevelDiagrams = 0;
  YungDiagram *nextLevelDiagrams = new YungDiagram[1];
  nextLevelDiagrams[0].SetMyProbability(1);
  nextLevelDiagrams[0].setColsNumber(1);
  nextLevelDiagrams[0].setCellsInCol(0, 1);

  size_t index = 0;
  levelSize = 0;
  size_t nextLevelSize = 1;
  size_t firstDiagramOnNextLevel = 0;
  size_t firstDiagramOnNextNextLevel  = 1;
  for (size_t level = 1; level < cellsNumber; level++)
  {
    currentLevelDiagrams = nextLevelDiagrams;
    levelSize = nextLevelSize;
    firstDiagramOnNextLevel = firstDiagramOnNextNextLevel;

    firstDiagramOnNextNextLevel = YungDiagramHandler::GetMaxNumberWithNCells(level + 1)._get_digit(0) - 1;
    nextLevelSize = firstDiagramOnNextNextLevel - firstDiagramOnNextLevel;
    nextLevelDiagrams = new YungDiagram[nextLevelSize];

    for (size_t i = 0; i < levelSize; ++i) 
    {
      currentLevelDiagrams[i].countAncestors();
      size_t ancestorsNumber = currentLevelDiagrams[i].getAncestorsNumber();
      size_t *ancestors = currentLevelDiagrams[i].getAncestors();
      long double deltaProb = currentLevelDiagrams[i].m_probability / (long double) ancestorsNumber;

      for (size_t j = 0; j < ancestorsNumber; j++) 
      {
        size_t idx = ancestors[j] - firstDiagramOnNextLevel;
        if (nextLevelDiagrams[idx].m_cellsNumber == 0)
        {
          nextLevelDiagrams[idx].m_cellsNumber = currentLevelDiagrams[i].m_cellsNumber + 1;
          nextLevelDiagrams[idx].m_cols = currentLevelDiagrams[i].m_cols;
          if (j == ancestorsNumber - 1)
            nextLevelDiagrams[idx].m_cols.push_back(1);
          else
            nextLevelDiagrams[idx].m_cols[currentLevelDiagrams[i].m_ancestorsColDifferent[j]]++;
        }
        nextLevelDiagrams[idx].m_probability += deltaProb;
      }
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

void YungDiagramHandler::CountAlphaProbabilities(size_t cellsNumber)
{
  YungDiagram *currentLevelDiagrams = 0;
  YungDiagram *nextLevelDiagrams = new YungDiagram[1];
  nextLevelDiagrams[0].SetMyProbability(1);
  nextLevelDiagrams[0].setColsNumber(1);
  nextLevelDiagrams[0].setCellsInCol(0, 1);

  size_t index = 0;
  levelSize = 0;
  size_t nextLevelSize = 1;
  size_t firstDiagramOnNextLevel = 0;
  size_t firstDiagramOnNextNextLevel  = 1;
  for (size_t level = 1; level < cellsNumber; level++)
  {
    currentLevelDiagrams = nextLevelDiagrams;
    levelSize = nextLevelSize;
    firstDiagramOnNextLevel = firstDiagramOnNextNextLevel;

    firstDiagramOnNextNextLevel = YungDiagramHandler::GetMaxNumberWithNCells(level + 1)._get_digit(0) - 1;
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
        if (nextLevelDiagrams[idx].m_cellsNumber == 0)
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
        numerator[j] = pow((x * x + y * y), s_alpha);
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
    fprintf(f, "%u %.16lf\n", numbers[i] + 1, probabilities[i]);

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
  boost::xint::integer number = 0;
  YungDiagram *currentDiagram = new YungDiagram(number);
  unsigned seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);

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

vector<size_t> YungDiagramHandler::getRandomWalkFrequencies(ProcessType processType, size_t cellsNumber, size_t bucketsNumber, size_t testsNumber)
{
  boost::xint::integer *leftBounds = new boost::xint::integer[bucketsNumber];
  boost::xint::integer maxNumberPlusOne = YungDiagramHandler::GetMaxNumberWithNCells(cellsNumber);
  boost::xint::integer minNumber = YungDiagramHandler::GetMaxNumberWithNCells(cellsNumber - 1);
  
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
  //  int p = 0;
  //  int q = bucketsNumber - 1;
   // int m = (p + q + 1) >> 1;
    //int num = number._get_digit(0);
   // while (p < q)
   // {
    //  number < leftBounds[m] ? q = (m - 1) : p = m;
   //   m = (p + q + 1) >> 1;
  //  }
    //buckets[p]++;
    delete d;
  }

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

  for (boost::xint::integer num : generatedDiagrams)
  {
    if (permutation.find(num) == permutation.end())
    {
      permutation.insert(std::pair<boost::xint::integer,boost::xint::integer>(num, *newNumIter));
      newNumIter++;
    }
  }
  
  for (size_t i = 0; i < testsNumber; i++)
  {
    boost::xint::integer number = generatedDiagrams[i];//permutation.at(generatedDiagrams[i]);
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
  delete[] leftBounds;

  cout << "Finished random walking.\n";
  return buckets;
}
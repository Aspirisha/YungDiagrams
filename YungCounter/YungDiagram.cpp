#include <fstream>
#include <iostream>
#include "YungDiagram.h"

using namespace std;
boost::xint::integer *YungDiagram::partitionsAmount = 0;
long long *YungDiagram::offsets = 0;
double *YungDiagram::probabilities = 0;
size_t *YungDiagram::numbers = 0;
size_t YungDiagram::levelSize = 0;

YungDiagram::YungDiagram() : m_cellsNumber(0), m_colsNumber(0), m_cols(0), m_numberIsCounted(false), 
  m_number(0), m_probability(0), m_ancestorsNumber(0), m_ancestors(0), m_ancestorsAreCounted(false)
{

}

YungDiagram::YungDiagram(const char *fileName) : m_cellsNumber(0), m_colsNumber(0), m_cols(0), m_numberIsCounted(false), 
  m_number(0), m_probability(0), m_ancestorsNumber(0), m_ancestors(0)
{
  std::ifstream in(fileName);
  in >> m_cellsNumber;
  in >> m_colsNumber;

  m_cols = new size_t[m_colsNumber];
  in >> m_cols[0];
  for (size_t i = 1; i < m_colsNumber; ++i)
  {
    in >> m_cols[i];
    if (m_cols[i] > m_cols[i - 1])
    {
      std::cout << "Input file contains incorrect diagram.\n";
    }
  }
}

YungDiagram::YungDiagram(const boost::xint::integer &number) :  m_cellsNumber(0), m_colsNumber(0), m_cols(0), m_numberIsCounted(true), m_number(0),
  m_probability(0), m_ancestorsNumber(0), m_ancestors(0)
{
  if (!partitionsAmount)
  {
    CountPartitionsAmount();
  }

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
  int *temp = new int[m_cellsNumber];
  int colSize = n;

  while (n)
  {
    if (colSize == 1)
      break;
    m_number += partitionsAmount[offsets[n - 1] + colSize - 2]; // with max k'th col less than current
    if (m_number > number)
    {
      m_number -= partitionsAmount[offsets[n - 1] + colSize - 2];
      colSize--;
    }
    else
    {
      temp[m_colsNumber] = colSize;
      n -= colSize;
      m_colsNumber++;
    }
  }

  m_cols = new size_t[m_colsNumber + n];
  memcpy(m_cols, temp, m_colsNumber * sizeof(int));
  for (int i = 0; i < n; i++, m_colsNumber++)
    m_cols[m_colsNumber] = 1;
  delete[] temp;
}

void YungDiagram::CountPartitionsAmount()
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
  if (!partitionsAmount)
  {
    CountPartitionsAmount();
  }

  for (size_t i = 0; i < m_cellsNumber - 1; ++i)
    m_number += partitionsAmount[offsets[i] + i];
  m_number += 1;

  int n = m_cellsNumber;
  for (size_t k = 0; k < m_colsNumber; k++)
  {
    int colSize = m_cols[k];
    if (colSize == 1)
      break;
    m_number += partitionsAmount[offsets[n - 1] + colSize - 2]; // with max k'th col less than current
    n -= colSize;
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
  out << m_colsNumber << std::endl;
  for (size_t i = 0; i < m_colsNumber; i++)
    out << m_cols[i] << " ";
}

void YungDiagram::countAncestors()
{
  m_ancestorsNumber = 2; // left column always can be increased
  for (size_t i = 1; i < m_colsNumber; i++) 
  {
    if (m_cols[i] < m_cols[i - 1])
      m_ancestorsNumber++;
  }

  m_ancestors = new size_t[m_ancestorsNumber];
  m_ancestorsColDifferent = new size_t[m_ancestorsNumber];
  m_cols[0]++;
  m_ancestors[0] = YungDiagram::GetSmallDiagramNumber(m_cellsNumber + 1, m_colsNumber, m_cols) - 1;
  m_ancestorsColDifferent[0] = 0;
  m_cols[0]--;

  int index = 1;
  for (size_t i = 1; i < m_colsNumber; i++) 
  {
    if (m_cols[i] < m_cols[i - 1])
    {
      m_cols[i]++;
      m_ancestors[index] = YungDiagram::GetSmallDiagramNumber(m_cellsNumber + 1, m_colsNumber, m_cols) - 1;
      m_ancestorsColDifferent[index++] = i;
      m_cols[i]--;
    }
  }

  m_ancestors[m_ancestorsNumber - 1] = YungDiagram::GetSmallDiagramNumber(m_cellsNumber + 1, m_colsNumber + 1, m_cols, true) - 1;
  m_ancestorsColDifferent[m_ancestorsNumber - 1] = m_colsNumber;

  m_ancestorsAreCounted = true;
}

size_t YungDiagram::GetSmallDiagramNumber(size_t cellsNumber, size_t colsNumber, const size_t *cols, bool withExtraCell)
{
  size_t number = 0;
  if (!partitionsAmount)
  {
    CountPartitionsAmount();
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

boost::xint::integer YungDiagram::GetMaxNumberWithNCells(size_t n)
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
    delete[] m_ancestorsColDifferent;
    m_probability = 0;
    m_ancestorsAreCounted = false;
  }
  delete m_cols;
  m_cols = new size_t[colsNumber];
  ZeroMemory(m_cols, sizeof(size_t) * colsNumber);
  m_colsNumber = colsNumber;
}

void YungDiagram::setCellsInCol(size_t colIndex, size_t cellsNumber)
{
  if (m_ancestorsAreCounted)
  {
    delete[] m_ancestors;
    delete[] m_ancestorsColDifferent;
    m_probability = 0;
    m_ancestorsAreCounted = false;
  }

  if (colIndex >= m_colsNumber)
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

  if (colIndex < m_colsNumber - 1)
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

void YungDiagram::CountProbabilities(size_t cellsNumber)
{
  size_t maxDiagramNumber = YungDiagram::GetMaxNumberWithNCells(cellsNumber)._get_digit(0);
  levelSize = maxDiagramNumber - YungDiagram::GetMaxNumberWithNCells(cellsNumber - 1)._get_digit(0);

  probabilities = new double[levelSize];
  numbers = new size_t[levelSize];

  YungDiagram *diagrams = new YungDiagram[maxDiagramNumber + 1];
  diagrams[0].SetMyProbabilityRichardson(1);
  diagrams[0].setColsNumber(1);
  diagrams[0].setCellsInCol(0, 1);

  size_t index = 0;
  for (size_t i = 0; i < maxDiagramNumber; ++i) 
  {
    if (diagrams[i].m_cellsNumber == cellsNumber)
    {
      probabilities[index] = diagrams[i].m_probability;
      numbers[index] = i;
      index++;
    }

    diagrams[i].countAncestors();
    size_t ancestorsNumber = diagrams[i].getAncestorsNumber();
    size_t *ancestors = diagrams[i].getAncestors();
    long double deltaProb = diagrams[i].m_probability / (long double) ancestorsNumber;

    for (size_t j = 0; j < ancestorsNumber; j++) 
    {
      size_t idx = ancestors[j];
      if (idx > maxDiagramNumber)
        continue;
      if (diagrams[idx].m_cellsNumber == 0)
      {
        diagrams[idx].m_cellsNumber = diagrams[i].m_cellsNumber + 1;
        diagrams[idx].m_colsNumber = diagrams[i].m_colsNumber;
        if (j == ancestorsNumber - 1)
        {
          diagrams[idx].m_colsNumber++;
          diagrams[idx].m_cols = new size_t[diagrams[idx].m_colsNumber];
          diagrams[idx].m_cols[diagrams[idx].m_colsNumber - 1] = 1;
        }
        else 
          diagrams[idx].m_cols = new size_t[diagrams[idx].m_colsNumber];

        for (size_t k = 0; k < diagrams[i].m_colsNumber; k++)
        {
          diagrams[idx].m_cols[k] = diagrams[i].m_cols[k] + (diagrams[i].m_ancestorsColDifferent[j] == k);
        }


      }
      diagrams[idx].m_probability += deltaProb;
    }

    delete[] ancestors;
  }
}

void YungDiagram::SortByProbability()
{
  sort(0, levelSize - 1);
}

void YungDiagram::sort(size_t l, size_t r)
{
  if (r <= l)
    return;
  double mp = probabilities[(r + l) / 2];
  int mi = (r + l) / 2;
  int i = l;
  int j = r;
  while (i < j)
  {
    while (probabilities[i] <= mp && i < mi)
      i++;
    while (probabilities[j] > mp && j > mi)
      j--;
    if (i == j)
      break;

    double tempProb = probabilities[i];
    probabilities[i] = probabilities[j];
    probabilities[j] = tempProb;

    size_t tempInd = numbers[i];
    numbers[i] = numbers[j];
    numbers[j] = tempInd;
    
    if (i == mi)
    {
      mi = j;
      i++;
    }
    else if (j == mi)
    {
      mi = i;
      j--;
    }
    else 
    {
      i++;
      j--;
    }

  }
  sort(l, mi - 1);
  sort(mi + 1, r);
}
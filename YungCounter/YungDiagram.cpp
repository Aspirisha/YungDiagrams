#include <fstream>
#include <iostream>
#include "YungDiagram.h"

boost::xint::integer *YungDiagram::partitionsAmount = 0;
long long *YungDiagram::offsets = 0;

YungDiagram::YungDiagram() : m_cellsNumber(0), m_colsNumber(0), m_cols(0), m_numberIsCounted(false), m_number(0)
{

}

YungDiagram::YungDiagram(const char *fileName) : m_cellsNumber(0), m_colsNumber(0), m_cols(0), m_numberIsCounted(false), m_number(0)
{
  std::ifstream in(fileName);
  in >> m_cellsNumber;
  in >> m_colsNumber;

  m_cols = new int[m_colsNumber];
  in >> m_cols[0];
  for (int i = 1; i < m_colsNumber; ++i)
  {
    in >> m_cols[i];
    if (m_cols[i] > m_cols[i - 1])
    {
      std::cout << "Input file contains incorrect diagram.\n";
    }
  }
}

YungDiagram::YungDiagram(const boost::xint::integer &number) :  m_cellsNumber(0), m_colsNumber(0), m_cols(0), m_numberIsCounted(true), m_number(0)
{
  if (!partitionsAmount)
  {
    CountPartitionsAmount();
  }

  for (int i = 0; i < maxCellsNumber; ++i)
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

  m_cols = new int[m_colsNumber + n];
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

  for (int i = 0; i < m_cellsNumber - 1; ++i)
    m_number += partitionsAmount[offsets[i] + i];
  m_number += 1;

  int n = m_cellsNumber;
  for (int k = 0; k < m_colsNumber; k++)
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

void YungDiagram::SaveToFile(const char *fileName) const
{
  std::ofstream out(fileName);

  out << m_cellsNumber << std::endl;
  out << m_colsNumber << std::endl;
  for (int i = 0; i < m_colsNumber; i++)
    out << m_cols[i] << " ";
}
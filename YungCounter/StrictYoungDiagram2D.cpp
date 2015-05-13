#include "StrictYoungDiagram2D.h"
#include <iostream>

using namespace std;

std::map<std::pair<size_t, size_t>, boost::xint::integer> StrictYoungDiagram2D::partitionsAmount; 
bool StrictYoungDiagram2D::partitionsCounted = false;

StrictYoungDiagram2D::StrictYoungDiagram2D() : m_cellsNumber(1), m_numberIsCounted(true), m_number(1)
{
  m_cols.push_back(1);
}

StrictYoungDiagram2D::~StrictYoungDiagram2D()
{

}

StrictYoungDiagram2D::StrictYoungDiagram2D(const string &fileName) : m_cellsNumber(0), m_numberIsCounted(false), m_number(0)
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

StrictYoungDiagram2D::StrictYoungDiagram2D(const boost::xint::integer &number) :  m_cellsNumber(0), m_numberIsCounted(true), m_number(0)
{
  if (!partitionsCounted)
    CountPartitionsAmount();

  for (size_t i = 0; i < maxCellsForPartitionsNumber; ++i)
  {
    m_number += partitionsAmount[make_pair(i, i)];
    if (m_number >= number)
    {
      m_number -= partitionsAmount[make_pair(i, i)];
      if (m_number == number)
        m_cellsNumber = i + 1;
      else
        m_cellsNumber = i;
      break;
    }
  }
  m_number += 1;

  int n = m_cellsNumber;
  int colSize = 1;

  while (n)
  {
    boost::xint::integer tmp = m_number + partitionsAmount[make_pair(n - colSize, min(colSize - 1, n - colSize))]; // with max k'th col less than current
    if (tmp >= number)
    {
      m_number += partitionsAmount[make_pair(n - colSize - 1, min(colSize - 2, n - colSize - 1))]; // O(n)
      n -= colSize;
      m_cols.push_back(colSize);
      colSize = 1;
    }
    else
    {
      colSize++;
    }
  }

  for (int i = 0; i < n; i++)
    m_cols.push_back(1);
}

void StrictYoungDiagram2D::CountPartitionsAmount()
{
  if (partitionsCounted)
    return;

  for (int i = 1; i < maxCellsForPartitionsNumber; i++)
    partitionsAmount[make_pair(i, 0)] = 0;
  for (int i = 0; i < maxCellsForPartitionsNumber; i++)
    partitionsAmount[make_pair(0, i)] = 1;

  for (int n = 1; n < maxCellsForPartitionsNumber; n++)
  {
    for (int k = 1; k <= n; k++)
    {
      partitionsAmount[make_pair(n, k)] = partitionsAmount[make_pair(n, k - 1)] + partitionsAmount[make_pair(n - k, min(k - 1, n - k))];
      cout << partitionsAmount[make_pair(n, k)] << " ";
    }
    cout << endl;
  }


  partitionsCounted = true;
}
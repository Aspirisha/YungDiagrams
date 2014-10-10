#pragma once
#include <boost\xint\xint.hpp>

class YungDiagram
{
public:
  YungDiagram();
  YungDiagram(const char *fileName);
  YungDiagram(const boost::xint::integer &number);
  void SaveToFile(const char *fileName) const;
  boost::xint::integer GetDiagramNumber() const;
private:
  static const int maxCellsNumber = 6000;
  static boost::xint::integer *partitionsAmount; 
  static long long *offsets;
  static void CountPartitionsAmount();

  int m_cellsNumber;
  int m_colsNumber;
  int *m_cols;
  mutable boost::xint::integer m_number;
  mutable bool m_numberIsCounted;
};
#pragma once
#include <boost\xint\xint.hpp>

class YungDiagram
{
public:
  static size_t GetSmallDiagramNumber(size_t cellsNumber, size_t colsNumber, const size_t *cols, bool withExtraCell = false);
  static boost::xint::integer GetMaxNumberWithNCells(size_t n);
  static void CountProbabilities(size_t cellsNumber);
  static void SortByProbability();
  static double *probabilities;
  static size_t *numbers;
  static size_t levelSize;

  YungDiagram();
  YungDiagram(const char *fileName);
  YungDiagram(const boost::xint::integer &number);
  void SaveToFile(const char *fileName) const;
  boost::xint::integer GetDiagramNumber() const;
  double GetMyProbabilityRichardson();
  void SetMyProbabilityRichardson(long double probability) { m_probability = probability; }
  void countAncestors();
  size_t getAncestorsNumber() { return m_ancestorsNumber; };
  size_t *getAncestors() {return m_ancestors; }
  void setColsNumber(size_t colsNumber);
  void setCellsInCol(size_t colIndex, size_t cellsNumber);
  size_t getColsNumber() { return m_colsNumber; }
  size_t getCellsNumber() { return m_cellsNumber; }

  long double m_probability; // to boost up it's public
private:
  static void sort(size_t l, size_t r); // returns position of pivot
  static const int maxCellsNumber = 200;
  static boost::xint::integer *partitionsAmount; 
  static long long *offsets;
  
  static void CountPartitionsAmount();
  
  size_t *m_ancestors; // numbers of ancestors diagrams
  size_t *m_ancestorsColDifferent; // number of column that differs in ancestor(i) and this diagram 
  size_t m_ancestorsNumber; 
  bool m_ancestorsAreCounted;

  size_t m_cellsNumber;
  size_t m_colsNumber;
  size_t *m_cols;
  mutable boost::xint::integer m_number;
  mutable bool m_numberIsCounted;
};
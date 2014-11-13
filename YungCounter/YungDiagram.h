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
  double GetMyProbabilityRichardson();
  void SetMyProbability(long double probability) { m_probability = probability; }
  void countAncestors();

  size_t getAncestorsNumber() { return m_ancestorsNumber; };
  size_t *getAncestors() {return m_ancestors; }
  void setColsNumber(size_t colsNumber);
  void setCellsInCol(size_t colIndex, size_t cellsNumber);
  size_t getColsNumber() { return m_colsNumber; }
  size_t getCellsNumber() { return m_cellsNumber; }
  void resetProbability() {m_probability = 0;}

  long double m_probability; // to boost up it's public
private:  
  friend class YungDiagramHandler;

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

class YungDiagramHandler
{
public:
  static boost::xint::integer GetMaxNumberWithNCells(size_t n);
  static void CountRichardsonProbabilities(size_t cellsNumber);
  static void CountAlphaProbabilities(size_t cellsNumber, double alpha);
  static void SortByProbability();
  static size_t GetSmallDiagramNumber(size_t cellsNumber, size_t colsNumber, const size_t *cols, bool withExtraCell = false);
  static bool isPartitionsAmountCounted() { return partitionsAmount != 0;}
  static void countPartitionsAmount();

  static const  boost::xint::integer *getPartitionsAmount() { return partitionsAmount; }
  static size_t getLevelSize() { return levelSize; }
  static const long long *getOffsets() { return offsets; }
  static void saveProbabilities(const char *fileName);
  
  static const int maxCellsNumber = 200;
private:
  static void sort(int l, int r); // sorts array probabilities and numbers
  
  static boost::xint::integer *partitionsAmount; 
  static long long *offsets;
  
  static size_t levelSize;
  static size_t *numbers;
  static double *probabilities;
};
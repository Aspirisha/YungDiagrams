#pragma once
#include <boost\xint\xint.hpp>
#include "globaldefinitions.h"

enum ProcessType 
{
  RICHARDSON,
  ALPHA,
  BETA,
  PLANSHEREL,
  PLANSHEREL_POWERED_GAMMA
};

class YungDiagram
{
public:
  YungDiagram();
  explicit YungDiagram(const std::string &fileName);
  explicit YungDiagram(const boost::xint::integer &number);
  ~YungDiagram();
  void SaveToFile(const char *fileName) const;
  boost::xint::integer GetDiagramNumber() const;
  double GetMyProbabilityRichardson();
  void SetMyProbability(long double probability) { m_probability = probability; }
  void countAncestors();
  void resetAncestors();

  size_t getAncestorsNumber();
  size_t *getAncestors() {return m_ancestors; }
  void setColsNumber(size_t colsNumber);
  void setCellsInCol(size_t colIndex, size_t cellsNumber);
  size_t getColsNumber() { return m_cols.size(); }
  size_t getCellsNumber() { return m_cellsNumber; }
  void resetProbability() {m_probability = 0;}

  bool addCell(size_t col);
  long double m_probability; // to boost up it's public
//private:  
  friend class YungDiagramHandler;

  size_t *m_ancestors; // numbers of ancestors diagrams
  std::vector<size_t> m_ancestorsColDifferent; // number of column that differs in ancestor(i) and this diagram 
  size_t m_ancestorsNumber; 
  bool m_ancestorsAreCounted;

  size_t m_cellsNumber;
  std::vector<size_t> m_cols;
  mutable boost::xint::integer m_number;
  mutable bool m_numberIsCounted;
};

class YungDiagramHandler
{
public:
  static boost::xint::integer GetFirstNumberWithNPlusOneCells(size_t n);
  static boost::xint::integer GetFirstNumberWithNCells(size_t n);
  static boost::xint::integer GetLastNumberWithNCells(size_t n);

  static void CountProbabilities(ProcessType processType, size_t cellsNumber);
  static void PrintPartitionsAmount(const std::string &fileName);

  static void SortByProbability();
  static size_t GetSmallDiagramNumber(size_t cellsNumber, const std::vector<size_t> &cols, bool withExtraCell = false);
  static size_t GetSmallDiagramNumber(size_t colsNumber, ...);

  static bool isPartitionsAmountCounted() { return partitionsAmount != 0;}
  static void countPartitionsAmount();

  static const  boost::xint::integer *getPartitionsAmount() { return partitionsAmount; }
  static size_t getLevelSize() { return levelSize; }
  static const long long *getOffsets() { return offsets; }
  static void saveProbabilities(const char *fileName);

  static YungDiagram *getRandomDiagram(ProcessType procType, size_t n);

  static std::vector<size_t> getRandomWalkFrequencies(ProcessType processType, size_t cellsNumber, size_t bucketsNumber, size_t testsNumber, bool needPermutation = false);
  static void setAlpha(double alpha) { s_alpha = alpha; }
  static size_t getMaxCellsNumber() { return maxCellsNumber; }

  static void saveColumnsToFile(const char* fileName, YungDiagram *d);

  static double countKantorovichDistance(YungDiagram &d1, YungDiagram &d2);
  static double countKantorovichDistance(size_t n1, size_t n2);

  static void printSmallDiagram(size_t n);
  static void printSmallDiagramsPair(size_t n1, size_t n2);

  static void getBall(size_t num, double r, std::vector<size_t> &ballDiagrams, std::vector<double> &ballDistances);
  static void getLinearCoefficientsEstimationKantorovich(size_t diagramNum, dVector &c, dMatrix &deltas, 
    size_t checkedNeighbours, std::vector<size_t> &diagrams, std::vector<double> &dists);
  static void getLinearCoefficientsEstimationKantorovich(size_t diagramNum, dVector &c, dMatrix &deltas, 
    size_t checkedNeighbours);

  static void getDistancesToSomeDiagrams(size_t num, size_t diagramsNumber, std::vector<size_t> &diagrams, std::vector<double> &distances);

private:
  static void sort(int l, int r); // sorts array probabilities and numbers
  static void countCoefficientsForKantorovichMetric(size_t diagramNumber, dVector &coefficients, std::vector<size_t> &predecessorsNums);
  
  static const int maxCellsNumber;
  static boost::xint::integer *partitionsAmount; 
  static long long *offsets;
  
  static size_t levelSize;
  static size_t *numbers;
  static double *probabilities;
  static double s_alpha;
  static double s_gamma;
};
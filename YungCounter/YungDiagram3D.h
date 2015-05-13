#pragma once

#include <boost\xint\xint.hpp>
#include <map>
#include <vector>
#include "globaldefinitions.h"

typedef std::pair<size_t, size_t> ind_pair;
typedef std::pair<ind_pair, size_t> col_element;

class YungDiagram3D
{
public:
  enum ProcessType3D
  {
    HOOKS = 0,
    HOOKS_POWERED,
    RICHARDSON
  };

  YungDiagram3D();
  explicit YungDiagram3D(const std::string &fileName);
  explicit YungDiagram3D(const boost::xint::integer &number);
  ~YungDiagram3D();

  inline bool isCornerCell(size_t x, size_t y) const;
  inline bool canRemoveCell(size_t x, size_t y) const;
  void getPredeccessors(std::map<boost::xint::integer, std::vector<boost::xint::integer> > &predeccessors,
    std::map<boost::xint::integer, std::vector<ind_pair> > &dif_cells);

  inline double countTransitiveProb(size_t x, size_t y) const;
  void saveToFile(const char *fileName, bool forMatlab = false) const;
  void printToConsole() const;
  void readFromFile(const char *fileName);
  boost::xint::integer GetDiagramNumber(bool forceRecount = false);

  // own data
  size_t m_cellsNumber;
  std::vector<size_t> m_rowsY; // size of it is number of rows in x
  mutable std::map<ind_pair, size_t> m_cols;

  mutable bool m_numberIsCounted;
  mutable boost::xint::integer m_number;
};

class YungDiagram3DHandler
{
public:
  static void countPrimes();
  static bool arePrimesCounted() { return primesAreCounted; }

  static YungDiagram3D *getRandomWalkDiagram(YungDiagram3D::ProcessType3D type, size_t cellsNumber);
  static YungDiagram3D *getRandomWalkDiagramFast(YungDiagram3D::ProcessType3D type, size_t cellsNumber, double power = 1.0);
  static double countHook(YungDiagram3D &d, size_t x, size_t y, size_t z);
//private:
  static const int maxCellsNumber;
  static std::vector<size_t> primes; // mb not needed
  static std::map<ind_pair, size_t> primesToCells; // each cell has it's own prime number
  static bool primesAreCounted;
  static double countDistance(YungDiagram3D &d1, YungDiagram3D &d2);
  static void countCoefficientsForKantorovichMetric(boost::xint::integer diagramNum, 
    std::map<boost::xint::integer, double> &probs, dVector &c, 
    std::map<boost::xint::integer, std::vector<boost::xint::integer> > &preds, std::map<boost::xint::integer, std::vector<ind_pair> > &dif_cells,
    std::map<boost::xint::integer, double> &newProbs);
};
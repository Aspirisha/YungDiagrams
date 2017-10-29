#pragma once

#include <map>
#include <unordered_map>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include "globaldefinitions.h"

typedef std::pair<size_t, size_t> ind_pair;
typedef std::pair<ind_pair, size_t> col_element;
typedef boost::multiprecision::cpp_int mpz_int;
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
  explicit YungDiagram3D(const mpz_int &number);
  ~YungDiagram3D();

  inline bool isCornerCell(size_t x, size_t y) const;
  inline bool canRemoveCell(size_t x, size_t y) const;
  void getPredeccessors(std::map<mpz_int, std::vector<mpz_int> > &predeccessors,
    std::map<mpz_int, std::vector<ind_pair> > &dif_cells);

  inline double countTransitiveProb(size_t x, size_t y) const;
  void saveToFile(const char *fileName, bool forMatlab = false) const;
  void printToConsole() const;
  void readFromFile(const char *fileName);
  mpz_int GetDiagramNumber(bool forceRecount = false);

  // own data
  size_t m_cellsNumber;
  std::vector<size_t> m_rowsY; // size of it is number of rows in x
  mutable std::unordered_map<ind_pair, size_t, std::function<size_t(ind_pair)>> m_cols;

  mutable bool m_numberIsCounted;
  mutable mpz_int m_number;
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
  static void countCoefficientsForKantorovichMetric(mpz_int diagramNum, 
    std::map<mpz_int, double> &probs, dVector &c, 
    std::map<mpz_int, std::vector<mpz_int> > &preds, std::map<mpz_int, std::vector<ind_pair> > &dif_cells,
    std::map<mpz_int, double> &newProbs);
};
#pragma once

#include <boost\xint\xint.hpp>
#include <map>

typedef std::pair<size_t, size_t> ind_pair;
typedef std::pair<ind_pair, size_t> col_element;

enum ProcessType3D
{
  HOOKS
};

class YungDiagram3D
{
public:
  YungDiagram3D();
  explicit YungDiagram3D(const std::string &fileName);
  explicit YungDiagram3D(const boost::xint::integer &number);
  ~YungDiagram3D();


  void saveToFile(const char *fileName) const;
  void printToConsole() const;
  void readFromFile(const char *fileName);
  boost::xint::integer GetDiagramNumber();
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

  static YungDiagram3D *getRandomWalkDiagram(ProcessType3D type, size_t cellsNumber);
  static YungDiagram3D *getRandomWalkDiagramFast(ProcessType3D type, size_t cellsNumber);
  static double countHook(YungDiagram3D &d, size_t x, size_t y, size_t z);
//private:
  static const int maxCellsNumber;
  static std::vector<size_t> primes; // mb not needed
  static std::map<ind_pair, size_t> primesToCells; // each cell has it's own prime number
  static bool primesAreCounted;
};
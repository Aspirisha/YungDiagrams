#pragma once

#include <boost\xint\xint.hpp>
#include <map>

typedef std::pair<size_t, size_t> ind_pair;
typedef std::pair<ind_pair, size_t> col_element;

class YungDiagram3D
{
public:
  YungDiagram3D();
  explicit YungDiagram3D(const std::string &fileName);
  explicit YungDiagram3D(const boost::xint::integer &number);
  ~YungDiagram3D();


  void SaveToFile(const char *fileName) const;
  boost::xint::integer GetDiagramNumber() const;
  // own data
  size_t m_cellsNumber;
  size_t m_rowsX;
  size_t m_rowsY;
  std::map<ind_pair, size_t> m_cols;

  bool m_numberIsCounted;
  boost::xint::integer m_number;
};

class YungDiagram3DHandler
{
public:
  static void countPrimes();
  static bool arePrimesCounted() { return primesAreCounted; }

//private:
  static const int maxCellsNumber;
  static std::vector<size_t> primes; // mb not needed
  static std::map<ind_pair, size_t> primesToCells; // each cell has it's own prime number
  static bool primesAreCounted;
};
#pragma once

#include <boost/xint/integer.hpp>
#include <string>
#include <map>

struct StrictYoungDiagram2D
{
public:
  StrictYoungDiagram2D();
  explicit StrictYoungDiagram2D(const std::string &fileName);
  explicit StrictYoungDiagram2D(const boost::xint::integer &number);
  ~StrictYoungDiagram2D();
  void SaveToFile(const char *fileName) const;
  boost::xint::integer GetDiagramNumber() const;

  //virtual std::vector<boost::xint::integer> getAncestors() = 0;

  size_t m_cellsNumber;
  std::vector<size_t> m_cols;
  mutable boost::xint::integer m_number;
  mutable bool m_numberIsCounted;
private:
  static std::map<std::pair<size_t, size_t>, boost::xint::integer> partitionsAmount; 
  static void CountPartitionsAmount();
  static bool partitionsCounted;
  static const size_t maxCellsForPartitionsNumber = 10;
};
#pragma once

#include <string>
#include <map>
#include <vector>

struct StrictYoungDiagram2D
{
public:
  StrictYoungDiagram2D();
  explicit StrictYoungDiagram2D(const std::string &fileName);
  explicit StrictYoungDiagram2D(const uint64_t &number);
  ~StrictYoungDiagram2D();
  void SaveToFile(const char *fileName) const;
  uint64_t GetDiagramNumber() const;

  //virtual std::vector<boost::xint::integer> getAncestors() = 0;

  size_t m_cellsNumber;
  std::vector<size_t> m_cols;
  mutable uint64_t m_number;
  mutable bool m_numberIsCounted;
private:
  static std::map<std::pair<size_t, size_t>, uint64_t> partitionsAmount;
  static void CountPartitionsAmount();
  static bool partitionsCounted;
  static const size_t maxCellsForPartitionsNumber = 10;
};
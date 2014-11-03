#include <iostream>
#include <boost\xint\xint.hpp>
#include "YungDiagram.h"

using namespace std;
#define _CRT_SECURE_NO_WARNINGS
int main(void)
{
  YungDiagram diagram("max_diagram_to_count_prob.txt");

  cout << diagram.GetDiagramNumber() << endl;

  boost::xint::integer number("28");
  YungDiagram d2(number);
  d2.SaveToFile("diagram5.txt");

  diagram.GetMyProbabilityRichardson();

  YungDiagram::CountProbabilities(50);
  YungDiagram::SortByProbability();

  FILE *f;
  fopen_s(&f, "probs.txt", "w");
  for (size_t i = 0; i < YungDiagram::levelSize; i++)
    fprintf(f, "%u %.16lf\n", YungDiagram::numbers[i], YungDiagram::probabilities[i]);

  fclose(f);
  return 0;
}
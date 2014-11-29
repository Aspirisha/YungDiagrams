#include <iostream>
#include <boost\xint\xint.hpp>
#include "YungDiagram.h"

using namespace std;
#define _CRT_SECURE_NO_WARNINGS

int main(void)
{
  //YungDiagram d2(0xFFFFFFFFFFFFFFFFULL);
  //d2.SaveToFile("Max_diagram_with_unsigned_number.txt");

  ProcessType procType = RICHARDSON;
  size_t cellsNumber = 5;
  double alpha = 0.16;
  const char *fileName = 0;

  /*switch (procType)
  {
  case RICHARDSON:
    YungDiagramHandler::CountRichardsonProbabilities(cellsNumber);
    fileName = "probs_richardson.txt";
    break;
  case ALPHA:
    YungDiagramHandler::CountAlphaProbabilities(cellsNumber, alpha);
    fileName = "probs_alpha.txt";
    break;
  }
  YungDiagramHandler::SortByProbability();
  YungDiagramHandler::saveProbabilities(fileName);*/


  size_t testsNumber = 1000000;
  cellsNumber = 400;
  size_t bucketsNumber = 1000;
  YungDiagramHandler::setAlpha(0.16);
  vector<size_t> buckets = YungDiagramHandler::getRandomWalkFrequencies(RICHARDSON, cellsNumber, bucketsNumber, testsNumber);
  
  ofstream out("Freqs_400_1000buckets_richardson_perm.txt");
  for (size_t i = 0; i < bucketsNumber; i++)
    out << buckets[i] << endl;
  out.close();


  int dummy = 0;
  cin >> dummy;
  return 0;
}
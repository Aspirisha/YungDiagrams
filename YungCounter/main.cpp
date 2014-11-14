#include <iostream>
#include <boost\xint\xint.hpp>
#include "YungDiagram.h"

using namespace std;
#define _CRT_SECURE_NO_WARNINGS

enum ProcessType 
{
  RICHARDSON,
  ALPHA
};

int main(void)
{
  YungDiagram d2(0xFFFFFFFFU);
  d2.SaveToFile("Max_diagram_with_unsigned_number.txt");

  ProcessType procType = RICHARDSON;
  size_t cellsNumber = 70;
  double alpha = 0.3;
  const char *fileName = 0;

  switch (procType)
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

  YungDiagramHandler::saveProbabilities(fileName);
  return 0;
}
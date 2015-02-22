#include <iostream>
#include <boost\xint\xint.hpp>
#include "YungDiagram.h"

using namespace std;
#define _CRT_SECURE_NO_WARNINGS

int main(void)
{
  //YungDiagram d2(0xFFFFFFFFFFFFFFFFULL);
  //d2.SaveToFile("Max_diagram_with_unsigned_number.txt");

  ProcessType procType = BETA;
  size_t cellsNumber = 50;
  double alpha = 0.16;
  const char *fileName = 0;
  YungDiagramHandler::setAlpha(alpha);

  /*********************precized probabilities***************************/
  /*YungDiagramHandler::CountProbabilities(procType, cellsNumber);
  switch (procType)
  {
  case RICHARDSON:
    fileName = "probs_richardson.txt";
    break;
  case ALPHA:
    fileName = "probs_alpha.txt";
    break;
  case BETA:
    fileName = "probs_beta.txt";
    break;
  }
  //YungDiagramHandler::SortByProbability();
  YungDiagramHandler::saveProbabilities(fileName);*/



  /***************** frequencies ******************/
 /* size_t testsNumber = 1000000;
  cellsNumber = 400;
  size_t bucketsNumber = 1000;
  YungDiagramHandler::setAlpha(0.16);
  vector<size_t> buckets = YungDiagramHandler::getRandomWalkFrequencies(RICHARDSON, cellsNumber, bucketsNumber, testsNumber);*/
  
 /* ofstream out("Freqs_400_1000buckets_richardson.txt");
  for (size_t i = 0; i < bucketsNumber; i++)
    out << buckets[i] << endl;
  out.close();*/

  //**********************asymptotic*******************************/
  //YungDiagramHandler::PrintPartitionsAmount("PartitionsAmount.txt");
  //YungDiagram *d = YungDiagramHandler::getRandomDiagram(PLANSHEREL_POWERED_GAMMA, 100000);
  //YungDiagramHandler::saveColumnsToFile("100000Cells_PLANSHEREL_GAMMA.txt", d);

  /*************************Kantorovich distance*************************/
  YungDiagram d1(75);
  YungDiagram d2(76);

  double dist = YungDiagramHandler::countKantorovichDistance(d1, d2);
  cout  << "dist(4, 5) = " << dist << endl;


  int dummy = 0;


  cout << "Insert anything to quit.\n";
  cin >> dummy;
  return 0;
}
#include <iostream>
#include <boost\xint\xint.hpp>
#include "YungDiagram.h"
#include "YungDiagram3D.h"
#include "BratelliDiagram.h"
#include "3DInterface.h"
#include "2DInterface.h"
#include "StrictYoungGraph.h"

using namespace std;
#define _CRT_SECURE_NO_WARNINGS


void readFromFileAndPrintDiagram3DNumber(const char *fileName)
{
  YungDiagram3D d;
  d.readFromFile(fileName);
  cout << d.GetDiagramNumber() << endl;
}

int main(void)
{
  //YungDiagram d2(0xFFFFFFFFFFFFFFFFULL);
  //d2.SaveToFile("Max_diagram_with_unsigned_number.txt");

  ProcessType procType = BETA;
  size_t cellsNumber = 50;
  double alpha = -0.16;
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
  //YungDiagram *d = YungDiagramHandler::getRandomDiagram(PLANSHEREL, 100000);
  //YungDiagramHandler::saveColumnsToFile("100000Cells_PLANSHEREL.txt", d);

  /*************************Kantorovich distance*************************/

  //countDistanceBetweenInputDiagrams();
  //countDistanceBetweenRandomDiagrams(RICHARDSON, 10);
  //countAllDistancesOnLevel(10, "distances.txt");
  //countDistanceBetweenUniformlyRandomDiagrams();
  //YungDiagram *d = YungDiagramHandler::getRandomDiagram(PLANSHEREL, 20);
  //writeKantorovichBallToFile("BallRandomPlansherel20_Radius_0_1.txt", d->GetDiagramNumber()._get_digit(0), 0.05);
  //writeKantorovichBallToFile("AllDistances20.txt", d->GetDiagramNumber()._get_digit(0), 1.1);
  //writeKantorovichBallToFile("AllDistances20FromLine.txt", YungDiagramHandler::GetFirstNumberWithNCells(20)._get_digit(0), 1.1);
  //printDiagramsInCycle();
  //writeKantorovichEstimationCoefs("EstimationCoefs.txt");
  //writeEstimationCoefficientsAndDistances("EstimationCoefsFabs_35.txt", "EstimationDistancesFabs_35.txt", 30);

  //****************3D diagrams**********************************************/
  //printDiagrams3DInCycle();
  //readFromFileAndPrintDiagram3DNumber("3D.txt");
  printRandomDiagram3DHooks();
  //printDiagramsWithDifferentPowers();
  //countTimeForRandom3D();
  
  /************************ Bratelli diagrams *****************************************/
 // PascalGraph pg;

  //int node1 = 0;
  //int node2 = 7;
  //int level = 10;
  //cout << "dist(" << node1 << ", " << node2 << ") = " << pg.countDistance(level, node1, node2);

  /************************** Strict diagram dists ***************************************/
  //countAllStrictDistancesOnLevel();

  /************************** 3D distances *************************************************/

  //size_t num1, num2;

  //cin >> num1 >> num2;
  //cout << YungDiagram3DHandler::countDistance(YungDiagram3D(num1), YungDiagram3D(num2)) << endl;
  
  /********************** for 9/pi^2 prove **********************************************/
  
  int dummy = 0;


  cout << "Insert anything to quit.\n";
  cin >> dummy;
  return 0;
}
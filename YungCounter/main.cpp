#include <iostream>
#include <boost\xint\xint.hpp>
#include <random>
#include <chrono>
#include <limits>
#include "YungDiagram.h"
#include "YungDiagram3D.h"

using namespace std;
#define _CRT_SECURE_NO_WARNINGS

#undef max

void countDistanceBetweenInputDiagrams() 
{
  size_t n1, n2;
  cout << "Insert numbers of diagrams to find distance between: \n";
  cout << "Diagram 1: ";
  cin >> n1;
  cout << "Diagram 2: ";
  cin >> n2;

  double dist = YungDiagramHandler::countKantorovichDistance(n1, n2);
  cout  << "dist(" << n1 << ", " << n2 << ") = "<< dist << endl;
}

void countAllDistancesOnLevel(size_t level, const char *fileName)
{
  size_t n1 = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level - 1)._get_digit(0) - 1;
  size_t n2 = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level)._get_digit(0) - 1;

 // n1++;
 // n2--;
  ofstream out(fileName);
  out.precision(15);
  for (int i = n1 + 1; i <= n2; i++)
  {
    for (int j = n1 + 1; j <= n2; j++)
    {
      out << fixed << YungDiagramHandler::countKantorovichDistance(i, j) << " ";
    }
    out << endl;
  }
  out.close();
}

void countDistanceBetweenRandomDiagrams(ProcessType type, size_t cellsNum)
{
  YungDiagram *d1 = YungDiagramHandler::getRandomDiagram(type, cellsNum);
  YungDiagram *d2 = YungDiagramHandler::getRandomDiagram(type, cellsNum);
  YungDiagram *d3 = YungDiagramHandler::getRandomDiagram(type, cellsNum);

  cout << "Random numbers are: " << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber()  << endl;
  //size_t n = YungDiagramHandler::GetSmallDiagramNumber(3, 5, 3, 2);
  double dist1 = YungDiagramHandler::countKantorovichDistance(*d1, *d2);
  double dist2 = YungDiagramHandler::countKantorovichDistance(*d1, *d3);
  double dist3 = YungDiagramHandler::countKantorovichDistance(*d3, *d2);
  //  double dist = YungDiagramHandler::countKantorovichDistance(2400, 2501);
  cout.precision(15);
  cout  << "dist(" << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ") = "<< fixed << dist1 << endl;
  cout  << "dist(" << d1->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = "<< fixed << dist2 << endl;
  cout  << "dist(" << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = "<< fixed << dist3 << endl;

  cout  << "dist(" << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ")  + " 
    << "dist(" << d1->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist1 + dist2 << endl;
  cout  << "dist(" << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ")  + " 
    << "dist(" << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist1 + dist3 << endl;
  cout  << "dist(" << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ")  + " 
    << "dist(" << d1->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist2 + dist3 << endl;
}

size_t getSmallUniformlyRandomDiagram(size_t cellsNum)
{
  size_t seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  size_t n1 = YungDiagramHandler::GetFirstNumberWithNCells(cellsNum)._get_digit(0);
  size_t n2 = YungDiagramHandler::GetLastNumberWithNCells(cellsNum)._get_digit(0);
  std::uniform_int_distribution<int> distribution(n1, n2);
  size_t d1 = distribution(generator);

  return d1;
}

void countDistanceBetweenUniformlyRandomDiagrams()
{
  size_t cellsNum;
  cout << "Insert number of cells in rndom diagrams to find distance between: \n";
  cin >> cellsNum;
  size_t seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);

  size_t n1 = YungDiagramHandler::GetFirstNumberWithNCells(cellsNum)._get_digit(0);
  size_t n2 = YungDiagramHandler::GetLastNumberWithNCells(cellsNum)._get_digit(0);

  std::uniform_int_distribution<int> distribution(n1, n2);
  size_t d1 = distribution(generator);
  size_t d2 = distribution(generator);

  cout << "Random numbers are: " << d1 << ", " << d2 << endl;
  cout.precision(15);

  YungDiagramHandler::printSmallDiagramsPair(d1, d2);
  double dist = YungDiagramHandler::countKantorovichDistance(d1, d2);
  cout  << "dist(" << d1 << ", " << d2 << ") = "<< fixed << dist << endl;
}

void writeKantorovichBallToFile(const char *fileName, size_t diagramNum, double r)
{
  vector<double> ballDists;
  vector<size_t> ballDiagrams;

  YungDiagramHandler::getBall(diagramNum, r, ballDiagrams, ballDists);

  ofstream out(fileName);
  for (size_t i = 0; i < ballDists.size(); i++)
  {
    out << ballDiagrams[i] << " " << ballDists[i] << endl;
  }
}

void writeKantorovichBallToFile(const char *fileName)
{
  size_t cellsNum;
  double r;

  cout << "Count Kantorovich ball.\n";
  cout << "Cells number: ";
  cin >> cellsNum;
  cout << "radius: ";
  cin >> r;

  writeKantorovichBallToFile(fileName, getSmallUniformlyRandomDiagram(cellsNum), r);
}

void printDiagramsInCycle()
{
  cout << "Insert diagrams numbers to print until you're done. To exit insert 0.\n";

  size_t num = 0;
  while (true)
  {
    cout << "Diagram number: ";
    cin >> num;

    if (num == 0)
      break;

    YungDiagramHandler::printSmallDiagram(num);
    YungDiagram d(num);
    cout << "Columns: \n";
    for (size_t i = 0; i < d.m_cols.size(); i++)
      cout << d.m_cols[i] << " ";
    cout << endl;
  }
}

void printDiagrams3DInCycle()
{
  cout << "Insert diagrams numbers to print until you're done. To exit insert 0.\n";

  boost::xint::integer num = 0;
  while (true)
  {
    cout << "Diagram number: ";

    cin >> num;
    while (cin.fail())
    {
      cin.clear();
      cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      cin >> num;
    }
    

    if (num == 0)
      break;

    YungDiagram3D d(num);
    d.printToConsole();
    cout << endl;
    cout.flush();
  }
}

void writeKantorovichEstimationCoefs(const char *fileName, size_t diagramNum)
{
  dVector c;
  dMatrix deltas;
  YungDiagramHandler::getLinearCoefficientsEstimationKantorovich(diagramNum, c, deltas);

  ofstream out(fileName);
  for (size_t i = 0; i < c.size(); i++)
    out << c[i] << " ";
  out << endl;
  out << endl;

  YungDiagram d(diagramNum);
  size_t cellsNum = d.m_cellsNumber;

  for (size_t i = 0; i < deltas.size1(); i++)
  {
    for (size_t j = 0; j < cellsNum; j++)
    {
      out << ceil(deltas(i, j)) << " ";
    }
    out << endl;
  }
}

void writeKantorovichEstimationCoefs(const char *fileName)
{
  cout << "Estimation for Kantorovich distance.\n";
  cout << "Cells number for random diagram: ";
  size_t cellsNum;
  cin >> cellsNum;

  writeKantorovichEstimationCoefs(fileName, getSmallUniformlyRandomDiagram(cellsNum));
}

void writeEstimationCoefficientsAndDistances(const char *fileNameCoefs, const char *fileNameDists, size_t diagramNumber)
{
  writeKantorovichEstimationCoefs(fileNameCoefs, diagramNumber);
  writeKantorovichBallToFile(fileNameDists, diagramNumber, 1.1);
}

void writeEstimationCoefficientsAndDistances(const char *fileNameCoefs, const char *fileNameDists)
{
  cout << "Estimation and distances for random diagrm.\n";
  cout << "Cells number for random diagram: ";
  size_t cellsNum;
  cin >> cellsNum;

  size_t diagramNumber = getSmallUniformlyRandomDiagram(cellsNum);
  cout << "Diagram number is " << diagramNumber << endl;
  writeEstimationCoefficientsAndDistances(fileNameCoefs, fileNameDists, diagramNumber);
}

void readFromFileAndPrintDiagram3DNumber(const char *fileName)
{
  YungDiagram3D d;
  d.readFromFile(fileName);
  cout << d.GetDiagramNumber() << endl;
}

void printRandomDiagram3DHooks()
{
  cout << "Random hooks process diagram generation.\n";
  cout << "Insert number of cells: \n";
  size_t n = 0;
  cin >> n;
  if (n == 0)
    return;
  YungDiagram3D *d = YungDiagram3DHandler::getRandomWalkDiagram(HOOKS, n);
  //d->printToConsole();
  d->saveToFile("3dHooksRandom.txt");
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
  //YungDiagram *d = YungDiagramHandler::getRandomDiagram(ALPHA, 1000000);
  //YungDiagramHandler::saveColumnsToFile("1000000Cells_ALPHA.txt", d);

  /*************************Kantorovich distance*************************/

  //countDistanceBetweenInputDiagrams();
  //countDistanceBetweenRandomDiagrams(RICHARDSON, 10);
  //countAllDistancesOnLevel(10, "distances.txt");
  //countDistanceBetweenUniformlyRandomDiagrams();
  //writeKantorovichBallToFile("AllDistances20.txt");
  //printDiagramsInCycle();
  //writeKantorovichEstimationCoefs("EstimationCoefs.txt");
  //writeEstimationCoefficientsAndDistances("EstimationCoefs.txt", "EstimationDistances.txt", 2223);

  //****************3D diagrams**********************************************/
  //printDiagrams3DInCycle();
  //readFromFileAndPrintDiagram3DNumber("3D.txt");
  printRandomDiagram3DHooks();
  
  int dummy = 0;


  cout << "Insert anything to quit.\n";
  cin >> dummy;
  return 0;
}
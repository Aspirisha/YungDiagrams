#include <iostream>
#include <boost\xint\xint.hpp>
#include "YungDiagram.h"

using namespace std;

int main(void)
{
  YungDiagram diagram("diagram3.txt");

  cout << diagram.GetDiagramNumber() << endl;

  boost::xint::integer number("28");
  YungDiagram d2(number);
  d2.SaveToFile("diagram5.txt");
  return 0;
}
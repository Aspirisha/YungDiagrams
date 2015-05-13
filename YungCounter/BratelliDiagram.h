#pragma once
#include <boost\xint\xint.hpp>
#include "globaldefinitions.h"
#include <vector>

using namespace std;

struct BratelliDiagram
{
  // first line contains n1...nk - number of vertices on each level
  // next k - 1 lines contain number of predessecors Si following Si indexes o predessecors
  BratelliDiagram(const char *fileName);

  BratelliDiagram() {}

  vector<vector<int> > prevVertices;
  vector<int> verticesNumber;
  
  double countDistance(int level, int node1, int node2);

  vector<vector<unsigned long long> > pathsNumber; // number of paths to given node; we use it to count co-probabilities
};

struct PascalGraph
{
  PascalGraph() {}
  double countDistance(int level, int node1, int node2);
};

#pragma once

#include "YungDiagram3D.h"

void printRandomDiagram3DHooks();
void countTimeForRandom3D();
void printDiagrams3DInCycle();
void printDiagramsWithDifferentPowers();
// not interactive
void printRandomDiagram3DHooksOffline(int cellsNum, YungDiagram3D::ProcessType3D processType, const char *fileName, double power = 1.0);
void printMany3DDiagrams();
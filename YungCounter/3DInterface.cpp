#include <iostream>
#include <sstream>
#include <chrono>
#include <string>
#include <fstream>
#include "3DInterface.h"

using namespace std;

#undef max
// allows to generate random 3D diagram and save it to file for further matlab 
// processing; Do not try to load this diagrams with this program - the format is for matlab!
void printRandomDiagram3DHooks() {
    cout << "Random hooks process diagram generation.\n";
    cout << "Insert number of cells: ";
    size_t n = 0;
    cin >> n;
    if (n == 0)
        return;

    int procType = 0;
    cout << "Insert process type (0 - HOOKS, 1 - POWERED_HOOKS, 2 - RICHARDSON): ";
    cin >> procType;
    if (procType != 0 && procType != 1 && procType != 2)
        return;

    double power = 1;
    if (procType == 1) {
        cout << "Insert power for process POWERED_HOOKS: ";
        cin >> power;
    }

    string fileName;
    cout << "Insert file name, where the diagram will be saved: ";
    cin >> fileName;

    cout << "Starting computation...\n";

    printRandomDiagram3DHooksOffline(n, static_cast<YungDiagram3D::ProcessType3D>(procType), fileName.c_str(), power);
}

void printRandomDiagram3DHooksOffline(int cellsNum, YungDiagram3D::ProcessType3D processType, const char *fileName, double power) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    YungDiagram3D *d = YungDiagram3DHandler::getRandomWalkDiagramFast(processType, cellsNum, power);
    end = std::chrono::system_clock::now();
    size_t elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    //d->printToConsole();
    cout << "elapsed time: " << elapsed_seconds << "s\n";
    d->saveToFile(fileName, true);
    delete d;
}

void printMany3DDiagrams() {
    string fileName = "RandomHooksDiagrams\\random_hooks_";

#pragma omp parallel
    {
#pragma omp for
        for (int i = 103000; i < 400000; i += 1000) {
            YungDiagram3D *d = YungDiagram3DHandler::getRandomWalkDiagramFast(YungDiagram3D::ProcessType3D::HOOKS, i);
            stringstream ss;
            ss << fileName << i << ".txt";
            d->saveToFile(ss.str().c_str(), true);
            delete d;
        }
    }
}

void printDiagramsWithDifferentPowers() {
    double step = 0.01;
    size_t cellsNumber = 100000;
    string fileName = "powers_";

#ifdef _OPENMP
    cout << "openmp enabled\n";
#endif
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < 1000; i++) {
            YungDiagram3D *d = YungDiagram3DHandler::getRandomWalkDiagramFast(YungDiagram3D::ProcessType3D::HOOKS_POWERED, cellsNumber, step * i);
            stringstream ss;
            ss << fileName << step * i << ".txt";
            d->saveToFile(ss.str().c_str(), true);
            delete d;
        }
    }
}

void countTimeForRandom3D() {
    cout << "Time estimation for 3D generation processes.\n";
    cout << "Insert max cells number: ";
    size_t maxCells = 0;
    cin >> maxCells;
    cout << "Insert cells step (every k*|step|-cells diagram is generated, for k*|step| <= max cells): ";
    size_t step = 0;
    cin >> step;
    int procType = 0;
    cout << "Insert process type (0 - HOOKS, 1 - POWERED_HOOKS): ";
    cin >> procType;
    if (procType != 0 && procType != 1)
        return;
    double power = 1;
    if (procType == 1) {
        cout << "Insert power for process POWERED_HOOKS: ";
        cin >> power;
    }
    string fileName;
    cout << "Insert file name, where statistics will be saved: ";
    cin >> fileName;

    ofstream out(fileName);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    for (size_t j = 1; j < maxCells; j += step) {
        cout << "processing " << j << " cells...\n";
        start = std::chrono::system_clock::now();
        YungDiagram3D *d = YungDiagram3DHandler::getRandomWalkDiagramFast(static_cast<YungDiagram3D::ProcessType3D>(procType), j, power);
        end = std::chrono::system_clock::now();

        int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>
            (end - start).count();

        out << j << " " << elapsed_seconds << endl;
        delete d;
    }
}


void printDiagrams3DInCycle() {
    cout << "Insert diagrams numbers to print until you're done. To exit insert 0.\n";

    mpz_int num(0);
    while (true) {
        cout << "Diagram number: ";

        cin >> num;
        while (cin.fail()) {
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            cin >> num;
        }


        if (num.is_zero())
            break;

        YungDiagram3D d(num);
        d.printToConsole();
        cout << endl;
        cout.flush();
    }
}
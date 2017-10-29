#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <string>
#include "YungDiagram.h"
#include "StrictYoungGraph.h"

using namespace std;

static void writeKantorovichBallToFile(const char *fileName, size_t diagramNum, double r);
static void writeEstimationCoefficientsAndDistances(const char *fileNameCoefs, const char *fileNameDists, size_t diagramNumber, size_t checkedNeighbours);



void countDistanceBetweenInputDiagrams() {
    size_t n1, n2;
    cout << "Insert numbers of diagrams to find distance between: \n";
    cout << "Diagram 1: ";
    cin >> n1;
    cout << "Diagram 2: ";
    cin >> n2;

    double dist = YungDiagramHandler::countKantorovichDistance(n1, n2);
    cout << "dist(" << n1 << ", " << n2 << ") = " << dist << endl;
}

size_t getSmallUniformlyRandomDiagram(size_t cellsNum) {
    size_t seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    size_t n1 = YungDiagramHandler::GetFirstNumberWithNCells(cellsNum);
    size_t n2 = YungDiagramHandler::GetLastNumberWithNCells(cellsNum);
    std::uniform_int_distribution<int> distribution(n1, n2);
    size_t d1 = distribution(generator);

    return d1;
}


void countDistanceBetweenUniformlyRandomDiagrams() {
    size_t cellsNum;
    cout << "Insert number of cells in rndom diagrams to find distance between: \n";
    cin >> cellsNum;
    size_t seed = (size_t)std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    size_t n1 = YungDiagramHandler::GetFirstNumberWithNCells(cellsNum);
    size_t n2 = YungDiagramHandler::GetLastNumberWithNCells(cellsNum);

    std::uniform_int_distribution<int> distribution(n1, n2);
    size_t d1 = distribution(generator);
    size_t d2 = distribution(generator);

    cout << "Random numbers are: " << d1 << ", " << d2 << endl;
    cout.precision(15);

    YungDiagramHandler::printSmallDiagramsPair(d1, d2);
    double dist = YungDiagramHandler::countKantorovichDistance(d1, d2);
    cout << "dist(" << d1 << ", " << d2 << ") = " << fixed << dist << endl;
}


void writeKantorovichBallToFile(const char *fileName, size_t diagramNum, double r) {
    vector<double> ballDists;
    vector<size_t> ballDiagrams;

    YungDiagramHandler::getBall(diagramNum, r, ballDiagrams, ballDists);

    ofstream out(fileName);
    for (size_t i = 0; i < ballDists.size(); i++) {
        out << ballDiagrams[i] << " " << ballDists[i] << endl;
    }
}

void writeKantorovichBallToFile() {
    size_t cellsNum;
    double r;
    string fileName;

    cout << "Writing Kantorovich ball to file.\n";

    cout << "Cells number: ";
    cin >> cellsNum;
    cout << "radius: ";
    cin >> r;
    cout << "file name: ";
    cin >> fileName;

    writeKantorovichBallToFile(fileName.c_str(), getSmallUniformlyRandomDiagram(cellsNum), r);
}

void printDiagramsInCycle() {
    cout << "Insert diagrams numbers to print until you're done. To exit insert 0.\n";

    size_t num = 0;
    while (true) {
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


void writeEstimationCoefficientsAndDistances(const char *fileNameCoefs, const char *fileNameDists, size_t diagramNumber, size_t checkedNeighbours) {
    dVector c;
    dMatrix deltas;
    vector<size_t> diagrams;
    vector<double> dists;
    YungDiagramHandler::getLinearCoefficientsEstimationKantorovich(diagramNumber, c, deltas, checkedNeighbours, diagrams, dists);

    ofstream out(fileNameCoefs);
    cout << "c.size() = " << c.size() << endl;
    for (size_t i = 0; i < c.size(); i++)
        out << c[i] << " ";
    out << endl;
    out << endl;

    YungDiagram d(diagramNumber);
    size_t cellsNum = d.m_cellsNumber;

    for (size_t i = 0; i < deltas.size1(); i++) {
        for (size_t j = 0; j < 2 * cellsNum; j++) {
            out << ceil(deltas(i, j)) << " ";
        }
        out << endl;
    }

    out.close();
    out.open(fileNameDists);
    for (size_t i = 0; i < dists.size(); i++) {
        out << diagrams[i] << " " << dists[i] << endl;
    }

}

void writeEstimationCoefficientsAndDistances() {
    cout << "Estimation and distances for random diagrm.\n";
    cout << "Cells number for random diagram: ";
    size_t cellsNum;
    cin >> cellsNum;

    string fileNameCoefs;
    string fileNameDists;
    size_t checkedNeighbours;
    cout << "Insert number of distances, based on which estimation is counted: ";
    cin >> checkedNeighbours;
    cout << "Insert file name where to store coefficients: ";
    cin >> fileNameCoefs;
    cout << "Insert file name where to store distances: ";
    cin >> fileNameDists;


    size_t diagramNumber = getSmallUniformlyRandomDiagram(cellsNum);
    cout << "Diagram number is " << diagramNumber << endl;
    writeEstimationCoefficientsAndDistances(fileNameCoefs.c_str(), fileNameDists.c_str(), diagramNumber, checkedNeighbours);
}

void countAllDistancesOnLevel() {
    size_t level;
    string fileName;

    cout << "Counting all distances on given Yung Graph level.\n";
    cout << "Insert level number: ";
    cin >> level;
    cout << "Insert file name where to store distances: ";
    cin >> fileName;

    size_t n1 = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level - 1) - 1;
    size_t n2 = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level) - 1;

    ofstream out(fileName);
    out.precision(15);
    for (size_t i = n1 + 1; i <= n2; i++) {
        YungDiagram d1(i);
        for (size_t j = n1 + 1; j <= n2; j++) {
            YungDiagram d2(j);
            out << fixed << YungDiagramHandler::countKantorovichDistance(d1, d2) << " ";
        }
        out << endl;
    }
    out.close();
}

void countAllStrictDistancesOnLevel() {
    size_t level;
    string fileName;

    cout << "Counting all distances on given Yung Graph level.\n";
    cout << "Insert level number: ";
    cin >> level;
    cout << "Insert file name where to store distances: ";
    cin >> fileName;

    size_t n1 = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level - 1) - 1;
    size_t n2 = YungDiagramHandler::GetFirstNumberWithNPlusOneCells(level) - 1;

    ofstream out(fileName);
    out.precision(15);
    for (size_t i = n1 + 1; i <= n2; i++) {
        YungDiagram d1(i);
        if (!d1.isStrict())
            continue;
        for (size_t j = n1 + 1; j <= n2; j++) {
            YungDiagram d2(j);
            if (!d2.isStrict())
                continue;

            out << fixed << StrictYoungGraph::countDistance(d1, d2) << " ";
        }
        out << endl;
    }
    out.close();
}

void countDistanceBetweenRandomDiagrams() {
    ProcessType type;
    size_t cellsNum;

    cout << "Counting distances between three random Yung diagrams.\n";
    cout << "Insert process type number. Available Processes are\n";
    cout << "0 : Richardson process\n";
    cout << "1 : alpha process\n";
    cout << "2 : beta process\n";
    cout << "3 : plansherel process\n";
    cout << "4 : plansherel powered process.\n";
    cout << "process number: ";
    int typeId;
    cin >> typeId;
    if (typeId < 0 || typeId > 4)
        return;
    type = static_cast<ProcessType>(typeId);
    cout << "Insert cells number: ";
    cin >> cellsNum;

    YungDiagram *d1 = YungDiagramHandler::getRandomDiagram(type, cellsNum);
    YungDiagram *d2 = YungDiagramHandler::getRandomDiagram(type, cellsNum);
    YungDiagram *d3 = YungDiagramHandler::getRandomDiagram(type, cellsNum);

    cout << "Random numbers are: " << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << endl;
    //size_t n = YungDiagramHandler::GetSmallDiagramNumber(3, 5, 3, 2);
    double dist1 = YungDiagramHandler::countKantorovichDistance(*d1, *d2);
    double dist2 = YungDiagramHandler::countKantorovichDistance(*d1, *d3);
    double dist3 = YungDiagramHandler::countKantorovichDistance(*d3, *d2);
    //  double dist = YungDiagramHandler::countKantorovichDistance(2400, 2501);
    cout.precision(15);
    cout << "dist(" << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ") = " << fixed << dist1 << endl;
    cout << "dist(" << d1->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist2 << endl;
    cout << "dist(" << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist3 << endl;

    cout << "dist(" << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ")  + "
        << "dist(" << d1->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist1 + dist2 << endl;
    cout << "dist(" << d1->GetDiagramNumber() << ", " << d2->GetDiagramNumber() << ")  + "
        << "dist(" << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist1 + dist3 << endl;
    cout << "dist(" << d2->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ")  + "
        << "dist(" << d1->GetDiagramNumber() << ", " << d3->GetDiagramNumber() << ") = " << fixed << dist2 + dist3 << endl;
}

#include "DistanceOnSimplexCounter.h"
#include "../YungCounter/LinearProblemSolver.h"
#include <iostream>

using namespace std;

typedef std::pair<size_t, size_t> num_pair;

double countDistanceOnSimplex(size_t vertNum, map<num_pair, double> vertDistances,
    std::vector<double> measure1, std::vector<double> measure2) {
    if (measure1.size() != vertNum || measure2.size() != vertNum || vertDistances.size() < vertNum * (vertNum - 1) / 2) {
        cerr << "Wrong operands in countDistanceOnSimplex.\n";
        return -1;
    }

    dMatrix x(vertNum, vertNum);
    dVector a(vertNum);
    dVector b(vertNum);
    dMatrix C(vertNum, vertNum);

    for (int i = 0; i < vertNum; i++) {
        a(i) = measure1[i];
        b(i) = measure2[i];
        for (int j = i + 1; j < vertNum; j++) {
            C(j, i) = C(i, j) = vertDistances[make_pair(i + 1, j + 1)];
        }
        C(i, i) = 0;
    }

    solveTransportationPotential(C, a, b, x);

    double dist = 0;
    for (int i = 0; i < vertNum; i++) {
        for (int j = 0; j < vertNum; j++) {
            dist += C(i, j) * x(i, j);
        }
    }

    return dist;
}
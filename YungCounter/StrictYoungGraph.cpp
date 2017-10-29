#include "globaldefinitions.h"
#include "StrictYoungGraph.h"
#include "LinearProblemSolver.h"
#include <vector>
#include <map>

using namespace std;

typedef std::pair<size_t, size_t> num_pair;

static void countCoefficientsForKantorovichMetric(size_t diagramNumber, dVector &coefficients, vector<size_t> &predecessorsNums) {
    YungDiagram d(diagramNumber);

    size_t colsNumber = d.m_cols.size();

    vector<double> a; // c[i] = a[i] / sum(a)
    vector<double> weights;
    vector<size_t> removedCells;

    size_t coefsIndex = 0;
    coefficients.resize(colsNumber);

    for (int s = 0; s < colsNumber; s++) {
        bool cellCanBeRemoved = true;

        if (s != colsNumber - 1) {
            if (d.m_cols[s] == d.m_cols[s + 1] + 1) // strict!
                cellCanBeRemoved = false;
        }

        if (cellCanBeRemoved) {
            d.m_cols[s]--;
            bool needResize = false;
            double c = 1.0;
            if (d.m_cols[s] != 0) {
                predecessorsNums.push_back(YungDiagramHandler::GetSmallDiagramNumber(d.m_cellsNumber - 1, d.m_cols));

                double col_s = d.m_cols[s];

                for (int j = 0; j < s; j++) {
                    c *= (d.m_cols[j] - col_s);
                    c /= (col_s + d.m_cols[j]);
                }
                for (int j = s + 1; j < colsNumber; j++) {
                    c *= (col_s - d.m_cols[j]);
                    c /= (col_s + d.m_cols[j]);
                }
            } else {
                d.m_cols.resize(colsNumber - 1);
                predecessorsNums.push_back(YungDiagramHandler::GetSmallDiagramNumber(d.m_cellsNumber - 1, d.m_cols));
                d.m_cols.resize(colsNumber);
            }

            weights.push_back(c);
            d.m_cols[s]++;
            removedCells.push_back(s);
        }
    }

    size_t predeccessorsNumber = weights.size();
    size_t power = colsNumber * (colsNumber - 1) / 2;

    for (int i = 0; i < predeccessorsNumber; ++i) {
        double sum = 0;
        for (int j = 0; j < predeccessorsNumber; ++j) {
            sum += pow((double)(d.m_cols[removedCells[j]]) / (double)(d.m_cols[removedCells[i]]), power) * weights[j];
        }

        coefficients(i) = weights[i] / sum;
    }

    coefficients.resize(predeccessorsNumber);
}


double StrictYoungGraph::countDistance(YungDiagram &d1, YungDiagram &d2) {
    vector<size_t> neededDiagrams;
    vector<size_t> numbersOfNeededDiagramsOnLevel;
    vector<char> ancestorFlags;
    // probably we don't need to count distance between diagrams that are predecessors of the same diagram from
    // te top level. If it's predecessor of first one, it has first bit to be 1, if of the second - then second bit to be 1.
    // BUT, IF it turns out that this diagram is predecessor for both of top level diagrams,
    // we should count it's flag to be 3 - bits for both

    neededDiagrams.push_back(d1.GetDiagramNumber());
    neededDiagrams.push_back(d2.GetDiagramNumber());
    ancestorFlags.push_back(1);
    ancestorFlags.push_back(2);

    numbersOfNeededDiagramsOnLevel.push_back(2);

    int n = d1.m_cellsNumber;
    int startIndexInVector = 0; // indicate start and end of current level diagrams in neededDiagrams
    int endIndexInVector = 1;

    for (int i = n - 1; i > 2; i--) // we assume that we have metric on the second floor
    {
        size_t neededDiagramsNumberOnLevel = 0;
        for (int j = startIndexInVector; j <= endIndexInVector; j++) {
            YungDiagram d(neededDiagrams[j]);
            size_t colsNumber = d.m_cols.size();
            d.m_cellsNumber--;
            for (int k = 0; k < colsNumber; k++) {
                bool cellCanBeRemoved = true;

                if (k != colsNumber - 1) {
                    if (d.m_cols[k] == d.m_cols[k + 1] + 1)
                        cellCanBeRemoved = false;
                }

                if (cellCanBeRemoved) {
                    d.m_cols[k]--;
                    bool needResize = false;
                    if (d.m_cols[k] == 0) {
                        d.m_cols.resize(colsNumber - 1);
                        needResize = true;
                    }
                    size_t num = YungDiagramHandler::GetSmallDiagramNumber(d.m_cellsNumber, d.m_cols);
                    // check if we haven't added thid diagram as needed already
                    vector<size_t>::iterator iter;
                    if ((iter = std::find(neededDiagrams.begin(), neededDiagrams.end(), num)) == neededDiagrams.end()) {
                        if (num > 3) {
                            neededDiagrams.push_back(num);
                            neededDiagramsNumberOnLevel++;
                            ancestorFlags.push_back(ancestorFlags[j]);
                        }
                    } else // maybe modify flags
                    {
                        size_t ind = iter - neededDiagrams.begin();
                        ancestorFlags[ind] |= ancestorFlags[j];
                    }

                    if (needResize)
                        d.m_cols.push_back(0);
                    d.m_cols[k]++;
                }
            }
            d.m_cellsNumber++; // for proper deletion
        }

        numbersOfNeededDiagramsOnLevel.push_back(neededDiagramsNumberOnLevel);
        startIndexInVector = endIndexInVector + 1;
        endIndexInVector = neededDiagrams.size() - 1;
    }

    map<num_pair, double> distances; // key is pair of diagram numbers, first number is less than second.

    distances.insert(pair<num_pair, double>(pair<size_t, size_t>(5, 6), 1.0));
    endIndexInVector = neededDiagrams.size() - 1;
    size_t solvedTranspProblems = 0;
    for (int i = 3, levelIndex = numbersOfNeededDiagramsOnLevel.size() - 1; i <= n; i++, levelIndex--) {
        startIndexInVector = endIndexInVector - numbersOfNeededDiagramsOnLevel[levelIndex] + 1;

        // last is not needed for it's distances will be counted during previous steps
        for (int j = startIndexInVector; j < endIndexInVector; j++) {
            YungDiagram d1(neededDiagrams[j]);

            size_t colsNumber1 = d1.m_cols.size();
            vector<size_t> nums1; // numbers of predecessors of d1
            dVector c1;
            countCoefficientsForKantorovichMetric(neededDiagrams[j], c1, nums1);

            for (int k = j + 1; k <= endIndexInVector; k++) {
                if (ancestorFlags[k] == ancestorFlags[j] && (ancestorFlags[k] ^ 3)) // boost!? yes, about 2 times
                    continue;

                dVector c2;
                vector<size_t> nums2; // numbers of predecessors of d2
                countCoefficientsForKantorovichMetric(neededDiagrams[k], c2, nums2);

                // now we need to solve transportation problem
                size_t s1 = c1.size();
                size_t s2 = c2.size();
                dMatrix costs(s1, s2);
                dMatrix x(s1, s2);
                for (size_t l = 0; l < s1; l++) {
                    for (size_t p = 0; p < s2; p++) {
                        size_t d1 = nums1[l];
                        size_t d2 = nums2[p];

                        if (d1 == d2)
                            costs(l, p) = 0;
                        else if (d1 < d2)
                            costs(l, p) = distances[num_pair(d1, d2)];
                        else
                            costs(l, p) = distances[num_pair(d2, d1)];
                    }
                }

                if (!solveTransportationPotential(costs, c1, c2, x))
                    cout << "Potential method error\n";
                else {
                    solvedTranspProblems++;
                    double dist = 0;
                    for (size_t l = 0; l < s1; l++) {
                        for (size_t p = 0; p < s2; p++) {
                            dist += costs(l, p) * x(l, p);
                        }
                    }
                    if (neededDiagrams[j] < neededDiagrams[k])
                        distances.insert(pair<num_pair, double>(num_pair(neededDiagrams[j], neededDiagrams[k]), dist));
                    else
                        distances.insert(pair<num_pair, double>(num_pair(neededDiagrams[k], neededDiagrams[j]), dist));
                }
            }
        }
        endIndexInVector = startIndexInVector - 1;
    }

    cout << "During distance finding were solved " << solvedTranspProblems << " transportation problems.\n";

    if (neededDiagrams[0] < neededDiagrams[1])
        return distances[num_pair(neededDiagrams[0], neededDiagrams[1])];
    return distances[num_pair(neededDiagrams[1], neededDiagrams[0])];
}
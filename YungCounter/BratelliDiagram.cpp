#include "BratelliDiagram.h"
#include <map>
#include "LinearProblemSolver.h"

// levels from 0
double PascalGraph::countDistance(int level, int node1, int node2) {
    vector<vector<int> > neededNodes1;
    vector<vector<int> > neededNodes2;
    vector<map<pair<int, int>, double> > dists;

    neededNodes1.resize(level + 1);
    neededNodes2.resize(level + 1);
    dists.resize(level + 1);

    neededNodes1[level].push_back(node1);
    neededNodes2[level].push_back(node2);

    for (int i = level - 1; i > 1; i--) {
        for (int j = 0; j < neededNodes1[i + 1].size(); j++) {
            int curNode1 = neededNodes1[i + 1][j];
            if (curNode1 > 0)
                neededNodes1[i].push_back(curNode1 - 1);
            if (curNode1 <= i)
                neededNodes1[i].push_back(curNode1);
        }

        for (int j = 0; j < neededNodes2[i + 1].size(); j++) {
            int curNode2 = neededNodes2[i + 1][j];
            if (curNode2 > 0)
                neededNodes2[i].push_back(curNode2 - 1);
            if (curNode2 <= i)
                neededNodes2[i].push_back(curNode2);
        }
    }

    dists[1][make_pair(0, 1)] = 1; // distance on 1 floor (2 verts)

    int solvedTranspProblems = 0;

    for (int i = 2; i <= level; i++) {
        for (int j = 0; j < neededNodes1[i].size(); j++) {
            int n1 = neededNodes1[i][j];
            for (int k = 0; k < neededNodes2[i].size(); k++) {
                int n2 = neededNodes2[i][k];

                pair<int, int> ind = make_pair(min(n1, n2), max(n1, n2));
                if (dists[i].find(ind) != dists[i].end())
                    continue; // already counted
                  // now we need to solve transportation problem
                vector<int> nums1;
                int ind1 = 0;
                int ind2 = 0;
                dVector c1(2);
                dVector c2(2);
                if (n1 > 0) {
                    c1(ind1++) = double(n1) / i;
                    nums1.push_back(n1 - 1);
                }
                if (n1 < i) {
                    c1(ind1++) = double(i - n1) / i;
                    nums1.push_back(n1);
                }
                c1.resize(ind1);
                vector<int> nums2;
                if (n2 > 0) {
                    c2(ind2++) = double(n2) / i;
                    nums2.push_back(n2 - 1);
                }
                if (n2 < i) {
                    c2(ind2++) = double(i - n2) / i;
                    nums2.push_back(n2);
                }
                c2.resize(ind2);

                size_t s1 = nums1.size();
                size_t s2 = nums2.size();
                dMatrix costs(s1, s2);
                dMatrix x(s1, s2);
                for (size_t l = 0; l < s1; l++) {
                    for (size_t p = 0; p < s2; p++) {
                        size_t d1 = nums1[l];
                        size_t d2 = nums2[p];

                        if (d1 == d2)
                            costs(l, p) = 0;
                        else if (d1 < d2)
                            costs(l, p) = dists[i - 1][make_pair(d1, d2)];
                        else
                            costs(l, p) = dists[i - 1][make_pair(d2, d1)];
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
                    dists[i][ind] = dist;
                }
            }
        }

    }
    cout << "During distance finding were solved " << solvedTranspProblems << " transportation problems.\n";
    return dists[level][make_pair(min(node1, node2), max(node1, node2))];
}
#include "DistanceOnSimplexCounter.h"
#include <fstream>

using namespace std;

typedef std::pair<size_t, size_t> num_pair;

int main() {
    map<num_pair, double> vertDists;
    const double d[] = { 1.0, 2.0, 2.0 };

    if (d[0] + d[1] < d[2] || d[0] + d[2] < d[1] || d[1] + d[2] < d[0]) {
        cerr << "Metric doesn't satisfy triangle inequivalence.\n";
        return 1;
    }

    vertDists.insert(pair<num_pair, double>(num_pair(1, 2), d[0]));
    vertDists.insert(pair<num_pair, double>(num_pair(1, 3), d[1]));
    vertDists.insert(pair<num_pair, double>(num_pair(2, 3), d[2]));

    static const double m1[] = { 0.1, 0.2, 0.7 };
    static const double m2[] = { 0.3, 0.1, 0.6 };
    vector<double> measure1(m1, m1 + sizeof(m1) / sizeof(double));
    vector<double> measure2(m2, m2 + sizeof(m2) / sizeof(double));

    double dist = countDistanceOnSimplex(3, vertDists, measure1, measure2);
    cout << dist << endl;

    double step = 0.01;
    double u = 0, v = 0;
    static const double m3[] = { 1.0 / 3, 1.0 / 3, 1.0 / 3 };
    measure1 = vector<double>(m3, m3 + sizeof(m3) / sizeof(double));
    double fact = sqrt(3.0) / 2.0;
    ofstream out("SimplexDistances.txt");
    out.precision(15);

    double x3 = d[0];
    double x2 = (d[0] * d[0] + d[1] * d[1] - d[2] * d[2]) / (2 * d[0]);
    double y2 = sqrt(d[1] * d[1] - x2 * x2);

    for (double u = 0; u < 1.0; u += step) {
        for (double v = 0; v < 1 - u; v += step) {
            double w = 1.0 - u - v;
            const double m3[] = { u, v, w };
            vector<double> measure2(m3, m3 + sizeof(m3) / sizeof(double));

            double x = v * x2 + w * x3;
            double y = v * y2;
            double dist = countDistanceOnSimplex(3, vertDists, measure1, measure2);

            out << x << " " << y << " " << fixed << dist << endl;
        }
    }


    int dummy;
    cin >> dummy;
    return 0;
}
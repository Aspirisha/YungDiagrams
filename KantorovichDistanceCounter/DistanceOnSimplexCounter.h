#pragma once

#include "../YungCounter/globaldefinitions.h"
#include <vector>
#include <map>

typedef std::pair<size_t, size_t> num_pair;

double countDistanceOnSimplex(size_t vertNum, std::map<num_pair, double> vertDistances, std::vector<double> measure1, std::vector<double> measure2);
#pragma once

#include "../YungCounter/globaldefinitions.h"
#include <vector>
#include <map>

double countDistanceOnSimplex(size_t vertNum, std::map<num_pair, double> vertDistances, std::vector<double> measure1, std::vector<double> measure2);
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include "svputils.hpp"
#include "vecops.hpp"

std::vector<std::vector<double>> lll_reduce(std::vector<std::vector<double>> basis);

std::vector<std::vector<double>> combogen(int r, int n);

std::vector<double> sumvec(std::vector<std::vector<double>> vecs, size_t ll, size_t ul);

double brute_force(std::vector<std::vector<double>> basis);

std::vector<double> se_enum(std::vector<std::vector<double>> basis, double bound);

double solve_svp(std::vector<std::vector<double>> basis);


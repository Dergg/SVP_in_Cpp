#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> gs(std::vector<std::vector<double>> basis);

bool check_lin_ind(std::vector<std::vector<double>> basis);

bool checkValid(std::vector<std::vector<double>> basis);

std::vector<std::vector<double>> build_basis(int argc, char* argv[]);
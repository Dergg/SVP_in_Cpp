#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>

double dotproduct(const std::vector<double>& v1, const std::vector<double>& v2);

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);

double norm(const std::vector<double>& v1);

std::vector<double> operator*(const std::vector<double>& v1, float scalar);
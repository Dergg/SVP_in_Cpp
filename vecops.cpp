#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include "vecops.hpp"

double dotproduct(const std::vector<double>& v1,
const std::vector<double>& v2) {
    double result = 0;
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

std::vector<double>
operator-(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result;
    for (int i = 0; i < v1.size(); i++) {
        result.push_back(v1[i] - v2[i]);
    }
    return result;
}

double norm(const std::vector<double>& v1) {
    double result = 0;
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v1[i];
    }
    return sqrt(result);
}

std::vector<double>
operator*(const std::vector<double>& v1, const float scalar) {
    std::vector<double> result;
    for (int i = 0; i < v1.size(); i++) {
        result.push_back(v1[i] * scalar);
    }
    return result;
}

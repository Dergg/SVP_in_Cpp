#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include "svputils.hpp"
#include "svp.hpp"

void runTests(std::vector<std::vector<double>> basis, double expected) {
    std::cout << std::fixed << std::setprecision(11);
    double actual = solve_svp(basis);
    std::cout << "For basis: " << std::endl;
    for (auto& vec : basis) {
        for (auto& val : vec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Expected: " << expected << std::endl;
    std::cout << "Actual: " << actual << std::endl;
    if (actual != expected) {
        std::cout << "We didn't get the expected result." << std::endl;
    }
    assert(actual == expected);
}

int main() {
    runTests({{1, 0}, {0, 1}}, 1);
    runTests({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1);
}


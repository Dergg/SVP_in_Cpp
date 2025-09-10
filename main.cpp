#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include "svputils.hpp"
#include "vecops.hpp"
#include "svp.hpp"

int main(int argc, char** argv) {
    std::cout << std::fixed << std::setprecision(11);
    std::vector<std::vector<double>> basis = build_basis(argc, argv);
    double shortest = solve_svp(basis);
    std::cout << "Shortest norm: " << shortest << std::endl;
    std::ofstream myfile;
    myfile.open("result.txt");
    myfile << std::fixed << std::setprecision(11) << shortest << std::endl;
    myfile.close();
}

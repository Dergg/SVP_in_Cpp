#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "vecops.hpp"
#include "svputils.hpp"

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
gs(std::vector<std::vector<double>> basis) {
    // Gram-Schmidt orthogonalisation
    std::vector<std::vector<double>> gsbasis = basis;
    int n = basis.size();
    std::vector<std::vector<double>> mu(n, std::vector<double>(n, 0));
    gsbasis[0] = basis[0];
    for (int i = 1; i < n; i++) {
        std::vector<double> v = basis[i];
        for (int j = i - 1; j > -1 ; j--) {
            mu[i][j] = dotproduct(basis[i], gsbasis[j])
            / dotproduct(gsbasis[j], gsbasis[j]);
            v = (v - (gsbasis[j] * mu[i][j]));
        }
        gsbasis[i] = v;
    }
    return std::make_pair(gsbasis, mu);
}

bool check_lin_ind(std::vector<std::vector<double>> basis) {
    // Checks if basis is linearly independent
    // Returns true if linearly independent, false otherwise
    std::vector<std::vector<double>> matrix = basis;
    int n = matrix.size();
    bool swapped;

    for (int i = 0; i < n; i++) {
        if (matrix[i][i] == 0) {
            for (int j = i + 1; j < n; j++) {
                if (matrix[j][i] != 0) {
                    std::swap(matrix[i], matrix[j]);
                    swapped = true;
                    break;
                }
            }
        }
        for (int j = i+1; j < n; j++) {
            double ratio = static_cast<double>(matrix[j][i]) / matrix[i][i];
            for (int k = 0; k < n; k++) {
                matrix[j][k] -= ratio * matrix[i][k];
            }
        }
    }
    return true;
}

bool checkValid(std::vector<std::vector<double>> basis) {
    if (basis.size() == 0) {
        std::cout << "Basis is empty." << std::endl;
        return false;  // Can't have a 0 length basis
    }
    for (int i = 0; i < basis.size(); i++) {
        if (basis[i].size() != basis.size()) {
            std::cout <<
            "Length of each vector not equal to number of vectors."
            << std::endl;
            return false;  // See above
        }
    }
    if (check_lin_ind(basis) == false) {
        std::cout << "Basis is not linearly independent." << std::endl;
        return false;  // Basis has to be linearly independent
    }

    return true;  // If we pass all those checks, we're all gravy, baby.
}

std::vector<std::vector<double>> build_basis(int argc, char* argv[]) {
    std::vector<std::vector<double>> basis;
    std::vector<double> sublist;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        arg.erase(std::remove(arg.begin(), arg.end(), ','), arg.end());
        if (arg.find('[') != std::string::npos) {
            arg.erase(std::remove(arg.begin(), arg.end(), '['), arg.end());
            sublist.push_back(std::stod(arg));
        } else if (arg.find(']') != std::string::npos) {
            arg.erase(std::remove(arg.begin(), arg.end(), ']'), arg.end());
            sublist.push_back(std::stod(arg));
            basis.push_back(sublist);
            sublist.clear();
        } else {
            sublist.push_back(std::stod(arg));
        }
    }

    return basis;
}

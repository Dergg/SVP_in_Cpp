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
#include "svp.hpp"

std::vector<std::vector<double>>
lll_reduce(std::vector<std::vector<double>> basis) {
    // LLL reduction, returning an LLL reduced basis
    // Based on pseudocode from https://www.math.auckland.ac.nz/~sgal018/crypto-book/main.pdf
    // Code was originally written in Python, then translated to C++
    int n = basis.size();
    auto gsify = gs(basis);
    std::vector<std::vector<double>> gsbasis = gsify.first;
    std::vector<std::vector<double>> mu = gsify.second;
    std::vector<double> bi;
    for (int i = 0; i < n; i++) {
        bi.push_back(dotproduct(gsbasis[i], gsbasis[i]));
    }
    int k = 1;
    while (k < n) {
        for (int j = k - 1; j > -1; j--) {
            double qj = abs(mu[k][j]);
            if (qj > 0.5) {
                basis[k] = (basis[k] - (basis[j] * qj));
                for (int j = 0; j < k; j++) {
                    mu[k][j] = dotproduct(basis[k], gsbasis[j])
                    / norm(gsbasis[j]);
                }
            }
        }
        if (bi[k] >= ((0.75 - pow(mu[k][k-1], 2)) * bi[k-1])) {
            k++;
        } else {
            std::swap(basis[k], basis[k-1]);
            std::vector<double> v1 = basis[k];
            std::vector<double> v2 = basis[k-1];
            for (int j = 0; j < k; j++) {
                mu[k-1][j] = dotproduct(basis[k-1], gsbasis[j])
                / dotproduct(gsbasis[j], gsbasis[j]);
                mu[k][j] = dotproduct(basis[k], gsbasis[j])
                / dotproduct(gsbasis[j], gsbasis[j]);
                v1 = v1 - (basis[j] * mu[k][j]);
                v2 = v2 - (basis[j] * mu[k-1][j]);
            }
            gsbasis[k] = v1;
            gsbasis[k-1] = v2;
            bi[k] = dotproduct(gsbasis[k], gsbasis[k]);
            bi[k-1] = dotproduct(gsbasis[k-1], gsbasis[k-1]);
            for (int i = k + 1; i < n; i++) {
                mu[i][k] = dotproduct(basis[i], gsbasis[k])
                / dotproduct(gsbasis[k], gsbasis[k]);
                mu[i][k-1] = dotproduct(basis[i], gsbasis[k-1])
                / dotproduct(gsbasis[k-1], gsbasis[k-1]);
            }
            k = std::max(k-1, 2);
        }
     }
    return basis;
}

std::vector<std::vector<double>> combogen(int r, int n) {
    if (r == 0) {
        return {{}};
    }

    std::vector<std::vector<double>> result;
    for (int i = -n; i <= n; i++) {
        auto subcomb = combogen(r - 1, n);
        for (auto& subcombos : subcomb) {
            subcombos.insert(subcombos.begin(), i);
            result.push_back(subcombos);
        }
    }
    return result;
}

std::vector<double>
sumvec(std::vector<std::vector<double>> vecs, size_t ll, size_t ul) {
    std::vector<double> result;
    for (size_t i = ll; i < ul; i++) {
        double sum = 0;
        for (size_t j = 0; j < vecs.size(); j++) {
            sum += vecs[j][i];
        }
        result.push_back(sum);
    }
    return result;
}

double brute_force(std::vector<std::vector<double>> basis) {
    int n = basis.size();
    float cutoff = 1.5;
    std::vector<double> shortest(n, 0);
    double snorm = INFINITY;
    double prevlen = 0;
    for (int r = 1; r <= n; r++) {
        for (const auto& combo : combogen(r, n)) {
            std::vector<double> v(n, 0);
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < n; j++) {
                    v[j] += basis[i][j] * combo[i];
                }
            }
            double len = norm(v);
            if (len > snorm * cutoff || len == prevlen) {
                continue;
            }
            if (len < snorm && len != 0) {
                shortest = v;
                snorm = len;
                std::cout << "Shortest norm found so far: "
                << snorm << std::endl;
            }
            prevlen = len;
        }
    }
    return snorm;
}

std::vector<double>
se_enum(std::vector<std::vector<double>> basis, double bound) {
    auto gsify = gs(basis);
    std::vector<std::vector<double>> gsbasis = gsify.first;
    std::vector<std::vector<double>> mu = gsify.second;
    double R2 = std::pow(bound, 2);
    int n = basis.size();
    std::vector<double> bstarnrm;
    for (int i = 0; i < n; i++) {
        bstarnrm.push_back(std::pow(norm(gsbasis[i]), 2));
    }
    std::vector<double> p(n+1, 0);
    std::vector<double> v(n, 0);
    v[0] = 1;
    std::vector<double> c(n, 0);
    std::vector<double> w(n, 0);
    std::vector<double> s(n, 0);
    int k = 0;
    int last_nz = 0;
    std::vector<std::vector<double>> temp;
    while (true) {
        p[k] = (p[k+1] + ((std::pow((v[k] - c[k]), 2)) * bstarnrm[k]));
        if (p[k] < R2) {
            if (k == 0) {
                R2 = p[k];
                for (int i = 0; i < n; i++) {
                    std::vector<double> prod = basis[i] * v[i];
                    temp.push_back(prod);
                    prod.clear();
                }
                s = sumvec(temp, 0, n);
                temp.clear();
            } else {
                k--;
                double sum = 0;
                for (int i = k + 1; i < n; i++) {
                    sum += mu[i][k] * v[i];
                }
                c[k] = -(std::abs(sum));
                v[k] = std::round(c[k]);
                w[k] = 1;
            }
        } else {
            k += 1;
            if (k == n) {
                return s;
            }
            if (k >= last_nz) {
                last_nz = k;
                v[k] += 1;
            } else {
                if (v[k] > c[k]) {
                    v[k] -= w[k];
                } else {
                    v[k] += w[k];
                }
                w[k] += 1;
            }
        }
        std::cout << "Shortest norm found so far: " << norm(s) << std::endl;
    }
}

double solve_svp(std::vector<std::vector<double>> basis) {
    if (checkValid(basis) == false) {
        return -1;
    } else {
        std::cout << "LLL reducing the basis." << std::endl;
        basis = lll_reduce(basis);
        std::cout << "LLL reduction done. Enumerating now." << std::endl;
        // double shortest = brute_force(basis);
        double shortest = norm(se_enum(basis, INFINITY));
        return shortest;
    }
}

// Poisson.hpp
#ifndef POISSON_HPP
#define POISSON_HPP

#include "header.hpp"

#include <vector>
#include "iostream"

class Poisson {
public:
    std::vector<double> x_left_idx;
    std::vector<double> x_right_idx;
    std::vector<double> y_left_idx;
    std::vector<double> y_right_idx;

    // 생성자
    Poisson();

    void poisson_solve(std::vector<double>& div, std::vector<double>& phi);
};

#endif
// TDMA_SYSTEM_HPP
#ifndef TDMA_SYSTEM_HPP
#define TDMA_SYSTEM_HPP

#include "header.hpp"

#include <vector>
#include "iostream"

class TDMASystem {
public:
    // u-momentum
    std::vector<double> rhs_x_u;
    std::vector<double> a_x_u, b_x_u, c_x_u, d_x_u;

    std::vector<double> rhs_y_u;
    std::vector<double> a_y_u, b_y_u, c_y_u, d_y_u;

    // v-momentum
    std::vector<double> rhs_y_v;
    std::vector<double> a_y_v, b_y_v, c_y_v, d_y_v;

    std::vector<double> rhs_x_v;
    std::vector<double> a_x_v, b_x_v, c_x_v, d_x_v;

    // 생성자
    TDMASystem();

    void update_rhs(std::vector<double>& rhs_1, std::vector<double>& rhs_2);

    void TDMA(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double>& d, int n1);

    void get_du(std::vector<double>& du);
    void get_dv(std::vector<double>& dv);
};




#endif // TDMA_SYSTEM_HPP
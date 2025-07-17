#include "TDMASystem.hpp"
#include "save.hpp"
#include "iostream" 

TDMASystem::TDMASystem() {
    // u-momentum
    rhs_x_u = std::vector<double>(Nx * (Ny + 1));
    a_x_u = std::vector<double>(Nx - 2, -coef_x);
    b_x_u = std::vector<double>(Nx - 2, 1.0 + 2.0 * coef_x);
    c_x_u = std::vector<double>(Nx - 2, -coef_x);
    d_x_u = std::vector<double>(Nx - 2, 0.0);

    rhs_y_u = std::vector<double>(Nx * (Ny + 1));
    a_y_u = std::vector<double>(Ny - 1, -coef_y);
    b_y_u = std::vector<double>(Ny - 1, 1.0 + 2.0 * coef_y);
    c_y_u = std::vector<double>(Ny - 1, -coef_y);
    d_y_u = std::vector<double>(Ny - 1, 0.0);
    a_y_u[0]     = -8.0 / 3.0 * coef_y;
    a_y_u[Ny-2]  = -4.0 / 3.0 * coef_y;
    b_y_u[0]     =  1.0 + 4.0 * coef_y;
    b_y_u[Ny-2]  =  1.0 + 4.0 * coef_y;
    c_y_u[0]     = -4.0 / 3.0 * coef_y;
    c_y_u[Ny-2]  = -8.0 / 3.0 * coef_y;

    // v-momentum
    rhs_y_v = std::vector<double>((Nx + 1) * Ny);
    a_y_v = std::vector<double>(Ny - 2, -coef_y);
    b_y_v = std::vector<double>(Ny - 2, 1.0 + 2.0 * coef_y);
    c_y_v = std::vector<double>(Ny - 2, -coef_y);
    d_y_v = std::vector<double>(Ny - 2, 0.0);

    rhs_x_v = std::vector<double>((Nx + 1) * Ny);
    a_x_v = std::vector<double>(Nx - 1, -coef_x);
    b_x_v = std::vector<double>(Nx - 1, 1.0 + 2.0 * coef_x);
    c_x_v = std::vector<double>(Nx - 1, -coef_x);
    d_x_v = std::vector<double>(Nx - 1, 0.0);
    a_x_v[0]     = -8.0 / 3.0 * coef_x;
    a_x_v[Nx-2]  = -4.0 / 3.0 * coef_x;
    b_x_v[0]     =  1.0 + 4.0 * coef_x;
    b_x_v[Nx-2]  =  1.0 + 4.0 * coef_x;
    c_x_v[0]     = -4.0 / 3.0 * coef_x;
    c_x_v[Nx-2]  = -8.0 / 3.0 * coef_x;

    std::vector<double> a_tmp;
    std::vector<double> b_tmp;
    std::vector<double> c_tmp;
}

void TDMASystem::update_rhs(std::vector<double>& rhs_1, std::vector<double>& rhs_2) {

    for (int j=1; j<Ny; ++j) {
        for (int i=1; i<Nx-1; ++i) {
            rhs_x_u[IDX(i, j, Nx)] = rhs_1[IDX(i, j, Nx)];
        }
    }

    for (int j=1; j<Ny-1; ++j) {
        for (int i=1; i<Nx; ++i) {
            rhs_y_v[IDX(i, j, Nx+1)] = rhs_2[IDX(i, j, Nx+1)];
        }
    }
}

void TDMASystem::TDMA(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double>& d, int n1) {
    
    int i;
    double r;

    d[0] = d[0]/b[0];
    c[0] = c[0]/b[0];

    for (i=1; i<n1; ++i) {
        r = 1.0/(b[i]-a[i]*c[i-1]);
        d[i] = r*(d[i]-a[i]*d[i-1]);
        c[i] = r*c[i];
    }

    for (i=n1-2; i>=0; --i) {
        d[i] = d[i]-c[i]*d[i+1];
    }
}

void TDMASystem::get_du(std::vector<double>& du) {
    // rhs_x_u (TDMA x) => rhs_y_u
    for (int j=1; j<Ny; ++j) {
        
        for (int i=1; i<Nx-1; ++i) {
            d_x_u[i-1] = rhs_x_u[IDX(i, j, Nx)];
        }
        
        TDMA(a_x_u, b_x_u, c_x_u, d_x_u, Nx-2);
        
        for (int i=1; i<Nx-1; ++i) {
            rhs_y_u[IDX(i, j, Nx)] = d_x_u[i-1];
        }
    }

    // rhs_y_u (TDMA y) => du
    for (int i=1; i<Nx-1; ++i) {

        for (int j=1; j<Ny; ++j) {
            d_y_u[j-1] = rhs_y_u[IDX(i, j, Nx)];
        }
        
        TDMA(a_y_u, b_y_u, c_y_u, d_y_u, Ny-1);

        for (int j=1; j<Ny; ++j) {
            du[IDX(i, j, Nx)] = d_y_u[j-1];
        }
    }
}

void TDMASystem::get_dv(std::vector<double>& dv) {
    // rhs_y_v (TDMA y) => rhs_x_v
    for (int i=1; i<Nx; ++i) {

        for (int j=1; j<Ny-1; ++j) {
            d_y_v[j-1] = rhs_y_v[IDX(i, j, Nx+1)];
        }

        TDMA(a_y_v, b_y_v, c_y_v, d_y_v, Ny-2);

        for (int j=1; j<Ny-1; ++j) {
            rhs_x_v[IDX(i, j, Nx+1)] = d_y_v[j-1];
        }
    }

    // rhs_x_v (TDMA x) => dv
    for (int j=1; j<Ny-1; ++j) {

        for (int i=1; i<Nx; ++i) {
            d_x_v[i-1] = rhs_x_v[IDX(i, j, Nx+1)];
        }

        TDMA(a_x_v, b_x_v, c_x_v, d_x_v, Nx-1);

        for (int i=1; i<Nx; ++i) {
            dv[IDX(i, j, Nx+1)] = d_x_v[i-1];
        }
    }
}
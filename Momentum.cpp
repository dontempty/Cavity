#include "Momentum.hpp"
#include "cmath"

void make_Momentum(std::vector<double>& u, std::vector<double>& v, std::vector<double>& u_1, std::vector<double>& v_1,
                    std::vector<double>& Momentum_x, std::vector<double>& Momentum_y) {
    
    int i, j;

    // x momentum
    double u2_ri;
    double u2_le;
    double uv_up;
    double uv_down;
    double u2_ri_1;
    double u2_le_1;
    double uv_up_1;
    double uv_down_1;
    for (j=1; j<Ny; ++j) {
        for (i=1; i<Nx-1; ++i) {

            // -dt/2 * (3H_{n} - H_{n-1}) + 2 * (A_x + A_y) * u_{n}
            u2_ri   = ( pow(  u[IDX(i+1, j, Nx)], 2) + pow(  u[IDX(i, j, Nx)], 2) ) / 2.0;
            u2_ri_1 = ( pow(u_1[IDX(i+1, j, Nx)], 2) + pow(u_1[IDX(i, j, Nx)], 2) ) / 2.0;
            u2_le   = ( pow(  u[IDX(i-1, j, Nx)], 2) + pow(  u[IDX(i, j, Nx)], 2) ) / 2.0;
            u2_le_1 = ( pow(u_1[IDX(i-1, j, Nx)], 2) + pow(u_1[IDX(i, j, Nx)], 2) ) / 2.0;

            if (j==Ny-1) {
                uv_up   = (  u[IDX(i, j+1, Nx)]) * (  v[IDX(i+1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) / 2.0;
                uv_up_1 = (u_1[IDX(i, j+1, Nx)]) * (v_1[IDX(i+1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) / 2.0;

                uv_down   = (  u[IDX(i, j-1, Nx)] +   u[IDX(i, j, Nx)]) * (  v[IDX(i+1, j-1, Nx+1)] +   v[IDX(i, j-1, Nx+1)]) / 4.0;
                uv_down_1 = (u_1[IDX(i, j-1, Nx)] + u_1[IDX(i, j, Nx)]) * (v_1[IDX(i+1, j-1, Nx+1)] + v_1[IDX(i, j-1, Nx+1)]) / 4.0;

                Momentum_x[IDX(i, j, Nx)] = 2.0 * coef_y * (4.0/3.0) * (2.0*u[IDX(i, j+1, Nx)] - 3.0*u[IDX(i, j, Nx)] + u[IDX(i, j-1, Nx)]);
            }
            else if (j==1) {
                uv_up     = (  u[IDX(i, j+1, Nx)] +   u[IDX(i, j, Nx)]) * (  v[IDX(i+1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) / 4.0;
                uv_up_1   = (u_1[IDX(i, j+1, Nx)] + u_1[IDX(i, j, Nx)]) * (v_1[IDX(i+1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) / 4.0;

                uv_down   = (  u[IDX(i, j-1, Nx)]) * (  v[IDX(i+1, j-1, Nx+1)] +   v[IDX(i, j-1, Nx+1)]) / 2.0;
                uv_down_1 = (u_1[IDX(i, j-1, Nx)]) * (v_1[IDX(i+1, j-1, Nx+1)] + v_1[IDX(i, j-1, Nx+1)]) / 2.0;

                Momentum_x[IDX(i, j, Nx)] = 2.0 * coef_y * (4.0/3.0) * (u[IDX(i, j+1, Nx)] - 3.0*u[IDX(i, j, Nx)] + 2.0*u[IDX(i, j-1, Nx)]);
            }
            else {
                uv_up     = (  u[IDX(i, j+1, Nx)] +   u[IDX(i, j, Nx)]) * (  v[IDX(i+1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) / 4.0;
                uv_up_1   = (u_1[IDX(i, j+1, Nx)] + u_1[IDX(i, j, Nx)]) * (v_1[IDX(i+1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) / 4.0;

                uv_down   = (  u[IDX(i, j-1, Nx)] +   u[IDX(i, j, Nx)]) * (  v[IDX(i+1, j-1, Nx+1)] +   v[IDX(i, j-1, Nx+1)]) / 4.0;
                uv_down_1 = (u_1[IDX(i, j-1, Nx)] + u_1[IDX(i, j, Nx)]) * (v_1[IDX(i+1, j-1, Nx+1)] + v_1[IDX(i, j-1, Nx+1)]) / 4.0;

                Momentum_x[IDX(i, j, Nx)] = 2.0 * coef_y * (u[IDX(i, j+1, Nx)] - 2.0*u[IDX(i, j, Nx)] + u[IDX(i, j-1, Nx)]);
            }

            Momentum_x[IDX(i, j, Nx)] += -1.5 * dt * ( (u2_ri - u2_le)/dx + (uv_up - uv_down)/dy );
            Momentum_x[IDX(i, j, Nx)] +=  0.5 * dt * ( (u2_ri_1 - u2_le_1)/dx + (uv_up_1 - uv_down_1)/dy );
            Momentum_x[IDX(i, j, Nx)] +=  2.0 * coef_x * (u[IDX(i+1, j, Nx)] - 2.0*u[IDX(i, j, Nx)] + u[IDX(i-1, j, Nx)]);
        }
    }

    // y momentum
    double vu_ri = 0;
    double vu_le = 0;
    double v2_up = 0;
    double v2_down = 0;
    double vu_ri_1 = 0;
    double vu_le_1 = 0;
    double v2_up_1 = 0;
    double v2_down_1 = 0;
    
    for (j=1; j<Ny-1; ++j) {
        for (i=1; i<Nx; ++i) {

            // -dt/2 * (3H_{n} - H_{n-1}) + 2 * (A_x + A_y) * v_{n}
            v2_up     = ( pow(  v[IDX(i, j+1, Nx+1)], 2) + pow(  v[IDX(i, j, Nx+1)], 2) ) / 2.0;
            v2_up_1   = ( pow(v_1[IDX(i, j+1, Nx+1)], 2) + pow(v_1[IDX(i, j, Nx+1)], 2) ) / 2.0;
            v2_down   = ( pow(  v[IDX(i, j-1, Nx+1)], 2) + pow(  v[IDX(i, j, Nx+1)], 2) ) / 2.0;
            v2_down_1 = ( pow(v_1[IDX(i, j-1, Nx+1)], 2) + pow(v_1[IDX(i, j, Nx+1)], 2) ) / 2.0;   

            if (i==Nx-1) {
                vu_ri   = (  v[IDX(i+1, j, Nx+1)]) * (  u[IDX(i, j+1, Nx)] +   u[IDX(i, j, Nx)]) / 2.0;
                vu_ri_1 = (v_1[IDX(i+1, j, Nx+1)]) * (u_1[IDX(i, j+1, Nx)] + u_1[IDX(i, j, Nx)]) / 2.0;

                vu_le   = (  v[IDX(i-1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) * (  u[IDX(i-1, j+1, Nx)] +   u[IDX(i-1, j, Nx)]) / 4.0;
                vu_le_1 = (v_1[IDX(i-1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) * (u_1[IDX(i-1, j+1, Nx)] + u_1[IDX(i-1, j, Nx)]) / 4.0;

                Momentum_y[IDX(i, j, Nx+1)] = 2.0 * coef_x * (4.0/3.0) * (2.0*v[IDX(i+1, j, Nx+1)] - 3.0*v[IDX(i, j, Nx+1)] + v[IDX(i-1, j, Nx+1)]);
            }
            else if (i==1) {
                vu_ri   = (  v[IDX(i+1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) * (  u[IDX(i, j+1, Nx)] +   u[IDX(i, j, Nx)]) / 4.0;
                vu_ri_1 = (v_1[IDX(i+1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) * (u_1[IDX(i, j+1, Nx)] + u_1[IDX(i, j, Nx)]) / 4.0;

                vu_le   = (  v[IDX(i-1, j, Nx+1)]) * (  u[IDX(i-1, j+1, Nx)] +   u[IDX(i-1, j, Nx)]) / 2.0;
                vu_le_1 = (v_1[IDX(i-1, j, Nx+1)]) * (u_1[IDX(i-1, j+1, Nx)] + u_1[IDX(i-1, j, Nx)]) / 2.0;
                
                Momentum_y[IDX(i, j, Nx+1)] = 2.0 * coef_x * (4.0/3.0) * (v[IDX(i+1, j, Nx+1)] - 3.0*v[IDX(i, j, Nx+1)] + 2.0*v[IDX(i-1, j, Nx+1)]);
            }
            else {
                vu_ri   = (  v[IDX(i+1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) * (  u[IDX(i, j+1, Nx)] +   u[IDX(i, j, Nx)]) / 4.0;
                vu_ri_1 = (v_1[IDX(i+1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) * (u_1[IDX(i, j+1, Nx)] + u_1[IDX(i, j, Nx)]) / 4.0;

                vu_le   = (  v[IDX(i-1, j, Nx+1)] +   v[IDX(i, j, Nx+1)]) * (  u[IDX(i-1, j+1, Nx)] +   u[IDX(i-1, j, Nx)]) / 4.0;
                vu_le_1 = (v_1[IDX(i-1, j, Nx+1)] + v_1[IDX(i, j, Nx+1)]) * (u_1[IDX(i-1, j+1, Nx)] + u_1[IDX(i-1, j, Nx)]) / 4.0;

                Momentum_y[IDX(i, j, Nx+1)] = 2.0 * coef_x * (v[IDX(i+1, j, Nx+1)] - 2.0*v[IDX(i, j, Nx+1)] + v[IDX(i-1, j, Nx+1)]);
            }

            Momentum_y[IDX(i, j, Nx+1)] += -1.5 * dt * ( (v2_up - v2_down)/dy + (vu_ri - vu_le)/dx );
            Momentum_y[IDX(i, j, Nx+1)] +=  0.5 * dt * ( (v2_up_1 - v2_down_1)/dy + (vu_ri_1 - vu_le_1)/dx );
            Momentum_y[IDX(i, j, Nx+1)] +=  2.0 * coef_y * ( v[IDX(i, j+1, Nx+1)] - 2.0*v[IDX(i, j, Nx+1)] + v[IDX(i, j-1, Nx+1)] );
        }
    }
}


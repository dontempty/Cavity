#include "Poisson.hpp"
#include <cmath>

Poisson::Poisson() {
    x_left_idx = std::vector<double>(Nx+1, 0.0);
    x_right_idx = std::vector<double>(Nx+1, 0.0);
    y_left_idx = std::vector<double>(Ny+1, 0.0);
    y_right_idx = std::vector<double>(Ny+1, 0.0);

    x_left_idx[1] = 1.0;
    x_right_idx[Nx-1] = 1.0;
    y_left_idx[1] = 1.0;
    y_right_idx[Ny-1] = 1.0;
}

void Poisson::poisson_solve(std::vector<double>& div, std::vector<double>& phi) {

    int i, j;

    double pi = 3.14159265358979323846;
    double omega = 2.0/(1.0 + pi/(My));

    int maxiter = 10000;
    double rsme;
    double dx2 = dx*dx;
    double dy2 = dy*dy;

    std::vector<double> phi_old((Nx+1)*(Ny+1), 0.0);

    for (int iter = 0; iter < maxiter; ++iter) {
        rsme = 0;
        
        // dx = dy 라고 가정
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx; ++i) {
                phi_old[IDX(i, j, Nx+1)] = phi[IDX(i, j, Nx+1)];

                phi[IDX(i, j, Nx+1)] = (1 - omega) * (phi[IDX(i, j, Nx+1)])
                                     + (    omega) * (phi[IDX(i+1, j, Nx+1)] * (1.0 - x_right_idx[i])
                                                    + phi[IDX(i-1, j, Nx+1)] * (1.0 - x_left_idx[i])
                                                    + phi[IDX(i, j+1, Nx+1)] * (1.0 - y_right_idx[j])
                                                    + phi[IDX(i, j-1, Nx+1)] * (1.0 - y_left_idx[j])
                                                    - dx2 * div[IDX(i, j, Nx+1)]) / ( 4.0 - (x_right_idx[i] + x_left_idx[i] + y_right_idx[j] + y_left_idx[j]) );

                rsme = std::max(rsme, std::fabs(phi_old[IDX(i, j, Nx+1)] - phi[IDX(i, j, Nx+1)]));
            }
        }

        // // vertex
        // phi_old[IDX(1, 1, Nx+1)] = phi[IDX(1, 1, Nx+1)];
        // phi[IDX(1, 1, Nx+1)] = (1-omega) * phi[IDX(1, 1, Nx+1)] + omega * ( phi[IDX(2, 1, Nx+1)] + phi[IDX(1, 2, Nx+1)]  - dx2 * div[IDX(1, 1, Nx+1)]) / 2.0;

        // phi_old[IDX(1, Ny-1, Nx+1)] = phi[IDX(1, Ny-1, Nx+1)];
        // phi[IDX(1, Ny-1, Nx+1)] = (1-omega) * phi[IDX(1, Ny-1, Nx+1)] + omega * ( phi[IDX(2, Ny-1, Nx+1)] + phi[IDX(1, Ny-2, Nx+1)]  - dx2 * div[IDX(1, Ny-1, Nx+1)]) / 2.0;

        // phi_old[IDX(Nx-1, 1, Nx+1)] = phi[IDX(Nx-1, 1, Nx+1)];
        // phi[IDX(Nx-1, 1, Nx+1)] = (1-omega) * phi[IDX(Nx-1, 1, Nx+1)] + omega * ( phi[IDX(Nx-1, 2, Nx+1)] + phi[IDX(Nx-2, 1, Nx+1)]  - dx2 * div[IDX(Nx-1, 1, Nx+1)]) / 2.0;

        // phi_old[IDX(Nx-1, Ny-1, Nx+1)] = phi[IDX(Nx-1, Ny-1, Nx+1)];
        // phi[IDX(Nx-1, Ny-1, Nx+1)] = (1-omega) * phi[IDX(Nx-1, Ny-1, Nx+1)] + omega * ( phi[IDX(Nx-1, Ny-2, Nx+1)] + phi[IDX(Nx-2, Ny-1, Nx+1)]  - dx2 * div[IDX(Nx-1, Ny-1, Nx+1)]) / 2.0;

        // // side
        // for (j=2; j<Ny-1; ++j) {
        //     phi_old[IDX(1, j, Nx+1)] = phi[IDX(1, j, Nx+1)];
        //     phi[IDX(1, j, Nx+1)] = (1-omega) * phi[IDX(1, j, Nx+1)] + omega * ( phi[IDX(2, j, Nx+1)] + phi[IDX(1, j+1, Nx+1)] + phi[IDX(1, j-1, Nx+1)] - dx2 * div[IDX(1, j, Nx+1)] ) / 3.0;

        //     phi_old[IDX(Nx-1, j, Nx+1)] = phi[IDX(Nx-1, j, Nx+1)];
        //     phi[IDX(Nx-1, j, Nx+1)] = (1-omega) * phi[IDX(Nx-1, j, Nx+1)] + omega * ( phi[IDX(Nx-2, j, Nx+1)] + phi[IDX(Nx-1, j+1, Nx+1)] + phi[IDX(Nx-1, j-1, Nx+1)] - dx2 * div[IDX(Nx-1, j, Nx+1)] ) / 3.0;
        // }
        // for (i=2; i<Nx-1; ++i) {
        //     phi_old[IDX(i, 1, Nx+1)] = phi[IDX(i, 1, Nx+1)];
        //     phi[IDX(i, 1, Nx+1)] = (1-omega) * phi[IDX(i, 1, Nx+1)] + omega * ( phi[IDX(i, 2, Nx+1)] + phi[IDX(i+1, 1, Nx+1)] + phi[IDX(i-1, 1, Nx+1)] - dx2 * div[IDX(i, 1, Nx+1)] ) / 3.0;
            
        //     phi_old[IDX(i, Ny-1, Nx+1)] = phi[IDX(i, Ny-1, Nx+1)];
        //     phi[IDX(i, Ny-1, Nx+1)] = (1-omega) * phi[IDX(i, Ny-1, Nx+1)] + omega * ( phi[IDX(i, Ny-2, Nx+1)] + phi[IDX(i+1, Ny-1, Nx+1)] + phi[IDX(i-1, Ny-1, Nx+1)] - dx2 * div[IDX(i, Ny-1, Nx+1)] ) / 3.0;
        // }

        // // inside
        // for (j=2; j<Ny-1; ++j) {
        //     for (i=2; i<Nx-1; ++i) {
        //         phi_old[IDX(i, j, Nx+1)] = phi[IDX(i, j, Nx+1)];
        //         phi[IDX(i, j, Nx+1)] = (1-omega) * phi[IDX(i, j, Nx+1)] + omega * (
        //                                         phi[IDX(i+1, j, Nx+1)] +
        //                                         phi[IDX(i-1, j, Nx+1)] +
        //                                         phi[IDX(i, j+1, Nx+1)] +
        //                                         phi[IDX(i, j-1, Nx+1)]
        //                                         - dx2 * div[IDX(i, j, Nx+1)]
        //                                     ) / 4.0;    
        //     }
        // }

        // for (j=1; j<Ny; ++j) {
        //     for (i=1; i<Nx; ++i) {
        //         rsme = std::max(rsme, std::fabs(phi_old[IDX(i, j, Nx+1)] - phi[IDX(i, j, Nx+1)]));
        //     }
        // }

        if (rsme < tol) {
            break;
        }
    }
}
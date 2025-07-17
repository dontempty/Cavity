#include "save.hpp"

#include <vector>
#include "iostream"  
#include "cmath"

void poisson_solve(std::vector<double>& div, std::vector<double>& phi,
                std::vector<double>& x_left_idx,
                std::vector<double>& x_right_idx,
                std::vector<double>& y_left_idx,
                std::vector<double>& y_right_idx,
                int Nx, int Ny, double dx, double dy) {

    double pi = 3.14159265358979323846;
    double omega = 2.0/(1.0 + 2.0*pi/(Ny + Nx));

    int maxiter = 10000;
    double error;
    double rsme;
    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double tol = 1e-12;
    int idx, idx_im, idx_ip, idx_jm, idx_jp;

    std::vector<double> phi_old((Nx+1)*(Ny+1), 0.0);

    for (int iter = 0; iter < maxiter; ++iter) {
        rsme = 0;
        
        // dx = dy 라고 가정
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx; ++i) {
                
                idx = j * (Nx+1) + i;
                idx_im = (j) * (Nx+1) + (i-1);
                idx_ip = (j) * (Nx+1) + (i+1);
                idx_jm = (j-1) * (Nx+1) + (i);
                idx_jp = (j+1) * (Nx+1) + (i);

                phi_old[idx] = phi[idx];

                // phi[idx] = (1 - omega) * (phi[idx])
                //          + (    omega) * (phi[idx_ip] * (1.0 - x_right_idx[i])
                //                         + phi[idx_im] * (1.0 - x_left_idx[i])
                //                         + phi[idx_jp] * (1.0 - y_right_idx[j])
                //                         + phi[idx_jm] * (1.0 - y_left_idx[j])
                //                         - dx2 * div[idx])
                // / (4.0 - x_right_idx[i] + x_left_idx[i] + y_right_idx[j] + y_left_idx[j]);

                if (i==1 && j==1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_ip] + phi[idx_jp] - dx2*div[idx]) / 2.0;
                }
                else if (i==1 && j==Ny-1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_ip] + phi[idx_jm] - dx2*div[idx]) / 2.0;
                }
                else if (i==Nx-1 && j==1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_im] + phi[idx_jp] - dx2*div[idx]) / 2.0;
                }
                else if (i==Nx-1 && j==Ny-1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_im] + phi[idx_jm] - dx2*div[idx]) / 2.0;
                }
                else if (i==1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_ip] + phi[idx_jp] + phi[idx_jm] - dx2*div[idx]) / 3.0;
                }
                else if (i==Nx-1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_im] + phi[idx_jp] + phi[idx_jm] - dx2*div[idx]) / 3.0;
                }
                else if (j==1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_ip] + phi[idx_im] + phi[idx_jp] - dx2*div[idx]) / 3.0;
                }
                else if (j==Ny-1) {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_ip] + phi[idx_im] + phi[idx_jm] - dx2*div[idx]) / 3.0;
                }
                else {
                    phi[idx] = (1 - omega) *phi[idx] + (omega) * (phi[idx_ip] + phi[idx_im] + phi[idx_jp] + phi[idx_jm] - dx2*div[idx]) / 4.0;
                }
            }
        }
        
        // cal error
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx; ++i) {
                idx = j * (Nx+1) + i;
                error = phi_old[idx] - phi[idx];
                rsme += error*error;
            }
        }
        rsme  = sqrt(rsme / (Nx+1) / (Ny+1));

        if (rsme < tol) {
            std::cout << "poisson iteration:" << iter << "| Error = " << rsme << std::endl;
            break;
        }
    }
}

int main() {
    double Pi = 3.14159265358979323846;

    int Mx = 32;
    int My = 32;
    int Nx = Mx+1;
    int Ny = My+1;
    double x0 = -0.5;
    double xN =  0.5;
    double y0 = -0.5;
    double yN =  0.5;
    double lx = (xN-x0);
    double ly = (yN-y0);
    double dx = lx/Mx;
    double dy = ly/My;

    int idx;


    std::vector<double> X(Nx+1);
    for (int i=0; i<Nx+1; ++i) {
        if (i==0) {
            X[i] = x0;
        }
        else if (i==Nx) {
            X[i] = xN;
        }
        else {
            X[i] = x0 + dx/2 + dx*(i-1);
        }
    }
    std::vector<double> Y(Ny+1);
    for (int i=0; i<Ny+1; ++i) {
        if (i==0) {
            Y[i] = y0;
        }
        else if (i==Ny) {
            Y[i] = yN;
        }
        else {
            Y[i] = y0 + dy/2 + dy*(i-1);
        }
    } 

    std::vector<double> exact_sol((Nx+1)*(Ny+1), 0.0);
    for (int j=0; j<Ny+1; ++j) {
        for (int i=0; i<Nx+1; ++i) {
            idx = j * (Nx+1) + i;

            exact_sol[idx] = sin(Pi * X[i]) * sin(Pi * Y[j]);

        }
    }


    std::vector<double> sol((Nx+1)*(Ny+1), 0.0);
    std::vector<double> f((Nx+1)*(Ny+1), 0.0);
    for (int j=0; j<Ny+1; ++j) {
        for (int i=0; i<Nx+1; ++i) {
            idx = (j) * (Nx+1) + (i);
            f[idx] = -2 * Pi*Pi * sin(Pi*X[i]) * sin(Pi*Y[j]);
        }
    }

    std::vector<double> x_left_idx(Nx+1, 0.0);
    std::vector<double> x_right_idx(Nx+1, 0.0);
    std::vector<double> y_left_idx(Ny+1, 0.0);
    std::vector<double> y_right_idx(Ny+1, 0.0);
    x_left_idx[1] = 1.0;
    x_right_idx[Nx-1] = 1.0;
    y_left_idx[1] = 1.0;
    y_right_idx[Ny-1] = 1.0;

    poisson_solve(f, sol, x_left_idx, x_right_idx, y_left_idx, y_right_idx, Nx, Ny, dx, dy);

    save_rhs_to_csv(f, Nx+1, Ny+1, "haha", "f.csv");
    save_rhs_to_csv(sol, Nx+1, Ny+1, "haha", "sol.csv");
    save_rhs_to_csv(exact_sol, Nx+1, Ny+1, "haha", "exact_sol.csv");

    return 0;
}
#include "header.hpp"
#include "Momentum.hpp"
#include "TDMASystem.hpp"
#include "Poisson.hpp"
#include "save.hpp"
#include "load.hpp"

#include <vector>
#include "iostream"   
#include <algorithm>  // std::max
#include <cmath> // std::fabs
#include <chrono>

int main() {

    // 저장공간 init
    
    // u, v, u_1, v_1 (staggered grid 사용)
    std::vector<double> u((Nx)*(Ny+1));
    std::vector<double> v((Nx+1)*(Ny));
    std::vector<double> u_1((Nx)*(Ny+1));
    std::vector<double> v_1((Nx+1)*(Ny));
    for (int i=0; i<Nx; ++i) {
        u[IDX(i, Ny, Nx)] = U;
        u[IDX(i, 0, Nx)] = 0.0;
        u_1[IDX(i, Ny, Nx)] = U;
        u_1[IDX(i, 0, Nx)] = 0.0;
    }
    for (int j=0; j<Ny; ++j) {
        v[IDX(Nx, j, Nx+1)] = 0.0;
        v[IDX(0, j, Nx+1)] = 0.0;
        v_1[IDX(Nx, j, Nx+1)] = 0.0;
        v_1[IDX(0, j, Nx+1)] = 0.0;
    }

    // TDMA
    std::vector<double> rhs_x_u((Nx)*(Ny+1));
    std::vector<double> rhs_y_v((Nx+1)*(Ny));

    TDMASystem TDMASystem;
    std::vector<double> du((Nx)*(Ny+1));
    std::vector<double> dv((Nx+1)*(Ny));
    

    // div, phi(pseudo pressure)
    Poisson Poisson_solver;
    std::vector<double> div((Nx+1)*(Ny+1));
    std::vector<double> phi((Nx+1)*(Ny+1));

    // time 시작
    int t;
    int vervose = 0;
    int check_num = 100;
    std::chrono::system_clock::time_point start_step1, start_step2; // 시작 시간
    std::chrono::duration<double> end_sec1, end_sec2; // 걸린 시간간

    start_step1 = std::chrono::system_clock::now();
    for (t=1; t<Nt+1; ++t) {

        if (t%check_num==0) {
            vervose = 1;
            std::cout << "Re=" << Re <<", iter:"<< t << "------------------------"<< std::endl; 
            start_step2 = std::chrono::system_clock::now();
        }
        else {
            vervose = 0;
        }

        // get momentum 
        make_Momentum(u, v, u_1, v_1, rhs_x_u, rhs_y_v);

        TDMASystem.update_rhs(rhs_x_u, rhs_y_v);

        TDMASystem.get_du(du);
        TDMASystem.get_dv(dv);
        
        // velocity update -----------------------------------------------

        // u_1 = u, u += du
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx-1; ++i) {
                u_1[IDX(i, j, Nx)]  =  u[IDX(i, j, Nx)];
                  u[IDX(i, j, Nx)] += du[IDX(i, j, Nx)];
            }
        }

        // v_1 = v, v += dv
        for (int j=1; j<Ny-1; ++j) {
            for (int i=1; i<Nx; ++i) {
                v_1[IDX(i, j, (Nx+1))]  =  v[IDX(i, j, (Nx+1))];
                  v[IDX(i, j, (Nx+1))] += dv[IDX(i, j, (Nx+1))];
            }
        }

        // Solve pseudo-pressure equation ------------------------------------------------------

        // cal divergence
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx; ++i) {
                div[IDX(i, j, Nx+1)]  = ( u[IDX(i, j, Nx)]   - u[IDX(i-1, j, Nx)]   ) / (dx * dt);
                div[IDX(i, j, Nx+1)] += ( v[IDX(i, j, Nx+1)] - v[IDX(i, j-1, Nx+1)] ) / (dy * dt);
            }
        }

        Poisson_solver.poisson_solve(div, phi);

        // velocity correction
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx-1; ++i) {
                u[IDX(i, j, Nx)] -= dt * ( phi[IDX(i+1, j, Nx+1)] - phi[IDX(i, j, Nx+1)] ) / dx;
            }
        }

        for (int j=1; j<Ny-1; ++j) {
            for (int i=1; i<Nx; ++i) {
                v[IDX(i, j, Nx+1)] -= dt * ( phi[IDX(i, j+1, Nx+1)] - phi[IDX(i, j, Nx+1)] ) / dy;
            }
        }

        // cal error
        double u_error = 0;
        double v_error = 0;
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx-1; ++i) {
                u_error = std::max(u_error, std::fabs(u[IDX(i, j, Nx)] - u_1[IDX(i, j, Nx)]));
            }
        }
        for (int j=1; j<Ny; ++j) {
            for (int i=1; i<Nx-1; ++i) {
                v_error = std::max(v_error, std::fabs(v[IDX(i, j, Nx+1)] - v_1[IDX(i, j, Nx+1)]));
            }
        }
        
        // cal divergence
        // std::vector<double> divergence((Nx+1)*(Ny+1));
        // for (int j=1; j<Ny; ++j) {
        //     for (int i=1; i<Nx; ++i) {
        //         divergence[IDX(i, j, Nx+1)]  = ( u[IDX(i, j, Nx)] -   u[IDX(i-1, j, Nx)] ) / (dx);
        //         divergence[IDX(i, j, Nx+1)] += ( v[IDX(i, j, Nx+1)] - v[IDX(i, j-1, Nx+1)] ) / (dy);
        //     }
        // }
        // save_rhs_to_csv(divergence, Nx+1, Ny+1, "RE_" + std::to_string(static_cast<int>(Re)) + "_" + std::to_string(Mx), "div_" + std::to_string(t) + ".csv", 15);
        
        
        // break
        if (u_error+v_error < tol) {
            std::cout << "Converge at " << t << std::endl;
            save_rhs_to_csv(u, Nx, Ny+1, "RE_" + std::to_string(static_cast<int>(Re)) + "_" + std::to_string(Mx), "u_" + std::to_string(t) + ".csv", 15);
            save_rhs_to_csv(v, Nx+1, Ny, "RE_" + std::to_string(static_cast<int>(Re)) + "_" + std::to_string(Mx), "v_" + std::to_string(t) + ".csv", 15);
            break;
        }

        // vervose
        if (vervose) {
            std::cout << "t = " << t << " | Error = " << u_error+v_error << std::endl;
            end_sec2 = std::chrono::system_clock::now() - start_step2;
            std::cout << "1 iter time: " << end_sec2.count() <<" seconds"<< std::endl;
            save_rhs_to_csv(u, Nx, Ny+1, "RE_" + std::to_string(static_cast<int>(Re)) + "_" + std::to_string(Mx), "u_" + std::to_string(t) + ".csv", 15);
            save_rhs_to_csv(v, Nx+1, Ny, "RE_" + std::to_string(static_cast<int>(Re)) + "_" + std::to_string(Mx), "v_" + std::to_string(t) + ".csv", 15);
        }
    }

    end_sec1 = std::chrono::system_clock::now() - start_step1;
    std::cout << "1 iter time: " << end_sec1.count() <<" seconds"<< std::endl;
}
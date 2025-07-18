#ifndef HEADER_H
#define HEADER_H

#define IDX(i, j, nx) ((j) * (nx) + (i))

// Mesh
#define Mx 128
#define My 128
#define Nx (Mx+1)
#define Ny (My+1)
#define x0 0.0
#define xN 1.0
#define y0 0.0
#define yN 1.0
#define lx (xN - x0)
#define ly (yN - y0)
#define dx ((lx) / (Mx))
#define dy ((ly) / (My))

// Parameter
#define tol 1e-12
#define rho 1.0
#define mu 1e-3
#define U 1.0
#define Re ( rho * U * lx / mu )
#define cfl 0.7

// Time
#define dt ( cfl * dx / U )
#define t0 0.0
#define Nt 30000
#define tN (dt * Nt) 

#define coef_x ( dt / (2.0 * Re * dx*dx) )
#define coef_y ( dt / (2.0 * Re * dy*dy) )

#endif // HEADER_H
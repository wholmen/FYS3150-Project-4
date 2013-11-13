#include <iostream>
#include <armadillo>
#include <fstream>
#include "functions.h"

using namespace arma;
using namespace std;

int main()
{
    // Declaring variables and conditions on diffusal equation
    double L = 1, T = 1;
    // Setting integration points and length of dimension
    int Nx = 10, Nt = 200;
    double dx = L / Nx, dt = T / Nt;

    // Declaring vector v. Given as v = u + us, where u is true solution.
    vec v = zeros<vec>(Nx); // zeros call set boundary conditions.

    // Setting the initial conditions.
    v(0) = 0; v(Nx-1) = 0;
    for (int i=1;i<Nx-1;i++){
        v(i) = i*dx - 1;
    }

    ForwardEuler (&v,dt,dx,Nx,Nt); // Changing v into spacial solution for t=T

    v(0) = 0; v(Nx-1) = 0;
    for (int i=1;i<Nx-1;i++){
        v(i) = i*dx - 1;
    }

    BackwardEuler(&v,dt,dx,Nx,Nt);

    v(0) = 0; v(Nx-1) = 0;
    for (int i=1;i<Nx-1;i++){
        v(i) = i*dx - 1;
    }

    CrankNicolson(&v,dt,dx,Nx,Nt);
    v = Analytical(dt,dx,Nx,Nt);

    return 0;
}



#include <iostream>
#include <armadillo>
#include <fstream>

using namespace arma;
using namespace std;

void ForwardEuler(vec *v, double dt, double dx, int Nx, int Nt);
void BackwardEuler(vec *v, double dt, double dx, int Nx, int Nt);
void CrankNicolson();
void Analytical();
vec tridiagonal(vec a1, vec a2, vec a3, vec b,int N);

int main()
{
    // Declaring variables and conditions on diffusal equation
    double v0 = 0, vL = 0, L = 1, T = 10;
    // Setting integration points and length of dimension
    int Nx = 1e2, Nt = 1e2;
    double dx = L / Nx, dt = T / Nt;

    // Declaring vector v. Given as v = u + us, where u is true solution.
    vec v = zeros<vec>(Nx); // zeros call set boundary conditions.

    // Setting initial conditions. Vector v is the spacial solution for t=0

    ForwardEuler(&v, dt, dx, Nx, Nt); // Changing v into spacial solution for t=T
    BackwardEuler(&v,dt,dx,Nx,Nt);


    return 0;
}


void ForwardEuler(vec *v, double dt, double dx, int Nx, int Nt){
    double a, a2; vec vnew = zeros<vec>(Nx), V;
    a = dt/dx/dx; a2 = 1-2*a; V = *v; ofstream myfile;

    myfile.open("ForwardEuler.txt"); // Writing results to file. Each line in file is
                                     // the spacial solution for time t.
    // algorithm
    for (int j=0; j<=Nt; j++){
        myfile << V(0) << " ";
        for (int i=1; i<Nx-1; i++){
            vnew(i) = a * V(i-1) + a2 * V(i) + a * V(i+1);
            myfile << vnew(i) << " "; }

        myfile << V(Nx-1) << endl;
        *v = vnew;
    }
    myfile.close();
}

void BackwardEuler(vec *v, double dt, double dx, int Nx, int Nt){
    double a; vec A1 = zeros<vec>(Nx), A2 = zeros<vec>(Nx), A3 = zeros<vec>(Nx);
    a = dt / dx / dx;
    for (int i=0; i<Nx; i++){
        A1(i) = -a; A2(i) = 1 + 2*a; A3(i) = -a;
    }
    for (int j=0; j<Nt; j++){
        *v = tridiagonal(A1,A2,A3,*v,Nx);
    }
}

void CrankNicolson(){

}

vec tridiagonal(vec a1, vec a2, vec a3, vec b,int N){
    // Forward Substitution. Row reducing the matrix equation
    vec v;
    for (int i = 1; i <= N-1; i++){
        float factor = a1[i]/a2[i-1];
        a2[i] = a2[i] - a3[i-1]*factor;
        b[i] = b[i] - b[i-1]*factor;
    }

    // Backward Substitution. Solving the equation for vector v.
    v[N-1] = b[N-1]/a2[N-1];
    for (int k = N-1; k >= 0; k--){
        v[k] = (b[k] - a3[k] * v[k+1])/a2[k];
    }
    return v;
}


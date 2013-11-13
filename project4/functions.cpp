#include "functions.h"


void ForwardEuler(vec *v, double dt, double dx, int Nx, int Nt){
    double a, a2, ui; vec vnew = zeros<vec>(Nx), V;
    a = dt/dx/dx; a2 = 1-2*a; V = *v; ofstream myfile;

    myfile.open("ForwardEuler.txt"); // Writing results to file. Each line in file is
                                     // the spacial solution for time t.
    // algorithm
    for (int j=1; j<=Nt; j++){
        myfile << V(0) + 1 << " "; // Writing in boundary conditions to file.

        for (int i=1; i<Nx-1; i++){
            vnew(i) = a * V(i-1) + a2 * V(i) + a * V(i+1);
            ui = vnew(i) + 1 - i*dx;
            myfile << ui << " "; } // Writing in results to file.

        myfile << V(Nx-1) << endl; // Writing in boundary conditions;
        V = vnew;
    }
    myfile.close();
}

void BackwardEuler(vec *v, double dt, double dx, int Nx, int Nt){
    double a; vec A1 = zeros<vec>(Nx), A2 = zeros<vec>(Nx), A3 = zeros<vec>(Nx), V;
    a = dt / dx / dx; ofstream myfile;
    myfile.open("BackwardEuler.txt");

    for (int i=0; i<Nx; i++){
        A1(i) = -a; A2(i) = 1 + 2*a; A3(i) = -a;
    }
    // Time loop
    V = *v;
    for (int j=1; j<Nt; j++){
        for (int i=0; i<Nx; i++){myfile << V(i) + 1 - i*dx << " ";}; myfile << endl;
        V = tridiagonal(A1,A2,A3,V,Nx);
    }
    myfile.close();
}

void CrankNicolson(vec *v, double dt, double dx, int Nx, int Nt){
    double a,a2,a3; vec A1 = zeros<vec>(Nx), A2 = zeros<vec>(Nx), A3 = zeros<vec>(Nx), vtilde, vnew;
    a = dt / dx / dx; a2 = 2 - 2*a; a3 = 2 + 2*a; ofstream myfile;

    myfile.open("CrankNicolson.txt");
    // Setting the diagonals on of the LHS-matrix.
    for (int i=0; i<Nx; i++){
        A1(i) = -a; A2(i) = a3; A3(i) = -a;
    }

    vnew = *v; vtilde = vnew;
    for (int j=1; j<Nt; j++){
        // Chaning the RHS vector v_old into vtilde.
        for (int i=0; i<Nx; i++){myfile << vnew(i) + 1 - i*dx << " ";}; myfile << endl;
        for (int i=1; i<Nx-1; i++){
            vtilde(i) = a*vnew(i-1) + a2*vnew(i) + a*vnew(i+1);
        }

        vnew = tridiagonal(A1,A2,A3,vtilde,Nx);
    }
    myfile.close();
}

vec tridiagonal(vec a1, vec a2, vec a3, vec b,int Nx){
    // Forward Substitution. Row reducing the matrix equation
    vec v = zeros<vec>(Nx);
    for (int i = 1; i < Nx; i++){
        float factor = a1(i)/a2(i-1);
        a2(i) = a2(i) - a3(i-1)*factor;
        b(i) = b(i) - b(i-1)*factor;
    }

    // Backward Substitution. Solving the equation for vector v.
    v(Nx-1) = b(Nx-1)/a2(Nx-1);
    for (int k = Nx-2; k >= 0; k--){
        v(k) = (b(k) - a3(k) * v(k+1))/a2(k);
    }
    return v;
}

vec Analytical(double dt, double dx, int Nx, int Nt){
    double pi = 3.14159265359; vec u; ofstream myfile;
    myfile.open("Analytical.txt");

    for (int j=0; j <= Nt; j++){
        u = zeros<vec>(Nx);
        for (int i=0; i<Nx; i++){
            for (int n=1; n < 15; n++){
                u(i) += -2/(n*pi)*sin(n*pi*i*dx) * exp(-(n*n*pi*pi*j*dt));
            }; myfile << u(i) + 1 - i*dx << " ";}; myfile << endl; }; myfile.close();
    return u;
}


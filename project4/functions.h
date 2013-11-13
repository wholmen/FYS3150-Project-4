#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#endif // FUNCTIONS_H

#include <iostream>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

void ForwardEuler(vec *v, double dt, double dx, int Nx, int Nt);
void BackwardEuler(vec *v, double dt, double dx, int Nx, int Nt);
void CrankNicolson(vec *v, double dt, double dx, int Nx, int Nt);
vec Analytical(double dt, double dx, int Nx, int Nt);
vec tridiagonal(vec a1, vec a2, vec a3, vec b,int Nx);

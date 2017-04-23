#ifndef VMC_H
#define VMC_H
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
using namespace arma;
using namespace std;





class VMC
{
public:
    VMC(int n, int cycles, double step, double w);
    double wavefunction (mat &r);
    double localEnergyder( mat &r);
    double localEnergyana(mat &r);
    void MCH();
    double Green(mat &r1, mat &r2);
    void test();
    void MonteCarlo();
    mat m_rold;
    mat m_rnew;
    vec SteepestDescend(mat A, vec b, vec x0);
    int m_cycles;
    int m_nelectrons;
    double m_w;
    double m_step;
    double m_alpha;
    double m_beta;
    double m_dt;
};

#endif // VMC_H

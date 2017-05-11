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
    void findoptParameter();//7777777777777777777777777777777777777777777777777777777777777777777777
    double EnergyCorrection();
    vec Quantumforce( mat &r);
    void MCH(int cycles);
    double Green(mat &r1, mat &r2);
    void test();
    void MonteCarlo();
    mat m_rold;
    mat m_rnew;
    void SteepestDescend(mat A, vec b, vec x0);
    void resetvalues();
    int m_cycles;
    int m_nelectrons;
    double m_w;
    double m_step;
    double m_alpha;
    double m_a;
    double m_beta;
    double m_dt;
    bool m_insteepestdesc;
    double m_localE;
    double m_localE2;
    vec m_dE;
};

#endif // VMC_H

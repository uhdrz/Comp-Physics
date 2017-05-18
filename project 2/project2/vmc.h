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
    double spwf (double x, double y, int nx, int ny);
    double wavefunction1 (mat &r);
    double SlaterDet(mat &r);
    double hermite(int n, double x);
    double localEnergy(mat &r);

    vec GradJastrow( mat &r); //
    vec GradSP0(mat &r);//
    vec GradSP1(mat &r);//
    vec GradSP2(mat &r);//
    double LapJastrow( mat &r);//
    double LapSP0(mat &r);//
    double LapSP1(mat &r);//
    double LapSP2(mat &r);//
    double LapDet0(mat &r);//
    vec postonum(int i);


    double localEnergyder( mat &r);
    double localEnergyana(mat &r);
    void findoptParameter();
    double relDis(mat &r, int i, int j);
    double pos2(mat &r, int i);
    mat Quantumforce( mat &r);
    void MCH();
    void test();
    void MonteCarlo();
    mat m_rold;
    mat m_rnew;
    int m_cycles;
    int m_nelectrons;
    double m_w;
    double m_step;
    double m_a;
    double m_dt;
    vec m_localEn;
    vec m_varpar;
};

#endif // VMC_H


inline double VMC::pos2(mat &r, int i){
    double r2 = 0;
    for (int j = 0; j<2; j++) {
        r2+=r(i,j)*r(i,j);
    }
    return r2;


}


inline double VMC::relDis(mat &r, int i, int j){
    double r_ij = 0;
    for (int d = 0; d <2; d++){
        r_ij += (r(i,d)-r(j,d))*(r(i,d)-r(j,d));
          }
    return sqrt(r_ij);

}

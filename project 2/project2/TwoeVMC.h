#ifndef TwoeVMC_H
#define TwoeVMC_H
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
using namespace arma;
using namespace std;





class TwoeVMC
{
public:
    TwoeVMC(int n, int cycles, double step, double w);
    double wavefunction (mat &r);
    double localEnergyder( mat &r);
    double localEnergyana(mat &r);
    void findoptParameter();
    double relDis(mat &r, int i, int j);
    double pos2(mat &r, int i);
    vec Quantumforce(int i, mat &r);
    void MCH(int c,char **v);
    void test();
    void MCbrute();
    //mat m_rold;
    //mat m_rnew;
    int m_cycles;
    int m_nelectrons;
    double m_w;
    double m_step;
    double m_a;
    double m_dt;
    vec m_localEn;
    vec m_varpar;
};

#endif // TwoeVMC_H


inline double TwoeVMC::pos2(mat &r, int i){
    double r2 = 0;
    for (int j = 0; j<2; j++) {
        r2+=r(i,j)*r(i,j);
    }
    return r2;


}


inline double TwoeVMC::relDis(mat &r, int i, int j){
    double r_ij = 0;
    for (int d = 0; d <2; d++){
        r_ij += (r(i,d)-r(j,d))*(r(i,d)-r(j,d));
          }
    return sqrt(r_ij);

}

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
    double wavefunction (mat &r,mat &Slaterup, mat &Slaterdown);
    double spwf (double x, double y, int nx, int ny);
    double localEnergy(mat &r, mat &InvUp, mat &InvDown);
    double locKin(mat &r, mat &InvUp, mat &InvDown);
    double locPot(mat &r);


    vec GradJastrow(int k, mat &r); //
    vec GradSP(int i, mat &r, mat &InvUp, mat &InvDown);
    vec DerivSP(int i, int nx, int ny, mat &r);
    double Deriv2SP(int i, int nx, int ny, mat &r);
    double DerivSPA(int i,int nx,int ny, mat &r);

    double LapJastrow(int k, mat &r);//
    double LapSP(int i , mat &r, mat &InvUp, mat &InvDown);
    vec postonum(int i);
    mat Slatermatrixdown(mat &r);
    mat Slatermatrixup(mat &r);
    mat SlaterUpInv(mat &Slaterup);
    mat SlaterDownInv(mat &Slaterdown);
    double derivJastrow(int i, int j, double rij);
    double deriv2Jastrow(int i, int j, double rij);
    void MCbrute();
    double DerivAlpha( mat &r, mat &InvUp, mat &InvDown);
    double DerivBeta( mat &r);



    void findoptParameter();
    double relDis(mat &r, int i, int j);
    double pos2(mat &r, int i);
    vec Quantumforce(int i, mat &r, mat &InvUp, mat &InvDown);
    void MCH(int c, char **v);
    int m_cycles;
    int m_nelectrons;
    double m_w;
    double m_step;
    vec m_a;
    double m_dt;
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

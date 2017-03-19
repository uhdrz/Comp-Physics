#ifndef HFCALC_H
#define HFCALC_H
#include "hfcalc.h"
#include "Coulomb_Functions.hpp"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

class HFcalc
{
public:
    HFcalc();
    HFcalc(int n, int level, double w);
    int m_natoms; //number of electrons
    int m_levels;
    int m_numofpossstates;
    vec m_obelements;
    vec m_hfenergy;
    double m_w;
    mat m_C;
    mat m_rho;
    mat m_HF;
    int degeneracy(int n);
    void compDensity();
    void compFock();
    double compHFenergy(int n, bool test);
    void diagonalize();
    vector<double> codec(int i);
    void codectest();

};

#endif // HFCALC_H

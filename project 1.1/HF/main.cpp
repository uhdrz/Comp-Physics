#include "hfcalc.h"
#include "Coulomb_Functions.hpp"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

int main(/*int argc, char *argv[]*/)
{

    HFcalc a(6,4,1);

    /*
    for(int i=0;i<12;i++){
    for(int j=0;j<3;j++){
        cout << a.codec(i)[j];
        cout << endl;

        }
     cout << "******" << endl;
    }
    */

    for(int i=0; i<10; i++) {
        a.compDensity();
        a.compFock();
        a.diagonalize();
    }

    cout << a.Energy() << endl;

    //a.m_HF.print();
    //a.m_hfenergy.print();

    //a.diagonalize();
    //double E=a.compHFenergy(100);

    //cout<<E << endl;
}


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




    double E=a.compHFenergy(100);

    cout<<E << endl;


}


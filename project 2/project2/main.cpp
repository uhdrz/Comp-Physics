#include "vmc.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
#include <fstream>
using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{






    VMC a(2,1000000,1.3,1);

    //A.print();






     //a.MCH();
   // a.MonteCarlo();
    a.findoptParameter();



    return 0;




}

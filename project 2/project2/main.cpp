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

//binary writing

   /* ofstream outFile;
    outFile.open("localEnergies.dat", ios::out | ios::binary);
    if (! outFile.is_open()) {
        cout << "Could not open file." << endl;
    }

    //for (MC cycle) ..
    for (int i = 0 ; i < 1000; i++) {
        outFile.write( (char*)&localE, sizeof(double));
    }

    fflush(stdout);
    outFile.close();*/




    VMC a(2,1000000,1.3,1);

    //A.print();






     a.MCH();
   // a.MonteCarlo();
    //a.findoptParameter();



    return 0;




}

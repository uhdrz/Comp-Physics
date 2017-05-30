#include "TwoeVMC.h"
#include "vmc.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
#include <fstream>
#include "mpi.h"
using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{       int c=argc;
        char** v=argv;



   TwoeVMC b(2,10000,1.3,1);
  // b.MCH(c, v);
  // b.MonteCarlo();
   VMC a(2,1e6,1.3,1);
   a.MCH();

    cout << endl;
    /*mat A(2, 2, fill::randu);
    VMC a(2,1000000,1.3,1);
    TwoeVMC b(2,10000,1.3,1);
    double wf1=a.wavefunction(A);
    double wf2=b.wavefunction(A);


    cout<< wf1<<" "<< wf2;*/
    //A.print();
     //b.MCH();
   // a.MonteCarlo();
    //b.findoptParameter();



    return 0;




}

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
#include <ctime>
using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{       int c=argc;
        char** v=argv;

    //clock_t t1,t2;
    //t1=clock();

   TwoeVMC b(2,1e6,1.3,0.01);

    //b.findoptParameter();
    b.MCH(c, v);
    //b.MonteCarlo();
   //VMC a(2,100000,1.3,1);
   //a.MCH();

    //cout << endl;
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

    //t2=clock();
    //float diff ((float)t2-(float)t1);
    //float seconds = diff / CLOCKS_PER_SEC;
     //cout<<seconds<<endl;




    return 0;




}

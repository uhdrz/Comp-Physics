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
{

    int c=argc;
    char **v=argv;



    //clock_t t1,t2;
    //t1=clock();



    //TwoeVMC b(2,1e6,1.3,0.5);
    //b.findoptParameter();
    //b.MCH(c,v);




    VMC a(12,500000,1.3,1);
    a.findoptParameter();
    //a.MCH(c, v);






    //t2=clock();
    //float diff ((float)t2-(float)t1);
    //float seconds = diff / CLOCKS_PER_SEC;
     //cout<<seconds<<endl;



    return 0;




}

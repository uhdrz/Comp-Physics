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







    //TwoeVMC b(2,1e6,1.3,0.5);
    //b.findoptParameter();
    //b.MCH(c,v);




    VMC a(6,1000000,1.3,1);
    //a.findoptParameter();
    a.MCH(c, v);









    return 0;




}

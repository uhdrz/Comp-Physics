#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <omp.h>
#include <hattree.h>

using namespace std;

int main(/*int argc, char *argv[]*/)
{
    hattree q(1,0,1,1.0,5);
    hattree p(1,0,1,1.0,5);
    hattree r(0,0,1,1.0,5);
    hattree s(0,0,1,1.0,5);
    double n=q.GaussHermiteIntegration(4);
    cout <<n;








}



//freq? var trafo???

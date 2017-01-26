#include "system.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>

using namespace std;

double system::wavefunc(double x, double w,double y, double A)
{
    double phi= A*hermite(m_nx,sqrt(w)*x)* hermite(m_ny,sqrt(w)*y)*exp(-w*0.5*(pow(x,2)+pow(y,2)));
    return phi;
}

system::system (int c)
{
    cut=c;
}

int system::numtopos()
{   int pos=0;
    int n = m_nx+m_ny;
    int temp=0;
    for (int i =0; i<n; i++)
    {
        temp+= 2*(i+1);
    }
    pos=temp+m_ny+m_s;
    return pos;
}


void system::postonum(int pos)
{
    m_s = pos%2;
    int temp=0;
    int n=0;
    for(int i=0; temp< pos; i++)
    {
      temp+= 2*(i+1);
      n=i;
    }

    m_ny=pos-m_s-temp;
    m_nx=n-m_ny;


    //return m, nx,ny;
}


double system::hermite(int n, double x)
{
    if (n==0)
    {return 1;}
    else if (n==1)
    {return x;}
    else
    {   double temp1=1;
        double temp2=x;
        double temp=0;
        for(int i=2; i<n; i++)
        {   temp=2*x*temp2-2*(i-1)*temp1;
            temp1=temp2;
            temp2=temp;
        }
        return temp;


    }


}


int system::cutoff(int n1, int n2, int cut_test)
{
    if((n1+n2<=cut_test))
    {
        return 1;
    }
    else {return 0;};
}


void system::get(int x, int y, int m)
{ int test= cutoff(x,y, cut);
    if(test==1){
    m_nx=x;
    m_ny=y;
    m_s=m;}
    else{cout<<"CUTOFF!!!";}
}



void system::genstate()
{
    int temp=0;
    for (int i =0; i<cut; i++)
    {
        temp+= 2*(i+1);
    }


}

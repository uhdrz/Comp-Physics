#include "system.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>

using namespace std;

double system::wavefunc(double x, double y,double w, double A)
{
    double phi= A*hermite(m_nx,sqrt(w)*x)* hermite(m_ny,sqrt(w)*y)*exp(-w*0.5*(pow(x,2)+pow(y,2)));
    return phi;
}

system::system (int c)
{
    cut=c;
    pos=0;
    m_nx=0;
    m_ny=0;
    m_s=0;
}

void system::numtopos()
{   int n = m_nx+m_ny;
    int temp=0;
    for (int i =0; i<n; i++)
    {
        temp+= 2*(i+1);
    }
    pos=temp+m_ny+m_s;
}


void system::postonum()
{
    m_s = pos%2;
    int temp1=0;
    int temp2=0;
    int n;
    for(n=0; temp1< pos; n++)
    { temp2=temp1;
      temp1+= 2*(n+1);
     }

    m_ny=pos-m_s-temp2;
    m_nx=n-1-m_ny;


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




void system::get(int x, int y, int m)
{    if((x+y<=cut))
    {
        m_nx=x;
        m_ny=y;
        m_s=m;
    }
     else{cout<<"CUTOFF!!!";}
}



void system::genstate()
{
    int temp=0;
    for (int i =0; i<cut; i++)
    {
        temp+= 2*(i+1);
    };
    state{new int[temp]};
    for(int i=0;i<temp; i++)
    {
        state[i]=0;
    };
    state[pos]=1;


}

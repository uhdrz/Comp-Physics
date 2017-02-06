#include "hattree.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace std;

hattree::hattree()
{cout<< "Hello";}

long hattree::fac(int n)
{
    if(n==1)return 1;
    else return n*fac(n-1);
}

double hattree:: energy()
{
    return m_w*(m_nx+m_ny+1);
}
double hattree::wavefunc(double x, double y)
{
    double phi=sqrt(m_w/M_PI)*1.0/(sqrt(pow(2,m_nx)*fac(m_nx))) *1.0/(sqrt(pow(2,m_ny)*fac(m_ny)))*hermite(m_nx,sqrt(m_w)*x)* hermite(m_ny,sqrt(m_w)*y)*exp(-m_w*0.5*(pow(x,2)+pow(y,2)));
    return phi;
}

hattree::hattree (int n_x, int n_y, int s,double w,int c)
{   cut=c;
    pos=0;
    m_w=w;
    m_nx=n_x;
    m_ny=n_y;
    m_s=s;
}

void hattree::numtopos()
{   int n = m_nx+m_ny;
    int temp=0;
    for (int i =0; i<n; i++)
    {
        temp+= 2*(i+1);
    }
    pos=temp+m_ny+m_s;
    cout<<pos;
}


void hattree::postonum()
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
    cout<<m_nx<<m_ny<<m_s;


    //return m, nx,ny;
}


double hattree::hermite(int n, double x)
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




void hattree::get(int x, int y, int m)
{    if((x+y<=cut))
    {
        m_nx=x;
        m_ny=y;
        m_s=m;
    }
     else{cout<<"CUTOFF!!!";}
}



void hattree::genstate()
{
    int temp=0;
    for (int i =0; i<cut; i++)
    {
        temp+= 2*(i+1);
    };
    state.reserve(temp);
    for(int i=0;i<temp; i++)
    {
        state[i]=0;
    };
    state[pos]=1;


}

void hattree::print()
{
    cout<<"Single Quantum State :"<<endl;
    cout<<"n_x="<<m_nx<<"  n_y="<<m_ny<<endl;
    cout<<"s_z ="<<m_s<< "Cut-off energy="<<cut<<endl;
    cout<<"Position"<<pos;
}


double hattree::potential (double x1, double x2, double y1, double y2)
{
    return 1.0/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


void hattree::GaussHermiteQuadrature(double *x, double *w, int n)
{
  int i,its,j,m;
  double p1,p2,p3,pp,z,z1;
  double Epsilon = 3.0e-14, PIM4 = 0.7511255444649425;
  int MaxIterations = 10;
  m=(n+1)/2;
  for (i=1;i<=m;i++) {
    if (i == 1) {
      z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if (i == 2) {
      z -= 1.14*pow((double)n,0.426)/z;
    } else if (i == 3) {
      z=1.86*z-0.86*x[0];
    } else if (i == 4) {
      z=1.91*z-0.91*x[1];
    } else {
      z=2.0*z-x[i-3];
    }
    for (its=1;its<=MaxIterations;its++) {
      p1=PIM4;
      p2=0.0;
      for (j=1;j<=n;j++) {
    p3=p2;
    p2=p1;
    p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
      }
      pp=sqrt((double)2*n)*p2;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= Epsilon) break;
    }
    if (its > MaxIterations) cout << "too many iterations in Hermite quadrature" << endl;
    x[i-1]=z;
    x[n-i] = -z;
    w[i-1]=2.0/(pp*pp);
    w[n-i]=w[i-1];
  }
} // Gaussian quadrature weights and integration points

double hattree::Integration(int n,hattree p, hattree q, hattree r, hattree s)
{   double temp= pow(1/M_PI,2);
    double *x=new double [n+1];
    double *w= new double [n];
    double sol =0;
    GaussHermiteQuadrature(x,w,n);
    for(int i=0; i<n; i++)
    {
        for(int j=0;j<n;j++)
        {
            for (int k=0;k<n;k++)
            {
                for(int l=0;l<n;l++)
                {
                    sol+=potential(x[i],x[j], x[k], x[l])*w[i]*w[j]*w[k]*w[l]*
                            p.hermite(p.m_nx,x[i])*p.hermite(p.m_ny,x[k])*q.hermite(q.m_nx,x[j])*q.hermite(q.m_ny,x[l])*r.hermite(r.m_nx,x[i])*r.hermite(r.m_ny,x[k])*s.hermite(s.m_nx,x[j])*s.hermite(s.m_ny,x[l])
                            *temp*sqrt(p.m_w)*sqrt(q.m_w)*sqrt(r.m_w)*sqrt(s.m_w)*1/(sqrt(pow(2,p.m_nx)*fac(p.m_nx))) *1/(sqrt(pow(2,p.m_ny)*fac(p.m_ny)))*1/(sqrt(pow(2,q.m_nx)*fac(q.m_nx))) *1/(sqrt(pow(2,q.m_ny)*fac(q.m_ny)))*1/(sqrt(pow(2,r.m_nx)*fac(r.m_nx))) *1/(sqrt(pow(2,r.m_ny)*fac(r.m_ny)))*1/(sqrt(pow(2,s.m_nx)*fac(s.m_nx))) *1/(sqrt(pow(2,s.m_ny)*fac(s.m_ny)));
                }
            }
        }
    }

/*freqency w??? normalisation?*/

    delete [] x;
    delete [] w;
    return sol;

}


double hattree::ASMatrixElement(int n, hattree p, hattree q, hattree r, hattree s)
{
    return Integration(n,p,q,r,s)-Integration(n,p,p,s,r);
}



#include "vmc.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
using namespace arma;
using namespace std;








VMC::VMC(int n, int cycles, double step, double w)
{
    m_nelectrons=n;
    m_cycles=cycles;
    m_step=step;
    m_w=w;
    m_alpha=1;
    m_beta=1;
    m_rold.zeros(m_nelectrons,2);
    m_rnew.zeros(m_nelectrons,2);
    m_dt=0.01;


}



double VMC::wavefunction(mat &r){
    //first index particle, second dimension
    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            sum+=r(i,j)*r(i,j);
        }

    }


    double phi= exp(-0.5*m_w*m_alpha*sum );

    return phi;


}




double VMC::localEnergyana(mat &r){


    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            sum+=r(i,j)*r(i,j);
        }

    }
    double abs=0;

    for(int i=0; i<m_nelectrons-1;i++){
        for(int j=i+1; j<m_nelectrons;j++){
            double temp=0;
            for(int d=0; d<2; d++){
                temp+=(r(i,d)-r(j,d))*(r(i,d)-r(j,d));
            }
             abs+=temp;
            }
        }

    double overr=sqrt(1/abs);







    double Elocal=0.5*m_w*m_w*(1-m_alpha*m_alpha)*sum+2*m_alpha*m_w+overr;
    return Elocal;

}




double VMC::localEnergyder(mat &r){
    double h=0.001;

    mat rp=zeros<mat>(m_nelectrons,2);
    mat rm=zeros<mat>(m_nelectrons,2);
    rp=r;
    rm=r;
    double wfm=0;
    double wfp=0;
    double wfc=wavefunction(r);
    double deriv=0;



    for (int i=0; i< m_nelectrons;i++){
        for(int j=0; j<2;j++){
            rp(i,j)+=h;
            rm(i,j)-=h;
            wfm=wavefunction(rm);
            wfp=wavefunction(rp);
            deriv+=(wfm+wfp-2*wfc);
            rp(i,j)=r(i,j);
            rm(i,j)=r(i,j);

        }}

    double Elocal=-0.5*deriv/(wfc*h*h);

    return Elocal;

}


void VMC::test(){
    mat r=zeros<mat>(m_nelectrons,2);
    double Elder=localEnergyder(r);
    double Elana=localEnergyana(r);
    cout<<Elder<<endl;
    cout<<Elana<<endl;



}


void VMC::SteepestDescend(mat A, vec b, vec x0){
    int IterMax=10000;
    int i=0;
    const double tolerance = 1.0e-14;
    double aold=0;
    double anew=0;
    double bold=0;
    double bnew=0;
    double diffa=0;
    double diffb=0;

    while(diffa>tolerance & diffb>tolerance||i<IterMax){ /////////////////////


        anew=aold+m_step*Ableitung(aold);
        bnew=bold+m_step*Ableitung(bold);

        diffa=anew-aold;
        diffb=bnew-bold;
        i++;
    }

    m_alpha=anew;
    m_beta=bnew;


}












/*void VMC::MCH(){
    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <>dis(0,1);
    dis(gen);


    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;
    double sum=0;
    double sumsquared=0;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            rold(i,j)=dis(gen)*sqrt(m_dt);
        }
    }

    rnew=rold;
    for(int n=0;n<m_cycles;n++){
        wfold=wavefunction(rold);
        vec F=2*Ableitung(rold); ///////////////////////////////////////////////////////////
        for(int i=0;i<m_nelectrons;i++){
                    for(int j=0;j<2;j++){
                        rnew(i,j)=rold(i,j)+0.5*F(j)*m_dt+dis(gen)*sqrt(m_dt);

                  }

            wfnew=wavefunction(rnew);
            double q=(Green(rold, rnew)*wfnew*wfnew)/(wfold*wfold*Green(rnew, rold));
            if(dis(gen)<=q){
                for(int j=0;j<2;j++){
                    rold(i,j)=rnew(i,j);
                    wfold=wfnew;
                }
            }
            else{
                for(int j=0; j<0;j++){
                    rnew(i,j)=rold(i,j);
                }
            }
            double temp=localEnergyana(rnew);
            sum+=temp;
            sumsquared+=temp*temp;

        }

    }

    double energy=sum/(m_cycles*m_nelectrons);
    double energysquared=sumsquared/(m_cycles*m_nelectrons);
    double variance=energysquared- energy*energy;




    cout<< "Energy"<< energy<<endl;
    cout<< "Variance"<< variance<<endl;



















}*/



/*double VMC::Green(mat &r1, mat &r2){

    vec F=2*Ableitung(rold);
    double constant =pow(2*M_PI*m_dt,1.5*m_nelectrons);
    double dr2=0;
    for(int i=0; i<m_nelectrons-1;i++){
            double temp=0;
            for(int d=0; d<2; d++){
                temp+=(r1(i,d)-r2(i,d)-0.5*m_dt*F(d))*(r1(i,d)-r2(i,d)-0.5*m_dt*F(d));
            }
             dr2+=temp;

        }



    double G=(1/constant)*exp(-(dr2/(2*m_dt)));
    return G;




}*/

void VMC::MonteCarlo(){


    random_device rnd;
    mt19937 gen(rnd());
    uniform_real_distribution <>dis(0,1);



    dis(gen);



    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;
    double sum=0;
    double sumsquared=0;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            rold(i,j)=m_step*dis(gen);
        }
    }

    rnew=rold;
    for(int n=0;n<m_cycles;n++){
        wfold=wavefunction(rold);
        for(int i=0;i<m_nelectrons;i++){
            for(int j=0;j<2;j++){
                rnew(i,j)=rold(i,j)+m_step*dis(gen);

          }

            wfnew=wavefunction(rnew);
            if(dis(gen)<=(wfnew*wfnew)/(wfold*wfold)){
                for(int j=0;j<2;j++){
                    rold(i,j)=rnew(i,j);
                    wfold=wfnew;
                }
            }
            else{
                for(int j=0; j<0;j++){
                    rnew(i,j)=rold(i,j);
                }
            }
            double temp=localEnergyana(rnew);
            sum+=temp;
            sumsquared+=temp*temp;

        }

    }

    double energy=sum/(m_cycles*m_nelectrons);
    double energysquared=sumsquared/(m_cycles*m_nelectrons);
    double variance=energysquared- energy*energy;




    cout<< "Energy"<< energy<<endl;
    cout<< "Variance"<< variance<<endl;



}















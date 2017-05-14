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
    m_rold.zeros(m_nelectrons,2);
    m_rnew.zeros(m_nelectrons,2);
    m_dt=0.01;
    m_a=1;
    m_localEn=zeros<vec>(2);
    m_varpar={1,0.3}; //0=alpha, 1 =beta




}


double VMC::wavefunction(mat &r){
    //first index particle, second dimension
    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            sum+=r(i,j)*r(i,j);
        }

    }


    double phi= exp(-0.5*m_w*m_varpar(0)*sum );

    return phi;


}





mat VMC::Quantumforce(mat &r){
    mat F=zeros(m_nelectrons,2);


    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            double spart=r(i,j)*m_varpar(0);

            double jpart=0;

            for(int k=0;k<m_nelectrons;k++){
                if(i!=k){

                    double r_ik=relDis(r,i, k);

                    jpart +=(r(i,j)-r(k,j))/(r_ik*pow(1+m_varpar(1)*r_ik,2));

               }
              }

            F(i,j)=2*(spart+jpart);


        }
    }



     return F;
}

double VMC::localEnergyana(mat &r){


    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        sum+=pos2(r,i);
    }

    double abs=0;
    for(int i=0; i<m_nelectrons-1;i++){
        for(int j=i+1; j<m_nelectrons;j++){
            abs+=relDis(r,i,j);
            }
        }

    double overr=1/abs; //careful for nele<2







    double Elocal=0.5*m_w*m_w*(1-m_varpar(0)*m_varpar(0))*sum+2*m_varpar(0)*m_w+overr;
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









void VMC::MCH(){
    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <double>norm(0,1);
    uniform_real_distribution <double>uni(0,1);




    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;
    double sum=0;
    double sumsquared=0;








    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            rold(i,j)=norm(gen)*sqrt(m_dt);
        }
    }

    rnew=rold;
    for(int n=0;n<m_cycles;n++){
        wfold=wavefunction(rold);
        mat Qforceold=Quantumforce(rold);
        for(int i=0;i<m_nelectrons;i++){
            for(int j=0;j<2;j++){
                rnew(i,j)=rold(i,j)+0.5*Qforceold(i,j)*m_dt+norm(gen)*sqrt(m_dt);

          }

            wfnew=wavefunction(rnew);
            mat Qforcenew=Quantumforce(rnew);
            double G=0;
            for(int j=0; j<2;j++){
                G += 0.5*(Qforcenew(i,j)+Qforceold(i,j))*(0.25*m_dt*(Qforceold(i,j)-Qforcenew(i,j))-rnew(i,j)+rold(i,j));

            }
            G=exp(G);

            if(uni(gen)<=G*(wfnew*wfnew)/(wfold*wfold)){

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

            double r12=relDis(rnew, 0,1);
            double den=1/(1+m_varpar(1)*r12);
            double l1=.5*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))*(1-m_varpar(0)*m_varpar(0))+2*m_varpar(0)*m_w+1/r12;
            double l2=2*m_varpar(0)*m_varpar(0)*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))-4*m_w*m_varpar(0)-2*m_a*m_varpar(0)*m_w*r12*den*den+2*m_a*den*den*(m_a*den*den+1/r12-2*m_varpar(1)*den);

            double temp=l1+l2;
            sum+=temp;
            sumsquared+=temp*temp;

        }

    }
    int samples=m_cycles*m_nelectrons;
    m_localEn(0)=sum/samples;
    m_localEn(1)=sumsquared/samples;

    double variance=m_localEn(1)-m_localEn(0)*m_localEn(0);
    m_localEn.print();
    cout<<variance;



}





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






void VMC::findoptParameter(){
    int updates=100;
    int cycles=10000;
    vec dE=zeros<vec>(2);
    double tolerance = 1.0e-14;
    vec parold={1,0};
    vec parnew=zeros<vec>(2);
    double diff=0;



    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <double>norm(0,1);
    uniform_real_distribution <double>uni(0,1);




    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;

    vec sumPsi=zeros<vec>(5); //0=en, 1,2=psi
    vec expVal=zeros<vec>(5); //0=en, 1,2=psi

    for( int m=0; m<updates; m++ ){

        for(int i=0;i<m_nelectrons;i++){
                for(int j=0;j<2;j++){
                    rold(i,j)=norm(gen)*sqrt(m_dt);
                }
            }

        rnew=rold;
        for(int n=0;n<m_cycles;n++){
            wfold=wavefunction(rold);
            mat Qforceold=Quantumforce(rold);
            for(int i=0;i<m_nelectrons;i++){
                for(int j=0;j<2;j++){
                    rnew(i,j)=rold(i,j)+0.5*Qforceold(i,j)*m_dt+norm(gen)*sqrt(m_dt);

              }

                wfnew=wavefunction(rnew);
                mat Qforcenew=Quantumforce(rnew);
                double G=0;
                for(int j=0; j<2;j++){
                    G+= 0.5*(Qforcenew(i,j)+Qforceold(i,j))*(0.25*m_dt*(Qforceold(i,j)-Qforcenew(i,j))-rnew(i,j)+rold(i,j));

                }
                G=exp(G);

                if(uni(gen)<=G*(wfnew*wfnew)/(wfold*wfold)){

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

                double r12=relDis(rnew, 0,1);
                double den=1/(1+m_varpar(1)*r12);
                double l1=.5*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))*(1-m_varpar(0)*m_varpar(0))+2*m_varpar(0)*m_w+1/r12;
                double l2=2*m_varpar(0)*m_varpar(0)*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))-4*m_w*m_varpar(0)-2*m_a*m_varpar(0)*m_w*r12*den*den+2*m_a*den*den*(m_a*den*den+1/r12-2*m_varpar(1)*den);

                double temp=l1+l2;
                double tempPsia=-0.5*m_w*(pos2(rnew,0)+pos2(rnew,1));
                double tempPsib=-m_a*r12*r12*den*den;
                sumPsi(0)+=temp;
                sumPsi(2)+=tempPsib;
                sumPsi(1)+=tempPsia;
                sumPsi(3)+=tempPsib*temp;
                sumPsi(4)+=tempPsia*temp;




            }





        }

        expVal=sumPsi/(cycles*m_nelectrons);



        dE(0)=2*(expVal(3)-expVal(1)*expVal(0));
        dE(1)=2*(expVal(4)-expVal(2)*expVal(0));
        parnew=parold-m_step*dE;


        m_varpar=parnew;
        diff= sqrt(dot(parnew-parold,parnew-parold));
        if(diff<tolerance){ break;}
        else{parold=parnew;}









}


    cout<<m_varpar(0)<<endl;
    cout<<m_varpar(1);




}












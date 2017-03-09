#include "hfcalc.h"
#include "Coulomb_Functions.hpp"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


HFcalc::HFcalc()
{

}



HFcalc::HFcalc(int n, int level, double w)
{   m_natoms=n;
    m_levels=level;
    m_w=w;
    int temp=0;
    int i=0;
    for (; temp<m_natoms; i++)
        {
            temp+=2*(i+1);
        };
    m_fermilevel=i-1;  //starts with 0, so i energy level


    m_numofpossstates=m_levels*(m_levels+1);


    m_obelements.zeros(m_numofpossstates);
    m_hfenergy.zeros(m_numofpossstates);


    for(int i=0;i<m_natoms; i++)
      {
          if(i<=1){
              m_obelements(i)=m_w;}
          else if(i<6 && i>=2){
              m_obelements(i)=2*m_w;}
          else if(i<12 && i>=6){
              m_obelements(i)=3*m_w;
          }
          else if(i<20 && i>=12){
              m_obelements(i)=4*m_w;
          }
          else  m_obelements(i)=0;

      };




    m_C.eye(m_numofpossstates,m_numofpossstates);
    m_rho.zeros(m_numofpossstates,m_numofpossstates);
    m_HF.zeros(m_numofpossstates,m_numofpossstates);








}


vector<double> HFcalc::codec(int i){
    int temp1=0;
    int temp2=0;
    int m=0;
    int si=i%2;           //0up, 1 down
    for (int n=0; temp1<=i; n++)
        {
        temp2=temp1;
        temp1+= 2*(n+1);
        m=n;

        };

    int leveli=m;               //count from 0
    int posinshell=i-temp2;
    int ni=0;

    int mi=-leveli;

    for(int o=1;o<=0.5*posinshell;o++){
        mi+=2;
        if((2*ni+abs(mi))<leveli){
            ni+=1;
        };
        if((2*ni+abs(mi))>leveli){
            ni-=1;
        };
    }

    vector<double> numbers;
    numbers.reserve(3);
    numbers[0]=ni;
    numbers[1]=mi;
    numbers[2]=si;

    return numbers;


}








double HFcalc::compRefenergy(){
    double E=0;
    double w=m_w;






    for(int i=0;i<m_natoms;i++){
        for(int j=0;j<m_natoms;j++){
            int ni=codec(i)[0];
            int mi=codec(i)[1];
            int nj=codec(j)[0];
            int mj=codec(j)[1];




           E+=0.5*(Coulomb_HO(w,ni,ni,nj,mj,ni,mi,nj,mj)-Coulomb_HO(w,ni,ni,nj,mj,nj,mj,ni,mi));
        }
        E+=m_obelements(i);
    };



return E;


        }





void HFcalc::compDensity(){
    //doppelte durch spin??? sieh ref energy....

    for(int g=0;g<m_numofpossstates;g++){
        for(int  d=0;d<m_numofpossstates;d++){

            double sum=0;
            for(int i=0; i<m_natoms;i++){
                sum+=m_C(i,g)*m_C(i,d);
            }
            m_rho(g,d)=sum;




}
}

}


void HFcalc::compFock(){
    double w=m_w;
    for(int a=0;a<m_numofpossstates;a++){
        int na=codec(a)[0];
        int ma=codec(a)[1];
        int sa=codec(a)[2];
        for(int b=0;b<m_numofpossstates;b++){
            int nb=codec(b)[0];
            int mb=codec(b)[1];
            int sb=codec(b)[2];
            double sum=0;

            for(int g=0;g<m_numofpossstates;g++){
                int ng=codec(g)[0];
                int mg=codec(g)[1];
                int sg=codec(g)[2];
                for(int d=0;d<m_numofpossstates;d++){

                    int nd=codec(d)[0];
                    int md=codec(d)[1];
                    int sd=codec(d)[2];





                    if(ma+mg==mb+md && sa==sb && sg==sd && sa+sg==sd+sb ){

                        sum +=m_rho(g,d)*Coulomb_HO(w,na,ma,ng,mg,nb,mb,nd,md);

                    }
                    if(ma+mg==mb+md && sa==sd && sg==sb && sa+sg==sd+sb  ){
                        sum -=m_rho(g,d)*Coulomb_HO(w,na,ma,ng,mg,nd,md,nb,mb);

                    }

                }//if(a==0 && b==8){cout<< sum <<' ';}

            }


            if(a==b){
                sum+=m_obelements(a);
            }
            m_HF(a,b)=sum;
            cout<<m_HF(a,b)<<' ';


        }
        cout<<endl;

    }

}




void HFcalc::diagonalize(){
    eig_sym(m_hfenergy,m_C,m_HF);
    /*for(int i=0; i<m_numofpossstates;i++){
    cout<<m_hfenergy(i)<<' ';}*/


}

double HFcalc::compHFenergy(int n){
    double lambda= pow(10,-8);
    int counter=0;
    double deltaE=100;

    vec temp=zeros<vec>(m_numofpossstates);


    while(deltaE/m_numofpossstates>lambda && counter<n){


        deltaE=0;
        compDensity();
        compFock();
        for(int i=0;i<m_numofpossstates;i++){
        temp(i)=m_hfenergy(i);};
        diagonalize();
        counter++;





       for(int i=0; i<m_numofpossstates;i++){
          deltaE+=m_hfenergy(i)-temp(i);
        };
    };

   /* for(int i=0; i<m_numofpossstates;i++){
        cout<<m_hfenergy(i)<<' ';
    }*/














    double E0hf=0;
    double w=m_w;
    for(int i=0; i<m_natoms;i++){
    for (int j=0; j<m_natoms;j++){

    for(int a=0;a<m_numofpossstates;a++){
        int na=codec(a)[0];
        int ma=codec(a)[1];
        int sa=codec(a)[2];
        for(int b=0;b<m_numofpossstates;b++){
            int nb=codec(b)[0];
            int mb=codec(b)[1];
            int sb=codec(b)[2];
            for(int g=0;g<m_numofpossstates;g++){
                int ng=codec(g)[0];
                int mg=codec(g)[1];
                int sg=codec(g)[2];
                for(int d=0;d<m_numofpossstates;d++){

                    int nd=codec(d)[0];
                    int md=codec(d)[1];
                    int sd=codec(d)[2];


                    if(ma+mb==mg+md && sa==sg && sb==sd){
                         E0hf-=0.5*m_C(i,a)*m_C(i,b)*m_C(j,g)*m_C(j,d)*Coulomb_HO(w,na,ma,nb,mb,ng,mg,nd,md);

                    }
                    if(ma+mb==mg+md && sa==sd && sb==sg){
                         E0hf+=0.5*m_C(i,a)*m_C(i,b)*m_C(j,g)*m_C(j,d)*Coulomb_HO(w,na,ma,nb,mb,nd,md,ng,mg);

                    }






                   }

            }}}

    }
    E0hf+=m_hfenergy(i);
    }








    return E0hf;}


























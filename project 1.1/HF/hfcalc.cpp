#include "hfcalc.h"
#include "Coulomb_Functions.hpp"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


HFcalc::HFcalc(int n, int level, double w)
{
    m_natoms=n;
    m_levels=level;
    m_w=w;



    m_numofpossstates=m_levels*(m_levels+1); //number of possible states


    m_obelements.zeros(m_numofpossstates);  // unperturbed energies

    m_hfenergy.zeros(m_numofpossstates);  //Hattree-Fock energy


    for(int i=0;i<m_levels; i++)         //generating unperturbed energies
    {       for(int j=i*(i+1);j<(i+1)*(i+2);j++){

            m_obelements(j)=(i+1.)*m_w;}


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












void HFcalc::compDensity(){
    m_rho.zeros();
    for(int g=0;g<m_numofpossstates;g++){
        for(int  d=0;d<m_numofpossstates;d++){
            double sum=0;
            for(int i=0; i<m_natoms;i++){
                sum+=m_C(g,i)*m_C(d,i);
            }
            m_rho(g,d)=sum;
        }
    }
}



void HFcalc::compFock(){
    double w=m_w;
    m_HF.zeros();
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

                    if( (ma+mg==mb+md) && (sa==sb) && (sg==sd) && (sa+sg==sd+sb) ){
                        sum +=m_rho(g,d)*Coulomb_HO(w,na,ma,ng,mg,nb,mb,nd,md);
                    }
                    if( (ma+mg==mb+md) && (sa==sd) && (sg==sb) && (sa+sg==sd+sb)  ){
                        sum -=m_rho(g,d)*Coulomb_HO(w,na,ma,ng,mg,nd,md,nb,mb);

                    }

                }

            }


            if(a==b){
                sum+=m_obelements(a);
            }

            m_HF(a,b) += sum;


        }


    }

}




void HFcalc::diagonalize(){
    eig_sym(m_hfenergy,m_C,m_HF);



}

double HFcalc::compHFenergy(int n){
    double lambda= pow(10,-8);
    int counter=1;

    vec temp=zeros<vec>(m_numofpossstates);

    compDensity();
    compFock();
    diagonalize();



    while( abs(m_hfenergy-temp).max() >lambda ){  //loop until it reaches small enough energy difference or number of steps && counter<n


        compDensity();
        compFock();
        temp  = m_hfenergy;
        diagonalize();
        counter++;





    };

    cout << counter << endl;



    double E0hf=0;
    double w=m_w;
    //m_hfenergy.print();
    for(int i=0; i<m_natoms;i++){
        E0hf+=m_hfenergy(i);
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


                           if( (ma+mg==mb+md) && (sa==sb) && (sg==sd) && (sa+sg==sd+sb) ){
                                E0hf-=0.5*m_C(a,i)*m_C(b,i)*m_C(g,j)*m_C(d,j)*Coulomb_HO(w,na,ma,ng,mg,nb,mb,nd,md);
                            }
                            if( (ma+mg==mb+md) && (sa==sd) && (sg==sb) && (sa+sg==sd+sb)  ){
                                E0hf+=0.5*m_C(a,i)*m_C(b,i)*m_C(g,j)*m_C(d,j)*Coulomb_HO(w,na,ma,ng,mg,nd,md,nb,mb);

                            }


                        }

                    }}}

        }

    }

    return E0hf;
}


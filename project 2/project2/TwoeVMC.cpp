#include "TwoeVMC.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
#include "mpi.h"
using namespace arma;
using namespace std;








TwoeVMC::TwoeVMC(int n, int cycles, double step, double w)
{
    m_nelectrons=n;
    m_cycles=cycles;
    m_step=step;
    m_w=w;
    //m_rold.zeros(m_nelectrons,2);
    //m_rnew.zeros(m_nelectrons,2);
    m_dt=0.01;
    m_a=1;
    m_localEn=zeros<vec>(2);
    m_varpar={1,0.4}; //0=alpha, 1 =beta




}

double TwoeVMC::wavefunction(mat &r){
    //first index particle, second dimension
    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        sum+=pos2(r,i);

    }





    double phi= exp(-0.5*m_w*m_varpar(0)*sum )*exp(m_a*relDis(r,0,1)/(1+m_varpar(1)*relDis(r,0,1)));

    return phi;


}

vec TwoeVMC::Quantumforce(int i,mat &r){
    vec F=zeros<vec>(2);




    for(int k=0;k<2;k++){
        double spart=-1.0*r(i,k)*m_varpar(0);

        double jpart=0;



        double r_12=relDis(r,0, 1);

        jpart =m_a*(r(0,k)-r(1,k))/(r_12*pow(1+m_varpar(1)*r_12,2));



        if(i==0){
            F(k)=2.0*(spart+jpart);
        }
        else{
            F(k)=2.0*(spart-jpart);

        }
    }

     return F;
}

double TwoeVMC::localEnergyana(mat &r){


    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        sum+=pos2(r,i);
    }

    //double abs=0;//relDis(r,0,1)

    /*for(int i=0; i<m_nelectrons-1;i++){
        for(int j=i+1; j<m_nelectrons;j++){
            abs+=relDis(r,i,j);
            }
        }
    */

   /* double rabs=sqrt((r(0,0) - r(1,0))*(r(0,0) - r(1,0)) + (r(0,1) - r(1,1))*(r(0,1) - r(1,1)));
    double overr=1.0/rabs;*/
    //cout<<((r(0,0) - r(1,0))*(r(0,0) - r(1,0)) + (r(0,1) - r(1,1))*(r(0,1) - r(1,1)));

    //careful for nele<2








    //double Elocal=0.5*m_w*m_w*(1.0-m_varpar(0)*m_varpar(0))*sum+2.0*m_varpar(0)*m_w+1.0/relDis(r,0,1);
    //overr;
    double r12=relDis(r,0,1);
    double den=1/(1+m_varpar(1)*r12);
    //double l1=0.5*m_w*m_w*(pos2(r,0)+pos2(r,1))*(1-m_varpar(0)*m_varpar(0))+2*m_varpar(0)*m_w+1/r12;
    double l2=-0.5*(m_w*m_w*m_varpar(0)*m_varpar(0)*(pos2(r,0)+pos2(r,1))-4*m_varpar(0)*m_w-2*m_a*m_varpar(0)*m_w*r12*den*den+2*m_a*den*den*(1/r12+m_a*den*den-2*m_varpar(1)*den))+0.5*m_w*m_w*(pos2(r,0)+pos2(r,1))+1/r12;
    double Elocal=l2;
    return Elocal;

}

double TwoeVMC::localEnergyder(mat &r){
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
    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        sum+=pos2(r,i);
    }

    double Elocal=-0.5*deriv/(wfc*h*h)+0.5*m_w*m_w*sum+1/relDis(r,0,1);

    return Elocal;

}

void TwoeVMC::test(){
    mat r=zeros<mat>(m_nelectrons,2);
    double Elder=localEnergyder(r);
    double Elana=localEnergyana(r);
    cout<<Elder<<endl;
    cout<<Elana<<endl;



}

void TwoeVMC::MCH(int c, char **v){
    int NumberProcesses, MyRank;
    //cout<<c;

    MPI_Init (&c, &v);
    MPI_Comm_size (MPI_COMM_WORLD, &NumberProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &MyRank);
    ofstream outFile;
    ofstream position;




    string filename = "mc"; // first command line argument after name of program
    string fileout = filename;
    string argument = to_string(MyRank);
    fileout.append(argument);
    fileout.append(".txt");
    outFile.open(fileout, ios::out);


    string filename2 = "pos"; // first command line argument after name of program
    string fileout2 = filename2;
    string argument2 = to_string(m_w);
    fileout2.append(argument2);
    fileout2.append(".txt");







    if ( MyRank == 0){position.open(fileout2,ios::out);}

    double TotalEnergy=0;
    double TotalEnergy2=0;
    double processEnergy=0;
    double processEnergy2=0;


    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <double>norm(0,1);
    uniform_real_distribution <double>uni(0,1);



    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;


    /*double sum=0;
    double sumsquared=0;*/
    int accept=0;



    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            rold(i,j)=norm(gen)*sqrt(m_dt);
        }
    }
    rnew=rold;
    wfold=wavefunction(rold);





    for(int n=0;n<m_cycles;n++){

        for(int i=0;i<m_nelectrons;i++){
            vec Qforceold=Quantumforce(i,rold);




            for(int j=0;j<2;j++){
                rnew(i,j)=rold(i,j)+0.5*Qforceold(j)*m_dt+norm(gen)*sqrt(m_dt);

          }





            wfnew=wavefunction(rnew);
            mat Qforcenew=Quantumforce(i,rnew);
            double G=0;
            double gexp = 0;
            for(int j=0; j<2;j++){
                double term1 = - (rold(i,j) - rnew(i,j)  -  0.5*m_dt*Qforcenew(j))*
                                        (rold(i,j) - rnew(i,j)  -  0.5*m_dt*Qforcenew(j));
                double term2 =   (-rold(i,j) + rnew(i,j) - 0.5*m_dt*Qforceold(j))*
                                        (-rold(i,j) + rnew(i,j) - 0.5*m_dt*Qforceold(j));

                gexp += term1+term2;




            }
            G=exp(gexp/(2.0*m_dt));//gexp/(2.0*m_dt)
          //cout<<(wfnew*wfnew)/(wfold*wfold)<<endl;
            if(uni(gen)<=G*(wfnew*wfnew)/(wfold*wfold)){
                accept++;
                for(int j=0;j<2;j++){
                    rold(i,j)=rnew(i,j);

                }
                wfold=wfnew;

            }
            else{
                for(int j=0; j<2;j++){
                    rnew(i,j)=rold(i,j);
                }
            }
            double l1=localEnergyana(rold);
            double temp=l1;
            processEnergy+=temp;
            processEnergy2+=temp*temp;
            //outFile.write( (char*)&temp, sizeof(double));
            outFile << l1 << endl;


        }



        if ( MyRank == 0){rnew.print(position);}





    }



    MPI_Reduce(&processEnergy, &TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&processEnergy2, &TotalEnergy2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    if ( MyRank == 0) {
        double Energy = TotalEnergy/( (double)NumberProcesses*m_cycles*m_nelectrons);
        double Energy2=TotalEnergy2/( (double)NumberProcesses*m_cycles*m_nelectrons);
        double variance=(Energy2-Energy*Energy)/( (double)NumberProcesses*m_cycles*m_nelectrons);


        cout << Energy << endl;
        cout << variance<<endl;

        cout<<(double)accept/(m_cycles*m_nelectrons)*100<<endl;
        cout<< m_varpar(0)<<endl;
        cout<<m_varpar(1)<<endl;

    }
    /*double Energy = processEnergy/( (double)m_cycles);
    double Variance = processEnergy2/( (double)m_cycles)-Energy*Energy;
    cout << Energy << endl;
    cout<< Variance <<  endl;*/







    fflush(stdout);
    outFile.close();

    MPI_Finalize ();




}

void TwoeVMC::MonteCarlo(){


    random_device rnd;
    mt19937 gen(rnd());
    uniform_real_distribution <double>dis(-1,1);
    uniform_real_distribution <double>dis1(0,1);







    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);

    double wfold=0.0;
    double wfnew=0.0;
    double sum=0.0;
    double sumsquared=0.0;
    int accept=0;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            rold(i,j)=m_step*dis(gen);
        }
    }
    rnew=rold;
    wfold=wavefunction(rold);



    for(int n=0;n<m_cycles;n++){


        for(int i=0;i<m_nelectrons;i++){
            for(int j=0;j<2;j++){
                rnew(i,j)=rold(i,j)+m_step*dis(gen);

          }



         wfnew=wavefunction(rnew);
          //cout<<(wfnew*wfnew)/(wfold*wfold)<<endl;
            if(dis1(gen)<=(wfnew*wfnew)/(wfold*wfold)){                       //exp(m_varpar(0)*m_w*(pos2(rold,0)+pos2(rold,1)-pos2(rnew,0)-pos2(rnew,1)))
                accept++;
                for(int j=0;j<2;j++){
                    rold(i,j)=rnew(i,j);
                    wfold=wfnew;
                }

            }
            else{
                for(int j=0; j<2;j++){
                    rnew(i,j)=rold(i,j);
                }
            }


        }
        double temp=localEnergyder(rnew);
        //cout<<temp<<endl;
        sum+=temp;
        sumsquared+=temp*temp;

    }

    double energy=(double)sum/m_cycles;
    double energysquared=(double)sumsquared/m_cycles;
    double variance=(double)(energysquared - energy*energy)/m_cycles;




    cout<< "Energy"<< energy<<endl;
    cout<< "Variance"<< variance<<endl;
    cout<<(double)accept/(m_nelectrons*m_cycles)*100<<endl;



}

void TwoeVMC::findoptParameter(){
    int updates=10000;
    int cycles=1000000;
    vec dE=zeros<vec>(2);
    double tolerance = 1.0e-14;
    vec parold={1,0.5};
    vec parnew=zeros<vec>(2);
    double diff=0;
    double step=1;

    vec expVal=zeros<vec>(5); //0=en, 1,2=psi







    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <double>norm(0,1);
    uniform_real_distribution <double>uni(0,1);




    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;




    m_varpar=parold;


    for( int m=0; m<updates; m++ ){
        vec sumPsi=zeros<vec>(5); //0=en, 1,2=psi
        for(int i=0;i<m_nelectrons;i++){
            for(int j=0;j<2;j++){
                rold(i,j)=norm(gen)*sqrt(m_dt);
            }
        }
        rnew=rold;



        for(int n=0;n<cycles;n++){
            wfold=wavefunction(rold);


            for(int i=0;i<m_nelectrons;i++){
                vec Qforceold=Quantumforce(i,rold);
                for(int j=0;j<2;j++){
                    rnew(i,j)=rold(i,j)+0.5*Qforceold(j)*m_dt+norm(gen)*sqrt(m_dt);

              }




                wfnew=wavefunction(rnew);
                vec Qforcenew=Quantumforce(i,rnew);
                double G=0;
                for(int j=0; j<2;j++){
                  G += 0.5*m_dt*(Qforceold(j)*Qforceold(j)-Qforcenew(j)*Qforcenew(j))+0.5*(rold(i,j)-rnew(i,j))*(Qforceold(j)+Qforcenew(j));
                }
                G=exp(G);

                if(uni(gen)<=G*(wfnew*wfnew)/(wfold*wfold)){
                    for(int j=0;j<2;j++){
                        rold(i,j)=rnew(i,j);
                        wfold=wfnew;
                    }
                    Qforceold=Qforcenew;

                }
                else{
                    for(int j=0; j<2;j++){
                        rnew(i,j)=rold(i,j);
                    }
                }


            }
            double r12=relDis(rnew, 0,1);
            double den=1/(1+m_varpar(1)*r12);
            double temp=localEnergyana(rnew);//0.5*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))*(1-m_varpar(0)*m_varpar(0))+2*m_varpar(0)*m_w+1/r12;




            double tempPsia=-0.5*m_w*(pos2(rnew,0)+pos2(rnew,1));
            double tempPsib=-m_a*r12*r12*den*den;
            sumPsi(0)+=temp;
            sumPsi(2)+=tempPsib;
            sumPsi(1)+=tempPsia;
            sumPsi(3)+=tempPsia*temp;
            sumPsi(4)+=tempPsib*temp;


        }










        expVal=sumPsi/(cycles);




        dE(0)=2*(expVal(3)-expVal(1)*expVal(0));
        dE(1)=2*(expVal(4)-expVal(2)*expVal(0));
        parnew=parold-step*dE;


        m_varpar=parnew;
        parold=parnew;



        diff= sqrt(dot(parnew-parold,parnew-parold));
        if(diff<tolerance){ break;}
        else{parold=parnew;}









}



    //cout<<m_varpar(0)<<endl;
    //cout<<m_varpar(1);




}
















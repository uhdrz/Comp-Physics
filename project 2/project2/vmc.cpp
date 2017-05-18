#include "vmc.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
//#include "mpi.h"
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
    m_dt=0.001;
    m_a=1;
    m_localEn=zeros<vec>(2);
    m_varpar={1,0.4}; //0=alpha, 1 =beta




}

double VMC::hermite(int n, double x){
    if (n==0)
    {
        return 1;}
    else if (n==1)
    {
        return 2*x;}
    else
    {   double temp1=1;
        double temp2=2*x;
        double temp=0;
        for(int i=2; i<n+1; i++)
        {   temp=2*x*temp2-2*(i-1)*temp1;
            temp1=temp2;
            temp2=temp;
            //cout<<temp;
        }
       // cout<<temp;
        return temp;


    }
}





double VMC::spwf(double x, double y, int nx, int ny){

    double phi=hermite(nx,sqrt(m_w)*x)* hermite(ny,sqrt(m_w)*y)*exp(-0.5*m_varpar(0)*m_w*(pow(x,2)+pow(y,2)));

    return phi;




}


double VMC::SlaterDet(mat &r){
    int n_el= m_nelectrons*0.5; //Spin separatly
    mat Slaterup = zeros(n_el,n_el);
    mat Slaterdown = zeros(n_el,n_el);
    for(int i=0;i<n_el;i++){
        for( int j=0;j<n_el;j++){
            vec n=postonum(j);
            Slaterup(i,j)=spwf(r(i,0),r(i,1),n(0),n(1));
        }

     }
    for(int i=n_el;i<m_nelectrons;i++){
        for( int j=n_el;j<m_nelectrons;j++){
            vec n=postonum(j);
            Slaterdown(i,j)=spwf(r(i,0),r(i,1),n(0),n(1));
        }

     }

  double SlaterDet=det(Slaterup)*det(Slaterdown);

  return SlaterDet;

}




double VMC::wavefunction1(mat &r){
    double Slaterdet=SlaterDet(r);
    double jastrow=0;
    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<i;j++){
            double rij=relDis(r,i,j);
            jastrow*=exp(m_a*rij/(1+m_varpar(1)*rij));

        }
    }

    return Slaterdet*jastrow;

}

double VMC::localEnergy(mat &r){
    double Pot=0;
    for(int i=0;i<m_nelectrons;i++){
        Pot+=0.5*m_w*m_w*pos2(r,i);
        for(int j=0;j<i;j++){
           Pot+=1/relDis(r,i,j);
        }
    }
    double Kin=0;
    for(int i=0;i<m_nelectrons;i++){
        Kin+=LapJastrow(r);
        if(i<2){
            Kin+=LapSP0(r);
            Kin+=2*dot(GradJastrow(r),GradSP0(r));
        }
        if(i<6){
            Kin+=LapSP1(r);
            Kin+=2*dot(GradJastrow(r),GradSP1(r));
        }
        if(i<12){
            Kin+=LapSP2(r);
            Kin+=2*dot(GradJastrow(r),GradSP2(r));
        }

     }


    double locEnergy=Kin+Pot;
    return locEnergy;


}






double VMC::wavefunction(mat &r){
    //first index particle, second dimension
    double sum=0;
    for(int i=0;i<m_nelectrons;i++){
        sum+=pos2(r,i);

    }





    double phi= exp(-0.5*m_w*m_varpar(0)*sum )*exp(m_a*relDis(r,0,1)/(1+m_varpar(1)*relDis(r,0,1)));

    return phi;


}





mat VMC::Quantumforce(mat &r){
    mat F=zeros(m_nelectrons,2);



    for(int i=0;i<m_nelectrons;i++){
        for(int k=0;k<2;k++){
            double spart=-1.0*r(i,k)*m_varpar(0);

            double jpart=0;

            for(int j=0;j<m_nelectrons;j++){
                if(i!=j){

                    double r_ij=relDis(r,i, j);

                    jpart +=m_a*(r(i,k)-r(j,k))/(r_ij*pow(1+m_varpar(1)*r_ij,2));

               }
              }

            F(i,k)=2.0*(spart+jpart);


        }
    }





     return F;
}




double VMC::localEnergyana(mat &r){


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








    //double Elocal=0.5*m_w*m_w*(1.0-m_varpar(0)*m_varpar(0))*sum+2.0*m_varpar(0)*m_w+1.0/relDis(r,0,1);//overr;
    double r12=relDis(r,0,1);
    double den=1/(1+m_varpar(1)*r12);
    //double l1=0.5*m_w*m_w*(pos2(r,0)+pos2(r,1))*(1-m_varpar(0)*m_varpar(0))+2*m_varpar(0)*m_w+1/r12;
    double l2=-0.5*(m_w*m_w*m_varpar(0)*m_varpar(0)*(pos2(r,0)+pos2(r,1))-4*m_varpar(0)*m_w-2*m_a*m_varpar(0)*m_w*r12*den*den+2*m_a*den*den*(1/r12+m_a*den*den-2*m_varpar(1)*den))+0.5*m_w*m_w*(pos2(r,0)+pos2(r,1))+1/r12;
    double Elocal=l2;
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


    /*MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &NumberProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &MyRank);*/

    // MPI_Bcast (&m_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);






    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <double>norm(0,1);
    uniform_real_distribution <double>uni(0,1);
    ofstream outFile;
    outFile.open("localEnergies.dat", ios::out | ios::binary);
    if (! outFile.is_open()) {
         cout << "Could not open file." << endl;}




    mat rold=zeros<mat>(m_nelectrons,2);
    mat rnew=zeros<mat>(m_nelectrons,2);
    double wfold=0;
    double wfnew=0;
    vec sum=zeros<vec>(2);
    /*double sum=0;
    double sumsquared=0;*/
    int accept=0;



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
              G += 0.5*m_dt*(Qforceold(i,j)*Qforceold(i,j)-Qforcenew(i,j)*Qforcenew(i,j))+0.5*(rold(i,j)-rnew(i,j))*(Qforceold(i,j)+Qforcenew(i,j));
            }
            G=exp(G);
          //cout<<(wfnew*wfnew)/(wfold*wfold)<<endl;
            if(uni(gen)<=G*(wfnew*wfnew)/(wfold*wfold)){                       //exp(m_varpar(0)*m_w*(pos2(rold,0)+pos2(rold,1)-pos2(rnew,0)-pos2(rnew,1)))
                accept++;
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
        //double r12=relDis(rnew, 0,1);
        //double den=1/(1+m_varpar(1)*r12);
        double l1=localEnergyana(rnew);//0.5*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))*(1-m_varpar(0)*m_varpar(0))+2*m_varpar(0)*m_w+1/r12;
       // double l2=2*m_varpar(0)*m_varpar(0)*m_w*m_w*(pos2(rnew,0)+pos2(rnew,1))-4*m_w*m_varpar(0)-2*m_a*m_varpar(0)*m_w*r12*den*den+2*m_a*den*den*(m_a*den*den+1/r12-2*m_varpar(1)*den);
        double temp=l1;


        outFile.write( (char*)&temp, sizeof(double));
        sum(0)+=temp;
        sum(1)+=temp*temp;

    }
    //MPI_Reduce(&LocalProcessEnergy, &TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    m_localEn=sum/m_cycles;

    double variance=(m_localEn(1)-m_localEn(0)*m_localEn(0))/m_cycles;
    m_localEn.print();
    cout<<variance<<endl;
    cout<<(double)accept/(m_cycles*m_nelectrons)*100;


    fflush(stdout);
    outFile.close();

    // MPI_Finalize ();



}





void VMC::MonteCarlo(){


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



    for(int n=0;n<m_cycles;n++){
        wfold=wavefunction(rold);

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
        double temp=localEnergyana(rnew);
        //cout<<temp<<endl;
        sum+=temp;
        sumsquared+=temp*temp;

    }

    double energy=(double)sum/m_cycles;
    double energysquared=(double)sumsquared/m_cycles;
    double variance=(double)(energysquared - energy*energy)/m_cycles;




    cout<< "Energy"<< energy<<endl;
    cout<< "Variance"<< variance<<endl;
    cout<<(double)accept/(m_nelectrons*m_cycles)*100;



}






void VMC::findoptParameter(){
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
            mat Qforceold=Quantumforce(rold);

            for(int i=0;i<m_nelectrons;i++){
                for(int j=0;j<2;j++){
                    rnew(i,j)=rold(i,j)+0.5*Qforceold(i,j)*m_dt+norm(gen)*sqrt(m_dt);

              }




                wfnew=wavefunction(rnew);
                mat Qforcenew=Quantumforce(rnew);
                double G=0;
                for(int j=0; j<2;j++){
                  G += 0.5*m_dt*(Qforceold(i,j)*Qforceold(i,j)-Qforcenew(i,j)*Qforcenew(i,j))+0.5*(rold(i,j)-rnew(i,j))*(Qforceold(i,j)+Qforcenew(i,j));
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



    cout<<m_varpar(0)<<endl;
    cout<<m_varpar(1);




}
















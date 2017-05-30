#include "vmc.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include "mpi.h"
using namespace arma;
using namespace std;

VMC::VMC(int n, int cycles, double step, double w)
{
    m_nelectrons=n;
    m_cycles=cycles;
    m_step=step;
    m_w=w;
    m_dt=0.001;
    m_a={1,1.0/3};  //0 for antipar
    m_localEn=zeros<vec>(2); //?????
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

vec VMC::postonum(int i){
    //position in slater matrix not number of electron, because of spin split

    int temp1=0;
    int temp2=0;
    vec n= zeros<vec>(2);
    int m;
    for(m=0; temp1<=i; m++)
    { temp2=temp1;
      temp1+= (m+1);
     }

    n(1)=i-temp2;
    n(0)=m-1-n(1);


    return n;

}

double VMC::spwf(double x, double y, int nx, int ny){

    double phi=hermite(nx,sqrt(m_w)*x)* hermite(ny,sqrt(m_w)*y)*exp(-0.5*m_varpar(0)*m_w*(pow(x,2)+pow(y,2)));

    return phi;




}

double VMC::SlaterDet(mat &r){
    int n_states= m_nelectrons/2; //Spin separatly

    mat Slaterup = zeros(n_states,n_states);
    mat Slaterdown = zeros(n_states,n_states);
    for(int i=0;i<n_states;i++){
        for( int j=0;j<n_states;j++){
            vec n=postonum(j);
            //n.print();
            Slaterup(i,j)=spwf(r(i,0),r(i,1),n(0),n(1));

        }

     }

    for(int i=0;i<n_states;i++){
        for( int j=0;j<n_states;j++){
            int id=i+n_states;
            vec n=postonum(j);
            Slaterdown(i,j)=spwf(r(id,0),r(id,1),n(0),n(1));

        }

     }

  double SlaterDet=det(Slaterup)*det(Slaterdown);




  return SlaterDet;

}

double VMC::wavefunction(mat &r){
    double Slaterdet=SlaterDet(r);
    double jastrow=1;
    int n_states=m_nelectrons/2;
    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<i;j++){
            double rij=relDis(r,i,j);
            if((i < n_states && j< n_states)||(i >= n_states && j>= n_states)){
                jastrow*=exp(m_a(1)*rij/(1+m_varpar(1)*rij));}
            else{
                jastrow*=exp(m_a(0)*rij/(1+m_varpar(1)*rij));
            }

        }
    }


    return Slaterdet*jastrow;

}

vec VMC::GradJastrow(int k, mat &r){

    vec deriv=zeros<vec>(2);
    for(int i=0;i<k;i++){
        for(int j=0;j<2;j++){
            double rik=relDis(r,i,k);
            deriv(j)+=(r(k,j)-r(i,j))/(rik*pow(1.0+m_varpar(1)*rik,2));

        }
    }
    for(int i=k+1;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            double rki=relDis(r,k,i);
            deriv(j)-=(r(i,j)-r(k,j))/(rki*pow(1.0+m_varpar(1)*rki,2));

        }
    }

    return deriv;


}

double VMC::LapJastrow(int k, mat &r){
    double lap=0;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<m_nelectrons;j++){
            if( (i!=k) && (j!=k) ){
                for( int d=0; d<2;d++){
                    double rki=relDis(r,k,i);
                    double rkj=relDis(r,k,j);
                    lap+=(r(k,d)-r(i,d))*(r(k,d)-r(j,d))/(rki*pow(1.0+m_varpar(1)*rki,2)*rkj*pow(1.0+m_varpar(1)*rkj,2));

                }
            }
        }
        if(i!=k){
            double rki=relDis(r,k,i);
            lap+=2/((rki*pow(1.0+m_varpar(1)*rki,2)))-2*m_varpar(1)/pow(1.0+m_varpar(1)*rki,3);
        }
    }
    return lap;

}

vec VMC::GradSP(int i, mat &r){


    vec Grad=zeros<vec>(2);
    int n_states= m_nelectrons/2; //Spin separatly

    mat Slaterup = zeros(n_states,n_states);
    mat Slaterdown = zeros(n_states,n_states);
    for(int l=0;l<n_states;l++){
        for( int j=0;j<n_states;j++){
            vec n=postonum(j);
            Slaterup(l,j)=spwf(r(l,0),r(l,1),n(0),n(1));

        }

     }
    for(int l=0;l<n_states;l++){
        for( int j=0;j<n_states;j++){
            int ld=l+n_states;
            vec n=postonum(j);
            Slaterdown(l,j)=spwf(r(ld,0),r(ld,1),n(0),n(1));

        }

     }
    mat SlaterupInv=inv(Slaterup);
    mat SlaterdownInv=inv(Slaterdown);



    for(int j=0;j<n_states;j++){
        vec n=postonum(j);
        if(i<n_states){
            Grad+=DerivSP(i,n(0),n(1),r)*spwf(r(i,0),r(i,1),n(0),n(1))*SlaterupInv(j,i);
         }
        else{
            int temp=i-n_states;
            Grad+=DerivSP(i,n(0),n(1),r)*spwf(r(i,0),r(i,1),n(0),n(1))*SlaterdownInv(j,temp);
         }

    }

    return Grad;





}

vec VMC::DerivSP(int i, int nx, int ny, mat &r){
    vec deriv=zeros<vec>(2);
    if((nx==0) && (ny==0)){
        deriv(0)=-m_varpar(0)*m_w*r(i,0);
        deriv(1)=-m_varpar(0)*m_w*r(i,1);
    }
    else if((nx==1) && (ny==0)){
        deriv(0)=sqrt(m_w)*(1-m_varpar(0)*m_w*r(i,0)*r(i,0));
        deriv(1)=-sqrt(m_w)*m_varpar(0)*m_w*r(i,1)*r(i,0);
    }
    else if((nx==0) && (ny==1)){
        deriv(1)=sqrt(m_w)*(1-m_varpar(0)*m_w*r(i,1)*r(i,1));
        deriv(0)=-sqrt(m_w)*m_varpar(0)*m_w*r(i,1)*r(i,0);
    }
    else if((nx==1) && (ny==1)){
        deriv(0)=m_w*(r(i,1)-m_varpar(0)*m_w*r(i,0)*r(i,0)*r(i,1));
        deriv(1)=m_w*(r(i,0)-m_varpar(0)*m_w*r(i,1)*r(i,1)*r(i,0));
    }
    else if((nx==2) && (ny==0)){
        deriv(0)=m_w*r(i,0)*(m_varpar(0)*(m_w*r(i,0)*r(i,0)-1)-2);
        deriv(1)=-m_varpar(0)*m_w*r(i,1)*(m_w*r(i,0)*r(i,0)-1);
    }
    else if((nx==0) && (ny==2)){
        deriv(1)=m_w*r(i,1)*(m_varpar(0)*(m_w*r(i,1)*r(i,1)-1)-2);
        deriv(0)=-m_varpar(0)*m_w*r(i,0)*(m_w*r(i,1)*r(i,1)-1);
    }

    return deriv;

}

double VMC::Deriv2SP(int i, int nx, int ny, mat &r){
    double lap=0;
    if((nx==0) && (ny==0)){
        lap=m_varpar(0)*m_varpar(0)*m_w*m_w*pos2(r,i)-2*m_varpar(0)*m_w;
    }
    else if((nx==1) && (ny==0)){
        lap=sqrt(m_w)*(-4*m_w*r(i,0)+pow(m_varpar(0),2)*pow(m_w,2)*r(i,0)*pos2(r,i));
    }
    else if((nx==0) && (ny==1)){
        lap=sqrt(m_w)*(-4*m_w*r(i,1)+pow(m_varpar(0),2)*pow(m_w,2)*r(i,1)*pos2(r,i));
    }
    else if((nx==1) && (ny==1)){
        lap=m_w*(-6*m_varpar(0)*m_w*r(i,1)*r(i,0)+pow(m_varpar(0),2)*pow(m_w,2)*r(i,1)*r(i,0)*pos2(r,i));
    }
    else if((nx==2) && (ny==0)){
        lap=(m_varpar(0)*pow(m_w,2)*pow(r(i,0),2)-m_varpar(0)*m_w)*(m_varpar(0)*m_w*pow(r(i,1),2)-1)+m_w*(pow(m_varpar(0),2)*pow(r(i,0),2)*m_w*(m_w*pow(r(i,0),2)-1)-5*m_varpar(0)*m_w*pow(r(i,0),2)+2+m_varpar(0));
    }
    else if((nx==0) && (ny==2)){
        lap=(m_varpar(0)*pow(m_w,2)*pow(r(i,1),2)-m_varpar(0)*m_w)*(m_varpar(0)*m_w*pow(r(i,0),2)-1)+m_w*(pow(m_varpar(0),2)*pow(r(i,1),2)*m_w*(m_w*pow(r(i,1),2)-1)-5*m_varpar(0)*m_w*pow(r(i,1),2)+2+m_varpar(0));
    }
    return lap;

}

double VMC::LapSP(int i ,mat &r){

    double laplace=0;
    int n_states= m_nelectrons/2; //Spin separatly

    mat Slaterup = zeros(n_states,n_states);
    mat Slaterdown = zeros(n_states,n_states);
    for(int l=0;l<n_states;l++){
        for( int j=0;j<n_states;j++){
            vec n=postonum(j);
            Slaterup(l,j)=spwf(r(l,0),r(l,1),n(0),n(1));

        }

     }
    for(int l=0;l<n_states;l++){
        for( int j=0;j<n_states;j++){
            int ld=l+n_states;
            vec n=postonum(j);
            Slaterdown(l,j)=spwf(r(ld,0),r(ld,1),n(0),n(1));

        }

     }
    mat SlaterupInv=inv(Slaterup);
    mat SlaterdownInv=inv(Slaterdown);



    for(int j=0;j<n_states;j++){
        vec n=postonum(j);
        if(i<n_states){
            laplace += Deriv2SP(i,n(0),n(1),r)*spwf(r(i,0),r(i,1),n(0),n(1))*SlaterupInv(j,i);

        }
        else{
            int temp=i-n_states;
            laplace += Deriv2SP(i,n(0),n(1),r)*spwf(r(i,0),r(i,1),n(0),n(1))*SlaterdownInv(j,temp);

        }

    }

    return laplace;




}

double VMC::localEnergy(mat &r){

    double Pot=0;
    for(int i=0;i<m_nelectrons;i++){
        Pot+=0.5*m_w*m_w*pos2(r,i);
        /*for(int j=0;j<i;j++){
           Pot+=1/relDis(r,i,j);
        }*/
    }
    double Kin=0;
    for(int i=0;i<m_nelectrons;i++){
        Kin += -0.5*LapSP(i,r); //LapJastrow(i,r)+2*dot(GradSP(i,r),GradJastrow(i,r))+
    }


    double locEnergy=Kin+Pot;
    return locEnergy;


}

mat VMC::Quantumforce(mat &r){ //fix m_a
    mat F=zeros(m_nelectrons,2);
    int n_el=m_nelectrons/2;
    for(int i=0;i<m_nelectrons;i++){
        for(int k=0;k<2;k++){
            double spart=-1.0*r(i,k)*m_varpar(0);

            double jpart=0;

            for(int j=0;j<m_nelectrons;j++){
                if(i!=j){

                    double r_ij=relDis(r,i, j);
                    if((i < n_el && j< n_el)||(i >= n_el && j>= n_el)){
                        jpart +=m_a(1)*(r(i,k)-r(j,k))/(r_ij*pow(1+m_varpar(1)*r_ij,2));}
                    else{
                        jpart +=m_a(0)*(r(i,k)-r(j,k))/(r_ij*pow(1+m_varpar(1)*r_ij,2));
                    }

               }
              }

            F(i,k)=2.0*(spart+jpart);


        }
    }
     return F;
}

void VMC::MCH(){
    //int NumberProcesses, MyRank;
    //cout<<c;

    /*MPI_Init (&c, &v);
    MPI_Comm_size (MPI_COMM_WORLD, &NumberProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &MyRank);
    ofstream outFile;




    string filename = "mc"; // first command line argument after name of program
    string fileout = filename;
    string argument = to_string(MyRank);
    fileout.append(argument);
    fileout.append(".txt");
    outFile.open(fileout, ios::out);*/


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

            double l1=localEnergy(rnew);
                    double temp=l1;


                    //outFile.write( (char*)&temp, sizeof(double));
                    //outFile << l1 << endl;
                    processEnergy+=temp;
                    processEnergy2+=temp*temp;


        }


    }



    //MPI_Reduce(&processEnergy, &TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&processEnergy2, &TotalEnergy2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    /*if ( MyRank == 0) {
        double Energy = TotalEnergy/( (double)NumberProcesses*m_cycles);
        double Variance = TotalEnergy2/( (double)NumberProcesses*m_cycles)-Energy*Energy;
        double StandardDeviation = sqrt(Variance/((double)NumberProcesses*m_cycles)); // over optimistic error

        cout << Energy << endl;
        cout<< Variance <<  endl;
        cout<< StandardDeviation << endl;
        //

    }*/
    double Energy = processEnergy/( (double)m_cycles*m_nelectrons);
    double Variance = processEnergy2/( (double)m_cycles)-Energy*Energy;
    cout << Energy << endl;
    cout<< Variance <<  endl;
    cout<<(double)accept/(m_cycles*m_nelectrons)*100;






    /*fflush(stdout);
    outFile.close();

    MPI_Finalize ();*/




}







/*
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




}*/

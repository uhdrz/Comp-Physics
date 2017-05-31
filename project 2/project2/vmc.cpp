#include "vmc.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <hermitepolynomials.h>
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
    m_a={0.599,1.0/3};  //0 for antipar
    m_localEn=zeros<vec>(2); //?????
    m_varpar={1,1e6}; //0=alpha, 1 =beta
    m_invDown.zeros(m_nelectrons/2,m_nelectrons/2);
    m_invUp.zeros(m_nelectrons/2,m_nelectrons/2);


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
    const double Hnx = HermitePolynomials::evaluate(nx,x,m_w,m_varpar(0));
    const double Hny = HermitePolynomials::evaluate(ny,y,m_w,m_varpar(0));


    double phi=Hnx* Hny*exp(-0.5*m_varpar(0)*m_w*(pow(x,2)+pow(y,2)));

    return phi;




}

mat VMC::Slatermatrixup(mat &r){
    int n_states= m_nelectrons/2; //Spin separatly

    mat Slaterup = zeros(n_states,n_states);
    for(int i=0;i<n_states;i++){
        for( int j=0;j<n_states;j++){
            vec n=postonum(j);
            //n.print();
            Slaterup(i,j)=spwf(r(i,0),r(i,1),n(0),n(1));

        }

     }
    return Slaterup;
}

mat VMC::Slatermatrixdown(mat &r){
    int n_states= m_nelectrons/2; //Spin separatly

    mat Slaterdown = zeros(n_states,n_states);

    for(int i=0;i<n_states;i++){
        for( int j=0;j<n_states;j++){
            int id=i+n_states;
            vec n=postonum(j);
            Slaterdown(i,j)=spwf(r(id,0),r(id,1),n(0),n(1));

        }

   }
    return Slaterdown;

}

mat VMC::SlaterUpInv(mat &Slaterup){

    mat SlaterInv=Slaterup.i();

    return SlaterInv;

}

mat VMC::SlaterDownInv(mat &Slaterdown){

    mat SlaterInv=Slaterdown.i();

    return SlaterInv;
}

double VMC::wavefunction(mat &r, mat &Slaterup, mat &Slaterdown){

    double Slaterdet = det(Slaterup)*det(Slaterdown);
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
            deriv(j)+=(r(k,j)-r(i,j))*derivJastrow(i,k,rik)   /rik;

        }
    }
    for(int i=k+1;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            double rki=relDis(r,k,i);
            deriv(j)-=(r(i,j)-r(k,j))*derivJastrow(k,i,rki)  /rki;
        }
    }

    return deriv;


}

double VMC::derivJastrow(int i, int j, double rij){
    int n_states=m_nelectrons/2;
    double deriv=0;

    if((i < n_states && j< n_states)||(i >= n_states && j>= n_states)){
        deriv=m_a(1)/(rij*pow(1.0+m_varpar(1)*rij,2));
    }
    else{
        deriv=m_a(0)/(rij*pow(1.0+m_varpar(1)*rij,2));

    }

    return deriv;
}

 double VMC::deriv2Jastrow(int i, int j, double rij){
     int n_states=m_nelectrons/2;
     double deriv2=0;

     if((i < n_states && j< n_states)||(i >= n_states && j>= n_states)){
         deriv2=m_a(1)*m_varpar(1)/(rij*pow(1.0+m_varpar(1)*rij,3));
     }
     else{
         deriv2=m_a(0)*m_varpar(1)/(rij*pow(1.0+m_varpar(1)*rij,3));

     }

     return deriv2;



 }

double VMC::LapJastrow(int k, mat &r){
    double lap=0;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<m_nelectrons;j++){
            if( (i!=k) && (j!=k) ){
                for( int d=0; d<2;d++){
                    double rki=relDis(r,k,i);
                    double rkj=relDis(r,k,j);
                    lap+=(r(k,d)-r(i,d))*(r(k,d)-r(j,d))*derivJastrow(k,i,rki)*derivJastrow(k,j,rkj)/(rki*rkj);

                }
            }
        }
        if(i!=k){
            double rki=relDis(r,k,i);
            lap+=deriv2Jastrow(k,i,rki)+2*derivJastrow(k,i,rki)/rki;
        }
    }
    return lap;

}

vec VMC::GradSP(int i, mat &r, mat &InvUp, mat &InvDown){


    vec Grad=zeros<vec>(2);
    int n_states= m_nelectrons/2; //Spin separatly

    ;



    for(int j=0;j<n_states;j++){
        vec n=postonum(j);
        if(i<n_states){
            Grad+=DerivSP(i,n(0),n(1),r)*InvUp(j,i);
         }
        else{
            int temp=i-n_states;
            Grad+=DerivSP(i,n(0),n(1),r)*InvDown(j,temp);
         }

    }

    return Grad;





}

vec VMC::DerivSP(int i,int nx,int ny, mat &r) {

    const double x = r(i,0);
    const double y = r(i,1);
    vec deriv=zeros<vec>(2);



    const double Hnx = HermitePolynomials::evaluate(nx,x,m_w,m_varpar(0));
    const double Hny = HermitePolynomials::evaluate(ny,y,m_w,m_varpar(0));

    const double dHnx = HermitePolynomials::evaluateDerivative(nx,x,m_w,m_varpar(0));
    const double dHny = HermitePolynomials::evaluateDerivative(ny,y,m_w,m_varpar(0));

    double r2 = x*x+y*y;
    deriv(0)=exp(-0.5*m_w*m_varpar(0)*r2) * Hny *( dHnx - Hnx*m_w*m_varpar(0)*x);
    deriv(1)=exp(-0.5*m_w*m_varpar(0)*r2) * Hnx *( dHny - Hny*m_w*m_varpar(0)*y);
    return deriv;



}

double VMC::Deriv2SP(int i, int nx, int ny, mat &r) {

    double x = r(i,0);
    double y = r(i,1);

    const double Hnx = HermitePolynomials::evaluate(nx,x,m_w,m_varpar(0));
    const double Hny = HermitePolynomials::evaluate(ny,y,m_w,m_varpar(0));

    const double dHnx = HermitePolynomials::evaluateDerivative(nx,x,m_w,m_varpar(0));
    const double dHny = HermitePolynomials::evaluateDerivative(ny,y,m_w,m_varpar(0));

    const double ddHnx = HermitePolynomials::evaluateDoubleDerivative(nx,x,m_w,m_varpar(0));
    const double ddHny = HermitePolynomials::evaluateDoubleDerivative(ny,y,m_w,m_varpar(0));

    double r2 = x*x + y*y;

    double m_omegaAlpha = m_w*m_varpar(0);

    return exp(-0.5*m_omegaAlpha*r2) *
               ( - 2*m_omegaAlpha*x*Hny*dHnx
                 - 2*m_omegaAlpha*y*Hnx*dHny
                 + m_omegaAlpha*Hnx*Hny*(m_omegaAlpha*r2-2)
                 + Hny*ddHnx
                 + Hnx*ddHny );

}

double VMC::LapSP(int i ,mat &r, mat &InvUp, mat &InvDown){

    double laplace=0;
    int n_states= m_nelectrons/2; //Spin separatly



    for(int j=0;j<n_states;j++){
        vec n=postonum(j);
        //n.print();
        if(i<n_states){
            laplace += Deriv2SP(i,n(0),n(1),r)*InvUp(j,i);

        }

        else{

            int temp=i-n_states;
            laplace += Deriv2SP(i,n(0),n(1),r)*InvDown(j,temp);

        }

    }

    //exit(1);

    return laplace;




}

double VMC::localEnergy(mat &r, mat &InvUp, mat &InvDown){

    double Pot=0;
    for(int i=0;i<m_nelectrons;i++){
        Pot+=0.5*m_w*m_w*pos2(r,i);
        for(int j=0;j<i;j++){
           Pot+=1.0/relDis(r,i,j);
        }
    }
    double Kin=0;
    for(int i=0;i<m_nelectrons;i++){
        Kin += -0.5*(LapSP(i,r,InvUp,InvDown));//+LapJastrow(i,r)+2*dot(GradSP(i,r,InvUp,InvDown),GradJastrow(i,r)));
    }


    double locEnergy=Kin+Pot;
    return locEnergy;


}

mat VMC::Quantumforce(mat &r, mat &InvUp, mat &InvDown){ //fix m_a

    mat F=zeros(m_nelectrons,2);


    for(int i=0;i<m_nelectrons;i++){
        vec GradSpace=GradSP(i,r,InvUp,InvDown);
        //vec jpart=GradJastrow(i,r);
        for(int k=0;k<2;k++){
          F(i,k) = 2.0*(GradSpace(k));//+jpart(k));
         }

    }
    return F;
}

void VMC::MCH(){
    //int NumberProcesses, MyRank;
    //cout<<c;

    /*MPI_Init (&c, &v);
    MPI_Comm_size (MPI_COMM_WORLD, &NumberProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &MyRank);*/
    ofstream outFile;




    string filename = "VMC"; // first command line argument after name of program
    string fileout = filename;
   // string argument = to_string(MyRank);
    //fileout.append(argument);
    fileout.append(".txt");
    outFile.open(fileout, ios::out);


    double TotalEnergy=0;

    double processEnergy=0;



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
        mat SlaterUpOld=Slatermatrixup(rold);
        mat SlaterDownOld=Slatermatrixdown(rold);
        mat SlaterInvUpOld=SlaterUpInv(SlaterUpOld);
        mat SlaterInvDownOld=SlaterUpInv(SlaterDownOld);


        wfold=wavefunction(rold,SlaterUpOld,SlaterDownOld);
        mat Qforceold=Quantumforce(rold,SlaterInvUpOld,SlaterInvDownOld);





        for(int i=0;i<m_nelectrons;i++){

            for(int j=0;j<2;j++) {
                rnew(i,j)=rold(i,j)+0.5*Qforceold(i,j)*m_dt+norm(gen)*sqrt(m_dt);

            }


            mat SlaterUpNew=Slatermatrixup(rnew);
            mat SlaterDownNew=Slatermatrixdown(rnew);
            mat SlaterInvUpNew=SlaterUpInv(SlaterUpNew);
            mat SlaterInvDownNew=SlaterUpInv(SlaterDownNew);


            wfnew=wavefunction(rnew,SlaterUpNew,SlaterDownNew);
            mat Qforcenew=Quantumforce(rnew,SlaterInvUpNew,SlaterInvDownNew);
            double G=0;
            for(int j=0; j<2;j++){
              G += 0.5*m_dt*(Qforceold(i,j)*Qforceold(i,j)-Qforcenew(i,j)*Qforcenew(i,j))+0.5*(rold(i,j)-rnew(i,j))*(Qforceold(i,j)+Qforcenew(i,j));
            }
            G=exp(G);

            if(uni(gen)<=G*(wfnew*wfnew)/(wfold*wfold)){
                accept++;
                for(int j=0;j<2;j++){
                    rold(i,j)=rnew(i,j);

                }
                Qforceold=Qforcenew;
                wfold=wfnew;
                SlaterUpOld=SlaterUpNew;
                SlaterDownOld=SlaterDownNew;
                SlaterInvUpOld=SlaterInvUpNew;
                SlaterInvDownOld=SlaterInvDownNew;


            }
            else{
                for(int j=0; j<2;j++){
                    rnew(i,j)=rold(i,j);
                }
            }

            double l1=localEnergy(rold,SlaterInvUpOld,SlaterInvDownOld);
            double temp=l1;


                    //outFile.write( (char*)&temp, sizeof(double));
                    outFile << l1 << endl;
            processEnergy+=temp;



        }


    }



    //MPI_Reduce(&processEnergy, &TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



    /*if ( MyRank == 0) {
        double Energy = TotalEnergy/( (double)NumberProcesses*m_cycles);

        cout << Energy << endl;

        //

    }*/
    double Energy = processEnergy/( (double)m_cycles*m_nelectrons);
    cout << Energy << endl;
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

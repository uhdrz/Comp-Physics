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
    m_dt=0.01;
    m_a={1,1.0/3};  //0 for antipar
    m_varpar={1,0.4}; //0=alpha, 1 =beta



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

double VMC::spwf(double x, double y, int nx, int ny){  //set up for harmonic wf

    const double Hnx = HermitePolynomials::evaluate(nx,x,m_w,m_varpar(0));
    const double Hny = HermitePolynomials::evaluate(ny,y,m_w,m_varpar(0));


    double phi=Hnx* Hny*exp(-0.5*m_varpar(0)*m_w*(pow(x,2)+pow(y,2)));

    return phi;
}

mat VMC::Slatermatrixup(mat &r){ //set up for the Slater matrix for Spin up
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

double VMC::wavefunction(mat &r, mat &Slaterup, mat &Slaterdown){ // set up for the entire wf

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

vec VMC::GradJastrow(int k, mat &r){ // calculation for the gradient of the Jastrow factor

    vec deriv=zeros<vec>(2);

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<2;j++){
            if(i!=k){
                double rik=relDis(r,i,k);
                deriv(j)+=(r(k,j)-r(i,j))*derivJastrow(i,k,rik)/rik;
            }

        }
    }


    return deriv;


}

double VMC::derivJastrow(int i, int j, double rij){  // first derivative of the exponent of jastrow
    int n_states=m_nelectrons/2;
    double deriv=0;

    if((i < n_states && j< n_states)||(i >= n_states && j>= n_states)){
        deriv = m_a(1)/(pow(1.0+m_varpar(1)*rij,2));
    }
    else{
        deriv = m_a(0)/(pow(1.0+m_varpar(1)*rij,2));
    }

    return deriv;
}

double VMC::deriv2Jastrow(int i, int j, double rij){ // second derivative of the exponent of jastrow

    int n_states=m_nelectrons/2;
    double deriv2=0;

    if((i < n_states && j< n_states)||(i >= n_states && j>= n_states)){
        deriv2= -2.0*m_a(1)*m_varpar(1)/(pow(1.0+m_varpar(1)*rij,3));
    }
    else{
        deriv2= -2.0*m_a(0)*m_varpar(1)/(pow(1.0+m_varpar(1)*rij,3));

    }

    return deriv2;



}

double VMC::LapJastrow(int k, mat &r){
    double lap=0;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<m_nelectrons;j++){
            if( i!=k){
                if(j!=k){
                    double rkj=relDis(r,k,j);
                    double rki=relDis(r,k,i);
                    for( int d=0; d<2;d++){

                        lap+=(r(k,d)-r(i,d))*(r(k,d)-r(j,d))*derivJastrow(k,i,rki)*derivJastrow(k,j,rkj)/(rki*rkj);

                    }


                 }
            }
        }
        if(i!=k){
            double rki=relDis(r,k,i);
            lap+=deriv2Jastrow(k,i,rki)+1*derivJastrow(k,i,rki)/rki;
        }
    }
    return lap;

}

vec VMC::GradSP(int i, mat &r, mat &InvUp, mat &InvDown){ // Gradient for the space part


    vec Grad=zeros<vec>(2);
    int n_states= m_nelectrons/2; //Spin separatly


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

vec VMC::DerivSP(int i,int nx,int ny, mat &r) { // first derivative for the harmonic wf

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

double VMC::DerivSPA(int i,int nx,int ny, mat &r){ // derivative of harmonic wf after alpha
    double deriv=0;
    const double x = r(i,0);
    const double y = r(i,1);
    const double Hnx = HermitePolynomials::evaluate(nx,x,m_w,m_varpar(0));
    const double Hny = HermitePolynomials::evaluate(ny,y,m_w,m_varpar(0));

    const double dHnx = HermitePolynomials::evaluateAlphaDerivative(nx,x,m_w,m_varpar(0));
    const double dHny = HermitePolynomials::evaluateAlphaDerivative(ny,y,m_w,m_varpar(0));
    double r2 = x*x+y*y;

    deriv=exp(-0.5*m_w*m_varpar(0)*r2)*(Hnx*dHny+Hny*dHnx+Hnx*Hny*(-0.5)*m_w*r2);

    return deriv;


}

double VMC::DerivAlpha( mat &r, mat &InvUp, mat &InvDown){ // whole derivative after alpha
    double deriv=0;
    int n_states= m_nelectrons/2;
    for(int i=0; i<m_nelectrons;i++){

        for(int j=0;j<n_states;j++){
            vec n=postonum(j);
            if(i<n_states){
                deriv+=DerivSPA(i,n(0),n(1),r)*InvUp(j,i);

            }
            else{
                int temp=i-n_states;
                deriv+=DerivSPA(i,n(0),n(1),r)*InvDown(j,temp);


            }
        }
    }


    return deriv;
}

double VMC::DerivBeta( mat &r){   // whole derivative after beta
    double jastrowb=0;
    int n_states= m_nelectrons/2;

    for(int i=0;i<m_nelectrons;i++){
        for(int j=0;j<i;j++){
            double rij=relDis(r,i,j);
            if((i < n_states && j< n_states)||(i >= n_states && j>= n_states)){
                jastrowb-=m_a(1)*rij*rij/(pow(1.0+m_varpar(1)*rij,2));
            }
            else{
                jastrowb-=m_a(0)*rij*rij/(pow(1.0+m_varpar(1)*rij,2));
            }

        }
    }


    return jastrowb;




}

double VMC::locPot(mat &r){  //Potential part of the local energy
    double Pot=0;
    for(int i=0;i<m_nelectrons;i++){
        Pot+=0.5*m_w*m_w*pos2(r,i);
        for(int j=0;j<i;j++){
            Pot+=1.0/relDis(r,i,j);
        }
    }
    return Pot;
}

double VMC::locKin(mat &r, mat &InvUp, mat &InvDown){ // kinetic part of the loca lenergy
    double Kin  = 0;
    double KinJ = 0;

    for(int i=0;i<m_nelectrons;i++){
        Kin  += -0.5*LapSP(i,r,InvUp,InvDown);
        KinJ += -0.5*(LapJastrow(i,r)+2*dot(GradSP(i,r,InvUp,InvDown),GradJastrow(i,r)));//
    }
    return Kin+KinJ;

}

vec VMC::Quantumforce(int i, mat &r, mat &InvUp, mat &InvDown){ //fix m_a

   vec F=zeros<vec>(2);



    vec GradSpace=GradSP(i,r,InvUp,InvDown);
    vec jpart=GradJastrow(i,r);

     for(int k=0;k<2;k++) {
         F(k) = 2.0*(GradSpace(k)+jpart(k));

     }


    return F;

}

void VMC::MCH(int c, char **v){
    int NumberProcesses, MyRank;
    MPI_Init (&c, &v);
    MPI_Comm_size (MPI_COMM_WORLD, &NumberProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &MyRank);
    double StartTime = MPI_Wtime();
    double StartTime1 = MPI_Wtime();


    ofstream outFile;
    ofstream position;




    string filename = "mc"; //
    string fileout = filename;
    string argument = to_string(MyRank);
    string argument3=to_string(m_w);
    fileout.append(argument3);
    fileout.append(argument);
    fileout.append(".txt");
    outFile.open(fileout, ios::out);


    string filename2 = "pos"; //
    string fileout2 = filename2;
    string argument2 = to_string(m_w);
    fileout2.append(argument2);
    fileout2.append(".txt");







    if ( MyRank == 0){position.open(fileout2,ios::out);}


    double TotalEnergy=0;
    double processEnergy=0;
    double TotalKin=0;
    double processKin=0;
    double TotalPot=0;
    double processPot=0;




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


    mat SlaterUpOld=Slatermatrixup(rold);
    mat SlaterDownOld=Slatermatrixdown(rold);
    mat SlaterInvUpOld=SlaterUpInv(SlaterUpOld);
    mat SlaterInvDownOld=SlaterUpInv(SlaterDownOld);
    wfold=wavefunction(rold,SlaterUpOld,SlaterDownOld);

    for(int n=0;n<m_cycles;n++){






        for(int i=0; i < m_nelectrons; i++) {
        vec Qforceold=Quantumforce(i,rold,SlaterInvUpOld,SlaterInvDownOld);






            for(int j=0;j<2;j++) {
                rnew(i,j)=rold(i,j)+0.5*Qforceold(j)*m_dt+norm(gen)*sqrt(m_dt);

            }


            mat SlaterUpNew=Slatermatrixup(rnew);
            mat SlaterDownNew=Slatermatrixdown(rnew);
            mat SlaterInvUpNew=SlaterUpInv(SlaterUpNew);
            mat SlaterInvDownNew=SlaterDownInv(SlaterDownNew);

            wfnew=wavefunction(rnew,SlaterUpNew,SlaterDownNew);
            vec Qforcenew=Quantumforce(i,rnew,SlaterInvUpNew,SlaterInvDownNew);

            double G=0;
            double gexp = 0;
            for(int j=0; j<2;j++){

                double term1 = - (rold(i,j) - rnew(i,j)  -  0.5*m_dt*Qforcenew(j))*
                        (rold(i,j) - rnew(i,j)  -  0.5*m_dt*Qforcenew(j));
                double term2 =   (-rold(i,j) + rnew(i,j) - 0.5*m_dt*Qforceold(j))*
                        (-rold(i,j) + rnew(i,j) - 0.5*m_dt*Qforceold(j));

                gexp += term1+term2;


            }

            G=exp(gexp/(2.0*m_dt));

            if(uni(gen) <= G*(wfnew*wfnew)/(wfold*wfold)){

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

                SlaterUpOld=Slatermatrixup(rold);
                SlaterDownOld=Slatermatrixdown(rold);
                SlaterInvUpOld=SlaterUpInv(SlaterUpOld);
                SlaterInvDownOld=SlaterUpInv(SlaterDownOld);

            }

            double temp=locPot(rold)+locKin(rold,SlaterInvUpOld,SlaterInvDownOld);
            double temp1=locPot(rold);
            double temp2=locKin(rold,SlaterInvUpOld,SlaterInvDownOld);


            processEnergy+=temp;
            processKin+=temp2;
            processPot+=temp1;

            outFile << temp << endl;
        }
        if ( MyRank == 0){rnew.print(position);}

    }
    double EndTime1 = MPI_Wtime();



    MPI_Reduce(&processEnergy, &TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&processPot, &TotalPot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&processKin, &TotalKin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



    if ( MyRank == 0) {
        double Energy = TotalEnergy/( (double)NumberProcesses*m_cycles*m_nelectrons);
        double Kin = TotalKin/( (double)NumberProcesses*m_cycles*m_nelectrons);
        double Pot = TotalPot/( (double)NumberProcesses*m_cycles*m_nelectrons);

        cout << Energy << endl;
        cout << Kin << endl;
        cout << Pot << endl;
        cout<<(double)accept/(m_cycles*m_nelectrons)*100<<endl;
        cout<< m_varpar(0)<<endl;
        cout<< m_varpar(1)<<endl;

        //

    }



    double EndTime = MPI_Wtime();
    double TotalTime = EndTime-StartTime;

    double TotalTime1= EndTime1-StartTime1;
    if ( MyRank == 0 ){

        cout << " total Time = " << TotalTime<<endl;
        cout << "process Time = " << TotalTime1<<endl;

    }








    fflush(stdout);
    outFile.close();

    MPI_Finalize ();




}

void VMC::MCbrute(){


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
    mat SlaterUpOld=Slatermatrixup(rold);
    mat SlaterDownOld=Slatermatrixdown(rold);

    wfold=wavefunction(rold,SlaterUpOld,SlaterDownOld);




    for(int n=0;n<m_cycles;n++){



        for(int i=0;i<m_nelectrons;i++){



            for(int j=0;j<2;j++){
                rnew(i,j)=rold(i,j)+m_step*dis(gen);

          }

         mat SlaterUpNew=Slatermatrixup(rnew);
         mat SlaterDownNew=Slatermatrixdown(rnew);




         wfnew=wavefunction(rnew,SlaterUpNew,SlaterDownNew);

            if(dis1(gen)<=(wfnew*wfnew)/(wfold*wfold)){
                accept++;
                for(int j=0;j<2;j++){
                    rold(i,j)=rnew(i,j);

                }
                SlaterUpOld=SlaterUpNew;
                SlaterDownOld=SlaterDownNew;
                wfold=wfnew;


            }
            else{
                for(int j=0; j<2;j++){
                    rnew(i,j)=rold(i,j);
                }
            }



            mat SlaterInvUpOld=SlaterUpInv(SlaterUpOld);
            mat SlaterInvDownOld=SlaterUpInv(SlaterDownOld);
            double temp=locKin(rold,SlaterInvUpOld,SlaterInvDownOld)+locPot(rold);

            sum+=temp;
            sumsquared+=temp*temp;


        }


    }

    double energy=(double)sum/(m_cycles*m_nelectrons);
    double energysquared=(double)sumsquared/(m_cycles*m_nelectrons);
    double variance=(double)(energysquared - energy*energy)/(m_cycles*m_nelectrons);




    cout<< "Energy"<< energy<<endl;
    cout<< "Variance"<< variance<<endl;
    cout<<(double)accept/(m_nelectrons*m_cycles)*100<<endl;



}

void VMC::findoptParameter(){



    int updates=1000;
    int cycles=10000;
    vec dE=zeros<vec>(2);
    double tolerance = 1.0e-14;
    vec parold=m_varpar;
    vec parnew=zeros<vec>(2);
    double diff=0;
    double step=0.01;






    random_device rnd;
    mt19937 gen(rnd());
    normal_distribution <double>norm(0,1);
    uniform_real_distribution <double>uni(0,1);






    for(int u=0;u<updates;u++){

        double processEnergy=0;
        double processPsia=0;
        double processPsib=0;
        double processPsiEa=0;
        double processPsiEb=0;
        mat rold=zeros<mat>(m_nelectrons,2);
        mat rnew=zeros<mat>(m_nelectrons,2);
        double wfold=0;
        double wfnew=0;

        for(int i=0;i<m_nelectrons;i++){
            for(int j=0;j<2;j++){
                rold(i,j)=norm(gen)*sqrt(m_dt);
            }
        }
        rnew=rold;
        mat SlaterUpOld=Slatermatrixup(rold);
        mat SlaterDownOld=Slatermatrixdown(rold);
        mat SlaterInvUpOld=SlaterUpInv(SlaterUpOld);
        mat SlaterInvDownOld=SlaterUpInv(SlaterDownOld);
        wfold=wavefunction(rold,SlaterUpOld,SlaterDownOld);



        for(int n=0;n<cycles;n++){

            for(int i=0; i < m_nelectrons; i++) {
            vec Qforceold=Quantumforce(i,rold,SlaterInvUpOld,SlaterInvDownOld);





                for(int j=0;j<2;j++) {
                    rnew(i,j)=rold(i,j)+0.5*Qforceold(j)*m_dt+norm(gen)*sqrt(m_dt);

                }


                mat SlaterUpNew=Slatermatrixup(rnew);
                mat SlaterDownNew=Slatermatrixdown(rnew);
                mat SlaterInvUpNew=SlaterUpInv(SlaterUpNew);
                mat SlaterInvDownNew=SlaterDownInv(SlaterDownNew);

                wfnew=wavefunction(rnew,SlaterUpNew,SlaterDownNew);
                vec Qforcenew=Quantumforce(i,rnew,SlaterInvUpNew,SlaterInvDownNew);

                double G=0;
                double gexp = 0;
                for(int j=0; j<2;j++){

                    double term1 = - (rold(i,j) - rnew(i,j)  -  0.5*m_dt*Qforcenew(j))*
                            (rold(i,j) - rnew(i,j)  -  0.5*m_dt*Qforcenew(j));
                    double term2 =   (-rold(i,j) + rnew(i,j) - 0.5*m_dt*Qforceold(j))*
                            (-rold(i,j) + rnew(i,j) - 0.5*m_dt*Qforceold(j));

                    gexp += term1+term2;


                }

                G=exp(gexp/(2.0*m_dt));

                if(uni(gen) <= G*(wfnew*wfnew)/(wfold*wfold)){


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

                    SlaterUpOld=Slatermatrixup(rold);
                    SlaterDownOld=Slatermatrixdown(rold);
                    SlaterInvUpOld=SlaterUpInv(SlaterUpOld);
                    SlaterInvDownOld=SlaterUpInv(SlaterDownOld);

                }
                double temp=locKin(rold,SlaterInvUpOld,SlaterInvDownOld)+locPot(rold);
                double tempPsia=DerivAlpha(rold,SlaterInvUpOld,SlaterInvDownOld);
                double tempPsib=DerivBeta(rold);
                processEnergy+=temp;
                processPsia+=tempPsia;
                processPsib+=tempPsib;
                processPsiEa+=tempPsia*temp;
                processPsiEb+=tempPsib*temp;


            }


        }











        double Energy = processEnergy/( (double)cycles*m_nelectrons);
        double Psia = processPsia/( (double)cycles*m_nelectrons);
        double Psib = processPsib/( (double)cycles*m_nelectrons);
        double PsiEa = processPsiEa/( (double)cycles*m_nelectrons);
        double PsiEb = processPsiEb/( (double)cycles*m_nelectrons);



         dE(0)=2*(PsiEa-Psia*Energy);
         dE(1)=2*(PsiEb-Psib*Energy);
         parnew=parold-step*dE;



          m_varpar=parnew;
          parold=parnew;
          vec difference=parnew-parold;



          diff= dot(difference,difference);









        if(diff<tolerance){ break;}

          cout <<u<<endl;


    }











}

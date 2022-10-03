#include <iostream>
#include <fstream>
#include "cuba.h"
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
/*****************************begin const****************************************/
#define M_PI 3.14159265358979323846
double sq(double a){return pow(a,2);}
double tao,Q2,x,y,maxr,factor,alphaem,Nc,sumeq,CA,CF,A,Sperp,x0,lambda,Qs2,qp,forL,forT,xB,Lambd;
double a[32400][3];
int nplot;
double taomin,taobin;
Vector2d k;

void dushujv(){
    ifstream inFile;
    inFile.open("/home/matata/CLionProjects/GTMD/rcBK_MVe_Heikki_table.txt");
    for(int i=0;i<32400;i++){
        for(int j=0;j<3;j++){
            inFile >> a[i][j];
        }
    }
}
void constants() {
    A=600.0;
    x0=2.47e-5;
    lambda=0.282;
    Lambd=0.241;
    Sperp=23.58;
    CA=3.0;
    CF=4.0/3;
    Nc=3.0;
    alphaem=1.0/137;
    Q2=4.0;
    x=1e-3;
    xB=x/(1-x);
    nplot=40;
    y= 0.8;
    maxr=1000.0;
    taomin=-1.0;
    taobin=-3.0/nplot;
    sumeq=1.0;
    Qs2=pow(x0/x,lambda)*pow(A,1.0/3);
    factor=Nc*alphaem*sumeq;
    forL=alphaem/xB/y/M_PI*(1.0-y);
    forT=alphaem/xB/y/M_PI*(1.0+sq(1.0-y))/2.0;
    k<<1.5,0;
}
/*****************************end const****************************************/

/*****************************hard function****************************************/
double S(double r){
    return a[int(400*log(0.01/x)/0.2)][2];
}
double HT(double xi){
    double epsilon=xi*(1-xi)*Q2;
    return xi*(1-xi)*(sq(xi)+sq(1-xi))*2*k.dot(k)*(sq(k.dot(k))+sq(epsilon))/pow(k.dot(k)+epsilon,4);
}
double HL(double xi){
    double epsilon=xi*(1-xi)*Q2;
    return pow(xi*(1-xi),3)*8*Q2*2*sq(k.dot(k))/pow(k.dot(k)+epsilon,4);

}
/*****************************end hard function****************************************/

/***************************** GBW ****************************************/
double xG(double q,double r){
    return 4.0/r*cyl_bessel_j(0,r*q)*(1-exp(-CA/4.0/CF*Qs2*sq(r)))*Sperp*(sq(Nc)-1)/pow(2*M_PI,3)/Nc;
}
double xh(double q,double r){
    return 0;
}
double xGmv(double q,double r){
    double mgj=(4 - 3*exp(1) *Lambd*r + 4 *sq(1 + exp(1) *Lambd*r)* log(exp(1) + 1.0/(Lambd*r)))/(sq(r*(1 + exp(1) *Lambd*r))*log(exp(1) + 1.0/(Lambd*r)));
    return r*cyl_bessel_j(0,r*q)*(1-exp(-CA/4.0/CF*Qs2*sq(r)*log(1.0/Lambd/r+exp(1))))*Sperp*(sq(Nc)-1)/pow(2*M_PI,3)/Nc*mgj;
}


/*****************************end GBW****************************************/
double integrandL(double xi,double q,double r){
    return HL(xi)*xGmv(q,r);
}

double integrandT(double xi,double q,double r){
    double epsilon=xi*(1-xi)*Q2;
    return HT(xi)*xGmv(q,r);
}

static int dsigmaL(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=maxr*xx[1],phi=2*M_PI*xx[2]
    ff[0]=2*M_PI*maxr*integrandL(xx[0],2*sqrt(tao*k.dot(k)),maxr*xx[1]);
    return 0;
}


static int dsigmaT(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    ff[0]=2*M_PI*maxr*integrandT(xx[0],2*sqrt(tao*k.dot(k)),maxr*xx[1]);
    return 0;
}
static int dsigmaLT(const int *ndim, const cubareal xx[],
                    const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=xx[1]*maxr
    ff[0]=4*sq(xx[0]*(1.0-xx[0]))*Q2*cyl_bessel_k(0.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))
          *cyl_bessel_k(0.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))*(1-exp(-CA/4.0/CF*Qs2*sq(xx[1]*maxr)*log(1.0/Lambd/maxr/xx[1]+exp(1))));
    return 0;
}
static int dsigmaTT(const int *ndim, const cubareal xx[],
                    const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=xx[1]*maxr
    ff[0]=(sq(xx[0])+sq(1.0-xx[0]))*Q2*xx[0]*(1.0-xx[0])*cyl_bessel_k(1.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))
          *cyl_bessel_k(1.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))*(1-exp(-CA/4.0/CF*Qs2*sq(xx[1]*maxr)*log(1.0/Lambd/maxr/xx[1]+exp(1))));
    return 0;
}


int main() {
    ofstream mycout("/home/matata/CLionProjects/GTMD/DijetEEC2.txt");
    constants();
    cout<<k;
    cout<<Qs2<<endl;
    double xT[nplot],yT[nplot];
    double xL[nplot],yL[nplot];
    double normT[nplot],normL[nplot];
    double sigmaT,sigmaL;
    int neval, fail;
    cubareal integral[1], error[1], prob[1];
    Vegas(2, 1, dsigmaTT, nullptr, 1,
          1e-2, 1e-12, 1, 0,
          1e2, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    sigmaT=integral[0];
    Vegas(2, 1, dsigmaLT, nullptr, 1,
          1e-2, 1e-12, 1, 0,
          1e2, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    sigmaL=integral[0];
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(3, 1, dsigmaT, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e5, 1e7, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xT[i] = factor*integral[0];
        yT[i] = factor*error[0];
        normT[i] = integral[0]/sigmaT;
    }
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(3, 1, dsigmaL, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e5, 1e7, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xL[i] = factor*integral[0];
        yL[i] = factor*error[0];
        normL[i] = integral[0]/sigmaL;
    }
    mycout<<"tao2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<pow(10,taomin+taobin*i)<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"tot2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<forL*xL[i]+forT*xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"nxT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<normT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"nxL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<normL[i]<<',';
    }
    mycout<<']'<<endl;

    mycout<<"ntot2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<(forL*xL[i]+forT*xT[i])/(forL*sigmaL+forT*sigmaT)<<',';
    }
    mycout<<']'<<endl;
    mycout.close();
    cout<<sigmaT<<endl;
    cout<<sigmaL<<endl;
    cout<<Qs2<<endl;
    return 0;
}
//
// Created by matata on 3/4/22.
//


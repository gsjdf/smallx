#include <iostream>
#include <fstream>
#include "cuba.h"
#include <cmath>
#include <Eigen/Dense>
#include "LHAPDF/PDF.h"

using namespace std;
using namespace Eigen;
using namespace LHAPDF;
/*****************************begin const****************************************/
#define M_PI 3.14159265358979323846
double sq(double a){return pow(a,2);}
double tao,Q2,x,y,maxr,factor,alphaem,Nc,sumeq,CA,CF,A,Sperp,x0,lambda,Qs2,qp,forL,forT,xB,pT2j,maxl;
double a[32400][3];
int nplot;
double taomin,taobin;
Vector2d k;
const PDF* FrF= mkPDF("JAM20-SIDIS_FF_hadron_nlo");
void constants() {
    A=79.0;
    x0=2.24e-4;
    lambda=0.27;
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
    maxr=100.0;
    taomin=-1.0;
    taobin=-3.0/nplot;
    sumeq=1.0;
    Qs2=pow(x0/x,0.27)*pow(A,1.0/3);
    maxl=Qs2;
    factor=Nc*alphaem;
    forL=alphaem/xB/y/M_PI*(1.0-y);
    forT=alphaem/xB/y/M_PI*(1.0+sq(1.0-y))/2.0;
    k<<1.5,0;
    ///////jet const////
    pT2j=0.12;
}
/*****************************end const****************************************/

/*****************************hard function****************************************/

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
    return 4/r*cyl_bessel_j(0,r*q)*(1-exp(-CA/4/CF*Qs2*sq(r)))*Sperp*(sq(Nc)-1)/pow(2*M_PI,3)/Nc;
}

/*****************************end GBW****************************************/
double integrandL(double xi,double q,double r){
    return HL(xi)*xG(q,r);
}

double integrandT(double xi,double q,double r){
    return HT(xi)*xG(q,r);
}
/*****************************jet function****************************************/
static int jet(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    //z1=xx[0], z2=xx[1]
    ff[0]=2.0/pT2j/(sq(xx[0]*(1-1e-2)+1e-2)+sq(xx[1]*(1-1e-2)+1e-2))*exp(-4*sq(xx[0]*(1-1e-2)+1e-2*xx[1]*(1-1e-2)+1e-2)/pT2j/(sq(xx[0]*(1-1e-2)+1e-2)+sq(xx[1]*(1-1e-2)+1e-2))*k.dot(k)*tao)
          *(((FrF->xfxQ2(1, xx[0]*(1-1e-2)+1e-2,Q2))+(FrF->xfxQ2(-1, xx[0]*(1-1e-2)+1e-2,Q2)))*((FrF->xfxQ2(1, xx[1]*(1-1e-2)+1e-2,Q2))+(FrF->xfxQ2(-1, xx[1]*(1-1e-2)+1e-2,Q2)))*8/9
            +((FrF->xfxQ2(2, xx[0]*(1-1e-2)+1e-2,Q2))+(FrF->xfxQ2(-2, xx[0]*(1-1e-2)+1e-2,Q2)))*((FrF->xfxQ2(2, xx[1]*(1-1e-2)+1e-2,Q2))+(FrF->xfxQ2(-2, xx[1]*(1-1e-2)+1e-2,Q2)))/9);
    return 0;
}
/*****************************end jet function****************************************/
static int dsigmaL(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=maxr*xx[1], l=xx[2]*maxl
    ff[0]=2*M_PI*maxr*integrandL(xx[0],1e-7,maxr*xx[1]);//2*M_PI*maxr*maxl*integrandL(xx[0],maxl*xx[2],maxr*xx[1]);
    return 0;
}


static int dsigmaT(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    ff[0]=2*M_PI*maxr*integrandT(xx[0],1e-7,maxr*xx[1]);//2*M_PI*maxr*maxl*integrandT(xx[0],maxl*xx[2],maxr*xx[1]);
    return 0;
}


int main() {
    ofstream mycout("/home/matata/CLionProjects/GTMD/DijetEEC2.txt");
    constants();
    double xT[nplot],yT[nplot];
    double xL[nplot],yL[nplot];
    int neval, fail;
    cubareal integral[1], error[1], prob[1];
    double factorL,factorT;

    Vegas(3, 1, dsigmaT, nullptr, 1,
          1e-1, 1e-12, 0, 0,
          1e6, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    factorT = integral[0];

    Vegas(3, 1, dsigmaL, nullptr, 1,
          1e-1, 1e-12, 0, 0,
          1e6, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    factorL = integral[0];
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(2, 1, jet, nullptr, 1,
              1e-1, 1e-12, 1, 0,
              1e6, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xT[i] = integral[0];
        xL[i] = integral[0];
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
    mycout<<"xL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xL[i]<<',';
    }
    mycout<<']'<<endl;

    mycout<<"tot2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<forL*xL[i]+forT*xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout.close();
    return 0;
}
//
// Created by matata on 3/4/22.
//


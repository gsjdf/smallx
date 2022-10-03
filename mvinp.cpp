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
double x,qmin,qmax,maxr,A,x0,lambda,Qs2,Lambd,q;
int nq;
void constants() {
    qmax=2.0;
    qmin=-6.0;
    A=600.0;
    x0=2.47e-5;
    lambda=0.282;
    Lambd=0.241;
    x=1e-3;
    nq=1e3;
    Qs2=pow(x0/x,lambda);
    maxr=1e4;
}


/*****************************end const****************************************/




static int dF(const int *ndim, const cubareal xx[],
                    const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=xx[1]*maxr
    double r=xx[0]*maxr;
    ff[0]=r*pow(exp(1)+1.0/r/Lambd,-Qs2*sq(r)/4)*cyl_bessel_j(0,r* q)/2/M_PI*maxr;
    return 0;
}

int main() {
    ofstream mycout("/home/matata/CLionProjects/GTMD/MVparton.txt");
    constants();
    cout<<1.0/Qs2/M_PI;
    double mv[nq];
    int neval, fail;
    cubareal integral[1], error[1], prob[1];
    for (int i=0;i<nq;i+=1) {
        q = pow(10,qmin+(qmax-qmin)/nq*i);
        Vegas(4, 1, dF, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1e3,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        mv[i]=integral[0];
    }
    for (int i=0;i<nq;i+=1){
        mycout<<pow(10,qmin+(qmax-qmin)/nq*i)<<" ";
    }
    mycout<<endl;

    for (int i=0;i<nq;i+=1){
        mycout<<mv[i]<<" ";
    }
    mycout<<endl;
    cout<<Qs2<<endl;
    return 0;
}
//
// Created by matata on 3/4/22.
//


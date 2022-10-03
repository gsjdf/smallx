#include <iostream>
#include <fstream>
#include "cuba.h"
#include <cmath>
#include <Eigen/Dense>
#include "spline.h"

using namespace std;
using namespace Eigen;
/*****************************begin const****************************************/
#define M_PI 3.14159265358979323846
double sq(double a){return pow(a,2);}
double tao,Q2,x,y,qmax,maxr,factor,alphaem,Nc,sumeq,CA,CF,A,Sperp,x0,lambda,Qs2,qp,forL,forT,xB,Lambd;
int nplot;
double taomin,taobin;
double a[2][1000];
tk::spline Fmv;
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
    qmax=100.0;
    maxr=100.0;
    taomin=-1;
    taobin=-3.0/nplot;
    sumeq=1.0;
    Qs2=pow(x0/x,lambda);
    qp=sqrt(Q2/2/x);
    factor=2*Nc*alphaem*sumeq*Sperp;
    forL=alphaem/xB/y/M_PI*(1.0-y);
    forT=alphaem/xB/y/M_PI*(1.0+sq(1.0-y))/2.0;
}

void dushujv(){
    ifstream inFile;
    inFile.open("/home/matata/CLionProjects/GTMD/MVparton.txt");
    for(int i=0;i<2;i++){
        for(int j=0;j<1000;j++){
            inFile >> a[i][j];
        }
    }
}
void Fmvcs(){
    vector<double> X(1000), Y(1000);
    for(int i=0;i<1000;i++)
    {
        X[i]=a[0][i];
        Y[i] = a[1][i];
        if (Y[i]<1e-6) {
            Y[i] = 0.0;
        }
    }
    Fmv.set_points(X,Y);
    Fmv.make_monotonic();
}
/*****************************end const****************************************/

/*****************************hard function****************************************/

double HT(Vector2d k, Vector2d q, double xi){
    Vector2d ans,m;
    m=k-q;
    double epsilon=xi*(1-xi)*Q2;
    ans=k/(k.dot(k)+epsilon)-m/(m.dot(m)+epsilon);
    return (sq(xi)+sq(1-xi))*ans.dot(ans);
}
double HL(Vector2d k, Vector2d q, double xi){
    Vector2d m;
    m=k-q;
    double epsilon=xi*(1-xi)*Q2,ans;
    ans=1/(k.dot(k)+epsilon)-1/(m.dot(m)+epsilon);
    return 4*sq(xi)*sq(1-xi)*sq(ans)*Q2;
}

/*****************************end hard function****************************************/

/*****************************general pdf****************************************/
double F1(Vector2d q){
    return exp(-q.dot(q)/Qs2)/M_PI/Qs2;
}

/*****************************end general pdf****************************************/
double integrandL(double xi,Vector2d k,Vector2d q){
    return pow(xi,3)*Fmv(sqrt(q.dot(q)))*HL(k,q,xi);
    //return pow(xi,3)*F1(q)*HL(k,q,xi);
}

double integrandT(double xi,Vector2d k,Vector2d q){
    return pow(xi,3)*Fmv(sqrt(q.dot(q)))*HT(k,q,xi);
    //return pow(xi,3)*F1(q)*HT(k,q,xi);
}

static int dsigmaL(const int *ndim, const cubareal xx[],
               const int *ncomp, cubareal ff[], void *userdata) {
    //k is k ,l is q
    Vector2d ktem,k,q,b,qtem;
    ktem<<1,0;
    k=ktem*xx[2]*qp*sqrt(tao);
    qtem<<cos(2*xx[1]*M_PI),sin(2*xx[1]*M_PI);
    q=qtem*qmax*xx[0];
    ff[0]=sq(2*M_PI)*sq(qmax)*xx[0]*integrandL(xx[2],k,q);
    return 0;
}


static int dsigmaT(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    Vector2d ktem,k,q,b,qtem;
    ktem<<1,0;
    k=ktem*xx[2]*qp*sqrt(tao);
    qtem<<cos(2*xx[1]*M_PI),sin(2*xx[1]*M_PI);
    q=qtem*qmax*xx[0];
    ff[0]=sq(2*M_PI)*sq(qmax)*xx[0]*integrandT(xx[2],k,q);
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
    ofstream mycout("/home/matata/CLionProjects/GTMD/SIDISEEC1.txt");
    constants();
    dushujv();
    Fmvcs();
    double xT[nplot],yT[nplot];
    double xL[nplot],yL[nplot];
    double normT[nplot],normL[nplot];
    double sigmaT,sigmaL;
    int neval, fail;
    cubareal integral[1], error[1], prob[1];
    Vegas(2, 1, dsigmaTT, nullptr, 1,
          1e-3, 1e-12, 1, 0,
          1e2, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    sigmaT=integral[0];
    Vegas(2, 1, dsigmaLT, nullptr, 1,
          1e-3, 1e-12, 1, 0,
          1e2, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    sigmaL=integral[0];
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(3, 1, dsigmaT, nullptr, 1,
              1e-3, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xT[i] = sq(qp)*factor*integral[0]/2;
        yT[i] = sq(qp)*factor*error[0]/2;
        normT[i] = integral[0]/sigmaT;
    }
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(3, 1, dsigmaL, nullptr, 1,
              1e-3, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xL[i] = sq(qp)*factor*integral[0]/2;
        yL[i] = sq(qp)*factor*error[0]/2;
        normL[i] = integral[0]/sigmaL;
    }
    mycout<<"tao=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<pow(10,taomin+taobin*i)<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xT=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yT=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xL=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yL=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"tot=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<forL*xL[i]+forT*xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"nxT=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<normT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"nxL=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<normL[i]<<',';
    }
    mycout<<']'<<endl;

    mycout<<"ntot=[";
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


#include <stdio.h>
#include <cmath>
#include <math.h>
#include <iostream>
//constantes 
#define m 800
#define h 1
#define w 60*2.6617E-6
#define G 60*60*6.67E-11
#define Mt 5.9736E24
#define Ml 0.07349E24
#define dtl 3.844E8
#define Rt 6.378160E6
#define Rl 1.7374E6
#define nabla G*Mt/(dtl*dtl*dtl)
#define nu Ml/Mt
#define N 10000

double cambiovelocidad(double v);
double calcrprima(double r, double phi, double t);
double cambiounidadr(double r);
double prinicial(double vr, double theta, double phi);
double mhiinicial(double r, double vr, double theta, double phi);
double cambiounidadmhi(double mhi);
double calck1(double pr);
double calck2(double r, double rprima, double phi, double mhi, double t);
double calck3(double r, double mhi);
double calck4(double r, double rprima, double phi, double t);
double calculoT(double pr, double mhi,double r);
double calpculoP(double r, double rprima);


int main(){
    double r, pr, phi, mhi, vr, theta, rprima, xl,yl,xc,yc, T, P, H; 
    double k11,k12,k13,k14;
    double k21, k22,k23,k24;
    double k31, k32,k33,k34;
    double k41, k42,k43,k44;
    double t=0;

    //inicializar constantes
    r=Rt;
    vr=11200;
    theta=15*M_PI/180;
    phi=0*M_PI/180;
    vr=cambiovelocidad(vr);
    r=cambiounidadr(r);
    pr=prinicial(vr,theta,phi);
    mhi=mhiinicial(r, vr, theta, phi);


    FILE*fposiciones;
    fposiciones = fopen("cohete.txt","w");

    FILE*fhamiltoniano;
    fhamiltoniano= fopen("hamiltoniano.txt", "w");

    //bucle principal
    
    for (int i = 0; i < N; i++)
    {
        xc=r*cos(phi);
        yc=r*sin(phi);
        xl=cos(w*t);
        yl=sin(w*t);


        rprima=calcrprima(r, phi, t);
        T=calculoT(pr,mhi,r);
        P=calpculoP(r, rprima);
        H=T+P;

        fprintf(fposiciones, "%lf, %lf\n",xc, yc);
        fprintf(fposiciones, "%lf, %lf\n\n", xl, yl);
        
        fprintf(fhamiltoniano, "%lf %lf \n", t, H);


        
    k11=h*calck1(pr);
    k21=h*calck2(r,rprima, phi, mhi,t);
    k31=h*calck3(r, mhi);
    k41=h*calck4(r, rprima, phi, t);

    rprima=calcrprima(r+k11/2.0, phi+k21/2.0, t+h/2);
    k12=h*calck1(pr+k21/2.0);
    k22=h*calck2(r+k11/2.0,rprima,phi+k21/2.0,mhi+k41/2.0, t+h/2);
    k32=h*calck3(r+k11/2.0,mhi+k41/2.0);
    k42=h*calck4(r+k11/2.0, rprima, phi+k31/2.0, t+h/2);

    rprima=calcrprima(r+k12/2.0, phi+k32/2.0, t+h/2);
    k13=h*calck1(pr+k22/2.0);
    k23=h*calck2(r+k12/2.0,rprima,phi+k32/2.0,mhi+k42/2.0, t+h/2);
    k33=h*calck3(r+k12/2.0,mhi+k42/2.0);
    k43=h*calck4(r+k12/2.0, rprima, phi+k32/2.0, t+h/2);

    rprima=calcrprima(r+k13, phi+k33, t+h);
    k14=h*calck1(pr+k23);
    k24=h*calck2(r+k13,rprima,phi+k33,mhi+k43, t+h);
    k34=h*calck3(r+k13,mhi+k43);
    k44=h*calck4(r+k13, rprima, phi+k33, t+h);

    r=r+(k11+2*k12+2*k13+k14)/6;
    pr=pr+(k21+2*k22+2*k23+k24)/6;
    phi=phi+(k31+2*k32+2*k33+k34)/6;
    mhi=mhi+(k41+2*k42+2*k43+k44)/6;



    t=t+h;
    }
    
    fclose(fposiciones);
    



   

    return 0;
}


double calcrprima(double r, double phi, double t)
{
    double rprima;
    rprima=sqrt(1+r*r-2*r*cos(phi-w*t));
    return rprima;
}

double cambiounidadr(double r){
    return r/dtl;
}
double cambiovelocidad(double v){
    return v*60/dtl;
}
double prinicial(double vr, double theta, double phi){
    double pr;
    pr=vr*cos(theta-phi);
    return pr; 
}
double mhiinicial(double r, double vr, double theta, double phi)
{
    double mhi;
    mhi=r*vr*sin(theta-phi);
    return mhi;
}
double cambiounidadmhi(double mhi)
{
    return mhi/(m*dtl*dtl);
}
//primer numero es la coordenada que corresponde, el segundo la K 1,2,3,4.
//1: r
//2: pr
//3: phi
//4: mhi
double calck1(double pr)
{
    return pr;
}

double calck2(double r, double rprima, double phi, double mhi, double t)
{
    double k21;
    k21=mhi*mhi/(pow(r,3))-nabla*(1/(r*r)+nu*(r-cos(phi-w*t))/(pow(rprima,3)));
    return k21;
}

double calck3(double r, double mhi)
{
    return mhi/(r*r);
}

double calck4(double r, double rprima, double phi, double t)
{
    double k41;
    k41=-1*nabla*nu*r*sin(phi-w*t)/pow(rprima,3);
    return k41;
}

double calculoT(double pr, double mhi,double r)
{
    //reconvierto las variables
    double prc, mhic, rc, T;
    prc=pr*m*dtl;
    mhic=mhi*m*dtl*dtl;
    rc=r*dtl;
    T=pow(prc,2)/(2*m)+pow(mhic,2)/(2*m*pow(rc,2));
    return T;
}

double calpculoP(double r, double rprima)
{
    double rc, rprimac, P;
    rc=r*dtl;
    rprimac=rprima*dtl;
    P=-G*m*Mt/rc-G*m*Ml/rprimac;
    return P;

}


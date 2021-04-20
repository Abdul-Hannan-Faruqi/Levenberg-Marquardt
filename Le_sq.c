#include "LM.h"
#include "LinAlg.h"

#define n 2
#define h 1e-5
#define nf 3

void Jacobian(double *x);
double ef[nf][1], Jt[n][nf];

double fx(double *x)
{
double f=0;
int i;

    ef[0][0] = pow(x[0],2)+2*pow(x[1],2)-0.3*cos(3*M_PI*x[0])-0.4*cos(4*M_PI*x[1])+0.7;
    ef[1][0] = pow(x[0],2)+2*pow(x[1],2)-0.3*cos(3*M_PI*x[0])*cos(4*M_PI*x[1])+0.3;
    ef[2][0] = pow(x[0],2)+2*pow(x[1],2)-0.3*cos(3*M_PI*x[0]+4*M_PI*x[1])+0.3;

    for(i=0;i<nf;i++)
        f = f + pow(ef[i][0],2);

return f/2;
}

void Jacobian(double *x){
int row, col;
double x_f[n] = {0}, x_b[n] = {0}, e_f, e_b;

for(col=0; col<nf; col++){
    for (row=0;row<n;row++){
        double shift[n] = {0};
        shift[row] = 1;
        vec_op(n,1,h,x, shift,x_f,0);
        vec_op(n,1,h,x, shift,x_b,1);
        fx(x_f);
        e_f = ef[col][0];
        fx(x_b);
        e_b = ef[col][0];
        Jt[row][col] = (e_f-e_b)/(2*h);
    }
}
}

void grad(double *x, double *g, int fac){
double temp[n][1];
int i;

    Jacobian(x);
    fx(x);
    mult(n,nf,Jt,1,ef,temp);

    for(i=0;i<n;i++)
        g[i] = temp[i][0];
}

void Hessian(double *x, double H[][n]){
double J[nf][n];
int i,j;

for(i=0; i<nf;i++){
   for(j=0; j<n; j++){
      J[i][j] = Jt[j][i];
      }
}

mult(n,n,Jt,n,J,H);
}

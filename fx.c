#include "LM.h"
#include "LinAlg.h"

#define n 20
#define d 1e-5

double fx(double *x)
{
int i;
double f, t1, t2 = 0;

for(i=0;i<n;i++){
if (i==0)
   t1 = pow(x[0]-1,2);
else{
   t1 = t1 + pow(x[i]-1,2);
   t2 = t2 + x[i]*x[i-1];
}
}
f = t1-t2;
return f;
}

void Hessian(double *x, double H[][n]){
int i, j;
double h = d/2;
double x_f[n] = {0}, x_b[n] = {0}, g_f[n] = {0}, g_b[n] = {0};

    for(i=0;i<n;i++){
        double shift[n] = {0};
        shift[i] = 1;
        vec_op(n,1,h,x, shift,x_f,0);
        vec_op(n,1,h,x, shift,x_b,1);
        grad(x_b,g_b,2);
        grad(x_f,g_f,2);
        for (j=0;j<n;j++)
            H[i][j] = (g_f[j]-g_b[j])/(2*h);
    }
}


void grad(double *x, double *g, int fac){
int i, j;
double x_f[n] = {0}, x_b[n] = {0};
double h = d/fac;

    for (i=0;i<n;i++){
        double shift[n] = {0};
        shift[i] = 1;
        vec_op(n,1,h,x, shift,x_f,0);
        vec_op(n,1,h,x, shift,x_b,1);
        g[i] = (fx(x_f)-fx(x_b))/(2*h);
    }
}

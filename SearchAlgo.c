#include "SearchAlgo.h"
#include "LinAlg.h"
#include "LM.h"

#define phi (1+sqrt(5))/2

int max(double *arr);

double Ex_LS(int n, double *x, double *g, double *d){
double xn[n], dp=10, alpha = 20, e = 1e-4, a, b, fa, fb, fval[2], lower, upper, s;
int i, flag = 0;

vec_op(n, 1, alpha, x, d, xn, 0);
grad(xn, g, 1);
dp = dot(n,d,g);
while(dp<0){
alpha += 5;
vec_op(n, 1, alpha, x, d, xn, 0);
grad(xn, g, 1);
dp = dot(n,d,g);
}

upper = alpha;
lower = 0;
s = upper - lower;
a = lower + (2-phi)*s;
b = lower + (phi-1)*s;
vec_op(n, 1, lower, x, d, xn, 0);
fval[0] = fx(xn);
vec_op(n, 1, a, x, d, xn, 0);
fa = fx(xn);
vec_op(n, 1, b, x, d, xn, 0);
fb = fx(xn);
vec_op(n, 1, upper, x, d, xn, 0);
fval[1] = fx(xn);

while(dp>e){
flag = 1;
i = max(fval);
if (i==0){
   lower = a;
   fval[0] = fa;
   a = b;
   fa = fb;
   s = upper - lower;
   b = lower + (phi-1)*s;
   vec_op(n, 1, b, x, d, xn, 0);
   fb = fx(xn);
}
else{
   upper = b;
   fval[1] = fb;
   b = a;
   fb = fa;
   s = upper-lower;
   a = lower + (2-phi)*s;
   vec_op(n, 1, a, x, d, xn, 0);
   fa = fx(xn);
}
grad(xn, g, 1);
dp = fabs(dot(n,d,g));
if(s<0.1)
break;
}
if(flag==1)
   alpha = (a+b)/2;

return alpha;
}

double Inex_LS(int n, double *x, double *g, double *d){
int i;
double alpha = 0.1, mu = 0.1, dp, diff=1, xn[n], gn[n], norm_g, ul;

norm_g = norm(n,g);
for(i=0;i<n;i++)
   gn[i] = g[i]/norm_g;

while(diff>0){
   alpha = alpha*1.2;
   ul = fx(x)+alpha*mu*dot(n,d,g);
   vec_op(n,1,alpha,x,d,xn,0);
   diff = ul-fx(xn);
}
return alpha/1.2;
}


int max(double *arr){
return arr[1]>arr[0];
}

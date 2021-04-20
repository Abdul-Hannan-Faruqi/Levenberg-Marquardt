#include <stdio.h>
#include <stdlib.h>
#include "LM.h"
#include "LinAlg.h"
#include "SearchAlgo.h"
#include <string.h>
#include <time.h>

#define d 1e-5            //difference used for central differencing
#define n 20

void Hessian(double *, double H[][n]);

int main(int argc, char *argv[]){

int i,j,iter;
int search = atoi(argv[1]);
double x[n] = {0}, xn[n], g[n], neg_g[n], H[n][n], dir[n];
double err = 1, e = 1e-7, alpha, del, gamma = atoi(argv[2]), mu1 = atoi(argv[3]), mu2 = atoi(argv[4]);
FILE *fp;
fp = fopen("log.dat", "a");     //Data logging
if(fp == NULL){
    printf("Couldn't open file\n");
    return 1;
}

iter = 0;

clock_t begin = clock();
while(err > e)
{
    iter++;
    printf("\niter = %d\n", iter);

/*    printf("x = \n");
    for(i=0;i<n;i++)
       printf("%lf\n",x[i]);
*/
    grad(x,g,1);
    Hessian(x, H);

    printf("f(x) = %lf\n", fx(x));
    printf("Norm(G) = %lf\n", norm(n,g));

//___________Modify Hessian____________
    for(i=0;i<n;i++)
        H[i][i] = H[i][i] + gamma;
//_____________________________________

    for(i=0;i<n;i++){
       neg_g[i] = -g[i];
    }

    GE(n, H, neg_g, dir);

    if(search==0)
       alpha = Ex_LS(n, x, g, dir);
    else
       alpha = Inex_LS(n, x, g, dir);

    vec_op(n, 1, alpha, x, dir, xn, 0);
    del = fx(xn)-fx(x);

    if(del<0){
        gamma = gamma/mu1;
        for(i=0; i<n; i++)
           x[i] = xn[i];
    }
    else
        gamma = mu2*gamma;

err = fabs(del);
printf("err = %lf\n", err);
}
clock_t end= clock();
double time=(double)(end-begin)/CLOCKS_PER_SEC;
    printf("x = \n");
    for(i=0;i<n;i++)
       printf("%lf\n",x[i]);
printf("\n\nCompleted in %g s\n", time);
fprintf(fp, "%d, %d, %g, %g, %g, %d, %g\n", search, atoi(argv[2]), mu1, mu2, fx(x), iter, time);
return 0;
}


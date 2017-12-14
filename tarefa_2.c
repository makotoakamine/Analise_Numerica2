#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "BasicFunc.h"
#define PI 3.14159265358979323846

double func(double x){
    double r = 0.000;
    r = sin(x);
    return r;
}

double Lk(double *X,int n,double x,int k){ //n e' o tamanho do vetor X
    int i = 0;
    double ret = 1;
    for(i=0;i<n;i++){
        if(i!=k){
            ret = (ret*(x-X[i]))/(X[k]-X[i]);
        }
    }
    return ret;
}

double Pol(double *X,double *Y,double x,int n){
    int i=0;
    double ret=0;
    for(i=0;i<n;i++){
        ret += Y[i]*Lk(X,n,x,i);
    }
    return ret;
}

int main(){
    
    int n = 0;
    for(n=2 ; n<=50 ; n++){
        double *Xeq = allocV(n+1);
        double *Yeq = allocV(n+1);
        double *Xcheb = allocV(n+1);
        double *Ycheb = allocV(n+1);
        double *x = allocV(10117);
        double max_eq = 0.0000;
        double max_cheb = 0.000;
        double temp = 0.0000;
        int i=0;
        
        //Particionar dominio
        for(i = 0;i<n+1; i++){
           Xeq[i] = -1.0000 + 2.0000*(double)i/(double)n;
           Yeq[i] = func(Xeq[i]);
           Xcheb[i] = -cos(PI*(double)i/(double)n);
           Ycheb[i] = func(Xcheb[i]);
        }
        
        //gerar dominio com 10117 pontos e ao mesmo testar
        x[0] = -1.0000;
        max_eq = abs(func(x[0]) - Pol(Xeq,Yeq,x[0],n+1));
        max_cheb = abs(func(x[0]) - Pol(Xcheb,Ycheb,x[0],n+1));
        for(i=1;i<10117;i++){
           x[i] = -1.0000 + (2.0000*((double)i))/(10116.0000);
           temp = fabs(func(x[i]) - Pol(Xeq,Yeq,x[i],n+1));
           if(temp > max_eq) {
               max_eq = temp;
           }
           
           temp = fabs(func(x[i]) - Pol(Xcheb,Ycheb,x[i],n+1));
           if(temp > max_cheb) {
               max_cheb = temp;
           }
        }
        printf("%d\t%.18lf\t%.18lf\n",n,log10(max_eq),log10(max_cheb));
    }
    return 0;
}
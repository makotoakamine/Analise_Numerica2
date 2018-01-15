#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "BasicFunc.h"
#include "Gauss.h"
#define PI 3.14159265358979323846

double* findM(double* Y, double x0, double h, double d0, double dn, int n);
int findFloorIntervIndex(double* X,int n,double x);
double func(double x){
    return 1.000/(1.0000+25.00*x*x);
}


int main(){
    int i;
    int j;
    int n = 10;
    double x0;
    double h;
    double h2;
    double d0;
    double dn;
    double Aj;
    double Bj;
    double* Y;
    double* X;
    double* X_aprox;
    double* Y_aprox;
    double* Y_real;
    double* VV;
    double* M;
    
    while(n<=200){
        x0 = -1.00;
        h = 2.00/(double)n;
        d0 = 0.00;//spline cúbico natural y''(0) = 0
        dn = 0.00;//spline cúbico natural y''(n) = 0
        M = allocV(n+1);
        
        Y = allocV(n+1);
        X = allocV(n+1);
        for(i=0;i<n+1;i++){
            X[i] = x0 + i*h; 
            Y[i] = func(X[i]);
        }
        M = findM(Y,x0,h,d0,dn,n);
        //printV(M,n+1); // Se for necessário visualizar o vetor M
        //intervalo de análise = [findFloorIntervIndex(x0,h,0.95); findFloorIntervIndex(x0,h,0.95) + 1]
        
        Y_aprox = allocV(10117);
        X_aprox = allocV(10117);
        Y_real = allocV(10117);
        h2 = 2.00/10116.00;
        for(i=0;i<10117;i++){
            X_aprox[i] = x0 + i*h2;
            j = findFloorIntervIndex(X,n+1,X_aprox[i]);
            Aj = (Y[j+1]-Y[j])/h - (M[j+1]-M[j])*(h/6.00);
            Bj = Y[j] - (M[j]*h*h)/6.00;
            Y_aprox[i] = (M[j]/(6.00*h))*pow(X[j+1]-X_aprox[i],3) + (M[j+1]/(6*h))*pow(X_aprox[i]-X[j],3) + Aj*(X_aprox[i]-X[j]) + Bj;
            Y_real[i] = func(X_aprox[i]);
        }
        VV = allocV(n+1);
        subVV(VV,Y_real,Y_aprox,n+1);
        printf("%.14lf;%.14lf\n",log10(n),log10(Ninf(VV,n+1)));
        //printf("%d\t%.14lf\n",n,Ninf(VV,n+1));
        //Free para o ponteiros 
        free(Y);
        free(X);
        free(X_aprox);
        free(Y_aprox);
        free(Y_real);
        free(VV);
        free(M);
        n+=10;
    }
    
    
        return 0;
}

double* findM(double* Y, double x0, double h, double d0, double dn, int n){
    //double h = 2.00/(double)n; //nos igualmente espacados
    double *M = allocV(n+1);
    double *d = allocV(n+1);
    double **A = allocM(n+1);
    int i=0;
    int j=0;
    //y''(x) = (50(75x²-1))/((25x²+1)³)
    d[0] = d0;
    d[n] = dn;
    
    for(i = 1;i<n; i++){
        d[i] = (3.00/h)*((Y[i+1] - 2*Y[i] + Y[i-1])/h);
    }
    for(i=0;i<n+1;i++){
        for(j=0;j<n+1;j++){
           if(i==j){
               A[i][i] = 2.00;
           }
           else if(j == i+1){
               A[i][j] = 0.50;
           }
           else if(i == j+1){
               A[i][j] = 0.50;
           }
           else{
               A[i][j] = 0.00;
           }
        }
    }
    A[0][1] = 0.00;
    A[n][n-1] = 0.00;
    //Eliminacao de gauss
    rowEl(A,d,n+1); 
    M = solveMXV(A,d,n+1);
    return M;
}

int findFloorIntervIndex(double* X,int n,double x){
    int i=0;
    int mid;
    if(x == X[n-1]){
        return n-2;
    }
    else{
        while(n-i>1){
            mid = fmax((n+i)/2,i);
            if(x>=X[mid]){
                i = mid;
            }
            else{
                n = mid;
            }
        }
    }
    return i;
}
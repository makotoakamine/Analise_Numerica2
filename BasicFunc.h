#ifndef BASICFUNC_H
#define BASICFUNC_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double** allocM(int n); // Aloca memoria para Matriz de dimensao n, nao inicializa
double* allocV(int n); // Aloca memoria para Matriz de dimensao n, nao inicializa

void fillM(double **M, int n); // Inicializa valores da Matriz M
void fillV(double *V, int n); // Inicializa valores da Matriz M
void zeroFyM(double **M, int n); // Inicializa os valores da Matriz com 0
void zeroFyV(double *V, int n); // Inicializa os valores do Vetor com 0
    
    
void printM(double **M, int n); // printa a Matriz M
void printV(double *V, int n); // print o vetor V


double** scalarMultM(double **lambdaM, double lambda, double **M, int n); // Multiplicacao Matriz por Escalar
double* scalarMultV(double *lambdaV, double lambda, double *V, int n); // Multiplicacao Vetor por Escalar
double* scalarMultMV(double *MV, double **M, double *V, int n); // Multiplicacao Matriz x Vetor

double** matrixMult(double **MN, double **M, double **N, int n); //Multiplicacao entre Matrizes nxn

double* sumVV(double *VV, double *V_1, double *V_2, int n); //Soma 2 vetores de tamanho n
double* subVV(double *VV, double *V_1, double *V_2, int n); //Subtrai V_2 de V_1

double** upperTriangularFy(double **M, int n); // Retorna uma Matriz triangular superior relacionado a Matriz M

void stripLU(double **M, double **L, double **R, int n); //Separa a matrix em L+U = M
void stripDR(double **M, double **D, double **R, int n); //Separa a diagonal da Matriz M na Matrix D e o restante na Matrix R
void sparse_stripDR(double **M, double *D, double **R, int n); //Separa a diagonal da Matriz M no vetor D e o restante na Matrix R
void stripLDU(double **M, double **L, double **D, double **U, int n); //Separa a diagonal da Matriz M na patrix D e o restante na Matrix R
double** diagInv(double **D, int n); // retorna a inversa da Matriz diagonal D

double Ninf(double *V, int n); // retorna a norma do vetor
double** allocM(int n){
    double **M;
    int i = 0;
    M = (double **)malloc(n*sizeof(double *));
    for(i=0;i<n;i++){
        M[i] = (double *)malloc(n*sizeof(double));
    }
    return M;
}

double* allocV(int n){
    double *v;
    int i = 0;
    v = (double *)malloc(n*sizeof(double));
    return v;
}

void fillM(double **M, int n){
    int i = 0;
    int j = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            scanf("%lf",&M[i][j]);
        }
    }
}

void fillV(double *V, int n){
    int i = 0;
    for(i=0;i<n;i++){
            scanf("%lf",&V[i]);
    }
}

void zeroFyM(double **M, int n){
    int i = 0;
    int j = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            M[i][j] = 0;
        }
    }
}


void zeroFyV(double *V, int n){
    int i = 0;
    for(i=0;i<n;i++){
            V[i] = 0;
        }
}

void identityFyM(double **M, int n){
    int i = 0;
    int j = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                M[i][i] = 1;
            }
            else{
                M[i][j] = 0;
            }
        }
    }
}

void printM(double **M, int n){
    int i = 0;
    int j = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%.10lf\t",M[i][j]);
        }
        printf("\n");
    }
}

void printV(double *V, int n){
    int i = 0;
    for(i=0;i<n;i++){
            printf("%.10lf\t",V[i]);
    }
    printf("\n");
}

double** scalarMultM(double **lambdaM, double lambda, double **M, int n){
   int i = 0;
   int j = 0;
   for(i=0;i<n;i++){
       for(j=0;j<n;j++){
           lambdaM[i][j] = lambda*M[i][j];
       }
   }
   return lambdaM;
}

double* scalarMultV(double *lambdaV, double lambda, double *V, int n){
    int i = 0;
    for(i=0;i<n;i++){
        lambdaV[i] = lambda*V[i];
    }
    return lambdaV;
}

double* scalarMultMV(double *MV, double **M, double *V, int n){
    int i;
    int j;
    for(i=0;i<n;i++){
        MV[i] = 0;
        for(j=0;j<n;j++){
            MV[i] += M[i][j] * V[j];
        }
    }
    return MV;
}
double** matrixMult(double **MN, double **M, double **N, int n){
    int i = 0;
    int j = 0;
    int k = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            MN[i][j] = 0;
            for(k=0;k<n;k++){
                MN[i][j] += M[i][k]*N[k][j];
            }
        }
    }
    return MN;
}

double* sumVV(double *VV, double *V_1, double *V_2, int n){
    int i;
    for(i=0;i<n;i++){
        VV[i] = V_1[i]+V_2[i];
    }
    return VV;
}

double* subVV(double *VV, double *V_1, double *V_2, int n){
    int i;
    for(i=0;i<n;i++){
        VV[i] = V_1[i]-V_2[i];
    }
    return VV;
}

double** upperTriangularFy(double **M, int n){
    int i = 0;
    int j = 0;
    int k = 0;
    double temp = 0.00000;
    for(i=0;i<n-1;i++){
        if(M[i][i] == 0){ //Se pivo = 0
            j=i+1;
            while(M[j][i] == 0 && j<n){
                j++;
            }
            if(j==n){
                printf("Erro , nao e' possivel triangularizar a matriz.\n");
            }
            else{ //Trocar linha 
                for(k=i;k<n;k++){
                    temp = M[i][k];
                    M[i][k] = M[j][k];
                    M[j][k] = temp;
                }
            }
        }
        //--------------------------------
        for(j=i+1;j<n;j++){
           if(M[j][i] != 0){
               temp = -(M[j][i]/M[i][i]);
               M[j][i] = 0;
               for(k=1+i;k<n;k++){
                   M[j][k] += temp*M[i][k];
               }
           }
        }
    }
    return M;
}

void stripLU(double **M, double **L, double **U, int n){
    int i;
    int j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i>=j){
                L[i][j] = M[i][j];
                U[i][j] = 0;
            }
            else{
                U[i][j] = M[i][j];
                L[i][j] = 0;
            }
        }
    }
    
}

void stripLDU(double **M, double **L, double **D, double **U, int n){
    int i;
    int j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                D[i][j] = M[i][j];
                U[i][j] = 0;
                L[i][j] = 0;
            }
            else if(i>j){
                D[i][j] = 0;
                L[i][j] = M[i][j];
                U[i][j] = 0;
            }
            else{
                D[i][j] = 0;
                U[i][j] = M[i][j];
                L[i][j] = 0;
            }
        }
    }
    
}

void stripDR(double **M, double **D, double **R, int n){
    int i;
    int j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                D[i][j] = M[i][j];
                R[i][j] = 0;
            }
            else{
                D[i][j] = 0;
                R[i][j] = M[i][j];
            }
        }
    }
    
}

void sparse_stripDR(double **M, double *D, double **R, int n){
    int i;
    int j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                D[j] = M[i][j];
                R[i][j] = 0;
            }
            else{
                R[i][j] = M[i][j];
            }
        }
    }
    
}

double** diagInv(double **D, int n){
    int i;
    int j;
    double **invD;
    invD = allocM(n);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                invD[i][i] = 1/D[i][i];
            }
            else{
                invD[i][j] = 0;
            }
        }
    }
    return invD;
}


double Ninf(double *V, int n){
    int i;
    int max=0;
    for(i=1;i<n;i++){
        if(fabs(V[i])>fabs(V[max])) max = i;
    }
    return fabs(V[max]);
}
#endif
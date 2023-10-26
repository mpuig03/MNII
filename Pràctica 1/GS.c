/*MARC PUIG CREIXELL*/
/*Aquesta funció serveix per aplicar el mètode de Gauss-Seidel per a resolució
de sistemes en el cas d'una matriu tridiagonal de diagional A[i][i]=[a1,...,an] i 
diagonals secundaries fixes A[i][i-1]=A[i][i+1]=b i terme independent v=[v1,...,vn]*/

#include<stdio.h>
#include<stdlib.h>
int GS(int n, double tol, double *a, double b, double *v, double *x){
    double d0, d1=0,e,*x0;
    int i;

    x0=(double*)malloc(n*sizeof(double));
    if(x0==NULL){
        printf("Error de reserva de memòria.\n");
        return -1;
    }
    do{
        for(i=0;i<n;i++){
            x0[i]=x[i];
        }

        /*Iteració Gauss Seidel*/
        /*Tenim que per components, es simplifiquen molts calculs per la forma de la matriu,
        així doncs, ens queda el següent*/

        x[0]=(1/a[0])*(v[0]-b*x[1]);

        for(i=1;i<n-1;i++){
            x[i]=(1/a[i])*(v[i]-b*x[i-1]-b*x0[i+1]);
        }

        x[n-1]=(1/a[n-1])*(v[n-1]-b*x[n-2]);

        /*Càlcul aproximació fita superior error*/
        d0=d1;
        d1=0;
        for(i=0;i<n;i++){
            d1+=(x[i]-x0[i])*(x[i]-x0[i]);
        }
        d1=sqrt(d1);
        e=d1*d1/(d0-d1);

    }while(fabs(e)>=tol); 

    free(x0);
    return 0;
}

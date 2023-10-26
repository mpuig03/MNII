/*MARC PUIG CREIXELL.*/
/*Funció principal. Reserva la memòria dinàmica, assigna els valors de la matriu i crida a les diferents funcions necessàries.*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "potencia.c"
#include "potencia_inversa.c"

int main(void){
    double L=5, V0=5, dx;
    double tol=1e-12;
    double *V, b;
    double *z, *aux, *y, smax, smin;
    int n=50, i;

    /*Reserva de memòria dinàmica*/
    V=(double*)malloc(n*sizeof(double));
    z=(double*)malloc(n*sizeof(double));
    aux=(double*)malloc(n*sizeof(double));
    y=(double*)malloc(n*sizeof(double));
    if(V==NULL||z==NULL||aux==NULL||y==NULL){
        printf("Error de reserva de memòria.\n");
        return -1;
    }

    /*ASSIGNACIÓ DE VALORS A LA MATRIU*/
    dx=2*L/n;
    /*Posem V a la diagonal, (Vi + 1/dx²)*/
    for(i=0;i<n;i++){
        V[i]=0.5*(-L+(i+1)*dx)*(-L+(i+1)*dx)+V0; /*Vi*/
        V[i]+=1/(dx*dx); /*+1/dx²*/
    }
    /*Posem el valor de les diagonals secundàries*/
    b=-1/(2*dx*dx);

    /*MÈTODE DE LA POTÈNCIA*/
    z[0]=1;
    for(i=1;i<n;i++){
        z[i]=0;
    }
    smax=potencia(n,tol,V,b,z,y,aux);
    printf("VAP màxim: %22.15e\n",smax);
    /*MÈTODE DE LA POTÈNCIA INVERSA*/
    z[0]=1;
    for(i=1;i<n;i++){
        z[i]=0;
    }
    smin=inv_pot(n,tol,V,b,z,y,aux);
    printf("VAP mínim: %22.15e\n",smin);


    free(z);
    free(V);
    free(y);
    free(aux);
    return 0;  
}

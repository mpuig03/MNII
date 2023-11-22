/*MARC PUIG CREIXELL.   NIUB:20576474*/
/*Funció principal. Reserva la memòria dinàmica, assigna els valors de la matriu i 
crida a les diferents funcions necessàries.*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double norma(int n, double *v);
void multmat(int n, double *a, double b,double *v1, double *v2);
int GS(int n, double tol, double *a, double b, double *v, double *x);
int potencia(int n, double tol, double *a, double b, double *z, double *y);
int inv_pot(int n, double tol, double *a, double b, double *z, double *y);

#define itermax 500000

int main(void){
    double L=5.000, V0=5.000, tol=1e-08;
    double dx;
    double *V, b;
    double *z, *y;
    int n=50, i;


    /*Reserva de memòria dinàmica*/
    V=(double*)malloc(n*sizeof(double));
    z=(double*)malloc(n*sizeof(double));
    y=(double*)malloc(n*sizeof(double));

    if(V==NULL||z==NULL||y==NULL){
        printf("Error de reserva de memòria.\n");
        return -1;
    }

    /*ASSIGNACIÓ DE VALORS A LA MATRIU*/
    dx=2*L/(n+1);
    printf("Step dx: %f\n", dx);

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
    if(potencia(n,tol,V,b,z,y)!=0){
        return -1;
    };
    /*MÈTODE DE LA POTÈNCIA INVERSA*/
    z[0]=1;
    for(i=1;i<n;i++){
        z[i]=0;
    }
    if(inv_pot(n,tol,V,b,z,y)!=0){
        return -1;
    };

    free(z);
    free(V);
    free(y);
    return 0;  
}
/*norma*/
double norm(int n,double *v){
    int i;
    double norm;
    norm=0;
    for(i=0;i<n;i++){
        norm+=v[i]*v[i];
    }
    return sqrt(norm);
}
/*MULTIPLICACIÓ DE MATRIUS*/
void multmat(int n, double *a, double b,double *v1, double *v2){
    int i;
    v2[0]=a[0]*v1[0]+b*v1[1];
    v2[n-1]=b*v1[n-2]+a[n-1]*v1[n-1];

    for(i=1;i<n-1;i++){
        v2[i]=a[i]*v1[i]+b*(v1[i-1]+v1[i+1]);
    }
}

/*GAUSS-SEIDEL*/
/*Aquesta funció calcula la solució d'un sistema amb matriu A com l'enunciat 
i terme independent v*/
int GS(int n, double tol, double *a, double b, double *v, double *x){
    double e, aux;
    int i, iter=0;

    do{
        e=0;
        /*Iteració Gauss Seidel*/
        /*Tenim que per components, es simplifiquen molts calculs per la forma de la matriu,
        així doncs, ens queda el següent*/

        aux=x[0];
        x[0]=(1/a[0])*(v[0]-b*x[1]);
        e+=(x[0]-aux)*(x[0]-aux);

        for(i=1;i<n-1;i++){
            aux=x[i];
            x[i]=(1/a[i])*(v[i]-b*x[i-1]-b*x[i+1]);
            e+=(x[i]-aux)*(x[i]-aux);
        }

        aux=x[n-1];
        x[n-1]=(1/a[n-1])*(v[n-1]-b*x[n-2]);
        e+=(x[n-1]-aux)*(x[n-1]-aux);
        e=sqrt(e);
        
        iter++;
    }while(e>tol && iter<itermax); 

    if(iter>=itermax){
        printf("El metode de Gauss-Seidel no convergeix.\n");
        return -1;
    }
    return 0;
}

/*MÈTODE DE LA POTÈNCIA*/
int potencia(int n, double tol, double *a, double b, double *z, double *y)
{
    double s, e, ny, aux;
    int i, iter = 0;
    do{
        /*Càlcul de y i z*/
        multmat(n, a, b, z, y);

        ny=norm(n,y);

        for (i = 0; i < n; i++){
            z[i] = y[i] / ny;
        }

        /*Càlcul de sigma i l'error*/
        /*Pel càlcul de s (sigma): fem primer y=Az i llavors s=<z,y>*/
        multmat(n, a, b, z, y);
        aux=s;
        s = 0;
        for (i = 0; i < n; i++){
            s += y[i] * z[i];
        }
        e=fabs(s-aux);
        iter++;
    } while (e >= tol && iter <= itermax);

    if (iter >= itermax)
    {
        printf("El mètode de la potència no convergeix.\n");
        return -1;
    }
  

    printf("Power method: lambda_max= %.8lf after %d iterations.\n", s, iter);
    return 0;
}

/*MÈTODE DE LA POTÈNCIA INVERSA*/
int inv_pot(int n, double tol, double *a, double b, double *z, double *y){
    double s, e, ny, aux;
    int i, iter=0;

    do{
        /*Càlcul de y*/
        if(GS(n,tol,a,b,z,y)!=0){
            printf("No s'ha pogut realitzar el mètode de la potència inversa.\n");
            return -1;
        }

        /*Càlcul de la norma de y i de z*/
        ny=norm(n,y);

        for(i=0;i<n;i++){
            z[i]=y[i]/ny;
        }


        /*Càlcul de sigma (s) i l'error*/
        /*Pel càlcul de s (sigma): fem primer y=Az i llavors s=<z,y>*/
        multmat(n,a,b,z,y);
        aux=s;
        s=0;
        for(i=0;i<n;i++){
            s+=y[i]*z[i];
        }
        e=fabs(s-aux);
        iter++;
    }while(fabs(e)>=tol && iter<=itermax);

    if(iter>=itermax){
        printf("El mètode de la potència inversa no convergeix.\n");
        return -1;
    }
    
    printf("Inverse power method: lambda_min= %.8lf after %d iterations.\n", s, iter);

    return 0;
}
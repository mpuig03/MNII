/*MARC PUIG CREIXELL*/
/*Aquesta funció aplica el mètode la potència inversa i retorna s, el VAP de
mòdul mínim de la matriu*/
#include "GS.c"

double inv_pot(int n, double tol, double *a, double b, double *z, double *y, double *aux){
    double s, d0, d1, e, ny;
    int i;

    do{
        for(i=0;i<n;i++){
            aux[i]=z[i];
        }
        /*Càlcul de y*/
        if(GS(n,tol,a,b,z,y)==-1){
            return -1;
        }

        /*Càlcul de la norma de y i de z*/
        ny=0;
        for(i=0;i<n;i++){
            ny+=y[i]*y[i];
        }
        ny=sqrt(ny);

        for(i=0;i<n;i++){
            z[i]=y[i]/ny;
        }


        /*Càlcul aproximació fita superior de l'error*/
        d0=d1;
        d1=0;
        for(i=0;i<n;i++){
            d1+=(z[i]-aux[i])*(z[i]-aux[i]);
        }
        d1=sqrt(d1);
        e=d1*d1/(d0-d1);
    }while(fabs(e)>=tol);

    /*Càlcul de s (sigma): fem primer y=Az i llavors s=<z,y>*/
        multmat(n,a,b,z,y);
        s=0;
        for(i=0;i<n;i++){
            s+=y[i]*z[i];
        }
    return s;
}

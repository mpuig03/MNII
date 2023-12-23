#include<stdio.h>
#include<math.h>
double f(double *p);
double fx(double *p);
double fy(double *p);
double newt(double y0, double imax, double prec, double tol);
int pred(double *p, double *v, double h, double tol);
int newt2 (double *p, double *pi, double h, int imax, double prec, double tol);

int main(void){
    FILE *f;
    int n=1000, imax=10, i;
    double prec=1e-8, tol=1e-12, h=0.005;
    double dist1, dist2, y0;
    double v1[2],v2[2], p[2], q[2], pi[2], qi[2];

    /*Obertura fitxer punts*/
    f=fopen("res.txt","w");

    if(f==NULL){
        printf("Error en l'obertura del fitxer.\n");
        return 1;
    }
    /*Càlcul punt inicial*/
    p[0]=0;
    q[0]=0;
    y0=newt(0,imax,prec,tol);
    p[1]=y0;
    q[1]=y0;
    fprintf(f,"%le %le\n",p[0],p[1]);
  
    /*Iteració en ambdós sentits, v1(+), v2(-). 
    Els inicialitzem a (1,0) i (-1,0) respectivament, ja que fy!=0, 
    i aleshores té component horitzontal i el producte amb aquests dos és no nul.*/
    v1[0]=1;
    v1[1]=0;
    v2[0]=-1;
    v2[1]=0;
    dist1=1;
    dist2=1;

    for(i=0;i<=n && dist1>h && dist2>h;i++){
        pi[0]=p[0];
        pi[1]=p[1];
        qi[0]=q[0];
        qi[1]=q[1];

        if(pred(p,v1,h,tol)!=0 || pred(q,v2,h,tol)!=0){
            printf("No s'ha pogut realitzar el programa.\n");
            return 1;
        };

        if(newt2(p,pi,h,imax,prec,tol)!=0 || newt2(q,qi,h,imax,prec,tol)!=0){
            printf("No s'ha pogut realitzar el programa.\n");
            return 1;
        }
        dist1=sqrt(p[0]*p[0]+(p[1]-y0)*(p[1]-y0));
        dist2=sqrt(q[0]*q[0]+(q[1]-y0)*(q[1]-y0));
        fprintf(f,"%le %le\n",p[0],p[1]);
        fprintf(f,"%le %le\n",q[0],q[1]);
    }

    return 0;
}

/*FUNCIONS AVALUAR: Avalua respectivament f, la seva parcial en x (fx) i la seva parcial en y (fy).*/
double f(double *p){
    double x=p[0], y=p[1];
    return x*x*x*x-x*x*x*y*y-4*x*x*x-x*x*y+6*x*x-3*x*y*y*y*y+x*y+exp(x*y)+x*exp(y)-4*x+2*y*y*y*y*y+y*y*y*y-8*y*y*y+24*y*y-32*y+9;
}

double fx(double *p){
    double x=p[0], y=p[1];
    return 4*x*x*x-3*x*x*y*y-12*x*x-2*x*y+12*x-3*y*y*y*y+y+y*exp(x*y)+exp(y)-4;
}

double fy(double *p){
    double x=p[0], y=p[1];
    return -2*x*x*x*y-x*x-12*x*y*y*y+x+x*exp(x*y)+x*exp(y)+10*y*y*y*y+4*y*y*y-24*y*y+48*y-32;
}

/*NEWTON RAPHSON D'UNA VARIABLE: Aplica el mètode de Newton-Raphson imposant x=0 
i en retorna l'arrel.*/
double newt(double y0, double imax, double prec, double tol){
/*Es calcula imposant x=0, i tractant-ho com una funció f:R-->R amb derivada fy*/
    int iter=0;
    double p[2];
    p[0]=0;
    p[1]=y0;
    while(fabs(f(p)>prec) && fabs(fy(p))>tol && iter<imax){
        y0=y0-f(p)/fy(p);
        p[1]=y0;
        iter++;
    }
    if(iter>=imax){
        printf("Iteracions màximes de Newton-Raphson (1-D).\n");
        return 0;
    }else if(fabs(fy(p))<=tol){
        printf("Divisió per zero. (NR-1D)\n");
        return -1;
    }
    return p[1];
}

/*FUNCIÓ PREDICTOR: Calcula el vector tangent a la corba en la direcció indicada, 
així com el nou punt de la successió. Retorna el punt i el vector a p i v, respectivament
i retorna 0 si s'ha pogut realitzar correctament, -1 altrament.*/

int pred(double *p, double *v, double h, double tol){
    /*Es pot veure que el vector (fy,-fx) és tangent a la corba en cada punt*/
    /*Es calcula la direcció v*/
    double px,py,norm;
    px=fx(p);
    py=fy(p);
    
    norm=sqrt(py*py+px*px);
    if(fabs(norm)<tol){
        printf("Divisió entre zero.(PRED)\n");
        return -1;
    }

    /*Es comprova que la direcció és la correcta i s'assigna la nova v*/
    if(py*v[0]-px*v[1]>0){    
        v[0]=py/norm;
        v[1]=-px/norm;
    }else{
        v[0]=-py/norm;
        v[1]=px/norm; 
    }

    /*Es calcula un punt a distància h de p en direcció v*/
    p[0]=p[0]+h*v[0];
    p[1]=p[1]+h*v[1];

    return 0;   
} 

/*NEWTON RAPHSON EN DUES VARIABLES: Passant el punt anterior pi, i un punt inicial de prova p, 
resol el sistema per trobar pi+1. Retorna 0 si tot ha anat bé, -1 altrament.*/
int newt2 (double *p, double *pi, double h, int imax, double prec, double tol){
    int iter=0;
    double DF[2][2],F[2], det;
    
    /*Definim F=(f(x,y), (x-xi)²+(y-yi)²-h²) i fem els càlculs necessaris*/
    /*DF=(fx fy // 2(x-xi)  2(y-yi))*/
    /*A⁻¹=1/|A|*(d -b//-c a)*/

    DF[0][0]=fx(p);
    DF[0][1]=fy(p);
    DF[1][0]=2*(p[0]-pi[0]);
    DF[1][1]=2*(p[1]-pi[1]);
    F[0]=f(p);
    F[1]=(p[0]-pi[0])*(p[0]-pi[0])+(p[1]-pi[1])*(p[1]-pi[1])-h*h;
    det=DF[0][0]*DF[1][1]-DF[0][1]*DF[1][0];

    while((fabs(F[0])>prec || fabs(F[1])>prec) && iter<imax && fabs(det)>tol){
        /*vector iteració*/
        p[0]=p[0]-(DF[1][1]*F[0]-DF[0][1]*F[1])/det;
        p[1]=p[1]-(-DF[1][0]*F[0]+DF[0][0]*F[1])/det;
        /*actualització iteració*/
        DF[0][0]=fx(p);
        DF[0][1]=fy(p);
        DF[1][0]=2*(p[0]-pi[0]);
        DF[1][1]=2*(p[1]-pi[1]);
        F[0]=f(p);
        F[1]=(p[0]-pi[0])*(p[0]-pi[0])+(p[1]-pi[1])*(p[1]-pi[1])-h*h;
        det=DF[0][0]*DF[1][1]-DF[0][1]*DF[1][0];
        iter++;
        
    }
    if(iter>=imax){
        printf("Iteracions màximes de Newton-Raphson (2-D).\n");
        return -1;
    }else if(fabs(det)<=tol){
        printf("Divisió per zero.(NR-2D)\n");
        return -1;
    }
    return 0;
}
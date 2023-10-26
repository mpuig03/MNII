/*MARC PUIG CREIXELL*/
/*Funció que multiplica una matriu tridiagonal de diagional A[i][i]=[a1,...,an] i 
diagonals secundaries fixes A[i][i-1]=A[i][i+1]=b per un vector v1=[v1,...,vn] 
i posa el resultat en un vector v2=[v1',...,vn']*/

void multmat(int n, double *a, double b,double *v1, double *v2){
    int i;
    v2[0]=a[0]*v1[0]+b*v1[1];
    v2[n-1]=b*v1[n-2]+a[n-1]*v1[n-1];

    for(i=1;i<n-1;i++){
        v2[i]=a[i]*v1[i]+b*(v1[i-1]+v1[i+1]);
    }
}

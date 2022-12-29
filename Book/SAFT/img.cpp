#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"saft.h"

IMG::IMG(){
	Xa[0]=0.0; Xa[1]=0.0;
	Wd[0]=1.0; Wd[1]=1.0;
};

void IMG::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void IMG::set_Wd(double W, double H){
	Wd[0]=W; Wd[1]=H;
	Xb[0]=Xa[0]+W;
	Xb[1]=Xa[1]+H;
};
void IMG::set_Ng(int nx, int ny){
	int i;
	Ng[0]=nx; Ng[1]=ny;
	for(i=0;i<2;i++){
		Ndiv[i]=Ng[i]-1;
		dx[i]=0.0;
		if(Ndiv[i]>0) dx[i]=Wd[i]/Ndiv[i];
	}
	IMG::mem_alloc();
	npnt=Ng[0]*Ng[1];

	for(i=0;i<Ng[0];i++) xcod[i]=Xa[0]+i*dx[0];
	for(i=0;i<Ng[1];i++) ycod[i]=Xa[1]+i*dx[1];
};
void IMG::mem_alloc(){
	double *pt=(double *)malloc(sizeof(double)*Ng[0]*Ng[1]);
	V=(double **)malloc(sizeof(double *)*Ng[0]);

	xcod=(double *)malloc(sizeof(double)*Ng[0]);
	ycod=(double *)malloc(sizeof(double)*Ng[1]);

	int i;
	for(i=0;i<Ng[0]*Ng[1];i++) pt[i]=0.0;
	for(i=0;i<Ng[0];i++) V[i]=pt+Ng[1]*i;


};

double IMG::get_xcod(int l){
	int i=l/Ng[1];
	return(Xa[0]+i*dx[0]);
};
double IMG::get_ycod(int l){
	int j=l%Ng[1];
	return(Xa[1]+j*dx[1]);
};
void IMG::set_V(int l,double val){
	int i=l/Ng[1];
	int j=l%Ng[1];
	V[i][j]=val;
};
void IMG::fwrite_V(char fname[128]){
	FILE *fp=fopen(fname,"w");
	int i,j;
	fprintf(fp,"# Xa[0]  Xa[1]\n");
	fprintf(fp,"%lf %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# Xb[0]  Xb[1]\n");
	fprintf(fp,"%lf %lf\n",Xb[0],Xb[1]);
	fprintf(fp,"# Ng[0]  Ng[1]\n");
	fprintf(fp,"%d %d\n",Ng[0],Ng[1]);
	fprintf(fp,"# Values\n");
	for(i=0;i<Ng[0];i++){
	for(j=0;j<Ng[1];j++){
		fprintf(fp," %lf\n",V[i][j]);
	}
	}
	fclose(fp);
};


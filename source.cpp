#include<stdio.h>
#include<stdlib.h>
#include"wveq2d.h"

void SOURCE::print(double *Xa, double *dx){
	printf("----------- Source Parameters -----------\n");
	printf("ng=%d\n",ng);
	double x,y;
	printf("wv ID=%d\n",iwv);

	double xofst=0.5, yofst=0.5;
	if(type==1) xofst=0.0;
	if(type==2) yofst=0.0;

	for(int i=0;i<ng;i++){
		printf("(i,j)=(%d, %d) ",isrc[i],jsrc[i]);
		x=Xa[0]+(isrc[i]+xofst)*dx[0];
		y=Xa[1]+(jsrc[i]+yofst)*dx[1];
		printf("(x,y)=(%lf, %lf)\n",x,y);
	};
};
void SOURCE::set_inc_ang(double th, double ct){
	double PI=4.0*atan(1.0);
	th=th/180.*PI;	
	Ctd=sin(th)/ct*dx[type-1];
};
void SOURCE::mem_alloc(int n){
	ng=n;
	isrc=(int *)malloc(sizeof(int)*ng);
	jsrc=(int *)malloc(sizeof(int)*ng);
	ksrc=(int *)malloc(sizeof(int)*ng);
	for(int i=0;i<ng;i++){
		isrc[i]=-1;
		jsrc[i]=-1;
		ksrc[i]=-1;
	};
};
void TRNSDCR::init_bwv(int nn, double dtau){
	Nt=nn;
	dt=dtau;
	double *p=(double *)malloc(sizeof(double)*Nt*ng);
	bwv=(double **)malloc(sizeof(double *)*ng);
	int i;
	for(int i=0; i<Nt*ng; i++) p[i]=0.0;
	for(i=0;i<ng;i++) bwv[i]=p+i*Nt;
};
void TRNSDCR::record(int jt, double **fld){
	int ig,i,j;
	for(ig=0; ig<ng; ig++){
		i=isrc[ig];
		j=jsrc[ig];
		if(type==1) i-=nml;
		if(type==2) j-=nml;
		bwv[ig][jt]=fld[i][j]; 
	};
};
void TRNSDCR::fwrite(){
	int i,j;
	char fname[128];
	sprintf(fname,"tr%d.out",ID);
	FILE *fp=fopen(fname,"w");
	fprintf(fp,"# Nt, dt\n");
	fprintf(fp,"%d, %lf\n",Nt,dt);
	double x,y;
	for(i=0;i<ng;i++){
		x=Xa[0]+dx[0]*(isrc[i]+0.5);
		y=Xa[1]+dx[1]*(jsrc[i]+0.5);
		fprintf(fp,"# %lf, %lf\n",x,y);
	for(j=0;j<Nt;j++){
		fprintf(fp,"%lf\n",bwv[i][j]);
	};
	};
	fclose(fp);
};
void ARRAY::init(int nn){
	nele=nn;
	actv=(int *)malloc(sizeof(int)*nele);
	a0=(double *)malloc(sizeof(double)*nele);
	tdly=(double *)malloc(sizeof(double)*nele);

	for(int i=0;i<nele;i++){
		actv[i]=0;
		a0[i]=1.0;
		tdly[i]=0.0;
	};
};
void ARRAY::print(){
	for(int i=0;i<nele;i++){
		printf("i=%d, actv=%d, a0=%lf, tdly=%lf\n",i,actv[i],a0[i], tdly[i]);
	};
};

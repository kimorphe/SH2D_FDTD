#include<stdio.h>
#include<stdlib.h>
#include"wveq2d.h"
void RECVR::mem_alloc(int n){
	ng=n;
	irec=(int *)malloc(sizeof(int)*ng);
	jrec=(int *)malloc(sizeof(int)*ng);
	for(int i=0;i<ng;i++){
		irec[i]=-1;
		jrec[i]=-1;
	};
};
void RECVR::print(){
	printf("----------- Receiver Setting -----------\n");
	printf("ng=%d\n",ng);
	double x,y;
	printf("type=%d, nml=%d, wd=%lf\n",type,nml,wd);
	for(int i=0;i<ng;i++){
		printf("(i,j)=(%d,%d)\n",irec[i],jrec[i]);
	};
};

void RECVR::init_bwv(int nn, double dtau){
	Nt=nn;
	dt=dtau;
	double *p=(double *)malloc(sizeof(double)*Nt*ng);
	bwv=(double **)malloc(sizeof(double *)*ng);
	int i;
	for(int i=0; i<Nt*ng; i++) p[i]=0.0;
	for(i=0;i<ng;i++) bwv[i]=p+i*Nt;
};

void RECVR::record(int jt, double **fld){
	int i;
	for(i=0;i<ng;i++) bwv[i][jt]=fld[irec[i]][jrec[i]]; 
	//printf("irec,jrec=%d, %d  ",irec[i-1],jrec[i-1]); 
	//printf("bwv=%lf\n",fld[irec[i-1]][jrec[i-1]]); 
};
void RECVR::set_cod(double *xa, double *dh){
	Xa=xa; 
	dx=dh;
};
void RECVR::fwrite_dump(int n_meas){
	int i,j;
	char fname[128];
	sprintf(fname,"T%d/bwv%d.out",n_meas,ID);
	FILE *fp=fopen(fname,"w");
	if(fp==NULL){
		printf("Error Canrot open %s\n",fname);
		exit(-1);
	};
	fprintf(fp,"# Nt, dt\n");
	fprintf(fp,"%d, %lf\n",Nt,dt);
	double x,y;
	double iofst=0.0, jofst=0.0;
	if(type==1) jofst=0.5;
	if(type==2) iofst=0.5;

	for(i=0;i<ng;i++){
		x=Xa[0]+dx[0]*(irec[i]+iofst);
		y=Xa[1]+dx[1]*(jrec[i]+jofst);
		fprintf(fp,"# %lf, %lf\n",x,y);
	for(j=0;j<Nt;j++){
		fprintf(fp,"%lf\n",bwv[i][j]);
	};
	};
	fclose(fp);
};
void RECVR::clear(){
	int i,j;
	for(i=0;i<ng;i++){
	for(j=0;j<Nt;j++){
		bwv[i][j]=0.0;
	}
	}
};
void RECVR::set_inc_ang(double th, double ct){
	double PI=4.0*atan(1.0);
	th=th/180.*PI;	
	Ctd=sin(th)/ct*dx[type-1];
};
double RECVR::val(int ig, double tt){
	if(tt<0.0) return(0.0);
	int it=int(tt/dt);
	if(it > Nt-1) return(0.0);
	if(it ==Nt-1) return(bwv[ig][Nt-1]);
	double xi=tt/dt-it;
	return(bwv[ig][it]*(1.0-xi)+bwv[ig][it+1]*xi);
};
double RECVR::amp_synth(int it){
	int ig;
	double amp=0.0;
	double tt=it*dt;
	for(ig=0; ig<ng; ig++){
		amp+=val(ig,tt-Ctd*ig);
	};
	return(amp);
};
void RECVR::fwrite(int n_meas){
	char fname[128];
	sprintf(fname,"T%d/rec%d.out",n_meas,ID);
	FILE *fp=fopen(fname,"w");
	if(fp==NULL){
		printf("Can't open file %s\n",fname);
		printf(" First create directory T* of appropriate numbers\n");
		exit(-1);
	};

	int i;
	double amp;
	fprintf(fp,"# Nt, dt\n");
	fprintf(fp,"%d, %lf\n",Nt,dt);

	double xsrc,ysrc,wdx,wdy;
	xsrc=Xa[0]+(irec[0]+irec[ng-1]+0.5)*dx[0];
	ysrc=Xa[1]+(jrec[0]+jrec[ng-1]+0.5)*dx[1];
	wdx=(irec[ng-1]-irec[0])*dx[0];	// grid-to-grid distance (not a cell extent)
	wdy=(jrec[ng-1]-jrec[0])*dx[1];	// grid-to-grid distance (not a cell extent)

	fprintf(fp,"# xsrc,ysrc\n");
	fprintf(fp,"%lf, %lf\n",xsrc,ysrc);
	fprintf(fp,"# wdx, wdy (grid-to-grid distance)\n");
	fprintf(fp,"%lf, %lf\n",wdx,wdy);
	fprintf(fp,"# amp\n");
	for(i=0; i<Nt; i++){
		amp=RECVR::amp_synth(i);
		fprintf(fp,"%lf\n",amp);
	};
	fclose(fp);
};	

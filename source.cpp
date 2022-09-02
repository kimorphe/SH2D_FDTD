#include<stdio.h>
#include<stdlib.h>
#include"wveq2d.h"

/*
void SOURCE::rng(double *dx, double *Xa){

	switch(type){
	case 1:
		xs1[0]=Xa[0]+dx[0]*(isrc[0]);
		xs2[0]=Xa[0]+dx[0]*(isrc[ng-1]);
		xs1[1]=Xa[1]+dx[1]*(jsrc[0]+0.5);
		xs2[1]=Xa[1]+dx[1]*(jsrc[ng-1]+0.5);
		break;
	case 2:
		xs1[0]=Xa[0]+dx[0]*(isrc[0]+0.5);
		xs2[0]=Xa[0]+dx[0]*(isrc[ng-1]+0.5);
		xs1[1]=Xa[1]+dx[1]*(jsrc[0]);
		xs2[1]=Xa[1]+dx[1]*(jsrc[ng-1]);
		break;
	}

};
*/
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

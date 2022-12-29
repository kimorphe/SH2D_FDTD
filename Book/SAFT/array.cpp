#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "saft.h"


void show_msg(char *fn){
	printf("Cannot find file %s \n",fn);
	printf(" --> abort process ...\n");
	exit(-1);

};
void Array::setup(){
	char fname[128]="../Book/array.inp";
	FILE *fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);

	char cbff[128];
	int i,j; 

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nele);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nscan);

	printf("nele=%d, nscan=%d\n",nele,nscan);

	Array::mem_alloc();

	int k=0,is_src;
	for(i=0; i<nscan; i++){
		fgets(cbff,128,fp);
	for(j=0; j<nele; j++){
		fscanf(fp,"%d, %lf, %lf\n",&is_src,A0+k,delay+k);
		if(is_src){
			actv[k]=true;
			isrc[i]=j;
		}

		k++;
	};
		printf("jsrc=%d\n",isrc[i]);
	};

	fclose(fp);
};

void Array::mem_alloc(){
	xcod=(double *)malloc(sizeof(double)*nele);
	ycod=(double *)malloc(sizeof(double)*nele);

	actv=(bool *)malloc(sizeof(bool)*nele*nscan);
	delay=(double *)malloc(sizeof(double)*nele*nscan);
	A0=(double *)malloc(sizeof(double)*nele*nscan);

	isrc=(int *)malloc(sizeof(int)*nscan);

	int i,j,k;
	k=0;
	for(i=0;i<nele;i++){
		xcod[i]=0.0;
		ycod[i]=0.0;
		
		for(j=0;j<nscan;j++){
			delay[k]=0.0;
			A0[k]=1.0;
			actv[k]=false;
			k++;
		}
	};
};

void Array::load(){
	char fname[128]="../Book/tr_elems.out";
	FILE *fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);

	char cbff[128];
	int i,j,ng,isrc,jsrc;
	double xx,yy;
	for(i=0;i<nele;i++){
		fgets(cbff,128,fp);
		fgets(cbff,128,fp);
		fgets(cbff,128,fp);
		fgets(cbff,128,fp);
		fgets(cbff,128,fp);
		fgets(cbff,128,fp);
		fscanf(fp,"%d\n",&ng);
		fgets(cbff,128,fp);
		xcod[i]=0.0;
		ycod[i]=0.0;
		for(j=0; j<ng; j++){
			fscanf(fp,"%d, %d, %lf, %lf\n",&isrc,&jsrc,&xx,&yy);
			xcod[i]+=xx;
			ycod[i]+=yy;
		};
		xcod[i]/=ng;
		ycod[i]/=ng;
		printf("%d, (%lf, %lf)\n",i,xcod[i],ycod[i]);
	};

	fclose(fp);
};
double* Array::set_transmitter(int iscan){
	xin[0]=xcod[isrc[iscan]];
	xin[1]=ycod[isrc[iscan]];
	return(xin);

};

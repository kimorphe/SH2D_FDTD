#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "saft.h"


void show_msg(char *fn){
	printf("Cannot find file %s \n",fn);
	printf(" --> abort process ...\n");
	exit(-1);

};
void Array::setup(char path_name[128]){
	sprintf(path,"%s",path_name);
	char fname[256]; //="../Book/array.inp";

	sprintf(fname,"%s/array.inp",path);
	printf("%s\n",fname);
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
	thij=(double *)malloc(sizeof(double)*nscan*nele);

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
	char fname[256]; //="../Book/tr_elems.out";
	sprintf(fname,"%s/tr_elems.out",path);
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
		//printf("%d, (%lf, %lf)\n",i,xcod[i],ycod[i]);
	};

	fclose(fp);
};
double* Array::set_transmitter(int iscan){
	xin[0]=xcod[isrc[iscan]];
	xin[1]=ycod[isrc[iscan]];
	return(xin);
};

void Array::set_angles(double xc, double yc){
	int i,j,k;
	double tin[2],tsc[2];
	double rin,rsc;
	double thin,thsc;
	double PI=4.0*atan(1.0);
	double *xi,xr[2];
	k=0;
	for(i=0;i<nscan;i++){
		xi=Array::set_transmitter(i);
		tin[0]=xi[0]-xc;
		tin[1]=xi[1]-yc;
		rin=sqrt(tin[0]*tin[0]+tin[1]*tin[1]);	
		tin[0]/=rin;
		tin[1]/=rin;
		thin=acos(tin[0]);
		if(tin[1]<0.0) thin=2.*PI-thin; 
	for(j=0;j<nele;j++){
		xr[0]=xcod[j];
		xr[1]=ycod[j];
		tsc[0]=xr[0]-xc;
		tsc[1]=xr[1]-yc;
		rsc=sqrt(tsc[0]*tsc[0]+tsc[1]*tsc[1]);
		tsc[0]/=rsc;
		tsc[1]/=rsc;
		thsc=acos(tsc[0]);
		if(tsc[1]<0.0) thsc=2.*PI-thsc; 
		thij[k]=(thin+thsc)*0.5;
		k++;
	}
	}
	th_min=thij[0];
	th_max=thij[0];
	for(k=0;k<nscan*nele;k++){
		if(th_min>thij[k]) th_min=thij[k];
		if(th_max<thij[k]) th_max=thij[k];
	}
	printf("th_min=%lf\n",th_min/PI*180.0);
	printf("th_max=%lf\n",th_max/PI*180.0);

	int kbin;
	nbin=101;
	hist=(int *)malloc(sizeof(int)*nbin);
	for(k=0;k<nbin;k++) hist[k]=0;

	dth=(th_max-th_min)/(nbin-1);
	for(k=0;k<nscan*nele;k++){
		kbin=int((thij[k]-th_min)/dth+0.5);
		hist[kbin]++;
	};
	//for(k=0;k<nbin;k++) printf("%lf, %d\n",(th_min+dth*k)/PI*180.,hist[k]);
};
int Array::get_hist_val(int ktr){
	// ktr: (T,R) combination
	int kbin=int((thij[ktr]-th_min)/dth+0.5);	// bin number
	return(hist[kbin]);  // count 
};

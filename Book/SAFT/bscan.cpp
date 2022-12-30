#include<stdio.h>
#include<stdlib.h>
#include"saft.h"

Bscan::Bscan(char path_name[128]){
	t0=0.0;
	alloc=false;
	sprintf(path,"%s",path_name);
};
void Bscan::load(int tnum,int wv_typ){
	char fname[256];
	//sprintf(fname,"../Book/T%d/ary.out",tnum);
	if(wv_typ==0) sprintf(fname,"%s/T%d/ary.out",path,tnum);
	if(wv_typ==1) sprintf(fname,"%s/T%d/ary_cor.out",path,tnum);
	FILE *fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);

	char cbff[128];

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&Ny);

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %lf\n",&Nt,&dt);

	if(!alloc) Bscan::mem_alloc();
	Td=t0+(Nt-1)*dt;

	int i,j;
	double tmp;
	for(i=0; i<Ny;i++){
		fgets(cbff,128,fp);
	for(j=0; j<Nt;j++){
		//fscanf(fp,"%lf\n",amp[i]+j);
		fscanf(fp,"%lf\n",&tmp);
		amp[i][j]=tmp;
	}
	}

	fclose(fp);
};

void Bscan::mem_alloc(){
	double *pt=(double *)malloc(sizeof(double)*Nt*Ny);
	amp=(double **)malloc(sizeof(double *)*Ny);
	int i;
	for(i=0;i<Nt*Ny;i++) pt[i]=0.0;
	for(i=0;i<Ny;i++) amp[i]=pt+Nt*i;

	//amp_sum=(double *)malloc(sizeof(double)*Nt);
	//for(i=0;i<Nt;i++) amp_sum[i]=0.0;

	alloc=true;

};
void Bscan::fwrite_Bscan(){
	FILE *fp=fopen("bscan.out","w");
	int i,j;

	fprintf(fp,"#\n");
	fprintf(fp,"%d\n",Ny);
	fprintf(fp,"#\n");
	fprintf(fp,"%lf, %lf, %d\n",t0,dt,Nt);
	fprintf(fp,"#\n");
	for(j=0;j<Ny;j++){
	for(i=0;i<Nt;i++){
		fprintf(fp,"%lf\n",amp[j][i]);
	}
	}
	fclose(fp);
};

double Bscan::get_amp(int iele, double time){
	if(time<t0) return(0.0);
	if(time>=Td) return(0.0);

	int it=int((time-t0)/dt);
	double xi=(time-t0)/dt-it;

	double a1=amp[iele][it];
	double a2=amp[iele][it+1];
	return(a1*(1.-xi)+a2*xi);
};

/*
double Bscan::stack_Ascans(double *tofs){
	int i,j,jd;
	double idly;
	tof=tofs[0];
	ntof=int((tof-t0)/dt);

	for(j=0;j<Nt;j++) amp_sum[j]=0.0;

	for(i=0;i<Ny;i++){
		idly=int((tofs[i]-tof)/dt);
		for(j=0;j<Nt;j++){
			jd=j+idly;
			if(jd<0) continue;
			if(jd>Nt-1) continue;
			amp_sum[j]+=amp[i][jd];
		}
//		amp_sum[i]/=Ny;
	}

	double amax=0.0;
	for(j=0;j<Nt;j++){
	     	if(amax < abs(amp_sum[j])) amax=abs(amp_sum[j]);
	}
	return(amax);
	
};
double Bscan::get_Asum(int inc_ntof){
	int mtof=inc_ntof+ntof;
	if(mtof<0) return(0);
	if(mtof>Nt-1) return(0);
	return(amp_sum[mtof]);
};
void Bscan::fwrite_Ascan(char fname[128]){
	FILE *fp=fopen(fname,"w");
	int i;
	for(i=0;i<Nt;i++){
		fprintf(fp,"%lf %lf\n",t0+dt*i,amp_sum[i]);
	}
	fprintf(fp,"\n");
	fprintf(fp,"%lf %lf\n",tof,5.0);
	fprintf(fp,"%lf %lf\n",tof,-5.0);
	fclose(fp);
};
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"wveq2d.h"

//--------------------------------------------------------------------------
void FIELD::fwrite_trim(int d_num, int num, int *NHa, int *NHb, double tout){
	char fname[128];
	if(type==0) sprintf(fname,"T%d/v3_%d.out",d_num,num);
	if(type==1) sprintf(fname,"T%d/q1_%d.out",d_num,num);
	if(type==2) sprintf(fname,"T%d/q2_%d.out",d_num,num);
	FILE *fp=fopen(fname,"w");
	if(fp==NULL) show_msg(fname);
	int i,j;
	double xll[2], wdt[2];
	int ngs[2];
	for(i=0;i<2;i++){
		xll[i]=Xa[i]+NHa[i]*dx[i];
		wdt[i]=Wd[i]-(NHa[i]+NHb[i])*dx[i];
		ngs[i]=Ng[i]-NHa[i]-NHb[i];
	}

	fprintf(fp,"# time\n");
	fprintf(fp,"%lf\n",tout);
	fprintf(fp,"# Xa[0:1]\n");
	fprintf(fp,"%lf, %lf\n",xll[0],xll[1]);

	fprintf(fp,"# Xb[0:1]\n");
	fprintf(fp,"%lf, %lf\n",xll[0]+wdt[0],xll[1]+wdt[1]);

	fprintf(fp,"# Ng[0:1]\n");
	fprintf(fp,"%d, %d\n",ngs[0],ngs[1]);
	fprintf(fp,"# field value\n");
	//for(i=NHa[0]; i<ngs[0]; i++){
	//for(j=NHa[1]; j<ngs[1]; j++){
	for(i=0; i<ngs[0]; i++){
	for(j=0; j<ngs[1]; j++){
		fprintf(fp,"%lf\n",F[i+NHa[0]][j+NHa[1]]);
	}
	}
	//printf("Xa=%lf, %lf\n",Xa[0],Xa[1]);
	//jprintf("dx=%lf, %lf\n",dx[0],dx[1]);
	//printf("NHa=%d, %d\n",NHa[0],NHa[1]);
	//printf("NHb=%d, %d\n",NHb[0],NHb[1]);
	fclose(fp);
}
void FIELD::fwrite(int d_num, int num){
	char fname[128];
	if(type==0) sprintf(fname,"T%d/v3_%d.out",d_num,num);
	if(type==1) sprintf(fname,"T%d/q1_%d.out",d_num,num);
	if(type==2) sprintf(fname,"T%d/q2_%d.out",d_num,num);
	FILE *fp=fopen(fname,"w");
	if(fp==NULL) show_msg(fname);
	int i,j;

	fprintf(fp,"# Xa[0:1]\n");
	fprintf(fp,"%lf, %lf\n",Xa[0],Xa[1]);

	fprintf(fp,"# Xb[0:1]\n");
	fprintf(fp,"%lf, %lf\n",Xa[0]+Wd[0],Xa[1]+Wd[1]);

	fprintf(fp,"# Ng[0:1]\n");
	fprintf(fp,"%d, %d\n",Ng[0],Ng[1]);
	fprintf(fp,"# field value\n");
	for(i=0; i<Ng[0]; i++){
	for(j=0; j<Ng[1]; j++){
		fprintf(fp,"%lf\n",F[i][j]);
	}
	}
	fclose(fp);
};
void FIELD::print_prms(){
	printf("---------FIELD SETTING-----------\n");
	printf("Xa=%lf, %lf\n",Xa[0],Xa[1]);
	printf("Wd=%lf, %lf\n",Wd[0],Wd[1]);
	printf("dx=%lf, %lf\n",dx[0],dx[1]);
	printf("Ndiv=%d, %d\n",Ndiv[0],Ndiv[1]);
	printf("Ng=%d, %d\n",Ng[0],Ng[1]);
	printf("type=%d\n",type);
	printf("ofst[2]=(%lf, %lf)\n",ofst[0],ofst[1]);
	printf("---------------------------------\n");
};
void FIELD::fwrite_prms(char *fn, char *mode, char *name){
	FILE *fp=fopen(fn,mode);
	if(fp==NULL) show_msg(fn);

	fprintf(fp,"--------------%s--------------\n",name);
	fprintf(fp,"type=%s (=%d)\n",stype,type);
	fprintf(fp,"Ndiv=%d, %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"Ng=%d, %d\n",Ng[0],Ng[1]);
	fprintf(fp,"ofst[2]=(%lf, %lf)\n",ofst[0],ofst[1]);
	fprintf(fp,"Nin=%d, Nbnd=%d, Nex=%d\n",Nin,Nbnd,Nex);

	fclose(fp);
};
void FIELD::set_IC(double xc, double yc, double sig, double f0){
	int i,j;
	double *xcod;
	double rx,ry,arg,val;

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		xcod=ij2xf(i,j);
		rx=(xcod[0]-xc)/sig;
		ry=(xcod[1]-yc)/sig;
		arg=sqrt(rx*rx+ry*ry);
		F[i][j]=f0*exp(-0.5*arg);
	}
	}
};
double* FIELD::ij2xf(int ii, int jj){
	xf[0]=Xa[0]+dx[0]*(ii+ofst[0]);
	xf[1]=Xa[1]+dx[1]*(jj+ofst[1]);
	return(xf);
};
int FIELD::ij2l(int i, int j){
	return(j+Ng[1]*i);
};
void FIELD::l2ij(int l, int *ix, int *jy){
	*ix=l/Ng[1];
	*jy=l%Ng[1];
};
void FIELD::print_F(){
	for(int i=0; i<Ng[0]; i++){
	for(int j=0; j<Ng[1]; j++){
		printf("%lf, ",F[i][j]);
	}
		printf("\n");
	}
};
void FIELD::init(int *ndiv, int ityp){
	type=ityp;
	Ndiv=ndiv;
	Ng[0]=Ndiv[0];
	Ng[1]=Ndiv[1];
	sprintf(stype,"v3");

	ofst[0]=0.5;
	ofst[1]=0.5;
	if(type==1){	// q1-grid
		ofst[0]=0.0;
		Ng[0]+=1;
		sprintf(stype,"q1(=s31)");
	};
	if(type==2){	// q2-grid
		ofst[1]=0.0;
		Ng[1]+=1;
		sprintf(stype,"q2(=s32)");
	}
	if(type==3){	// vertices
		ofst[0]=0.0;
		ofst[1]=0.0;
		Ng[0]+=1;
		Ng[1]+=1;
		sprintf(stype,"vertice(not used)");
	};

	ndat=Ng[0]*Ng[1];
	double *p=(double *)malloc(sizeof(double)*ndat);
	F=(double **)malloc(sizeof(double *)*Ng[0]);
	int k=0;
	for(int i=0; i<Ng[0]; i++){
		F[i]=p+i*Ng[1];
	for(int j=0; j<Ng[1]; j++){
		p[k++]=0.0; 
	}
	}
};
void FIELD::setup(double *xa, double *wdt, double *dh){
	Xa=xa;
	Wd=wdt;
	dx=dh;
};

//		GENERATE 1D INDEX FOR q1 
void FIELD::gen_indx1(int **kcell){
	int i,j;
	int il,ir;

	for(int m=0;m<2;m++){	//count grids at m=0, results stored at m=1; 
		Nin=0;
		Nbnd=0;
		Nex=0;
	for(i=0;i<Ng[0];i++){
	for(j=0;j<Ng[1];j++){
		il=1; ir=1;
		if(i>0) il=kcell[i-1][j];	
		if(i<Ndiv[0]) ir=kcell[i][j];	
		switch(il+ir){
		case 0:	// interior grid 
			if(m==1) kint[Nin]=i*Ng[1]+j;
			Nin++;
			break;
		case 1: // boundary grid
			if(m==1){
				kbnd[Nbnd]=i*Ng[1]+j;
				if(il==1) kbnd[Nbnd]*=-1;
			}
			Nbnd++; 
			break;
		case 2: // exterior grid
			Nex++;
			break;
		};
	}
	}

		if(m==0){
			kint=(int *)malloc(sizeof(int)*Nin);
			kbnd=(int *)malloc(sizeof(int)*Nbnd);
		}
	}
};
//		GENERATE 1D INDEX FOR q2 
void FIELD::gen_indx2(int **kcell){
	int i,j;
	int iu,id;

	for(int m=0;m<2;m++){
		Nin=0;
		Nbnd=0;
		Nex=0;
	for(i=0;i<Ng[0];i++){
	for(j=0;j<Ng[1];j++){
		id=1; iu=1;
		if(j>0) id=kcell[i][j-1];	
		if(j<Ndiv[1]) iu=kcell[i][j];	
		switch(iu+id){
		case 0: // interior grid
			if(m==1) kint[Nin]=i*Ng[1]+j;
			Nin++;
			break;
		case 1: // boundary grid
			if(m==1){
				kbnd[Nbnd]=i*Ng[1]+j;
				if(id==1) kbnd[Nbnd]*=-1; 
			}
			Nbnd++; 
			break;
		case 2: //exterior grid
			Nex++;
			break;
		};
	}
	}
		if(m==0){
			kint=(int *)malloc(sizeof(int)*Nin);
			kbnd=(int *)malloc(sizeof(int)*Nbnd);
		}
	}
};
//		GENERATE 1D INDEX FOR v3 
void FIELD::gen_indx0(int **kcell){
	int i,j;

	for(int m=0;m<2;m++){
		Nin=0;
		Nbnd=0;
		Nex=0;
	for(i=0;i<Ng[0];i++){
	for(j=0;j<Ng[1];j++){
		switch(kcell[i][j]){
		case 0: // interior grid
			if(m==1) kint[Nin]=i*Ng[1]+j;
			Nin++;
			break;
		case 1: // exterior grid
			Nex++;
			break;
		case 2: //exterior grid	(not to be used)
			if(m==1){
				kbnd[Nbnd]=i*Ng[1]+j;
			}
			Nbnd++; 
			break;
		};
	}
	}
		if(m==0){
			kint=(int *)malloc(sizeof(int)*Nin);
			kbnd=(int *)malloc(sizeof(int)*Nbnd);
		}
	}
};
void FIELD::clear(){
	int i,j;
	for(i=0;i<Ng[0];i++){ 
		for(j=0;j<Ng[1];j++) F[i][j]=0.0;
	};
};


void FIELD::fwrite_bnd(char *fname){
	FILE *fp=fopen(fname,"w");
	if(fp==NULL) show_msg(fname);

	int k,i,j;
	for(k=0; k<Nbnd; k++){
		l2ij(abs(kbnd[k]),&i,&j);
		ij2xf(i,j);
		fprintf(fp,"%lf, %lf\n",xf[0],xf[1]);
	}
	fclose(fp);
};

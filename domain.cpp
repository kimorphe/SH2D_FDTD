#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"wveq2d.h"

double* DOMAIN::ij2xf(int i, int j, int type){

	if(type==0){	// midpoint grid 
		xf[0]=Xa[0]+(i+0.5)*dx[0];
		xf[1]=Xa[1]+(j+0.5)*dx[1];
		return(xf);
	};

	if(type==1){ // v1-grid  (leftword shifted grid)
		xf[0]=Xa[0]+i*dx[0];
		xf[1]=Xa[1]+(j+0.5)*dx[1];
		return(xf);
	};

	if(type==2){ // v2-grid (upward shifted grid)
		xf[0]=Xa[0]+(i+0.5)*dx[0];
		xf[1]=Xa[1]+j*dx[1];
		return(xf);
	};

	if(type==3){ // vertices 
		xf[0]=Xa[0]+i*dx[0];
		xf[1]=Xa[1]+j*dx[1];
		return(xf);
	};
}	
void DOMAIN::setup(double *xll, double *wdt, double *dh){
	Xa=xll;
	Wd=wdt;
	dx=dh;
};
void DOMAIN::print_prms(){
	printf("---------- DOMAIN SETTING-------------\n");
	printf("Xa=%lf, %lf\n",Xa[0],Xa[1]);
	printf("Wd=%lf, %lf\n",Wd[0],Wd[1]);
	printf("Ndiv=%d, %d\n",Ndiv[0],Ndiv[1]);
	printf("dx=%lf, %lf\n",dx[0],dx[1]);
	printf("--------------------------------------\n");
};
void DOMAIN::init(int *nxy){
	Ndiv=nxy;
	ndat=Ndiv[0]*Ndiv[1];

	int *p=(int *)malloc(sizeof(int)*ndat);
	kcell=(int **)malloc(sizeof(int*)*Ndiv[0]);

	int k=0;
	for( int i=0; i<Ndiv[0];i++){
		kcell[i]=p+i*Ndiv[1];
	for( int j=0; j<Ndiv[1];j++){
		p[k++]=0;
	}
	}
};
void DOMAIN::print_kcell(){
	for( int i=0; i<Ndiv[0];i++){
	for( int j=0; j<Ndiv[1];j++){
		printf("%d ", kcell[i][j]);
	}
		printf("\n");
	}
}
//  ----------- DOMAIN PERFORATION  --------------
void DOMAIN::perfo_ellip(char *fname){
	int i,j;
	double xcod[2],xc[2],a,b;
	double xx,yy;

	FILE *fp;
	char cbff[8];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo_ellip !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,8,fp) !=NULL){
		if(strcmp(cbff,"##Ellip")==0){
			fscanf(fp,"%lf, %lf, %lf, %lf\n",xc,xc+1,&a,&b);
			printf("%lf %lf %lf %lf\n",xc[0],xc[1],a,b);

			for(i=0;i<Ndiv[0];i++){
				xcod[0]=Xa[0]+dx[0]*(i+0.5);	
			for(j=0;j<Ndiv[1];j++){
				xcod[1]=Xa[1]+dx[1]*(j+0.5);	
				xx=(xcod[0]-xc[0])/a;
				yy=(xcod[1]-xc[1])/b;
				if(xx*xx+yy*yy <1.0) kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
};
//-------------EXPORT GEOMETRY DATA  ------------------
void DOMAIN:: out_kcell(){

	FILE *fp=fopen("kcell.dat","w");
	int i,j;
	double Xb[2];
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];


	fprintf(fp,"# Xa[2], Xb[2] (computational domian)\n");
	fprintf(fp," %lf, %lf\n %lf, %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);
	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d, %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# kcell[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp, "%d\n",kcell[i][j]);
	}
	}

	fflush(fp);
};
void DOMAIN::topography(char *fname){
	FILE *fp=fopen(fname,"r");
	int typ=1;	//piecewise linear function
	double xx,yy;
	double *xknot,*yknot;
	int nseg;
	char cbff[9];

	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo_ellip !");
		puts(fname);
		exit(-1);
	}

	int i,j,jy;
	int isgn=1;
	while(fgets(cbff,9,fp) !=NULL){
	if(strcmp(cbff,"##Tpgrph")==0){
		fscanf(fp,"%d\n",&typ);
		printf("type=%d\n",typ);
		if(typ<0) isgn=-1;
		typ=abs(typ);

		if(typ<3){
			fscanf(fp,"%d\n",&nseg);
			printf("nseg=%d\n",nseg);
			xknot=(double *)malloc(sizeof(double)*nseg);
			yknot=(double *)malloc(sizeof(double)*nseg);
			for(i=0;i<nseg;i++) fscanf(fp,"%lf, %lf\n",xknot+i,yknot+i);

			for(i=0;i<Ndiv[0];i++){
				xx=Xa[0]+(i+0.5)*dx[0];
				yy=pwfun(typ,xx,xknot,yknot,nseg,1.0);
				jy=int((yy-Xa[1])/dx[1]);
				if(isgn==1){
					for(j=jy;j<Ndiv[1];j++) kcell[i][j]=1;
				}
				if(isgn==-1){
					for(j=0;j<jy;j++) kcell[i][j]=1;
				}
			};
		}	
	}
	};
	fclose(fp);

};
double pwfun(
	int typ, 	// piecewise constant(0), linear(1)
	double x, 	// x-coordinate
	double *xknot, double *yknot,  // discontinuous/broken  points
	int ndiv, // number of intervals
	double y0 // magnicifation
){

	int idiv=0;
	while(xknot[idiv]<=x) idiv++;

	for(idiv=0; idiv<ndiv; idiv++){
		if(xknot[idiv]>x) break;
	};
	if(idiv==0) return(yknot[0]);
	if(idiv==ndiv) return(yknot[ndiv-1]);

	double y;
	switch(typ){
	case 1: // piecsewise linear function
		double xl,xr,xi;
		xl=xknot[idiv-1];
		xr=xknot[idiv];
		xi=(x-xl)/(xr-xl);
		y=(1.0-xi)*yknot[idiv-1]+xi*yknot[idiv];
		break;
	case 2: //piecewise constant function
		return(yknot[idiv]*y0);
	};

	return(y*y0);
};
int DOMAIN::find_v1bnd(int nx, double y){

	int jy=int((y-Xa[1])/dx[1]);
	int i,ibnd,il,ir;
	switch(nx){ 
	case 1:
		ir=1;
		ibnd=-1;
		for(i=Ndiv[0]-1;i>=0;i--){
			il=kcell[i][jy];
			if((ir-il)==1){
				ibnd=i+1;
				break;
			};
			ir=il;
		};
		break;
	case -1:
		il=1;
		ibnd=-1;
		for(i=0;i<Ndiv[0];i++){
			ir=kcell[i][jy];
			if((il-ir)==1){
				ibnd=i;
				break;
			};
			il=ir;
		};
		break;
	};
	return(ibnd);
};
int DOMAIN::find_v2bnd(int ny, double x){
	int ix=int((x-Xa[0])/dx[0]);
	int j,jbnd,iu,id;

	switch(ny){
	case 1:
		iu=1;
		jbnd=-1;
		for(j=Ndiv[1]-1;j>=0;j--){
			id=kcell[ix][j];
			if((iu-id)==1){
				jbnd=j+1;
				break;
			};
			iu=id;
		};
		break;
	case -1:
		id=1;
		jbnd=-1;
		for(j=0; j<Ndiv[1]; j++){
			iu=kcell[ix][j];
			if((id-iu)==1){
				jbnd=j;
				break;
			};
			id=iu;
		};
		break;
	};
	return(jbnd);
};

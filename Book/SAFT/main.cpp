#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "saft.h"


double dist2d(double *x1, double *x2){
	double rx=x1[0]-x2[0];
	double ry=x1[1]-x2[1];

	return(sqrt(rx*rx+ry*ry));
};

int main(int argc, char *argv[]){

	IMG img;
	Array ary;

	//-------Read Imaging Conditions----------------
	FILE *finp=fopen("saft.inp","r");
	if(finp==NULL){
		printf("Cannt find saft.inp\n");
		exit(-1);
	};
	char cbff[128];
	char path1[128],path2[128];
	char fimg[128];

	fgets(cbff,128,finp);
	fscanf(finp,"%s\n",path1);

	fgets(cbff,128,finp);
	fscanf(finp,"%s\n",path2);

	printf("Performing saft...\n");
	printf(" <-- %s\n",path1);
	printf(" <-- %s\n",path2);

	fgets(cbff,128,finp);
	fscanf(finp,"%s\n",fimg);
	printf(" --> %s\n",fimg);

	fgets(cbff,128,finp);
	int wv_typ;
	fscanf(finp,"%d\n",&wv_typ);

	double xc,yc,Wx,Wy;
	int Nx,Ny;
	fgets(cbff,128,finp);
	fscanf(finp,"%lf, %lf\n",&xc, &yc);
	fscanf(finp,"%lf, %lf\n",&Wx, &Wy);

	fgets(cbff,128,finp);
	fscanf(finp,"%d, %d\n",&Nx, &Ny);
	
	double ct;	// wave speed
	fgets(cbff,128,finp);
	fscanf(finp,"%lf\n",&ct);
	int iwgt;	// weight(0: uniform, 1: point desity considered)
	fgets(cbff,128,finp);
	fscanf(finp,"%d\n",&iwgt);
	printf("ct=%lf, iwgt=%d\n",ct,iwgt);
	fclose(finp);
	//--------------------------------------

	img.set_Xa(xc-0.5*Wx, yc-0.5*Wy);
	img.set_Wd(Wx, Wy);
	img.set_Ng(Nx,Ny);
	ary.setup(path1);
	ary.load();
	ary.set_angles(xc,yc);

	Bscan bwv(path2);

	//bwv.fwrite_Bscan();
	int i,j;
	int kx,ky;
	double xr[2],xsc[2];
	double *xi;
	double tof,dI;


	int k,count;
	double rin,rsc;
	k=0;
	for(i=0;i<ary.nscan;i++){
		printf("iscan=%d\n",i);
		xi=ary.set_transmitter(i);
		bwv.load(i,wv_typ);
	for(j=0;j<ary.nele;j++){
		xr[0]=ary.xcod[j];
		xr[1]=ary.ycod[j];
		count=1;
		if(iwgt==1) count=ary.get_hist_val(k);
		for(kx=0;kx<img.Ng[0];kx++){
		for(ky=0;ky<img.Ng[1];ky++){
			xsc[0]=img.xcod[kx];
			xsc[1]=img.ycod[ky];

			rin=dist2d(xsc,xi);
			rsc=dist2d(xsc,xr);
			tof=(rin+rsc)/ct;
			dI=bwv.get_amp(j,tof)/count;
			//dI=bwv.get_amp(j,tof);
			img.V[kx][ky]+=dI;
		};
		};
		k++;
	};
	};
	img.fwrite_V(fimg);


	return(0);
};



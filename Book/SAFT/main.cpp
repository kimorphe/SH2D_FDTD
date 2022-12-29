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

	img.set_Xa(-10.0, 10.0);
	img.set_Wd( 20.0, 15.0);
	img.set_Ng(201, 151);

	char fimg[128]="saft.out";

	ary.setup();
	ary.load();

	Bscan bwv;

	//bwv.fwrite_Bscan();

	int i,j,kx,ky;
	double xr[2],xsc[2];
	double *xi;
	double rin,rsc;

	double ct=1.0,tof,dI;

	for(i=0;i<ary.nscan;i++){
		printf("iscan=%d\n",i);
		xi=ary.set_transmitter(i);
		bwv.load(i);
	for(j=0;j<ary.nele;j++){
			xr[0]=ary.xcod[j];
			xr[1]=ary.ycod[j];
		for(kx=0;kx<img.Ng[0];kx++){
		for(ky=0;ky<img.Ng[1];ky++){
			xsc[0]=img.xcod[kx];
			xsc[1]=img.ycod[ky];

			rin=dist2d(xsc,xi);
			rsc=dist2d(xsc,xr);
			tof=(rin+rsc)/ct;
			dI=bwv.get_amp(j,tof);
			img.V[kx][ky]+=dI;
		};
		};
	};
	};
	img.fwrite_V(fimg);


	return(0);
};



#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

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
	if(typ==0){
		return(yknot[idiv]*y0);
	}else{
		double xl,xr,xi;
		xl=xknot[idiv-1];
		xr=xknot[idiv];
		xi=(x-xl)/(xr-xl);
		y=(1.0-xi)*yknot[idiv-1]+xi*yknot[idiv];
	};

	return(y*y0);
};

int main(){

	int ndiv=5;
	double xknot[5],yknot[5];
	xknot[0]=-1.0; yknot[0]=1.0;
	xknot[1]= 1.0; yknot[1]=3.0;
	xknot[2]= 5.0; yknot[2]=-1.0;
	xknot[3]= 8.0; yknot[3]=1.0;
	xknot[4]=11.0; yknot[4]=2.0;
	double y0=1.0;

	double dx=0.2;
	double x,y;
	for(int i=0;i<101;i++){
		x=i*dx-2.0;	
		y=pwfun(0,x,xknot,yknot,ndiv,y0);
		printf("%lf, %lf\n",x,y);
	};
};

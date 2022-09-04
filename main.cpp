#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"wveq2d.h"

int main(int argc, char *argv[]){

	CNTRL ctr;
	char fgeom[128]="geom.inp";
	char ftset[128]="tset.inp";
	char finwv[128]="inwv0.dat";
	char fsrce[128]="src.inp";
	char frecs[128]="recs.inp";
	char farry[128]="array.inp";

	ctr.setup_domain(fgeom);
	ctr.time_setting(ftset);
		ctr.src_setting(fsrce);
	ctr.wvfm_setting(finwv);
	//ctr.rec_setting(frecs);
	ctr.array_setting(farry);
	
	int it,isum,m;
	
	for(m=0;m<ctr.ary.nmeas;m++){

	isum=0;
	ctr.v3.fwrite(isum++);
	for(it=0; it<ctr.Nt; it++){
		ctr.v2q(it); 
		//ctr.record(it);
		ctr.q2v(it);
		ctr.capture(it);
		if((it+1)%10==0){
			printf("it=%d,isum=%d\n",it,isum);
		       	ctr.v3.fwrite(isum++);
		};
	};
	//for(int i=0; i<ctr.nrec; i++) ctr.recs[i].fwrite();
	for(int i=0; i<ctr.nsrc; i++) ctr.srcs[i].fwrite();

	ctr.round++;
	ctr.clear();
	}

	return(0);

};

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
	ctr.setup_domain(fgeom);
	ctr.time_setting(ftset);
	ctr.src_setting(fsrce);
	ctr.wvfm_setting(finwv);
	ctr.rec_setting(frecs);

	
	int isum=0;
	ctr.pr.fwrite(isum++);

	int it;
	for(it=0; it<ctr.Nt; it++){
		ctr.p2v(it); 
		ctr.record(it);
		ctr.v2p(it);
		if((it+1)%10==0){
			printf("it=%d,isum=%d\n",it,isum);
		       	ctr.pr.fwrite(isum++);
		};
	};
	for(int i=0;i <ctr.nrec;i++) ctr.recs[i].fwrite();
	return(0);

};

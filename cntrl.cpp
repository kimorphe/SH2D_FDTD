#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"wveq2d.h"

void CNTRL::setup_domain(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	if(fp==NULL){
		printf("Cannot find file %s \n",fname);
		printf(" --> abort process ...\n");
		exit(-1);
	};

	// Load domain setting
	fgets(cbff,128,fp); 
	fscanf(fp,"%lf, %lf\n",Xa,Xa+1);	// lower-left corner
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",Wd,Wd+1);	// width & height

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",Ndiv,Ndiv+1);	// number of cells along x & y axes

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",&ct, &rho);
	amu=rho*ct*ct;

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",Ha,Ha+1);	// PML thickness 
	fscanf(fp,"%lf, %lf\n",Hb,Hb+1); 	// PML thickness

	fclose(fp);

	dx[0]=Wd[0]/Ndiv[0];
	dx[1]=Wd[1]/Ndiv[1];

	// domain extention 
	for(int i=0; i<2; i++){
		NHa[i]=ceil(Ha[i]/dx[i]);
		NHb[i]=ceil(Hb[i]/dx[i]);
		Ha[i]=dx[i]*NHa[i];
		Hb[i]=dx[i]*NHb[i];

		Xa[i]-=Ha[i];
		Wd[i]=Wd[i]+Ha[i]+Hb[i];
		Ndiv[i]=Ndiv[i]+NHa[i]+NHb[i];
	}

	dm.setup(Xa,Wd,dx);
	dm.ct=ct; dm.rho=rho; dm.amu=amu;
	dm.init(Ndiv);

	double gmm=1.e-04;	// expected decay 
	dm.Ha=Ha;
	dm.Hb=Hb;
	dm.PML_setup(gmm);

	int i;
	FILE *ftmp=fopen("log.txt","w");
	double tmp;
	for(i=0; i<Ndiv[0]; i++){
		tmp=Xa[0]+dx[0]*(i+0.5);
		fprintf(ftmp,"%lf, %lf\n",tmp,dm.PML_dcy(0,tmp));
	}
	fprintf(ftmp,"\n");
	for(i=0; i<Ndiv[1]; i++){
		tmp=Xa[1]+dx[1]*(i+0.5);
		fprintf(ftmp,"%lf, %lf\n",tmp,dm.PML_dcy(1,tmp));
	}
	fclose(ftmp);

	dm.perfo_ellip(fname);
	puts("ellip done");
	dm.topography(fname);
	puts("toppog done");
	dm.out_kcell();
	
	// Setup staggered grid system 
	v3.init(Ndiv,0);	
	v3x.init(Ndiv,0);	
	v3y.init(Ndiv,0);	
	q1.init(Ndiv,1);	
	q2.init(Ndiv,2);	

	v3.setup(Xa,Wd,dx);
	v3x.setup(Xa,Wd,dx);
	v3y.setup(Xa,Wd,dx);
	q1.setup(Xa,Wd,dx);
	q2.setup(Xa,Wd,dx);

	//dm.print_prms();
	//v3.print_prms();
	//q1.print_prms();
	//q2.print_prms();
	
	v3.gen_indx0(dm.kcell);
	v3x.gen_indx0(dm.kcell);
	v3y.gen_indx0(dm.kcell);
	q1.gen_indx1(dm.kcell);
	q2.gen_indx2(dm.kcell);

};

bool CNTRL::out_time(int it){
	if(it==iout){
		printf("it/Nt=%d/%d (Nout=%d)\n",it,Nt,iout);
		iout+=Ninc;
		return(true);
	}
	return(false);
};
void CNTRL::time_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %d\n",&Tf, &Nt);
	dt=Tf/(Nt-1);
	printf("dt=%lf, Nt=%d, Tf=%lf\n",dt,Nt,Tf);
	CNTRL::CFL();

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf, %d\n",&tout_s, &tout_e, &Nout);
	Ninc=(tout_e-tout_s)/((Nout-1)*dt);
	if(Ninc<1) Ninc=1;
	iout=floor(tout_s/dt);
	if(iout==0) iout+=Ninc;
	iout0=iout;
	printf("Ninc=%d, iout=%d\n",Ninc, iout);

	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	int ftyp;
	fscanf(fp,"%d\n",&ftyp);
	if(ftyp==1){
		double xs,ys,sig,f0;
		fgets(cbff,128,fp);
		fscanf(fp,"%lf, %lf, %lf, %lf\n",&xs,&ys,&sig,&f0);
		printf("xs=(%lf, %lf), sig=%lf, p0=%lf\n",xs,ys,sig,f0);
		if(ftyp==0) v3.set_IC(xs,ys,sig,f0);
		if(ftyp==1) q1.set_IC(xs,ys,sig,f0);
		if(ftyp==2) q2.set_IC(xs,ys,sig,f0);
	}

	fclose(fp);
};
void CNTRL::wvfm_setting(char *fname){
	wvs=(Wv1D *)malloc(sizeof(Wv1D)*nwv);

	for(int i=0;i<nwv;i++){
		sprintf(fname,"inwv%d.dat",i);
		wvs[i].gen_wv(fname);
	}
	
	//char fnout[128];
	//wv.FFT(1);
	//sprintf(fnout,"awvt.out");
	//wv.out_amp(fnout,',');
	//awv.Butterworth(16.0,10.0);
	//sprintf(fnout,"awvw.out");
	//wv.out_Amp(fnout,0);
};
int CNTRL::src_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	double xy,wdt,w2;
	int ityp;	// 1: v1-source, 2: v2-source
	int nml;	// 1,-1 (normal vector, nx/ny)
	int iwv;	// waveform No.
	double xsrc,ysrc,ain;
	int i0,i,j,i1,i2,j1,j2,ng;
	int isrc,jsrc,isum;

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nsrc);
	srcs=(TRNSDCR *)malloc(sizeof(TRNSDCR)*nsrc);
	fgets(cbff,128,fp);
	nwv=0;
	for(int k=0; k<nsrc; k++){	// k-th source element
		fscanf(fp,"%d, %lf, %lf, %d, %d, %lf\n",&ityp, &xy, &wdt, &nml, &iwv, &ain);
		if(iwv>nwv) nwv=iwv;
		srcs[k].ID=k;
		srcs[k].nml=nml;
		srcs[k].type=ityp;
		srcs[k].iwv=iwv;
		srcs[k].Xa=Xa;
		srcs[k].dx=dx;
		srcs[k].set_inc_ang(ain,ct);
		w2=wdt*0.5;
		isum=0;
		switch(ityp){
		case 1:	// v1-source
			ng=ceil(wdt/dx[1]);
			if(ng==0) ng=1;
			if(ng%2==0){
				i0=int((xy-Xa[1])/dx[1]+0.5);
				xy=Xa[1]+i0*dx[1];
			}else{
				i0=int((xy-Xa[1])/dx[1]);
				xy=Xa[1]+(i0+0.5)*dx[1];
			};
			printf("xy=%lf\n",xy);
			j1=int((xy-w2-Xa[1])/dx[1]);
			srcs[k].mem_alloc(ng);
			for(j=0; j<ng; j++){
				jsrc=j+j1;
				if(jsrc<0) continue;
				if(jsrc>=Ndiv[1]) break;
				ysrc=Xa[1]+(jsrc+0.5)*dx[1];
				isrc=dm.find_q1bnd(nml,ysrc);
				srcs[k].isrc[isum]=isrc;
				srcs[k].jsrc[isum]=jsrc;
				srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				isum++;
			}
			break;
		case 2: // v2-source
			ng=ceil(wdt/dx[0]);
			if(ng%2==0){
				i0=int((xy-Xa[0])/dx[0]+0.5);
				xy=Xa[0]+i0*dx[0];
			}else{
				i0=int((xy-Xa[0])/dx[0]);
				xy=Xa[0]+(i0+0.5)*dx[0];
			};
			printf("xy=%lf\n",xy);
			i1=int((xy-w2-Xa[0])/dx[0]);
			if(ng==0) ng=1;
			srcs[k].mem_alloc(ng);
			for(i=0; i<ng; i++){
				isrc=i+i1;
				if(isrc<0) continue;
				if(isrc>=Ndiv[0]) break;
				xsrc=Xa[0]+(0.5+isrc)*dx[0];
				jsrc=dm.find_q2bnd(nml,xsrc);
				srcs[k].isrc[isum]=isrc;
				srcs[k].jsrc[isum]=jsrc;
				srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				isum++;
			}
			break;
		};
		ng=isum;
		srcs[k].print();
		srcs[k].set_center();
		srcs[k].init_bwv(Nt,dt);
	}
	nwv++;
	printf("nwv=%d\n",nwv);
	fclose(fp);
	return(nwv);
};
void CNTRL::array_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	if(fp==NULL){
		printf("File %s cannot found !\n",fname);
		printf(" --> abort process..\n");
		exit(-1);
	};
	int nele,nmeas;
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nele);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nmeas);
	if(nele != nsrc){
		printf(" Error: nele(=%d) must be equal to nsrc(=%d) !\n",nele,nsrc);
		exit(-1);
	};
	printf("nmeas=%d\n",nmeas);	

	ary.init(nsrc,nmeas);
	int i,j,k=0;
	for(j=0;j<nmeas;j++){ 
		fgets(cbff,128,fp);
		for(i=0;i<nsrc;i++){
			fscanf(fp,"%d, %lf, %lf\n",ary.actv+k, ary.a0+k, ary.tdly+k);
			k++;
		};
	};
	round=0;
	ary.print();
	fclose(fp);
};
int CNTRL::rec_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	int ityp,nml,ng;
	int i,j,i0;
	int i1,i2,j1,j2,irec,jrec,isum;
	double xrec,yrec;
	double w2,wdt,xy,thr;
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nrec);
	fgets(cbff,128,fp);
	recs=(RECVR *)malloc(sizeof(RECVR)*nrec);
	for(int ir=0; ir<nrec; ir++){
		fscanf(fp,"%d, %lf, %lf, %d, %lf\n",&ityp, &xy, &wdt, &nml, &thr);
		recs[ir].type=ityp;
		recs[ir].nml=nml;
		recs[ir].wd=wdt;
		recs[ir].ID=ir;
		recs[ir].set_cod(Xa,dx);
		recs[ir].set_inc_ang(thr,ct);
		w2=0.5*wdt;
		isum=0;
		switch(ityp){
		case 1:	// v1-receiver
			ng=ceil(wdt/dx[1]);
			if(ng==0) ng=1;
			if(ng%2==0){
				i0=int((xy-Xa[1])/dx[1]+0.5);
				xy=Xa[1]+i0*dx[1];
			}else{
				i0=int((xy-Xa[1])/dx[1]);
				xy=Xa[1]+(i0+0.5)*dx[1];
			};

			j1=int((xy-w2-Xa[1])/dx[1]);
			recs[ir].mem_alloc(ng);
			for(j=0; j<ng; j++){
				jrec=j+j1;
				if(jrec<0) continue;
				if(jrec>=Ndiv[1]) break;
				yrec=Xa[1]+(jrec+0.5)*dx[1];
				irec=dm.find_q1bnd(nml,yrec);
				if(nml==1) irec--;
				recs[ir].irec[isum]=irec;
				recs[ir].jrec[isum]=jrec;
				isum++;
			}
			break;
		case 2: // v2-receiver
			ng=ceil(wdt/dx[1]);
			if(ng==0) ng=1;
			if(ng%2==0){
				i0=int((xy-Xa[0])/dx[0]+0.5);
				xy=Xa[0]+i0*dx[0];
			}else{
				i0=int((xy-Xa[0])/dx[0]);
				xy=Xa[0]+(i0+0.5)*dx[0];
			};

			i1=int((xy-w2-Xa[0])/dx[0]);
			recs[ir].mem_alloc(ng);
			for(i=0; i<ng; i++){
				irec=i+i1;
				if(irec<0) continue;
				if(irec>=Ndiv[0]) break;
				xrec=Xa[0]+(0.5+irec)*dx[0];
				jrec=dm.find_q2bnd(nml,xrec);
				if(nml==1) jrec--;
				recs[ir].irec[isum]=irec;
				recs[ir].jrec[isum]=jrec;
				isum++;
			};
			break;
		};
		ng=isum;
		recs[ir].print();
		recs[ir].init_bwv(Nt,dt);
	};

	return(nrec);
};

void CNTRL::record(int jt){
	int i,j,type;
	for(j=0;j<nrec;j++){
		type=recs[j].type;
		if(type==1) recs[j].record(jt,q1.F);
		if(type==2) recs[j].record(jt,q2.F);
	}
};

void CNTRL::capture(int jt){
	int i,j,type;
	for(j=0;j<nsrc;j++) srcs[j].record(jt,v3.F);
};
double CNTRL::CFL(){
	double dh=dx[0];
	double Crt;
	if(dh >dx[1]) dh=dx[1];
	dh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
	Crt=ct*dt/dh;
	printf("Crt=%lf\n",Crt);
	if(Crt>1.0){
		printf(" stability condition is not satisfied !!\n --> abort proces\n");
		exit(-1);
	};
	return(Crt);
};
void CNTRL::q2v(int itime){
	int i,j,k,l;
	int nx,ny;
	double **F1,**F2;
	double dFx, dFy;
	double beta,gmm, *xcod;

	nx=v3.Ng[0];
	ny=v3.Ng[1];
	F1=q1.F;
	F2=q2.F;

	for(k=0; k<v3.Nin; k++){
		l=v3.kint[k];
		v3.l2ij(l,&i,&j);
		xcod=dm.ij2xf(i,j,0);

		beta=0.5*dt*dm.PML_dcy(0,xcod[0]);
		gmm=dt/(rho*dx[0]);
		dFx=(F1[i+1][j]-F1[i][j])*gmm;
		v3x.F[i][j]=((1.-beta)*v3x.F[i][j]+dFx)/(1.+beta);

		beta=0.5*dt*dm.PML_dcy(1,xcod[1]);
		gmm=dt/(rho*dx[1]);
		dFy=(F2[i][j+1]-F2[i][j])*gmm;
		v3y.F[i][j]=((1.-beta)*v3y.F[i][j]+dFy)/(1.+beta);

		//v3.F[i][j]+=(dFx+dFy);
		v3.F[i][j]=v3x.F[i][j]+v3y.F[i][j];
	}
};
void CNTRL::v2q(int itime){
	int i,j,k,l;
	int nx,ny;
	double **F;
	double dFx,dFy,alph;

	F=v3.F;
	nx=q1.Ng[0]; ny=q1.Ng[1];

	for(k=0; k<q1.Nin; k++){
		l=q1.kint[k];
		q1.l2ij(l,&i,&j);
		alph=amu*dt/dx[0];
		//dFx=(F[i][j]-F[i-1][j])/dx[0]*dt/rho;
		dFx=alph*(F[i][j]-F[i-1][j]);
		q1.F[i][j]+=dFx;
	};

	for(k=0; k<q2.Nin; k++){
		l=q2.kint[k];
		q2.l2ij(l,&i,&j);
		alph=amu*dt/dx[1];
		//dFy=(F[i][j]-F[i][j-1])/dx[1]*dt/rho;
		dFy=alph*(F[i][j]-F[i][j-1]);
		q2.F[i][j]+=dFy;
	};

	//SOURCE src;
	TRNSDCR src;

	double bvl,tdly,tdly0,a0;
	int sgn,sft,iwv,i0;
	i0=nsrc*round;
	for(int isrc=0;isrc<nsrc;isrc++){
		if(ary.actv[isrc+i0]==0) continue;

		src=srcs[isrc]; // SOURCE
		iwv=src.iwv;

		a0=ary.a0[isrc+i0];
		tdly0=ary.tdly[isrc+i0];
		//bvl=wvs[iwv].amp[itime];
		sft=1;
		sgn=src.nml;
		if(sgn ==-1) sft=0;

		switch(src.type){
		case 1:	// v1 source
			for(k=0; k<src.ng; k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				tdly=src.Ctd*k+tdly0;
				bvl=wvs[iwv].val(itime*dt-tdly)*a0;
				//dFx=sgn*(bvl-F[i-sft][j])/dx[0]*dt/rho;
				//q1.F[i][j]+=dFx;
				q1.F[i][j]=bvl*sgn;
			}
			break;
		case 2: // v2 source
			for(k=0; k<src.ng; k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				tdly=src.Ctd*k+tdly0;
				bvl=wvs[iwv].val(itime*dt-tdly)*a0;
				//dFy=sgn*(bvl-F[i][j-sft])/dx[1]*dt/rho;
				//q2.F[i][j]+=dFy;
				q2.F[i][j]=bvl*sgn;
			}
			break;
		};
	}

};

void CNTRL::clear(){
	v3.clear();
	q1.clear();
	q2.clear();
	v3x.clear();
	v3y.clear();
	for(int i=0;i<nsrc;i++) srcs[i].clear();
	for(int i=0;i<nrec;i++) recs[i].clear();
	iout=iout0;
};
void CNTRL::fwrite_ary(){
	int j,k;
	char fname[128];
	FILE *fp;

	//	ARRAY WAVEFORMS 
	sprintf(fname,"T%d/ary.out",round);
	fp=fopen(fname,"w");
	fprintf(fp,"# nele, Nt\n");
	fprintf(fp,"%d, %d\n",nsrc,Nt);
	for(j=0;j<nsrc;j++){
		fprintf(fp,"# e=%d, %lf, %lf\n",j,srcs[j].x0,srcs[j].y0);
		for(k=0;k<Nt;k++){
			fprintf(fp,"%lf, %lf\n",dt*k, srcs[j].amp_synth(k));
		}
	};
	fclose(fp);
};
void CNTRL::snapshot(int n_meas, int isum){
	v3.fwrite_trim(n_meas,isum,NHa,NHb);
};

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"wveq2d.h"

void CNTRL::setup_domain(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	

	// Load domain setting
	fgets(cbff,128,fp); 
	fscanf(fp,"%lf, %lf\n",Xa,Xa+1);	// lower-left corner
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",Wd,Wd+1);	// width & height

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",Ndiv,Ndiv+1);	// number of cells along x & y axes

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",&cp, &rho);
	almb=rho*cp*cp;

	printf("cp=%lf, rho=%lf\n",cp,rho);
	fclose(fp);

	dx[0]=Wd[0]/Ndiv[0];
	dx[1]=Wd[1]/Ndiv[1];
	dm.setup(Xa,Wd,dx);
	dm.init(Ndiv);

	dm.perfo_ellip(fname);
	dm.topography(fname);
	dm.out_kcell();
	
	// Setup staggered grid system 
	pr.init(Ndiv,0);	
	v1.init(Ndiv,1);	
	v2.init(Ndiv,2);	
	pr.setup(Xa,Wd,dx);
	v1.setup(Xa,Wd,dx);
	v2.setup(Xa,Wd,dx);

	//dm.print_prms();
	//pr.print_prms();
	//v1.print_prms();
	//v2.print_prms();
	
	pr.gen_indx0(dm.kcell);
	v1.gen_indx1(dm.kcell);
	v2.gen_indx2(dm.kcell);

};

void CNTRL::time_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %d\n",&dt, &Nt);
	printf("dt=%lf, Nt=%d\n",dt,Nt);
//	dm.CFL(dt);
	CNTRL::CFL();

	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	int ftyp;
	fscanf(fp,"%d\n",&ftyp);
	double xs,ys,sig,f0;
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf, %lf, %lf\n",&xs,&ys,&sig,&f0);
	printf("xs=(%lf, %lf), sig=%lf, p0=%lf\n",xs,ys,sig,f0);

	if(ftyp==0) pr.set_IC(xs,ys,sig,f0);
	if(ftyp==1) v1.set_IC(xs,ys,sig,f0);
	if(ftyp==2) v2.set_IC(xs,ys,sig,f0);
	//pr.fwrite(0);

	fclose(fp);
};
void CNTRL::wvfm_setting(char *fname){
	wvs=(Wv1D *)malloc(sizeof(Wv1D)*nwv);

	for(int i=0;i<nwv;i++){
		sprintf(fname,"inwv%d.dat",i);
		printf("fname=%s\n",fname);
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
	double xsrc,ysrc;
	int i,j,i1,i2,j1,j2,ng;
	int isrc,jsrc,isum;

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nsrc);
	srcs=(SOURCE *)malloc(sizeof(SOURCE)*nsrc);
	fgets(cbff,128,fp);
	nwv=0;
	for(int k=0; k<nsrc; k++){	// k-th source element
		fscanf(fp,"%d, %lf, %lf, %d, %d\n",&ityp, &xy, &wdt, &nml, &iwv);
		if(iwv>nwv) nwv=iwv;
		srcs[k].nml=nml;
		srcs[k].type=ityp;
		srcs[k].iwv=iwv;
		w2=wdt*0.5;
		isum=0;
		switch(ityp){
		case 1:	// v1-source
			j1=int((xy-w2-Xa[1])/dx[1]);
			j2=int((xy+w2-Xa[1])/dx[1]);
			ng=j2-j1+1;
			srcs[k].mem_alloc(ng);
			for(j=0; j<ng; j++){
				jsrc=j+j1;
				if(jsrc<0) continue;
				if(jsrc>=Ndiv[1]) break;
				ysrc=Xa[1]+(jsrc+0.5)*dx[1];
				isrc=dm.find_v1bnd(nml,ysrc);
				srcs[k].isrc[isum]=isrc;
				srcs[k].jsrc[isum]=jsrc;
				srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				isum++;
			}
			break;
		case 2: // v2-source
			i1=int((xy-w2-Xa[0])/dx[0]);
			i2=int((xy+w2-Xa[0])/dx[0]);
			ng=i2-i1+1;
			srcs[k].mem_alloc(ng);
			for(i=0; i<ng; i++){
				isrc=i+i1;
				if(isrc<0) continue;
				if(isrc>=Ndiv[0]) break;
				xsrc=Xa[0]+(0.5+isrc)*dx[0];
				jsrc=dm.find_v2bnd(nml,xsrc);
				srcs[k].isrc[isum]=isrc;
				srcs[k].jsrc[isum]=jsrc;
				srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				isum++;
			}
			break;
		};
		ng=isum;
		srcs[k].print(Xa, dx);
	}
	nwv++;
	printf("nwv=%d\n",nwv);
	fclose(fp);
	return(nwv);
};
int CNTRL::rec_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	int ityp,nml,ng;
	int i,j;
	int i1,i2,j1,j2,irec,jrec,isum;
	double xrec,yrec;
	double w2,wdt,xyrec;
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nrec);
	fgets(cbff,128,fp);
	recs=(RECVR *)malloc(sizeof(RECVR)*nrec);
	for(int ir=0; ir<nrec; ir++){
		fscanf(fp,"%d, %lf, %lf, %d\n",&ityp, &xyrec, &wdt, &nml);
		recs[ir].type=ityp;
		recs[ir].nml=nml;
		recs[ir].wd=wdt;
		recs[ir].ID=ir;
		recs[ir].set_cod(Xa,dx);
		w2=0.5*wdt;
		isum=0;
		switch(ityp){
		case 1:	// v1-receiver
			j1=int((xyrec-w2-Xa[1])/dx[1]);
			j2=int((xyrec+w2-Xa[1])/dx[1]);
			ng=j2-j1+1;
			recs[ir].mem_alloc(ng);
			for(j=0; j<ng; j++){
				jrec=j+j1;
				if(jrec<0) continue;
				if(jrec>=Ndiv[1]) break;
				yrec=Xa[1]+(jrec+0.5)*dx[1];
				irec=dm.find_v1bnd(nml,yrec);
				recs[ir].irec[isum]=irec;
				recs[ir].jrec[isum]=jrec;
				isum++;
			}
			break;
		case 2: // v2-receiver
			i1=int((xyrec-w2-Xa[0])/dx[0]);
			i2=int((xyrec+w2-Xa[0])/dx[0]);
			ng=i2-i1+1;
			recs[ir].mem_alloc(ng);
			for(i=0; i<ng; i++){
				irec=i+i1;
				if(irec<0) continue;
				if(irec>=Ndiv[0]) break;
				xrec=Xa[0]+(0.5+irec)*dx[0];
				jrec=dm.find_v2bnd(nml,xrec);
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
		if(type==1) recs[j].record(jt,v1.F);
		if(type==2) recs[j].record(jt,v2.F);
	}
};
double CNTRL::CFL(){
	double dh=dx[0];
	double Crt;
	if(dh >dx[1]) dh=dx[1];
	dh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
	Crt=cp*dt/dh;
	printf("Crt=%lf\n",Crt);
	if(Crt>1.0){
		printf(" stability condition is not satisfied !!\n --> abort proces\n");
		exit(-1);
	};
	return(Crt);
};
void CNTRL::v2p(int itime){
	int i,j;
	int nx,ny;
	double **F1,**F2;

	nx=pr.Ng[0];
	ny=pr.Ng[1];
	F1=v1.F;
	F2=v2.F;
	double dFx, dFy;
	int k,l;
	for(k=0;k<pr.Nin;k++){
		l=pr.kint[k];
		pr.l2ij(l,&i,&j);
		dFx=(F1[i+1][j]-F1[i][j])/dx[0]*dt*almb;
		dFy=(F2[i][j+1]-F2[i][j])/dx[1]*dt*almb;
		pr.F[i][j]+=(dFx+dFy);
	}
};
void CNTRL::p2v(int itime){
	int i,j;
	int nx,ny;
	double **F;
	double dFx, dFy;
	F=pr.F;
	nx=v1.Ng[0]; ny=v1.Ng[1];

	int k,l;
	for(k=0;k<v1.Nin;k++){
		l=v1.kint[k];
		v1.l2ij(l,&i,&j);
		dFx=(F[i][j]-F[i-1][j])/dx[0]*dt/rho;
		v1.F[i][j]+=dFx;
	};

	for(k=0;k<v2.Nin;k++){
		l=v2.kint[k];
		v2.l2ij(l,&i,&j);
		dFy=(F[i][j]-F[i][j-1])/dx[1]*dt/rho;
		v2.F[i][j]+=dFy;
	};

	SOURCE src;

	double bvl;
	int sgn,sft,iwv;

	for(int isrc=0;isrc<nsrc;isrc++){
		//src=dm.srcs[isrc];
		src=srcs[isrc];
		iwv=src.iwv;
		bvl=wvs[iwv].amp[itime];
		sft=1;
		sgn=src.nml;
		if(sgn ==-1) sft=0;

		switch(src.type){
		case 1:	// v1 source
			for(k=0;k<src.ng;k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				dFx=sgn*(bvl-F[i-sft][j])/dx[0]*dt/rho;
				v1.F[i][j]+=dFx;
			}
			break;
		case 2: // v2 source
			for(k=0;k<src.ng;k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				dFy=sgn*(bvl-F[i][j-sft])/dx[1]*dt/rho;
				v2.F[i][j]+=dFy;
			}
			break;
		};
	}

};


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"wave1d.h"
//--------------------------------------------------------------------------
class FIELD{
	public:
		int Ng[2];
		int *Ndiv;
		int NHa[2], NHb[2]; 
		int ndat;
		double **F;
		void print_F();
		void init(int *Ndiv, int ityp);
		int type; 	// 0: pr, 1:v1, 2: v2, 3:vertices
		void setup(double *Xa, double *Wd, double *dx);
		double *Xa,*Wd,*dx;
		double xf[2];
		double *ij2xf(int i, int j);
		int ij2l(int i, int j);
		void l2ij(int l, int *ix, int *jy);
		double  ofst[2];
		void fwrite(int d_num, int num);
		void fwrite_trim(int d_num, int num, int *NHa, int *NHb);
		void set_IC(double xc, double yc, double sig, double f0);
		void print_prms();
		void gen_indx0(int **kcell);
		void gen_indx1(int **kcell);
		void gen_indx2(int **kcell);
		int *kbnd, *kint;
		int Nin,Nbnd,Nex;
		void clear();
	private:
};
class SOURCE{
	public:
		int ID;
		double wd;
		int ng;	// number of grids
		void mem_alloc(int n);
		int *isrc,*jsrc,*ksrc;
		int type;	// 1:v1, 2:v2
		int nml;	// normal vector (nx or ny)
		void set(double x, double y, double w);
		void print();
		int iwv;	// waveform number 
		double Ctd;	// delay time gradient
		bool Ton; // activated 
		double *dx, *Xa;
		//double dt;
		double x0,y0;
		void set_center();
		void set_inc_ang(double th, double ct);
	private:
};
class TRNSDCR: public SOURCE{
	public:
		bool Ron; //activated
		int Nt;
		double dt;
		void init_bwv(int n, double dtau);
		double **bwv;
		double val(int ig, double tt);
		double amp_synth(int it);
		void record(int it, double **fld);
		void fwrite();
		double mean_amp(int i);
		void clear();
	private:
};
class RECVR{
	public:
		int ID;		
		double wd;
		int ng;
		int *irec, *jrec;
		int type;	// 1:v1, 2:v2
		int nml;
		double Ctd;	// delay time gradient
		void print();
		double **bwv;
		void mem_alloc(int n);
		void init_bwv(int n, double dtau);
		int Nt;
		double dt;
		void record(int it, double **fld);
		double *Xa, *dx;
		void set_cod(double *Xa, double *dx);
		void set_inc_ang(double th, double ct);
		void fwrite_dump(int n_meas);
		void fwrite(int n_meas);
		double val(int ig, double tt);
		double amp_synth(int it);
		void clear();
	private:
};
class ARRAY{
	public:
		int nele;	// number of array elements
		int nmeas;	
		int *actv;	// activation
		double *a0;	// apodization
		double *tdly;	// time delay
		void init(int ne, int nm);
		int i0;
		void print();
//		TRNSDCR *trs;
		void fwrite();
	private:
};
class DOMAIN{
	public:
		double *Xa, *Wd, *dx; 	// domain location, width, cell size
		double *Ha, *Hb;	// PML thickness 
		int **kcell; 	// geometry ( binary image )
		int *Ndiv,ndat; // number of cells 
		void init(int *Ndiv); // allocate 2D kcell array
		void setup(double *Xa, double *Wd, double *dx); // set domain size
		double xf[2]; // grid coordinate
		double *ij2xf(int i, int j, int type); 
		void print_kcell(); //print kcell data
		void print_prms();// print domain related parameters
		double ct;	// phase velocity
		double rho,amu; // density, shear rigidity
		double dt;	// time step
		void perfo_ellip(char *fname); // domain perforation
		void out_kcell(); // write geometry(kcell) data
		void topography(char *fname); // set surface topography
		int find_q1bnd(int nx, double y); // find boundary v1-grid 
		int find_q2bnd(int ny, double x); // find boundary v2-grid
		int nsrc; // number of source elements
		double PML_dcy(int idir, double xy); // evaluate PML decay constant 
		void PML_setup(double gm);
		double A0[2],B0[2];
	private:
};
class CNTRL{
	public:
		double ct,rho,amu; 
		double Xa[2],Wd[2],dx[2];
		int Ndiv[2],Nt;
		double Ha[2], Hb[2];
		int NHa[2], NHb[2];
		double dt,Tf;
		double tout_s,tout_e;
		int Nout,Ninc,iout,iout0;
		bool out_time(int i);
		DOMAIN dm;
		FIELD v3,q1,q2; 
		FIELD v3x,v3y;
		TRNSDCR *srcs;		
		RECVR *recs;
		Wv1D *wvs;
		ARRAY ary;
		int round;
		int nwv;
		int nsrc;
		int nrec;
		void setup_domain(char *fname);
		void time_setting(char *fname);
		void wvfm_setting(char *fname);
		int src_setting(char *fname);
		int rec_setting(char *fname);
		void array_setting(char *fname);
		void v2q(int it);
		void q2v(int it);
		double CFL();
		void record(int ii);
		void capture(int ii);
		void clear();
		void fwrite_ary();
		void snapshot(int meas, int isum);
	private:
};
double pwfun(
	int typ, 	// piecewise constant(0), linear(1)
	double x, 	// x-coordinate
	double *xknot, double *yknot,  // discontinuous/broken  points
	int ndiv, // number of intervals
	double y0 // magnicifation
);


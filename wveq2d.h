#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"wave1d.h"
//--------------------------------------------------------------------------
class FIELD{
	public:
		int Ng[2];
		int *Ndiv;
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
		void fwrite(int num);
		void set_IC(double xc, double yc, double sig, double f0);
		void print_prms();
		void gen_indx0(int **kcell);
		void gen_indx1(int **kcell);
		void gen_indx2(int **kcell);
		int *kbnd, *kint;
		int Nin,Nbnd,Nex;
	private:
};
class SOURCE{
	public:
		double wd;
		int ng;	// number of grids
		void mem_alloc(int n);
		int *isrc,*jsrc,*ksrc;
		int type;	// 1:v1, 2:v2
		int nml;	// normal vector (nx or ny)
		void set(double x, double y, double w);
		void print(double *Xa, double *dx);
		int iwv;	// waveform number 
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
		void print();
		double **bwv;
		void mem_alloc(int n);
		void init_bwv(int n, double dtau);
		int Nt;
		double dt;
		void record(int it, double **fld);
		double *Xa, *dx;
		void set_cod(double *Xa, double *dx);
		void fwrite();
	private:
};
class DOMAIN{
	public:
		double *Xa, *Wd, *dx; 
		int **kcell;
		int *Ndiv,ndat;
		void init(int *Ndiv);
		void setup(double *Xa, double *Wd, double *dx);
		double xf[2];
		double *ij2xf(int i, int j, int type);
		void print_kcell();
		void print_prms();
		double cp;
		double dt;
		void perfo_ellip(char *fname);
		void out_kcell();
		void topography(char *fname);
		int find_v1bnd(int nx, double y);
		int find_v2bnd(int ny, double x);
		int nsrc;
	private:
};
class CNTRL{
	public:
		double cp,rho,almb;
		double Xa[2],Wd[2],dx[2];
		int Ndiv[2],Nt;
		double dt;
		DOMAIN dm;
		FIELD pr,v1,v2; 
		SOURCE *srcs;		
		RECVR *recs;
		Wv1D *wvs;
		int nwv;
		int nsrc;
		int nrec;
		void setup_domain(char *fname);
		void time_setting(char *fname);
		void wvfm_setting(char *fname);
		int src_setting(char *fname);
		int rec_setting(char *fname);
		void p2v(int it);
		void v2p(int it);
		double CFL();
		void record(int ii);
	private:
};
double pwfun(
	int typ, 	// piecewise constant(0), linear(1)
	double x, 	// x-coordinate
	double *xknot, double *yknot,  // discontinuous/broken  points
	int ndiv, // number of intervals
	double y0 // magnicifation
);


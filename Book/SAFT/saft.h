
void show_msg(char *fn);
/*
class Recs{
	public:
		int nrec; // number of receiver points
		double xr1[2]; // start
		double xr2[2]; // end --> equi-spaced nrec receivers will be generated
		void set_xr1(double x, double y);
		void set_xr2(double x, double y);
		void gen(int n);
		double *xr,*yr,*tf;
		double dxr[2];
		void print_cods();
		void set_tof(double val);
	private:
};
*/
class Array{
	public:
		int nele; // number of receiver points
		int nscan;
		double *xcod,*ycod;	// position
		double *delay,*A0;
		bool *actv;
		int *isrc;
		void load();
		void setup();
		double xin[2];
		double *set_transmitter(int iscan);
	private:
		void mem_alloc();
};

class Bscan{
	public:
		bool alloc;
		Bscan();
		char fname[128];
		int Ny;
		double y0,dy;
		int Nt;
		double t0,dt,Td;
		double **amp;
		void fwrite_Bscan();
		double get_amp(int iele, double tt);

		double *amp_sum;
		double stack_Ascans(double* delay);
		void fwrite_Ascan(char fname[128]);
		void load(int tnum);
		void mem_alloc();
		double tof;
		int ntof;
		double get_Asum(int inc_ntof);
	private:
};
/*
class Plate{
	public:
		double zb;	//bottom surface
		double zt;	//top surface
		double ht;	//plate thickness
		double zs[2];
		Plate(double z_bottom, double z_top);
		double mirror(double zf, int ud, int nref);
		double path(double xs[2],double xr[2], int ud, int nref);
		double TOF(double xs[2],double xr[2], int ud, int nref,double c);
		void TOFs(int ud, int nref,double c);
		void write_tx(char fn[128],char mode[2]);
		Recs rcv;	// receiver points
		double xs[2];	// source point
		void set_src(double x,double y);
	private:
};
*/

class IMG{
	public:
		double Xa[2],Xb[2];	// imaging area
		double dx[2];		// pixel size
		double Wd[2];		// physical image size
		int Ndiv[2],Ng[2]; 	// pixel size Ng, Ndiv=Ng-1
		int npnt;		// total number of the pixels
		double **V; 		// image intensity
		double *xcod,*ycod;
		IMG();			// contructor
		void set_Xa(double x, double y);	// set lower-left conrner
		void set_Wd(double W, double H);	// set width x height
		void set_Ng(int nx, int ny);		// set grid size 
		double get_xcod(int l);			// index --> x conversion
		double get_ycod(int l);			// index --> y conversion
		void set_V(int l,double val);		// set image intensity
		void fwrite_V(char fn[128]);		// export data to file
	private:
		void mem_alloc();		// memory allocation
};

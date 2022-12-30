
void show_msg(char *fn);
class Array{
	public:
		int nele; // number of receiver points
		int nscan;
		double *xcod,*ycod;	// position
		double *delay,*A0;
		bool *actv;
		int *isrc;
		void load();
		void setup(char path_name[128]);
		char path[128];
		double xin[2];
		double *set_transmitter(int iscan);
		void set_angles(double xc, double yc);
		double *thij;
		double th_min,th_max,dth;
		int *hist,nbin;
		int get_hist_val(int ktr);

	private:
		void mem_alloc();
};

class Bscan{
	public:
		bool alloc;
		Bscan(char path_name[128]);
		char path[128];
		char fname[128];
		int Ny;
		double y0,dy;
		int Nt;
		double t0,dt,Td;
		double **amp;
		void fwrite_Bscan();
		double get_amp(int iele, double tt);

		//double *amp_sum;
		//double stack_Ascans(double* delay);
		//void fwrite_Ascan(char fname[128]);
		void load(int tnum,int wv_typ);
		void mem_alloc();
		//double tof;
		//int ntof;
		//double get_Asum(int inc_ntof);
	private:
};

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

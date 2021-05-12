#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
//#include <math.h>
#ifndef __FFT__
	#define __FFT__
	#include "fft.h"
#endif

using namespace std;

class Wv1D{
	public:
		int Nt; // original data length
		int Np;	// FFT data length
		double *amp,*time,*tg;
		complex<double> *Amp;
		double *phi;
		double t1,t2,dt;
		double df;
		Wv1D();
		Wv1D(int Nt);
		Wv1D(char *fname);
		int load(char fname[128]);
		void print_info();
		bool mllc;
		char data_file[128];
		int FFT(int sign);
		void unwrap(double f0);
		void renew_phase();
		void out_Amp(char *fn);
		void out_Amp_polar(char *fn);
		void out_Amp(char *fn,int ofst);
		void out_amp(char *fn,char dlm);
		int fft_stat;
		double L2(double t1,double t2);
		double max(double t1,double t2);
		int arg_max(double t1,double t2);
		void Butterworth(double tb, double Tw_6dB);
		void Gauss(double tb, double sig);
		void Tri(double tb, double sig);
		void Sigmoid(double tb, double t90);
		double gdelay();
		int get_fnum(double ff);
		double get_f(int n);
	private:
	protected:
};
Wv1D corr(Wv1D wv1, Wv1D wv2, double *tmax, double *Amax);
class Grid{
	public:
		int Nx, Ny;
		double *Xcod, *Ycod;
		double *xcod, *ycod;
		double dx,dy;
		int *vals;
		int load(char *fname);
		char fname[128];
	private:
	protected:
};
class Array2D{
	public:
		double **A,*A2;
		int Nx,Ny,ndat;
		double dx[2];
		double Xa[2],Xb[2],Wd[2];
		Array2D();
		Array2D(int nx,int ny);
		void set_Xa(double x, double y);
		void set_Wd();
		void set_dx(double x, double y);
//		~Array2D();
		void out(char *fn);
	private:
	protected:
};
class Array3D{
	public:
		double ***A,**A2,*A3;
		int Nd[3];
		int Nx,Ny,Nz;
		int ndat;
		double dx[3];
		double Xa[3],Xb[3],Wd[3];
		Array3D(int nx,int ny, int nz);
		~Array3D();
		void set_Xa(double x, double y, double z);
		void set_Wd();
		void set_dx(double x, double y, double z);
		void set_val(double v);
		void print_dim();
		void load(char *dir_name);
		Wv1D awv;
		void CorrY();
		void CorrX();
		void Lp(int p); 	// p=0 for L_inf, p=2 for L2
		Array2D proj();
		void Butterworth(double cx, double cy);
		Array2D gdelay(double cy);
	private:
	protected:
};
//------------------------------------------------------------

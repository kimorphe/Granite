#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "waves.h"

using namespace std;

class Voigt{
	public:
		double f1,f2,df;
		double nf;	
		double alph,beta,s0;
		complex<double> *M;//Complex modulus
		double *s_mdl;// slowness (model)
		double *s_msd;	// slowness (measurment)
		void setup(double fs,double fe, int n);
		void eval0(double a, double b);
		double s0fit();
		double cost();
	private:
	protected:
};
void Voigt::setup(double fs, double fe, int n){
	f1=fs;
	f2=fe;
	nf=n;
	df=(f2-f1)/(nf-1);
	s_mdl= (double *)malloc(sizeof(double)*nf);
	s_msd= (double *)malloc(sizeof(double)*nf);
	M= (complex<double> *)malloc(sizeof(complex<double>)*nf);
};

double Voigt::s0fit(){
	int i;
	double xx=0.0,xy=0.0;
	for(i=0;i<nf;i++){
		xx+=(s_mdl[i]*s_mdl[i]);
		xy+=(s_mdl[i]*s_msd[i]);
	}
	s0=xy/xx;
	return(s0);
};

double Voigt::cost(){
	int i;
	double err;
	double L2=0.0;
	for(i=0;i<nf;i++){
		err=s_msd[i]-s_mdl[i]*s0;
		L2+=(err*err);
	};
	return(L2);
};

void Voigt::eval0(double a, double b){
	int i;
	double PI2=8.0*atan(1.0);
	complex<double> zi=complex<double>(0.0,1.0);
	double omg,Mb,dlt;
	alph=a;
	beta=b;
	for(i=0;i<nf;i++){
		omg=PI2*(f1+df*i);
		M[i]=1.0+pow((-zi*omg),alph)*beta;
		dlt=arg(M[i]);
		Mb=fabs(abs(M[i]));
		s_mdl[i]=cos(dlt*0.5)/Mb;
	}
};
int main(){
	Wv1D awv,arf;

	char fnref[128]="../1MHznew.csv";	// Reference signal
	char fname[128]="../Quartz/scope_423.csv"; // Transmitted waveform

	arf.load(fnref);	// load refrence signal
	arf.print_info();

	awv.load(fname);	// load transmitted waveform
	awv.print_info();


	char fout[128];
	double tb,sig;
	sprintf(fout,"ref0.dat");
	arf.out_amp(fout,' ');	// export raw data
	arf.Sigmoid(11.0,1.0);	// trucation by Sigmoid function
	arf.Gauss(11.8,0.5); 	// apply Gaussian window
	sprintf(fout,"ref1.dat");
	arf.out_amp(fout,' '); // export waveform data
	arf.FFT(1);		// perform FFT


	tb=12.5;sig=0.5;	// set Window parameter for the A-scan
	sprintf(fout,"wv0.dat");
	awv.out_amp(fout,' ');	// export raw data
	awv.Sigmoid(11.5,1.0);	// trunction by Sigmoid function
	awv.Gauss(tb,sig);	// applly Gaussian window
	sprintf(fout,"wv1.dat");
	awv.out_amp(fout,' ');	// export waveform data
	sprintf(fout,"wv1.fft");
	awv.FFT(1);	// perform FFT

	int i;
	for(i=0;i<awv.Nt;i++) awv.Amp[i]/=arf.Amp[i];	// Deconvolution
	awv.renew_phase();	// evaluate the phase spectrum 
	awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward
	awv.out_Amp_polar(fout);// export frequency spectrum

	FILE *fp=fopen("cp.dat","w");
	double omg,cp;
	double PI=4.0*atan(1.0);
	double ht=3.42;	// plate thickness = travel distance [mm]
	for(i=0;i<awv.Np;i++){
		omg=awv.df*PI*2*i;	// angular frequecy
		cp=omg*ht/awv.phi[i];	// phase velocity
		fprintf(fp,"%lf, %lf, %lf\n",i*awv.df,cp,1./cp);
	};
	fclose(fp);

	double f1,f2;	// frequency range
	int nf1,nf2,nf;
	f1=0.5;
	f2=2.0;
	nf1=awv.get_fnum(f1);
	nf2=awv.get_fnum(f2);
	f1=awv.get_f(nf1);
	f2=awv.get_f(nf2);
	nf=nf2-nf1+1;

	printf("f1,f2, df=%lf %lf %lf nf=%d\n",f1,f2,awv.df,nf);
	Voigt vgt;
	vgt.setup(f1,f2,nf);
	vgt.eval0(1.0,0.1);
	for(i=nf1;i<=nf2;i++){
		omg=awv.df*PI*2*i;	// angular frequecy
		vgt.s_msd[i-nf1]=awv.phi[i]/(omg*ht);
	};
	printf("s0=%lf\n",vgt.s0fit());
	printf("L2==%lf\n",vgt.cost());


	fp=fopen("sp.dat","w");
	for(i=0;i<nf;i++){
		fprintf(fp,"%lf, %lf, %lf\n",(i+nf1)*awv.df,vgt.s_msd[i],vgt.s_mdl[i]*vgt.s0);
	}
	fclose(fp);

	return(0);
};

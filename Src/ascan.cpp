#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "waves.h"
#include "voigt.h"

using namespace std;

int main(){
	Wv1D awv,arf;

	char fnref[128]="../1MHznew.csv";	// Reference signal
	char fname[128]="../Quartz/scope_336.csv"; // Transmitted waveform
	//char fname[128]="../Na_Feldspar/scope_13.csv"; // Transmitted waveform
	//char fname[128]="../Na_Feldspar/scope_123.csv"; // Transmitted waveform

	// load refrence signal
	arf.load(fnref);
	arf.print_info();
	// load transmitted waveform
	awv.load(fname);
	awv.print_info();

	char fout[128];
	double tb,sig;
	sprintf(fout,"ref0.dat");
	arf.out_amp(fout,' ');	// export raw data
	arf.Sigmoid(11.0,1.0);	// trucation by Sigmoid function
	//arf.Gauss(11.8,0.5); 	// apply Gaussian window
	arf.Tri(11.8,1.0); 	// apply Gaussian window
	sprintf(fout,"ref1.dat");
	arf.out_amp(fout,' '); // export windowed waveform 
	arf.FFT(1);		// perform FFT
	char fdbg[128]="Aref.out";
//	arf.out_Amp(fdbg,0);
	arf.out_Amp_polar(fdbg);

	tb=12.5;sig=0.5;	// set Window parameter for the A-scan
	sprintf(fout,"wv0.dat");
	awv.out_amp(fout,' ');	// export raw data
	awv.Sigmoid(11.5,1.0);	// trunction by Sigmoid function
	//awv.Gauss(tb,sig);	// applly Gaussian window
	awv.Tri(12.5,1.0); 	// apply Gaussian window
	sprintf(fout,"wv1.dat");
	awv.out_amp(fout,' ');	// export windowed waveform 
	sprintf(fout,"wv1.fft");
	awv.FFT(1);	// perform FFT

	int i,j;
	for(i=0;i<awv.Nt;i++) awv.Amp[i]/=arf.Amp[i];	// Deconvolution
	awv.renew_phase();	// evaluate the phase spectrum 
	awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward
	awv.out_Amp_polar(fout);// export frequency spectrum

	// Frequency Range
	double f1,f2,omg;	
	int nf1,nf2,nf;
	f1=0.6; f2=1.5; 
	nf1=awv.get_fnum(f1);
	nf2=awv.get_fnum(f2);
	f1=awv.get_f(nf1);
	f2=awv.get_f(nf2);
	nf=nf2-nf1+1;

	Voigt vgt;	// Voigt Complex Modulus (Class)
	vgt.setup(f1,f2,nf);	// set frequency range
	double PI2=8.0*atan(1.0);
	double ht=3.42;	// plate thickness = travel distance [mm]
	for(i=nf1;i<=nf2;i++){
		omg=awv.df*PI2*i;	// angular frequecy
		vgt.s_msd[i-nf1]=awv.phi[i]/(omg*ht); // set measured slowness(real part of s)
	};
	double a,b;
	vgt.linfit_cp(&a,&b);

	FILE *fp=fopen("sp.dat","w");
	double alph,beta,s1,s2,s3,cfit;
	double alph1,alph2,nalph,dalph;
	alph1=0.1;
	alph2=1.0;
	nalph=10;
	dalph=(alph2-alph1)/(nalph-1);
	for(i=0;i<nalph;i++){
		alph=alph1+i*dalph;
		printf("alph=%lf\n",alph);
		beta=vgt.argmin_beta(alph);
		vgt.eval0(alph,beta);
		for(j=0;j<nf;j++){
			s1=vgt.s_msd[j];	// Re[s] (measured)
			s2=vgt.s_mdl[j]*vgt.s0;	// Re[s] (slowness)
			s3=vgt.a_mdl[j]*vgt.s0;	// Im[s] (decay)
			cfit=a*(vgt.f1+vgt.df*j)+b;
			fprintf(fp,"%lf, %lf, %lf, %lf, %lf, %lf %lf\n",(j+nf1)*awv.df,s1,s2,1./s1,1./s2,s3,cfit);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);


	return(0);
};

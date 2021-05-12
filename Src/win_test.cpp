#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "waves.h"

using namespace std;

int main(){
	Wv1D awv,arf;

	char fnref[128]="../1MHznew.csv";	// Reference signal
	//char fname[128]="../Quartz/scope_300.csv"; // Transmitted waveform
	//char fname[128]="../Na_Feldspar/scope_13.csv"; // Transmitted waveform
	char fname[128]="../Na_Feldspar/scope_153.csv"; // Transmitted waveform


	double sig,dsig,sig1;
	double freq,time;
	int nsig;
	int i,j,imax;

	sig1=0.5;
	dsig=0.1;
	nsig=6;
	double tb_ref=12.0;
	double tb=12.8;
	FILE *fp=fopen("wv_win.dat","w");
	FILE *fp2=fopen("wv_win.fft","w");
	for(j=0;j<nsig;j++){	
		arf.load(fnref);
		//arf.print_info();

		awv.load(fname);
		//awv.print_info();

		if(j==0){
		imax=arf.arg_max(0.0,30.0);
		printf("tmax,amax=%lf, %lf\n",arf.dt*imax+arf.t1,arf.amp[imax]);
		imax=awv.arg_max(0.0,30.0);
		printf("tmax,amax=%lf, %lf\n",awv.dt*imax+awv.t1,awv.amp[imax]);
		}

		//arf.Sigmoid(11.0,1.0);	// trucation by Sigmoid function
		//arf.Gauss(11.8,0.5); 	// apply Gaussian window

		sig=sig1+dsig*j;
		arf.Tri(tb_ref,sig); 	// apply Gaussian window
		arf.FFT(1);		// perform FFT


		//awv.Sigmoid(11.5,1.0);	// trunction by Sigmoid function
		//awv.Gauss(tb,sig);	// applly Gaussian window
		awv.Tri(tb,sig); 	// apply Gaussian window
		awv.FFT(1);	// perform FFT
		for(i=0;i<arf.Nt;i++){
			time=arf.t1+arf.dt*i;
			fprintf(fp,"%lf, %lf, %lf\n",time,arf.amp[i],awv.amp[i]);
		};
		fprintf(fp,"\n");

		for(i=0;i<awv.Nt;i++) awv.Amp[i]/=arf.Amp[i];	// Deconvolution
		awv.renew_phase();	// evaluate the phase spectrum 
		awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward
		for(int i=0;i<awv.Np;i++){
			freq=i*awv.df;
			fprintf(fp2,"%le %le %le\n",freq,abs(awv.Amp[i]),awv.phi[i]);
		};
		fprintf(fp2,"\n");

	}
	fclose(fp);
	fclose(fp2);
	return(0);
};

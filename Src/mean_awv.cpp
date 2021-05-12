#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "waves.h"

int main(){
	Wv1D awv1,awv2;
	char fname[128],fout[128];

	sprintf(fname,"amp_mean_Qt.dat");
	awv1.load(fname);	// Mean waveform

	sprintf(fname,"../Quartz/scope_33.csv");
	awv2.load(fname);	// raw measured waveform

	sprintf(fout,"wv2.t");
	awv2.out_amp(fout,' ');


	//double tb=12.5;	// Window parameter 1 
	double sig=0.6; // Window parameter 2
	awv1.Gauss(12.5,sig);	// Apply Gaussian window 
	awv2.Gauss(12.5,sig);	// 

	sprintf(fout,"wv1.t");
	awv1.out_amp(fout,' ');	// export A-scan
	sprintf(fout,"wv2.t");
	awv2.out_amp(fout,' ');	// export A-scan
	
	// Fourier Transform & Unwraping

	awv1.FFT(1);		// perform FFT
	//awv1.unwrap(0.5);	// unwrap phase from 0.5MHz up and downward

	awv2.FFT(1);
	int i;
	for(i=0;i<awv1.Nt;i++) awv2.Amp[i]/=awv1.Amp[i];
	awv2.renew_phase();
	awv2.unwrap(0.5);

	sprintf(fout,"wv1.f");
	awv1.out_Amp_polar(fout);
	sprintf(fout,"wv2.f");
	awv2.out_Amp_polar(fout);

	return(0);
};

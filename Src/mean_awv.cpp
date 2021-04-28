#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "waves.h"

int main(){
	Wv1D awv,awv2;
	char fname[128],fout[128];

	sprintf(fname,"amp_mean_Qt.dat");
	awv.load(fname);

	sprintf(fname,"../Quartz/scope_320.csv");
	awv2.load(fname);

	sprintf(fout,"wv2.t");
	awv2.out_amp(fout,' ');


	double tb=12.5;
	double sig=5.0;
	awv.Gauss(tb,sig);

	sprintf(fout,"wv.t");
	awv.out_amp(fout,' ');	// export a-scan
	awv.FFT(1);		// perform FFT
	sprintf(fout,"wv.f");
	awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward
	awv.out_Amp_polar(fout);
	return(0);
};

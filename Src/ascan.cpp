#include <stdio.h>
#include <stdlib.h>
#include "waves.h"

using namespace std;
int main(){
	Wv1D awv,arf;

	char fnref[128]="../1MHznew.csv";
	char fname[128]="../Qt_scope_33.csv";

	arf.load(fnref);
	arf.print_info();

	awv.load(fname);
	awv.print_info();


	char fout[128];
	double tb,sig;

	sprintf(fout,"ref0.dat");
	arf.out_amp(fout,' ');
	arf.Gauss(11.8,0.5);

	sprintf(fout,"ref1.dat");
	arf.out_amp(fout,' ');
	arf.FFT(1);


	tb=12.5;sig=0.5;

	sprintf(fout,"wv0.dat");
	awv.out_amp(fout,' ');
	awv.Gauss(tb,sig);
	sprintf(fout,"wv1.dat");
	awv.out_amp(fout,' ');

	sprintf(fout,"wv1.fft");
	awv.FFT(1);

	int i;
	for(i=0;i<awv.Nt;i++){
		awv.Amp[i]/=arf.Amp[i];
	};
	awv.renew_phase();
	awv.unwrap(0.5);
	awv.out_Amp_polar(fout);

	FILE *fp=fopen("cp.dat","w");
	double omg,cp;
	double PI=4.0*atan(1.0);
	double ht=3.42;	//[mm]
	for(i=0;i<awv.Np;i++){
		omg=awv.df*PI*2*i;
		cp=omg*ht/awv.phi[i];
		fprintf(fp,"%lf, %lf\n",i*awv.df,cp);
	};
	fclose(fp);

	return(0);
};

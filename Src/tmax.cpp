#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "waves.h"

using namespace std;

class data_dir{
	public:
		char M[3]; 
		int nfile;
		char dir_name[128];
		int nchar;
		void set_dir(char M[3]);
		void set_file_name(int ii);
		double rho;
		char fname[128];
	private:
	protected:
};
void data_dir::set_dir(char M[3]){

	if(strcmp(M,"Qt")==0){
		nfile=589;
		rho=2.65;
		sprintf(dir_name,"%s","../Quartz/");
		nchar=strlen(dir_name);
	}else if(strcmp(M,"K")==0){
		nfile=824;
		rho=2.56;
		sprintf(dir_name,"%s","../K_Feldspar/");
		nchar=strlen(dir_name);
	}else if(strcmp(M,"Na")==0){
		nfile=548;
		rho=2.62;
		sprintf(dir_name,"%s","../Na_Feldspar/");
		nchar=strlen(dir_name);
	}else{
		printf("Data folder for %s cannot be found !\n",M);
		printf(" --> process terminated.\n");
		exit(-1);
	};
}
void data_dir::set_file_name(int num){
	char tmp[128];
	strcpy(fname,dir_name);
	sprintf(tmp,"scope_%d.csv",num);
	strcat(fname,tmp);
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(){

	double PI2=8.0*atan(1.0);
	int i,j,imax;
	double freq,time;
	double sp,cp,omg;

	char fname[128];

	char M[3]="K";	// {Qt,Na,K}
	char fnref[128]="../1MHznew.csv";	// Reference signal

	double tb_ref=11.8;
	double tb=12.45;
	double sig=1.0;
	sig=1.0;

	double ht=3.42;	// plate thickness [mm]

	double cmin=3.0;
	double cmax=8.0;
	int ibin,nbin=50;
	double dc=(cmax-cmin)/nbin;
	int *cp_hist=(int *)malloc(sizeof(int)*nbin);
	for(i=0;i<nbin;i++) cp_hist[i]=0;

	sprintf(fname,"tmax_%s.out",M);
	FILE *fout=fopen(fname,"w");

	sprintf(fname,"phi_%s.out",M);
	FILE *fp2=fopen(fname,"w");

	data_dir ddr;
	Wv1D awv,arf;
	ddr.set_dir(M);

	// Reference Signal
	arf.load(fnref);
	arf.Tri(tb_ref,sig); 	// apply Gaussian window
	arf.FFT(1);		// perform FFT

	double *cp_mean=(double *)malloc(sizeof(double)*arf.Np);
	for(i=0;i<arf.Np;i++) cp_mean[i]=0.0;

	fprintf(fout,"# %d\n",ddr.nfile);
	int count=0;
	int nf1=arf.get_fnum(1.0);
	int nf2=arf.get_fnum(2.0);

	for(j=0;j<ddr.nfile;j++){
		ddr.set_file_name(j);
		awv.load(ddr.fname);
		imax=awv.arg_max(0.0,20.0);
		time=awv.dt*imax+awv.t1;
		fprintf(fout,"%d, %lf, %lf\n",j,time,awv.amp[imax]);

		if(awv.amp[imax]>0.0) continue;
		if(time>13.5) continue;
		count++;

		awv.Tri(tb,sig); 	// apply Gaussian window
		//awv.Gauss(tb,sig*0.5); 	// apply Gaussian window
		awv.FFT(1);	// perform FFT

		for(i=0;i<awv.Nt;i++) awv.Amp[i]/=arf.Amp[i];	// Deconvolution
		awv.renew_phase();	// evaluate the phase spectrum 
		awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward

		for(int i=0;i<awv.Np;i++){
			freq=i*awv.df;
			omg=freq*PI2;
			sp=awv.phi[i]/(omg*ht);
			cp=1./sp;
			fprintf(fp2,"%le %le %le, %le\n",freq,abs(awv.Amp[i]),awv.phi[i],cp);
			cp_mean[i]+=cp;

			if(i <nf1) continue;
			if(i >nf2) continue;
			ibin=int((cp-cmin)/dc+0.5);
			if(ibin<0) continue;
			if(ibin>=nbin) continue;
			cp_hist[ibin]++;
		};
		fprintf(fp2,"\n");
	}

	
	fclose(fp2);
	sprintf(fname,"cp_mean_%s.out",M);
	fp2=fopen(fname,"w");
	for(i=0;i<arf.Np;i++){
	       	freq=i*awv.df;
		fprintf(fp2,"%le, %le\n",freq,cp_mean[i]/=count);
	}
	fclose(fp2);

	sprintf(fname,"cp_hist_%s.out",M);
	fp2=fopen(fname,"w");
	double cbin;
	for(i=0;i<nbin;i++){
		cbin=cmin+dc*(i+0.5);
		fprintf(fp2,"%lf, %d\n",cbin,cp_hist[i]);
	};
	fclose(fp2);


	return(0);
};

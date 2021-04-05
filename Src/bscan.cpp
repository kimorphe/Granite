#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "waves.h"

using namespace std;

class Wv2D{
	public:
		int nfile;	// number of files
		int Nt;		// number of time steps
		int Np;		// number of FFT data points
		double df,dt,t1,t2;
		Wv1D awv;
		Wv1D arf;
		void set_data_dir(char M[3]);
		void print_setup();
		char dir_name[128];
		int nchar;
		double rho;	//mass density
		void load(char M[3]);
		double **amp;
		complex<double> **Amp;
		double **phi;
		void Gauss(double tb, double sig);
		void FFT_all();
		void write_Amp(int i, char fname[128]);
		void write_amp(int i, char fname[128]);
		void write_Amp_polar(int i, char fname[128]);
		void set_refwv(char fn[128],double tb, double sig);
		bool refwv_ready;
		Wv2D();
		void dcnv_all();
		void unwrap(double f0);
	private:
	protected:
};
Wv2D::Wv2D(){
	refwv_ready=false;
};
void Wv2D::set_data_dir(char M[3]){

	puts(M);
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
void Wv2D::load(char M[3]){

	int i;
	char fname[128],tmp[128];
	double *pt;
	complex<double> *zpt;

	Wv2D::set_data_dir(M);
	Wv2D::print_setup();

	sprintf(tmp,"scope_%d.csv",0);
	strcpy(fname,dir_name);
	strcat(fname,tmp);

	awv.load(fname);
	awv.FFT(1);
	Nt=awv.Nt;
	dt=awv.dt;
	t1=awv.t1;
	t2=awv.t2;

	df=awv.df;
	Np=awv.Np;
	
	amp=(double **)malloc(sizeof(double*)*nfile);
	pt=(double *)malloc(sizeof(double)*nfile*Nt);
	for(i=0;i<nfile;i++) amp[i]=pt+Nt*i;

	phi=(double **)malloc(sizeof(double*)*nfile);
	pt=(double *)malloc(sizeof(double)*nfile*Np);
	for(i=0;i<nfile;i++) phi[i]=pt+Np*i;

	Amp=(complex<double> **)malloc(sizeof(complex<double>*)*nfile);
	zpt=(complex<double> *)malloc(sizeof(complex<double>)*nfile*Np);
	for(i=0;i<nfile;i++) Amp[i]=zpt+Np*i;

	printf(" Reading waveform data ....");
	for(i=0;i<nfile;i++){
		awv.amp=amp[i];
		strcpy(fname,dir_name);
		sprintf(tmp,"scope_%d.csv",i);
		strcat(fname,tmp);
		awv.load(fname);
	};
	puts(" Done !");
};

void Wv2D::print_setup(){
	printf("Data Directory: %s\n",dir_name);
	printf("Number of files: %d\n",nfile);
	printf("Mass density[g/cm3]: %lf\n",rho);
};
void Wv2D::FFT_all(){
	int i;
	printf(" Performign FFT ....");
	for(i=0;i<nfile;i++){
		awv.amp=amp[i];
		awv.phi=phi[i];
		awv.Amp=Amp[i];
		awv.fft_stat=1;	// avoid multiply allocating arrays
		awv.FFT(1);
	};
	puts(" Done !");
};
void Wv2D::Gauss(double tb,double sig){
	int i;
	for(i=0;i<nfile;i++){
	       	awv.amp=amp[i];
		awv.Gauss(tb,sig);
	};
}
void Wv2D::write_Amp(int i, char fname[128]){
	FILE *fp=fopen(fname,"w");
	complex<double> *Z;
	Z=Amp[i];
	for(int j=0;j<Np;j++){
		fprintf(fp,"%le, %le, %le, %le\n",j*df,Z[j].real(),Z[j].imag(),abs(Z[j]));
	}
	fclose(fp);
};
void Wv2D::write_Amp_polar(int i, char fname[128]){
	FILE *fp=fopen(fname,"w");
	complex<double> *Z;
	Z=Amp[i];
	for(int j=0;j<Np;j++){
		fprintf(fp,"%le, %le, %le\n",j*df,abs(Z[j]),phi[i][j]);
	}
	fclose(fp);
};
void Wv2D::write_amp(int i, char fname[128]){
	FILE *fp=fopen(fname,"w");
	for(int j=0;j<Nt;j++) fprintf(fp,"%lf, %lf\n",t1+dt*j,amp[i][j]);
	fclose(fp);
};
void Wv2D::set_refwv(char fname[128],double tb,double sig){
	arf.load(fname);
	arf.Gauss(tb,sig);
	arf.FFT(1);
	refwv_ready=true;
};
void Wv2D::dcnv_all(){
	if(!refwv_ready){
		printf(" Reference signal has not been specified!\n");
		printf(" --> process terminated\n");
		exit(-1);
	};
	int i,j;
	for(i=0;i<nfile;i++){
	for(j=0;j<Np;j++){
		Amp[i][j]/=arf.Amp[j];
		phi[i][j]=arg(Amp[i][j]);
	}
	}
};
void Wv2D::unwrap(double f0){
	int i;
	for(i=0;i<nfile;i++){
		awv.Amp=Amp[i];
		awv.phi=phi[i];
		awv.unwrap(f0);
	}
};

int main(){

	Wv2D bwv;
	double tb,sig;
	char fnref[128]="../1MHznew.csv";
	char M[3]; 

	sprintf(M,"%s","K");
	sprintf(M,"%s","Na");
	sprintf(M,"%s","Qt");

	bwv.load(M);
	tb=11.8, sig=0.5;
	bwv.set_refwv(fnref,tb,sig);
	bwv.dcnv_all();

	char fout[128]="awv.dat";
	bwv.write_amp(33,fout);

	tb=12.5, sig=0.5;
	bwv.Gauss(tb,sig);
	sprintf(fout,"awv_win.dat");
	bwv.write_amp(33,fout);

	bwv.FFT_all();
	bwv.dcnv_all();
	bwv.unwrap(0.5);
	sprintf(fout,"awv.fft");
	bwv.write_Amp_polar(33,fout);

	int i,j;
	double omg,freq;
	double ht=3.42,cp;
	double PI=4.0*atan(1.0);

	freq=1.0;
	j=int(freq/bwv.df);
	omg=bwv.df*j*2.*PI;

	for(i=0;i<bwv.nfile;i++){
		cp=omg*ht/bwv.phi[i][j];
		printf("%d %lf\n",i,cp);
	};
	printf("\n");

	freq=1.2;
	j=int(freq/bwv.df);
	omg=bwv.df*j*2.*PI;
	for(i=0;i<bwv.nfile;i++){
		cp=omg*ht/bwv.phi[i][j];
		printf("%d %lf\n",i,cp);
	};

};

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "waves.h"

using namespace std;


class Hist2D{
	public:
		int **Ht;
		double range[2];
		int nbin;
		int ndat[2];
		int axis;
		void init(int dsize[2], int ax, int nb);
		void gen(double **data, double y1, double y2);
		int Nx,Ny;
	private:
	protected:
};

void Hist2D::init(int dsize[2], int ax, int nb){
	axis=ax;
	nbin=nb;
	ndat[0]=dsize[0];
	ndat[1]=dsize[1];

	if(axis==0){
		Nx=nbin;
		Ny=ndat[1];
	}
	if(axis==1){
		Nx=ndat[0];
		Ny=nbin;
	}

	Ht=(int **)malloc(sizeof(int *)*Nx);
	int *pt=(int *)malloc(sizeof(int)*Nx*Ny);
	for(int i=0;i<Nx;i++) Ht[i]=pt+i*Ny;
};

void Hist2D::gen(double **data, double y1, double y2){
	range[0]=y1;
	range[1]=y2;
	double dy=(y2-y1)/nbin;
	int i,j;

	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		Ht[i][j]=0;
	}
	}

	int l,ibin,nsum;
	if(axis==0){
		for(j=0;j<Ny;j++){
			nsum=0;
		for(l=0;l<ndat[0];l++){
			ibin=int((data[l][j]-y1)/dy);
			if(ibin<0) continue;
			if(ibin>=nbin) continue;
			Ht[ibin][j]++;
			nsum++;
		}
		}
	};

	if(axis==1){
		for(i=0;i<Nx;i++){
		for(l=0;l<ndat[1];l++){
			ibin=int((data[i][l]-y1)/dy);
			if(ibin<0) continue;
			if(ibin>=nbin) continue;
			Ht[i][ibin]++;
		}
		}
	};
};

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
		double **cp;
		void Gauss(double tb, double sig);
		void Sigmoid(double tb,double t90);
		void FFT_all();
		void write_Amp(int i, char fname[128]);
		void write_amp(int i, char fname[128]);
		void write_Amp_polar(int i, char fname[128]);
		void set_refwv(char fn[128],double tb, double sig);
		bool refwv_ready;
		Wv2D();
		void dcnv_all();
		void unwrap(double f0);
		double ht;
		void PhaseVel(double h);
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

	int Np2=Np/2;
	cp=(double **)malloc(sizeof(double*)*nfile);
	pt=(double *)malloc(sizeof(double)*nfile*(Np2));
	for(i=0;i<nfile;i++) cp[i]=pt+Np2*i;

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
void Wv2D::Sigmoid(double tb,double t90){
	int i;
	for(i=0;i<nfile;i++){
	       	awv.amp=amp[i];
		awv.Sigmoid(tb,t90);
	};
}
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
	arf.Sigmoid(11.0,1.0);	// trucation by Sigmoid function
	arf.Gauss(12.5,3.0);
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
void Wv2D::PhaseVel(double h){
	ht=h;
	int i,j,Np2=Np/2;
	double PI2=8.0*atan(1.0);
	double freq,omg;
	
	for(i=0;i<nfile;i++){
	for(j=0;j<Np2;j++){
		freq=df*j;
		omg=freq*PI2;
		cp[i][j]=omg*ht/phi[i][j];
	}
	}
	
};
//-----------------------------------------------------------------
int main(){

	Wv2D bwv;
	double tb,sig;
	char fnref[128]="../1MHznew.csv";
	char M[3]; 

	sprintf(M,"%s","Qt");	// chose mineral type
	sprintf(M,"%s","Na");
	sprintf(M,"%s","K");

	bwv.load(M);			// Load waveform data (B-scan)
	tb=11.8, sig=0.5;		// Gaussian window parameter
	bwv.set_refwv(fnref,tb,sig);	// load reference data and apply Gaussian window
//	bwv.dcnv_all();

//	char fout[128]="awv.dat";
//	bwv.write_amp(33,fout);

	//bwv.Sigmoid(11.5,1.0);		// Truncation by Sigmoid function
	bwv.Sigmoid(11.0,1.0);		// Truncation by Sigmoid function
	bwv.Gauss(12.5,3.0);		// Apply Guassian window
	tb=12.5, sig=0.5;
	bwv.Gauss(tb,sig);		// Apply Guassian window
//	sprintf(fout,"awv_win.dat");
//	bwv.write_amp(33,fout);

	bwv.FFT_all();	// Peform FFT 
	bwv.dcnv_all();	// Deconvolution
	double f0=0.5;	// [MHz]
	bwv.unwrap(0.5);	// Unwrap phase spectrum from f0 [MHz] 
//	sprintf(fout,"awv.fft");
//	bwv.write_Amp_polar(33,fout);

	int i,j;
	double omg,freq;
	double ht=3.42;	// sample thicnkess [mm]
		bwv.PhaseVel(ht);
	double cp;
	double PI=4.0*atan(1.0);
	double f1=0.5;	// [MHz]
	double f2=1.5;	// [MHz]
	int j1=int(f1/bwv.df);
	int j2=int(f2/bwv.df);
	FILE *fp=fopen("cp.dat","w");
	fprintf(fp,"#f1, f2, nf\n");
	f1=bwv.df*j1;
	f2=bwv.df*j2;
	fprintf(fp,"%lf, %lf, %d\n",f1,f2,j2-j1+1);
	fprintf(fp,"#number of waveforms\n");
	fprintf(fp,"%d\n",bwv.nfile);
	fprintf(fp,"# cp(w)[km/s] as a function of frequency\n");
	for(i=0;i<bwv.nfile;i++){
	for(j=j1;j<=j2;j++){
		freq=bwv.df*j;
		omg=freq*2.*PI;
		//fprintf(fp,"%lf, %lf\n",freq,bwv.cp[i][j]);
		fprintf(fp,"%lf\n",bwv.cp[i][j]);
	}
	}


	// Histogram
	Hist2D ht2;
	double cp1=3.0;
	double cp2=8.0;
	int nbin=30;
	double dcp=(cp2-cp1)/nbin;
	int dsize[2];
	dsize[0]=bwv.nfile;
	dsize[1]=int(bwv.Np/2);
	ht2.init(dsize,0,nbin);
	ht2.gen(bwv.cp,cp1,cp2);


	FILE *fh=fopen("cp.hist","w");
	double x;
	fprintf(fh,"#f1, f2, nf\n");
	fprintf(fh,"%lf, %lf, %d\n",f1,f2,j2-j1+1);
	fprintf(fh,"#cp1, cp2, nbin\n");
	fprintf(fh,"%lf, %lf, %d\n",cp1,cp2,nbin);
	fprintf(fh,"# count\n");
	for(j=j1;j<=j2;j++){
	for(i=0;i<nbin;i++){
		fprintf(fh,"%d\n",ht2.Ht[i][j]);
	}
	}
	fclose(fh);
};

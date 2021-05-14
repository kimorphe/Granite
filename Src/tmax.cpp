#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "waves.h"

using namespace std;

class data_dir{
	public:
		char mnrl[3]; 	// {Qt,Na,K}
		int mid;	// 0=Qt, 1=Na, 2=K
		int nfile;
		char dir_name[128];
		int nchar;
		void set_dir(char M[3]);
		void set_dir(int MID);
		void set_file_name(int ii);
		double rho;
		char fname[128];
		void show();
	private:
	protected:
};
void data_dir::show(){
	puts("---------------------------------");
	printf("mnrl=%s, mid=%d\n",mnrl,mid);
	printf("nfile=%d\n",nfile);
	printf("dir_name=%s\n",dir_name);
};
void data_dir::set_dir(int MID){
	mid=MID;
	if(mid==0) sprintf(mnrl,"%s","Qt");
	if(mid==1) sprintf(mnrl,"%s","Na");
	if(mid==2) sprintf(mnrl,"%s","K");
	data_dir::set_dir(mnrl);

};
void data_dir::set_dir(char M[3]){

	strcpy(mnrl,M);
	//if(strcmp(M,"Qt")==0){
	if(strcmp(mnrl,"Qt")==0){
		nfile=589;
		rho=2.65;
		sprintf(dir_name,"%s","../Quartz/");
		nchar=strlen(dir_name);
		mid=0;
	}else if(strcmp(mnrl,"K")==0){
		nfile=824;
		rho=2.56;
		sprintf(dir_name,"%s","../K_Feldspar/");
		nchar=strlen(dir_name);
		mid=2;
	}else if(strcmp(mnrl,"Na")==0){
		nfile=548;
		rho=2.62;
		sprintf(dir_name,"%s","../Na_Feldspar/");
		nchar=strlen(dir_name);
		mid=1;
	}else{
		printf("Data folder for %s cannot be found !\n",M);
		printf(" --> process terminated.\n");
		exit(-1);
		mid=-1;
	};
}
void data_dir::set_file_name(int num){
	char tmp[128];
	strcpy(fname,dir_name);
	sprintf(tmp,"scope_%d.csv",num);
	strcat(fname,tmp);
};

class Hist1D{
	public:
		double cmin;
		double cmax;
		int nbin;
		double dc;
		double cbin(int k);
		int *hist;
		void setup(double c1, double c2, int n);
		int count(double data);
		void clear();
	private:
	protected:
};
void Hist1D::setup(double c1,double c2, int n){
	cmin=c1;
	cmax=c2;
	nbin=n;
	dc=(cmax-cmin)/(nbin-1);
	hist=(int *)malloc(sizeof(int)*nbin);
	for(int i=0;i<nbin;i++) hist[i]=0;
};
void Hist1D::clear(){
	for(int i=0;i<nbin;i++) hist[i]=0;
};
double Hist1D::cbin(int k){
	return(cmin+dc*k);
};
int Hist1D::count(double data){
	int ibin=int((data-cmin)/dc+0.5);
	if(ibin<0) return(0);
	if(ibin>=nbin) return(0);
	hist[ibin]++;
	return(1);

};

class Hist2D{
	public:
		double xmin,ymin;
		double xmax,ymax;
		int nxbin,nybin,nbin;
		double dx,dy;
		double xbin(int k);
		double ybin(int k);
		int **hist;
		void set_xaxis(double x1, double x2, int nx);
		void set_yaxis(double y1, double y2, int yx);
		void mem_alloc();
		int count(double X, double Y);
		void show_axes();
		void clear();
		void print_header(FILE *fp);
};
void Hist2D::set_xaxis(double x1, double x2, int nx){
	xmin=x1;
	xmax=x2;
	nxbin=nx;
	dx=(xmax-xmin)/(nxbin-1);
};
void Hist2D::set_yaxis(double y1, double y2, int ny){
	ymin=y1;
	ymax=y2;
	nybin=ny;
	dy=(ymax-ymin)/(nybin-1);

};
double Hist2D::xbin(int k){
	return(xmin+(k)*dx);
};
double Hist2D::ybin(int k){
	return(ymin+(k)*dy);
};
void Hist2D::show_axes(){
	printf("x1,x2,nxbin,dx=%lf %lf %d %lf\n",xmin,xmax,nxbin,dx);
	printf("y1,y2,nybin,dy=%lf %lf %d %lf\n",ymin,ymax,nybin,dy);
};
void Hist2D::print_header(FILE *fp){
	fprintf(fp,"#f1, f2 [MHz], nf\n");
	fprintf(fp,"#%lf, %lf, %d\n",xmin,xmax,nxbin);
	fprintf(fp,"#c1, c2 [km/s], nf\n");
	fprintf(fp,"#%lf, %lf, %d\n",ymin,ymax,nybin);
	fprintf(fp,"#count\n");
};
void Hist2D::mem_alloc(){
	nbin=nxbin*nybin;
	int *p=(int *)malloc(sizeof(int)*nbin);
	hist=(int **)malloc(sizeof(int*)*nxbin);

	int i;
	for(i=0;i<nbin;i++) p[i]=0;
	for(i=0;i<nxbin;i++) hist[i]=p+nybin*i;
};
void Hist2D::clear(){
	int i,j;
	for(i=0;i<nxbin;i++){
	for(j=0;j<nybin;j++){
	     	hist[i][j]=0;
	}
	}
};
int Hist2D::count(double X, double Y){
	int ibin=int((X-xmin)/dx+0.5);
	int jbin=int((Y-ymin)/dy+0.5);

	if(ibin<0) return(0);
	if(jbin<0) return(0);
	if(ibin>=nxbin) return(0);
	if(jbin>=nybin) return(0);

	hist[ibin][jbin]++;
	return(1);
};

void linfit(double *x, double *y, int n,double a[2]){
	double x1,x2,xy,y1;
	x1=0.0; 
	x2=0.0;
	y1=0.0; 
	xy=0.0;

	int i;
	for(i=0;i<n;i++){
		x1+=x[i];
		y1+=y[i];
		x2+=(x[i]*x[i]);
		xy+=(x[i]*y[i]);
	};
	x1/=n;
	x2/=n;
	y1/=n;
	xy/=n;
	double D=x2-x1*x1;
	a[0]=(xy-x1*y1)/D;
	a[1]=(x2*y1-x1*xy)/D;
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(){

	double PI2=8.0*atan(1.0);
	int i,j,imax;
	double freq,time;
	double sp,cp,omg;

	char fname[128];

	char M[3];	// {Qt,Na,K}
	double tb_ref=11.8;
	double tb=12.45;
	double sig=0.6;

	double dtb;

	char fnref[128]="../1MHznew.csv"; // Reference signal

	Wv1D awv,arf;
	// Reference Signal
	arf.load(fnref);
	//arf.Tri(tb_ref,sig); 	// apply Gaussian window
	arf.Hann(tb_ref,sig); 	// apply Gaussian window
	arf.FFT(1);		// perform FFT
	double *cp_mean=(double *)malloc(sizeof(double)*arf.Np);
	int *cp_count=(int *)malloc(sizeof(int)*arf.Np);

	double ht=3.42;	// plate thickness [mm]

	double cmin=3.0;
	double cmax=8.0;
	int nbin=50;
	double fmin=1.0;	// [MHz]
	double fmax=2.0;	// [MHz]
	int nf1=arf.get_fnum(fmin);
	int nf2=arf.get_fnum(fmax);
	int nfbin=nf2-nf1+1;
	fmin=arf.get_f(nf1);
	fmax=arf.get_f(nf2);

	double *cpw=(double *)malloc(sizeof(double)*arf.Np);

	Hist1D hst1,hst_tmax;
	Hist2D hst2;


	int ntbin=35;
	hst_tmax.setup(12.0,13.5,ntbin);

	hst1.setup(cmin,cmax,nbin);
	hst2.set_yaxis(cmin,cmax,nbin);
	hst2.set_xaxis(fmin,fmax,nfbin);
	hst2.mem_alloc();

	FILE *fout, *fp2;

	data_dir ddr;

	double P[2];
	int count,isum;

	//ddr.set_dir(M);
	for(int idir=0; idir<3;idir++){	// Data Directories
		ddr.set_dir(idir);
		ddr.show();

		for(i=0;i<arf.Np;i++){
		       	cp_mean[i]=0.0;
		       	cp_count[i]=0;
		}

		sprintf(fname,"tmax_%s.out",ddr.mnrl);
		fout=fopen(fname,"w");
		fprintf(fout,"# %d\n",ddr.nfile);

		sprintf(fname,"phi_%s.out",ddr.mnrl);
		fp2=fopen(fname,"w");

		double tmax_ave=0.0;
		double amax_ave=0.0;
		count=0;
		hst_tmax.clear();
		for(j=0;j<ddr.nfile;j++){
			ddr.set_file_name(j);
			awv.load(ddr.fname);
			//imax=awv.arg_max(0.0,20.0);
			imax=awv.arg_max(0.0,13.5);
			time=awv.dt*imax+awv.t1;
			tmax_ave+=time;
			amax_ave+=awv.amp[imax];
			count++;
			hst_tmax.count(time);
		};
		tmax_ave/=count;
		amax_ave/=count;
		printf("tmax=%lf, amax=%lf\n",tmax_ave,amax_ave);

		count=0;
		hst1.clear();
		hst2.clear();
		for(j=0;j<ddr.nfile;j++){
			ddr.set_file_name(j);
			awv.load(ddr.fname);
			//imax=awv.arg_max(0.0,20.0);
			imax=awv.arg_max(0.0,13.5);
			time=awv.dt*imax+awv.t1;
			fprintf(fout,"%d, %lf, %lf\n",j,time,awv.amp[imax]);

			if(awv.amp[imax]>0.0) continue;
			if(time>13.5) continue;
			count++;

			dtb=time-tmax_ave;
			//awv.Tri(tb+dtb,sig); 	// apply Gaussian window
			awv.Hann(tb+dtb,sig); 	// apply Gaussian window
			//awv.Gauss(tb,sig*0.5); 	// apply Gaussian window
			awv.FFT(1);	// perform FFT

			for(i=0;i<awv.Nt;i++) awv.Amp[i]/=arf.Amp[i];	// Deconvolution
			awv.renew_phase();	// evaluate the phase spectrum 
			awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward

			isum=0;
			for(int i=0;i<awv.Np;i++){
				freq=i*awv.df;
				omg=freq*PI2;
				sp=awv.phi[i]/(omg*ht);
				cp=1./sp;
				fprintf(fp2,"%le %le %le, %le\n",freq,abs(awv.Amp[i]),awv.phi[i],cp);
				if(cp >0.0 && cp <10.0){
					cp_mean[i]+=cp;
					cp_count[i]++;
				}
				cpw[i]=cp;
	
				if(i <nf1) continue;
				if(i >nf2) continue;
				hst1.count(cp);
				hst2.count(freq,cp);
			};
			fprintf(fp2,"\n");
//			linfit(awv.freq+nf1,cpw+nf1,nfbin,P);
//			printf("P=%lf %lf\n",P[0],P[1]);
		}

	
		fclose(fout);
		fclose(fp2);
		sprintf(fname,"cp_mean_%s.out",ddr.mnrl);
		fp2=fopen(fname,"w");
		//for(i=0;i<arf.Np;i++){
		for(i=nf1;i<=nf2;i++){
		       	freq=i*awv.df;
			fprintf(fp2,"%le, %le\n",freq,cp_mean[i]/=cp_count[i]);
		}
		fclose(fp2);

		sprintf(fname,"cp_hist_%s.out",ddr.mnrl);
		fp2=fopen(fname,"w");
		for(i=0;i<nbin;i++){
			fprintf(fp2,"%lf, %d\n",hst1.cbin(i),hst1.hist[i]);
		};
		fclose(fp2);

		sprintf(fname,"tmax_hist_%s.out",ddr.mnrl);
		fp2=fopen(fname,"w");
		for(i=0;i<ntbin;i++){
			fprintf(fp2,"%lf, %d\n",hst_tmax.cbin(i),hst_tmax.hist[i]);
		};
		fclose(fp2);

		sprintf(fname,"cp_hist2d_%s.out",ddr.mnrl);
		fp2=fopen(fname,"w");
		hst2.print_header(fp2);
		for(j=0;j<nfbin;j++){
		for(i=0;i<nbin;i++){
			fprintf(fp2,"%lf, %d\n",hst2.ybin(i),hst2.hist[j][i]);
		}
			fprintf(fp2,"#\n");
		}

	}
	return(0);
};

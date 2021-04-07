#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#ifndef __WVS__
	#define __WVS__
	#include "waves.h"
#endif

using namespace std;

//------------------------------------------------------------
Wv1D::Wv1D(){
	amp=0;
	time=0;
	mllc=false;
	fft_stat=0;
	Nt=0;
	dt=0.0;
	t1=0.0;t2=1.0;
	sprintf(data_file,"not_specified");
};
Wv1D::Wv1D(int ndat){
	Nt=ndat;
	amp=(double *)malloc(sizeof(double)*Nt);
	time=(double *)malloc(sizeof(double)*Nt);
	mllc=true;
	fft_stat=0;
	strcpy(data_file,"");
	t1=0.0;
	t2=0.0;
	dt=0.0;
};
Wv1D::Wv1D(char *fname){
	amp=0;
	time=0;
	mllc=false;
	fft_stat=0;
	Wv1D::load(fname);
};
void Wv1D::print_info(){
	printf("------wave data parameters-------\n");
	printf("File: %s\n",data_file);
	printf("(t1,t2)=%lf,%lf\n",t1,t2);
	printf("dt=%lf\n",dt);
	printf("Nt=%d\n",Nt);
	printf("---------------------------------\n");
}
int Wv1D::load(char fname[128]){

	FILE *fp=fopen(fname,"r");
	char cbff[128];

	if(fp==NULL){	// check if sucessfully opened
		printf("File %s not found!\n",fname);
		printf(" --> process terminated.");
		exit(-1);
	}
	Nt=0;
	while(fgets(cbff,128,fp)!=NULL){
		Nt++;	// count number of lines
	}	
	//printf("Nt=%d\n",Nt);
	fclose(fp);

	if(!mllc){	// allocate memory
		amp=(double *)malloc(sizeof(double)*Nt);
		time=(double *)malloc(sizeof(double)*Nt);
		mllc=true;
	}

	fp=fopen(fname,"r");
	strcpy(data_file,fname);
	double sum=0.0;
	for(int i=0;i<Nt;i++){
		fscanf(fp,"%lf, %lf\n",time+i,amp+i);
		time[i]*=1.e06;	// to micro sec
		sum+=amp[i];
	};
	fclose(fp);
	sum/=Nt;
	for(int i=0;i<Nt;i++) amp[i]-=sum;
	t1=time[0];
	t2=time[Nt-1];
	dt=time[1]-time[0];
}
int Wv1D::FFT(int isgn){
	int p=ceil(log2(Nt));
	Np=pow(2,p);

	if(fft_stat==0){
		Amp=(complex<double> *)malloc(sizeof(complex<double>)*Np);
		phi=(double *)malloc(sizeof(double)*Np);
	}
	if(isgn==1){
		for(int i=0;i<Nt;i++) Amp[i]=complex<double>(amp[i],0.0);
		for(int i=Nt;i<Np;i++) Amp[i]=complex<double>(0.0,0.0);
	}

	fft(Amp,Np,isgn);
	fft_stat=isgn;

	int i;
	for(i=0;i<Np;i++) phi[i]=arg(Amp[i]); 
	df=1./dt/Np;
	return(Np);
};
void Wv1D::renew_phase(){
	int i;
	for(i=0;i<Np;i++) phi[i]=arg(Amp[i]);
};
void Wv1D::unwrap(double f0){	// unwraping from f=f0
	int i,nf=int(f0/df),nwrap;
	double dphi;
	double phi2,phi1;
	double PI2=8.0*atan(1.0);
	double aPI2=0.9*PI2;

	phi1=phi[nf];
	nwrap=0;
	for(i=nf+1;i<Np;i++){
		phi2=phi[i];
		dphi=phi2-phi1;
		if(dphi>aPI2) nwrap--;
		if(dphi<-aPI2) nwrap++;
		phi[i]+=(PI2*nwrap);
		phi1=phi2;
	};

	phi2=phi[nf];
	nwrap=0;
	for(i=nf-1;i>=0;i--){
		phi1=phi[i];
		dphi=phi2-phi1;
		if(dphi>aPI2) nwrap++;
		if(dphi<-aPI2) nwrap--;
		phi[i]+=(PI2*nwrap);
		phi2=phi1;
	}
};
void Wv1D::out_Amp(char *fn){
	FILE *fp=fopen(fn,"w");
	df=1./dt/Np;
	double freq;
	for(int i=0;i<Np;i++){
		freq=i*df;
		fprintf(fp,"%le %le %le\n",freq,Amp[i].real(),Amp[i].imag());
	}
	fclose(fp);
};
void Wv1D::out_Amp_polar(char *fn){
	FILE *fp=fopen(fn,"w");
	df=1./dt/Np;
	double freq;
	for(int i=0;i<Np;i++){
		freq=i*df;
		//fprintf(fp,"%le %le %le\n",freq,abs(Amp[i]),arg(Amp[i]));
		fprintf(fp,"%le %le %le\n",freq,abs(Amp[i]),phi[i]);
	}
	fclose(fp);
};

double Wv1D::L2(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double sum=0.0;
	for(i=i1;i<=i2;i++) sum+=(amp[i]*amp[i]);
	return(sqrt(sum)/(i2-i1+1));
};
double Wv1D::max(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double amax=fabs(amp[i2]);
	for(i=i1;i<i2;i++){
	       	if(amax < fabs(amp[i])) amax=fabs(amp[i]);
	}
	return(amax);
};
void Wv1D::out_Amp(char *fn,int ofst){
	FILE *fp=fopen(fn,"w");
	int j;
	double xx,x1,dx;
	double df;
	if(fft_stat==1){
		df=1./dt/Np;
		dx=df;
		x1=0.0;
	}else{
		x1=t1;
		dx=dt;
	}
	for(int i=0;i<Np;i++){
		j=(ofst+i)%Np;
		xx=(i-ofst)*dx;
		fprintf(fp,"%le %le %le %le\n",xx,Amp[j].real(),Amp[j].imag(),abs(Amp[j]));
	}
	fclose(fp);
};
void Wv1D::out_amp(char *fn,char dlm){
	FILE *fp=fopen(fn,"w");
	for(int i=0;i<Nt;i++) fprintf(fp,"%le%c %le\n",i*dt+t1,dlm,amp[i]);

	fclose(fp);
};
void Wv1D::Sigmoid(double tb, double t90){
	double a=log(0.9/(1.-0.9))/t90;
	double Wt,tt;
	for(int i=0;i<Nt;i++){
		tt=time[i]-tb;
		Wt=1./(1.+exp(-a*tt));
		amp[i]*=Wt;
	}
};
void Wv1D::Gauss(double tb, double sig){
	double arg;
	int i;
	for(i=0;i<Nt;i++){
		arg=time[i]-tb;
		arg*=arg;
		amp[i]*=exp(-arg*0.5);
	}
};
void Wv1D::Butterworth(double tb, double Tw_6dB){

	int p=4;
	double tt;
	double t0=Tw_6dB*0.5;
	double arg;
	for(int i=0; i<Nt;i++){
		tt=t1+dt*i;
		arg=(tt-tb)/t0;
		arg=pow(arg,p);
		amp[i]/=(1.+arg);
	};
};
double Wv1D::gdelay(){
	if(fft_stat==0) FFT(1);
	tg=(double *)malloc(sizeof(double)*(Nt-1));
	double *phi=(double *)malloc(sizeof(double)*Nt);
	double *tmp=(double *)malloc(sizeof(double)*Nt);

	int i;
	double phi0=arg(Amp[0]),phi1;
	double PI2=8.0*atan(1.0);
	double tol=0.60;
	tol*=PI2;
	FILE *fp=fopen("tg.out","w");

	// phase spectrum
	double cost,sint,zlen;
	for(i=0;i<Nt;i++){
		zlen=abs(Amp[i]);
		cost=real(Amp[i])/zlen;	
		sint=imag(Amp[i])/zlen;	
		tmp[i]=acos(cost);
		if( sint < 0.0) tmp[i]=PI2-tmp[i];
		//phi[i]=arg(Amp[i]);
	};

	// unwrapping
	double dphi,ofst=0.0;
	//fprintf(fp,"%lf\n",phi[0]);
	phi[0]=tmp[0];
	for(i=1;i<Nt;i++){
		dphi=tmp[i]-tmp[i-1];
		phi[i]=tmp[i];
		if(dphi>tol){
			ofst-=PI2;
		}
		if(dphi<-tol){
			ofst+=PI2;
		}
		phi[i]+=ofst;
		//fprintf(fp,"%lf\n",phi[i]);
	};

	// group delay 
	double df=1./dt/Np;
	double Cf=PI2*df;
	double f1=0.05, f2=1.0;
	int nf1=int(f1/df), nf2=int(f2/df);
	for(i=0;i<Nt-1;i++){
		tg[i]=(phi[i+1]-phi[i])/Cf;
//		fprintf(fp,"%lf %lf %lf %lf\n",df*i,tg[i],1./tg[i],phi[i]);
	};
	fclose(fp);
	double tgb=0.0;
	for(i=nf1;i<=nf2;i++) tgb+=tg[i];
	return(tgb/(nf2-nf1+1));
};
//------------------------------------------------------------
Wv1D corr(Wv1D wv1, Wv1D wv2, double *tmax, double *Amax){
	wv1.FFT(1);
	wv2.FFT(1);

	Wv1D wv3(wv1.Np);
	wv3.FFT(1);
	double A1=0.0,A2=0.0;
	complex <double > Z1,Z2;
	Z1=complex<double>(0.0,0.0);
	Z2=complex<double>(0.0,0.0);
	wv3=wv1;
	for(int i=0;i<wv1.Np;i++){
		Z1=wv1.Amp[i];
		Z2=conj(wv2.Amp[i]);
		wv3.Amp[i]=Z1*Z2;
		A1+=abs(Z1*Z1);
		A2+=abs(Z2*Z2);
	};
	wv3.FFT(-1);
	A1=sqrt(A1*A2);
	int imax=0;
	double A0=abs(wv3.Amp[0].real()/A1);
	for(int i=0;i<wv1.Np;i++){
		wv3.Amp[i]/=A1;
		if(A0 < abs(wv3.Amp[i].real())){
			A0=abs(wv3.Amp[i].real());
			imax=i;
		}
	}
	double t0=imax*wv1.dt;
	if(imax>wv1.Np/2) t0=-(wv1.Np-imax)*wv1.dt;
	(*tmax)=t0;
	(*Amax)=A0;
	//printf("%d %le %le\n",imax,t0,A0);
	return(wv3);
};
int Wv1D::get_fnum(double ff){
	if(ff<0.0) return(0);
	int nf=int(ff/df+0.5);
	if(nf>=Np) return(Np-1);
	return(nf);
};
double Wv1D::get_f(int i){
	if(i<0) return(0.0);
	if(i>=Np) i=Np-1;
	return(i*df);
};
//---------------------------------------------------------
#if DB == 11 
int main(){
	char fname[128];
	sprintf(fname,"../W20H30_fine/scope_%d.csv",80);
	Wv1D awv(fname);
	//awv.Butterworth(16.0,10.0);
	awv.FFT(1);

	// Output
	sprintf(fname,"awvt.out");
	awv.out_amp(fname);
	sprintf(fname,"awvw.out");
	awv.out_Amp(fname,0);
	printf("tg=%lf\n",awv.gdelay());
	return(0);
};
#endif

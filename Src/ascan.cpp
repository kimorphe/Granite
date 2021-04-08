#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "waves.h"

using namespace std;

class Voigt{
	public:
		double f1,f2,df;
		double nf;	
		double alph,beta,s0;
		complex<double> *M;//Complex modulus
		double *s_mdl;// slowness (model, real part)
		double *a_mdl;// slowness (model, imaginary part)
		double *s_msd;	// slowness (measurment)
		void setup(double fs,double fe, int n);
		void eval0(double a, double b);
		double s0fit();
		double cost();
		double argmin_beta(double alpha);
	private:
	protected:
};
void Voigt::setup(double fs, double fe, int n){
	f1=fs;
	f2=fe;
	nf=n;
	df=(f2-f1)/(nf-1);
	s_mdl= (double *)malloc(sizeof(double)*nf);
	a_mdl= (double *)malloc(sizeof(double)*nf);
	s_msd= (double *)malloc(sizeof(double)*nf);
	M= (complex<double> *)malloc(sizeof(complex<double>)*nf);
};

double Voigt::s0fit(){
	int i;
	double xx=0.0,xy=0.0;
	for(i=0;i<nf;i++){
		xx+=(s_mdl[i]*s_mdl[i]);
		xy+=(s_mdl[i]*s_msd[i]);
	}
	s0=xy/xx;
	return(s0);
};

double Voigt::cost(){
	int i;
	double err;
	double L2=0.0;
	Voigt::s0fit();
	for(i=0;i<nf;i++){
		err=s_msd[i]-s_mdl[i]*s0;
		L2+=(err*err);
	};
	return(L2);
};

void Voigt::eval0(double a, double b){
	int i;
	double PI2=8.0*atan(1.0);
	complex<double> zi=complex<double>(0.0,1.0);
	double omg,Mb,dlt;
	alph=a;
	beta=b;
	for(i=0;i<nf;i++){
		omg=PI2*(f1+df*i);
		M[i]=1.0+pow((-zi*omg),alph)*beta;
		dlt=arg(M[i]);
		Mb=fabs(abs(M[i]));
		s_mdl[i]= cos(dlt*0.5)/Mb;
		a_mdl[i]=-sin(dlt*0.5)/Mb;
	}
};

double Voigt::argmin_beta(double alpha){
	double a,b,c,x;
	double fa,fb,fc,fx;
	double W=(3.-sqrt(5.0))*0.5;

	a=0.0;
	c=1.0;

	Voigt::eval0(alpha,a); fa=Voigt::cost();
	Voigt::eval0(alpha,c); fc=Voigt::cost();

	// Bracketing
	int i,imax=10,itr;
	double rt=0.1;
	b=c;
	itr=0;
	for(i=0;i<imax;i++){
		itr++;
		b*=rt;
		Voigt::eval0(alpha,b);
		fb=Voigt::cost();
		if(fb>=fa) continue;
		if(fb>=fc) continue;
		break;
	};
	printf("iteration=%d\n",itr);
	double bmin=a;
	double fmin=fa;
	if(itr==imax){
		printf("Braketing failed !");
		if(fc < fa){
		       	bmin=c;
			fmin=fc;
		}
		printf("b=%lf, fb=%lf\n",bmin,fmin);
		return(bmin);
	};

	// Refine approximation
	double xm;
	double a2,b2,c2;

	for(i=0;i<20;i++){
	xm=0.5*(a+c);
	if( b< xm){	// a < b < x < c
		x=c-W*(c-a);
		Voigt::eval0(alpha,x);
		fx=Voigt::cost();
		if(fx< fb){ // (b,x,c) --> (a,b,c)
			a=b; fa=fb;
			b=x; fb=fx;
		}else{ // (a,b,x) --> (a,b,c)
			c=x; fc=fx;
		}
	}else{	// a < x < b < c
		x=a+W*(c-a);
		Voigt::eval0(alpha,x);
		fx=Voigt::cost();
		if(fx < fb){ // (a,x,b) --> (a,b,c)
			c=b; fc=fb;
			b=x; fb=fx;
		}else{ // (x,b,c) --> (a,b,c)
			a=x; fa=fx;
		}
	};
		printf("itr=%d, b=%lf, fb=%lf",i,b,fb);
		a2=a*a, b2=b*b, c2=c*c;
		x=fa*(b2-c2)+fb*(c2-a2)+fc*(a2-b2);
		x=x/(fa*(b-c)+fb*(c-a)+fc*(a-b));
		x*=0.5;
		Voigt::eval0(alpha,x);
		fx=Voigt::cost();
		printf(" (x=%lf, fx=%lf)\n",x,fx);

	};
	bmin=x;
	if(fb<fx) bmin=b;
	return(bmin);

};
int main(){
	Wv1D awv,arf;

	char fnref[128]="../1MHznew.csv";	// Reference signal
	//char fname[128]="../Quartz/scope_223.csv"; // Transmitted waveform
	//char fname[128]="../Na_Feldspar/scope_223.csv"; // Transmitted waveform
	char fname[128]="../K_Feldspar/scope_323.csv"; // Transmitted waveform

	arf.load(fnref);	// load refrence signal
	arf.print_info();

	awv.load(fname);	// load transmitted waveform
	awv.print_info();


	char fout[128];
	double tb,sig;
	sprintf(fout,"ref0.dat");
	arf.out_amp(fout,' ');	// export raw data
	arf.Sigmoid(11.0,1.0);	// trucation by Sigmoid function
	arf.Gauss(11.8,0.5); 	// apply Gaussian window
	sprintf(fout,"ref1.dat");
	arf.out_amp(fout,' '); // export waveform data
	arf.FFT(1);		// perform FFT


	tb=12.5;sig=0.5;	// set Window parameter for the A-scan
	sprintf(fout,"wv0.dat");
	awv.out_amp(fout,' ');	// export raw data
	awv.Sigmoid(11.5,1.0);	// trunction by Sigmoid function
	awv.Gauss(tb,sig);	// applly Gaussian window
	sprintf(fout,"wv1.dat");
	awv.out_amp(fout,' ');	// export waveform data
	sprintf(fout,"wv1.fft");
	awv.FFT(1);	// perform FFT

	int i,j;
	for(i=0;i<awv.Nt;i++) awv.Amp[i]/=arf.Amp[i];	// Deconvolution
	awv.renew_phase();	// evaluate the phase spectrum 
	awv.unwrap(0.5);	// unwrap the phase from 0.5MHz up and down ward
	awv.out_Amp_polar(fout);// export frequency spectrum

	FILE *fp=fopen("cp.dat","w");
	double omg,cp;
	double PI=4.0*atan(1.0);
	double ht=3.42;	// plate thickness = travel distance [mm]
	for(i=0;i<awv.Np;i++){
		omg=awv.df*PI*2*i;	// angular frequecy
		cp=omg*ht/awv.phi[i];	// phase velocity
		fprintf(fp,"%lf, %lf, %lf\n",i*awv.df,cp,1./cp);
	};
	fclose(fp);

	double f1,f2;	// frequency range
	int nf1,nf2,nf;
	f1=0.6;
	f2=1.5;
	nf1=awv.get_fnum(f1);
	nf2=awv.get_fnum(f2);
	f1=awv.get_f(nf1);
	f2=awv.get_f(nf2);
	nf=nf2-nf1+1;

	printf("f1,f2, df=%lf %lf %lf nf=%d\n",f1,f2,awv.df,nf);
	Voigt vgt;
	vgt.setup(f1,f2,nf);
	vgt.eval0(1.0,0.1);
	for(i=nf1;i<=nf2;i++){
		omg=awv.df*PI*2*i;	// angular frequecy
		vgt.s_msd[i-nf1]=awv.phi[i]/(omg*ht);
	};
	printf("L2==%lf\n",vgt.cost());


	fp=fopen("sp.dat","w");

	double alph,beta,s1,s2,s3;
	for(i=0;i<6;i++){
		alph=0.5+i*0.1;
		printf("alph=%lf\n",alph);
		beta=vgt.argmin_beta(alph);
		vgt.eval0(alph,beta);
		for(j=0;j<nf;j++){
			s1=vgt.s_msd[j];
			s2=vgt.s_mdl[j]*vgt.s0;
			s3=vgt.a_mdl[j]*vgt.s0;
			fprintf(fp,"%lf, %lf, %lf, %lf, %lf, %lf\n",(j+nf1)*awv.df,s1,s2,1./s1,1./s2,s3);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	return(0);
};

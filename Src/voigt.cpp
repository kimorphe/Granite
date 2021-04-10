#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "voigt.h" 

using namespace std;

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

void Voigt::linfit_cp(double *a, double *b){
	int i;
	double x,y;
	double  a11,a12,a21,a22, b1,b2;
	
	a11=0.0;
	a12=0.0;
	a22=0.0;
	b1=0.0;
	b2=0.0;
	for(i=0;i<nf;i++){
		x=f1+df*i;
		y=1./s_msd[i];
		a11+=(x*x);
		a12+=x;
		b1+=(x*y);
		b2+=y;
	};

	a11/=nf;
	a12/=nf;
	a21=a12;
	a22=1.0;
	b1/=nf;
	b2/=nf;

	double D=a11*a22-a12*a21;
	*a=( a22*b1-a12*b2)/D;
	*b=(-a21*b1+a11*b2)/D;
	//double yb=0.5*(*a)*(f1+f2)+(*b);
	//printf("a=%lf, b=%lf mean=%lf %lf\n",*a,*b,b2,yb);
};

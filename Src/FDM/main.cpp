#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdm2d.h"
#include <random>

using namespace std;

double cp_gen(int num){
	static std::mt19937 eg1(-1);
	static std::mt19937 eg2(-1);
	static std::mt19937 eg3(-1);

	static std::normal_distribution<> Gss1(2.0,1.50);
	static std::normal_distribution<> Gss2(4.0,1.50);
	static std::normal_distribution<> Gss3(6.5,1.50);

	double cp;
	
	cp=0.0;
	if(num==0) cp=Gss1(eg1);
	if(num==1) cp=Gss2(eg2);
	if(num==2) cp=Gss3(eg3);

	return(cp);
};
int main(){
	std::random_device seed;
	std::mt19937 engine(-1);
	std::uniform_real_distribution<> urnd(0,1.0);
	std::normal_distribution<> Gss(1.0,1.0);


	char fgeom[124]="fdm.inp";
	char fout[124]="kcell.dat";
	Dom2D dom(fgeom);

	double xc[2],a,b;
	int i,j;
	a=0.3;
	b=0.3;
	dom.perfo_ellip(fgeom);
	int Np=400;
	for(i=0;i<Np;i++){
		xc[0]=urnd(engine)*dom.Wd[0];
		xc[1]=urnd(engine)*dom.Wd[1];
		//dom.perfo_ellip(xc,a,b,i%3);
		dom.perfo_ellip(xc,a,b,i);
	};

	double *cp0;

	cp0=(double *)malloc(sizeof(double)*Np);
	for(i=0;i<Np;i++) cp0[i]=cp_gen(i%3);


	int typ;
	for(i=0;i<dom.Ndiv[0];i++){
	for(j=0;j<dom.Ndiv[1];j++){
		typ=dom.kcell[i][j];
		dom.cp[i][j]=cp0[typ];
		dom.kcell[i][j]=typ%3;
	}
	}
	dom.out_cp();
	dom.out_kcell();

	return(0);
};

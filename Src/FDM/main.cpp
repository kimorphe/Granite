#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdm2d.h"
#include <random>

using namespace std;
int main(){
	std::random_device seed;
	std::mt19937 engine(-1);
	std::uniform_real_distribution<> urnd(0,1.0);
	std::normal_distribution<> Gss(1.0,1.0);


	char fgeom[124]="fdm.inp";
	char fout[124]="kcell.dat";
	Dom2D dom(fgeom);

	double xc[2],a,b;
	int i;
	a=0.6;
	b=0.4;
	dom.perfo_ellip(fgeom);
	for(i=0;i<400;i++){
		xc[0]=urnd(engine)*dom.Wd[0];
		xc[1]=urnd(engine)*dom.Wd[1];
		dom.perfo_ellip(xc,a,b,i%4);
	};
	dom.out_kcell();
	return(0);
};

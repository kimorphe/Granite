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

	static std::normal_distribution<> Gss1(3.0,0.5);
	static std::normal_distribution<> Gss2(4.0,0.5);
	static std::normal_distribution<> Gss3(5.0,0.5);

	double cp;
	
	cp=0.0;
	if(num==0) cp=Gss1(eg1);
	if(num==1) cp=Gss2(eg2);
	if(num==2) cp=Gss3(eg3);

	return(cp);
};
int main(){
	//std::random_device seed;
	std::mt19937 engine(-1);
	std::uniform_real_distribution<> urnd(0,1.0);
	std::normal_distribution<> Gss(1.0,1.0);


	char finwv[128]="inwv0.inp"; // incident waveform 
	char fgeom[128]="fdm.inp";   // general input
	char fout[128]="kcell.dat";  // domain data output
	char frec[128]="recs.inp";  // receiver array


	Dom2D dom(fgeom);	// load general input & set computational domain 
	dom.set_wvfm(finwv);	// set excitation(incident) waveform
//	dom.set_recs();
	dom.set_rec_array(frec);

/*
	double xc[2],a,b;
	int i,j;

	// create simple heterogeneous medium
	a=3.0; b=3.0;
	int Np=400;
	for(i=0;i<Np;i++){
		xc[0]=urnd(engine)*dom.Wd[0];
		xc[1]=urnd(engine)*dom.Wd[1];
		dom.perfo_ellip(xc,a,b,i); // inclution No.i
	};

	double *cp0;
	cp0=(double *)malloc(sizeof(double)*Np); // phase velocity data
	for(i=0;i<Np;i++) cp0[i]=cp_gen(i%3);

	int typ;
	for(i=0;i<dom.Ndiv[0];i++){
	for(j=0;j<dom.Ndiv[1];j++){
		typ=dom.kcell[i][j];
		dom.cp[i][j]=cp0[typ]; // assign phase velocity to material grid
		dom.kcell[i][j]=typ%3; 
	}
	}
*/
	dom.out_cp();
	dom.out_kcell();
	dom.set_cplim();
	//dom.fd2.out(1);
	dom.CFL();

	char fnwv[128]="bwv.out";
	double amp;
	int i,idat=0;
	printf("Nt=%d\n",dom.Nt);
	for(i=0;i<dom.Nt;i++){
		dom.fd2.s2v();	// stress --> velocity
		dom.fd2.v2s();  // velocity --> stress
		if(dom.bc==1) dom.fd2.periodicBC(); // periodic BC
		amp=dom.inwv.get_amp(i); // get incident wave amplitude
		dom.fd2.apply_Bcon(-amp); // stress BC
		if(i%100==0){
			printf("i=%d (%d)\n",i,idat);
			dom.fd2.out(0,idat);
			dom.fd2.out(1,idat);
			idat++;
		}
//	 	 dom.fd2.rec();	// record waveform data
		dom.fd2.record();
	};
//	dom.fd2.write_bwvs(fnwv);
	dom.fd2.write_bwv_array();
	return(0);
};

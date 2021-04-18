#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

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
		double linfit_cp(double *a, double *b);
	private:
	protected:
};

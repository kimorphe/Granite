#include<stdio.h>
#include<stdlib.h>
#include<math.h>


double fun(double x){
	double PI=4.0*atan(1.0);
	return(sin(2.*PI*x));
};
int main(){
	double x1,x2,x3,xs;
	double y1,y2,y3;

	x1=0.0;
	x2=0.75;

	x3=0.5*(x1+x2);

	y1=fun(x1);
	y2=fun(x2);
	y3=fun(x3);



	return(0);
};

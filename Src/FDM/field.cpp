#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdm2d.h"

void Field::init(int n[2]){
	Ndiv[0]=n[0];
	Ndiv[1]=n[1];
	Nv[0]=Ndiv[0]-1;
	Nv[1]=Ndiv[1]-1;
	Field::mem_alloc();
};
void Field::s2v(){
	int i,j;
	double rh=rho*dh;
	double dSx,dSy;
	for(i=0;i<Nv[0];i++){
	for(j=0;j<Nv[1];j++){
		dSx=(s[i+1][j]-s[i][j]+s[i+1][j+1]-s[i][j+1])/rh*dt;
		dSy=(s[i][j+1]-s[i][j]+s[i+1][j+1]-s[i+1][j])/rh*dt;
		v1[i][j]+=dSx;
		v2[i][j]+=dSy;
	}
	}
};
void Field::v2s(){
	int i,j;
	double dvx,dvy;
	double lmb;
	for(i=1;i<Nv[0];i++){
	for(j=1;j<Nv[1];j++){
		lmb=cp[i][j];
		lmb*=lmb;
		lmb*=(dt/dh);
		dvx=(v1[i][j-1]-v1[i-1][j-1]+v1[i][j]-v1[i-1][j])*lmb;
		dvy=(v2[i-1][j]-v2[i-1][j-1]+v2[i][j]-v2[i][j-1])*lmb;
		s[i][j]+=(dvx+dvy);
	}
	}
};
void Field::apply_Bcon(double amp){
	int j;
	for(j=0;j<Ndiv[1];j++){
		s[0][j]=amp;
	};
};

void Field::periodicBC(){
	int i,Ny;
	double dvx,dvy;
	double lmb;
	Ny=Nv[1]-1;
	for(i=1;i<Nv[0];i++){
		lmb=cp[i][0];
		lmb*=lmb;
		lmb*=(dt/dh);
		dvx=(v1[i][Ny]-v1[i-1][Ny]+v1[i][0]-v1[i-1][0])*lmb;
		dvy=(v2[i-1][0]-v2[i-1][Ny]+v2[i][0]-v2[i][Ny])*lmb;
		s[i][0]+=(dvx+dvy);
		s[i][Nv[1]]=s[i][0];
	};
};
void Field::mem_alloc(){
	int ndat=Ndiv[0]*Ndiv[1];
	double *ptmp;
	int i;

	ptmp=(double *)malloc(sizeof(double)*ndat);
	s=(double **)malloc(sizeof(double *)*Ndiv[0]);
	for(i=0;i<ndat;i++) ptmp[i]=0.0;
	for(i=0;i<Ndiv[0];i++) s[i]=ptmp+i*Ndiv[1];

	ndat=Nv[0]*Nv[1];
	ptmp=(double *)malloc(sizeof(double)*ndat);
	v1=(double **)malloc(sizeof(double *)*Nv[0]);
	for(i=0;i<ndat;i++) ptmp[i]=0.0;
	for(i=0;i<Nv[0];i++) v1[i]=ptmp+i*Nv[1];

	ptmp=(double *)malloc(sizeof(double)*ndat);
	v2=(double **)malloc(sizeof(double *)*Nv[0]);
	for(i=0;i<ndat;i++) ptmp[i]=0.0;
	for(i=0;i<Nv[0];i++) v2[i]=ptmp+i*Nv[1];

};

void Field::out(int num){
	char fname[128];
	sprintf(fname,"v%d.out",num);
	FILE *fp=fopen(fname,"w");
	int i,j;
	fprintf(fp,"# Xa[2], Xb[2]\n");
	fprintf(fp,"%lf, %lf, %lf, %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);
	fprintf(fp,"# Ndiv[2]\n");
	fprintf(fp,"%d, %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# vx, vy\n");
	for(i=0; i<Nv[0];i++){
	for(j=0; j<Nv[1];j++){
		fprintf(fp,"%lf, %lf\n",v1[i][j],v2[i][j]);
	}
	}
	fclose(fp);
};

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
	double rh=2.*rho*dh;
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
	double dvx,dvy,rdt2h=rho*dt*0.5/dh;
	double lmb;
	for(i=1;i<Nv[0];i++){
	for(j=1;j<Nv[1];j++){
		lmb=cp[i][j];
		lmb*=lmb;
		lmb*=rdt2h;
		dvx=(v1[i][j-1]-v1[i-1][j-1]+v1[i][j]-v1[i-1][j]);
		dvy=(v2[i-1][j]-v2[i-1][j-1]+v2[i][j]-v2[i][j-1]);
		s[i][j]+=(lmb*(dvx+dvy));
	}
	}
//	printf("cp=%lf, dh=%lf, dt=%lf\n",cp[i-1][j-1],dh,dt);
};
void Field::apply_Bcon(double amp){
	int j;
	for(j=0;j<Ndiv[1];j++){
		s[0][j]=amp;
	};
};

void Field::periodicBC(){
	int i,Ny;
	double dvx,dvy,rdt2h=rho*dt*0.5/dh;
	double lmb;
	Ny=Nv[1]-1;
	for(i=1;i<Nv[0];i++){
		lmb=cp[i][0];
		lmb*=lmb;
		lmb*=rdt2h;
		dvx=(v1[i][Ny]-v1[i-1][Ny]+v1[i][0]-v1[i-1][0]);
		dvy=(v2[i-1][0]-v2[i-1][Ny]+v2[i][0]-v2[i][Ny]);
		s[i][0]+=(lmb*(dvx+dvy));
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
void Field::mem_alloc_wvs(int nrx,int nry, int nt){
	int ndat=nrx*nry*nt;
	swv=(double *)malloc(sizeof(double)*ndat);
	v1wv=(double *)malloc(sizeof(double)*ndat);
	v2wv=(double *)malloc(sizeof(double)*ndat);
	Nt=nt;

};
void Field::set_xrec(double xrec[2], int nr){
	nrx=nr;
	irecx=(int *)malloc(sizeof(int)*nrx);
	int i,ir;
	double xr,dX=0.0;
	if(nrx>1) dX=(xrec[1]-xrec[0])/(nrx-1);
	for(i=0;i<nrx;i++){
		xr=xrec[0]+dX*i;
		ir=int((xr-Xa[0])/dh);
		irecx[i]=ir;
	};
};
void Field::set_yrec(double yrec[2], int nr){
	nry=nr;
	irecy=(int *)malloc(sizeof(int)*nry);
	int i,ir;
	double yr,dY=0.0;
	if(nry>1) dY=(yrec[1]-yrec[0])/(nry-1);
	for(i=0;i<nry;i++){
		yr=yrec[0]+dY*i;
		ir=int((yr-Xa[1])/dh);
		irecy[i]=ir;
	};
};
void Field::rec(){
	int i,j,ir,jr;
	for(ir=0; ir<nrx; ir++){
		i=irecx[ir];
	for(jr=0; jr<nry; jr++){
		j=irecy[jr];
		v1wv[iwv]=v1[i][j];
		v2wv[iwv]=v2[i][j];
		swv[iwv]=s[i][j];
		iwv++;
	}
	}
};

void Field::write_bwvs(char *fn){
	FILE *fp=fopen(fn,"w");
	int i,j,k,nrec=nrx*nry;
	int iwv0=0;
	double x1,x2,y1,y2;
	x1=Xa[0]+(irecx[0]+0.5)*dh;
	x2=Xa[0]+(irecx[nrx-1]+0.5)*dh;

	y1=Xa[1]+(irecy[0]+0.5)*dh;
	y2=Xa[1]+(irecy[nry-1]+0.5)*dh;

	fprintf(fp,"# xrec[2], nrx\n");
	fprintf(fp,"%lf, %lf, %d\n",x1,x2,nrx);
	fprintf(fp,"# yrec[2], nry\n");
	fprintf(fp,"%lf, %lf, %d\n",y1,y2,nry);
	fprintf(fp,"# Nt, dt \n");
	fprintf(fp,"%d, %lf\n",Nt,dt);
	fprintf(fp,"# v1, v2, s (for ix< nrx (for iy < nry (for it <Nt) )))\n");
	for(i=0;i<nrx;i++){
	for(j=0;j<nry;j++){
		iwv=iwv0;
		for(k=0;k<Nt;k++){
			fprintf(fp,"%lf, %lf, %lf\n",v1wv[iwv],v2wv[iwv],swv[iwv]);
			iwv+=nrec;
		}
		iwv0++;
	}
	}
	fclose(fp);
};
void Field::init_rec_array(int nn){
	nary=nn;
	rary=(rec_array *)malloc(sizeof(rec_array)*nary);
};
void Field::record(){
	int i,j;
	int iu,iv,ix,iy,jdat;
//	for(i=0;i<nary;i++) printf("npnt=%d\n",rary[i].npnt);

	for(i=0;i<nary;i++){
		jdat=rary[i].cntr;
	for(j=0;j<rary[i].npnt;j++){
//		printf("jdat=%d\n",jdat);
		iu=rary[i].irecs[j];
		iv=rary[i].jrec;
		ix=iu;
		iy=iv;
		if(rary[i].idir==1){
			ix=iv;
			iy=iu;
		};
		rary[i].v1[jdat]=v1[ix][iy];
		rary[i].v2[jdat]=v2[ix][iy];
		rary[i].s[jdat]=s[ix][iy];
		jdat++;
	}
	rary[i].cntr=jdat;
	};
};
void Field::write_bwv_array(){
	int i;
	char fname[128];
	printf("Field::write_bwv_array::nary=%d\n",nary);
	for(i=0;i<nary;i++){
		sprintf(fname,"bwv%d.out",i);
		printf("%s\n",fname);
	       	rary[i].write2(fname);
	}
};

void Field::out(int type, int num){
	char fname[128];
	FILE *fp;
	int i,j;
	if(type==0){
		sprintf(fname,"s%d.out",num);
		fp=fopen(fname,"w");
		fprintf(fp,"# Xa[2], Xb[2]\n");
		fprintf(fp,"%lf, %lf, %lf, %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);
		fprintf(fp,"# Ndiv[2]\n");
		fprintf(fp,"%d, %d\n",Ndiv[0],Ndiv[1]);
		fprintf(fp,"# s\n");
		for(i=0; i<Ndiv[0];i++){
		for(j=0; j<Ndiv[1];j++){
			fprintf(fp,"%lf\n",s[i][j]);
		}
		}
		fclose(fp);
	}
	if(type==1){
		sprintf(fname,"v%d.out",num);
		fp=fopen(fname,"w");
		fprintf(fp,"# Xa[2], Xb[2]\n");
		fprintf(fp,"%lf, %lf, %lf, %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);
		fprintf(fp,"# Ndiv[2]\n");
		fprintf(fp,"%d, %d\n",Nv[0],Nv[1]);
		fprintf(fp,"# vx, vy\n");
		for(i=0; i<Nv[0];i++){
		for(j=0; j<Nv[1];j++){
			fprintf(fp,"%lf, %lf\n",v1[i][j],v2[i][j]);
		}
		}
		fclose(fp);
	}
};

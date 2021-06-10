#include<iostream> 
#include<stdio.h>
#include<math.h> 
#include<string.h>
#include<stdlib.h>
#include "fdm2d.h"
using namespace std;


//  ---------CONSTRUCTOR -----------
Dom2D::Dom2D(char *fname){ //Contructor

	int i,ndim=2;
	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");	
	if(fp==NULL){
		printf("File cannot found !\n");
		printf(" --> abort process\n");
		exit(-1);
	};

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf",Xa,Xa+1);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);

	dx[0]=(Xb[0]-Xa[0])/Ndiv[0];
	dx[1]=(Xb[1]-Xa[1])/Ndiv[1];

	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];

	fd2.init(Ndiv);
	fd2.dh=dx[0];
	fd2.rho=1.0;
	fd2.Xa=Xa;
	fd2.Xb=Xb;

	fclose(fp);

	mem_alloc();	// allocate & initialize kcell[i][j]

	fd2.cp=cp;
};
void Dom2D::set_wvfm(char fname[128]){
	inwv.setup(fname);
	Nt=inwv.Nt;
	dt=inwv.dt;
	fd2.dt=dt;
};

//  ----------- MEMORY ALLOCATION --------------
void Dom2D::mem_alloc(){
	int i,j,*ptmp;
	ptmp=(int *)malloc(sizeof(int)*Ndiv[0]*Ndiv[1]);
	kcell=(int **)malloc(sizeof(int*)*Ndiv[0]);

	for(i=0; i<Ndiv[0];i++){
		kcell[i]=(ptmp+(i*Ndiv[1]));
	}

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		kcell[i][j]=0;
	}}

	double *tmp=(double *)malloc(sizeof(double)*Ndiv[0]*Ndiv[1]);
	cp=(double **)malloc(sizeof(double*)*Ndiv[0]);
	for(i=0;i<Ndiv[0]*Ndiv[1];i++) tmp[i]=0.0;
	for(i=0; i<Ndiv[0];i++) cp[i]=tmp+(i*Ndiv[1]);
};

//  ----------- DOMAIN PERFORATION  --------------
void Dom2D::perfo_ellip(double xc[2], double a, double b,int num){
	int i,j;
	double xcod[2];
	double xx,yy;

	for(i=0;i<Ndiv[0];i++){
		xcod[0]=Xa[0]+dx[0]*(i+0.5);	
		for(j=0;j<Ndiv[1];j++){
			xcod[1]=Xa[1]+dx[1]*(j+0.5);	
			xx=(xcod[0]-xc[0])/a;
			yy=(xcod[1]-xc[1])/b;
			if(xx*xx+yy*yy <1.0) kcell[i][j]=num;
		}
	}
}
void Dom2D::perfo_ellip(char *fname){
	int i,j;
	double xcod[2],xc[2],a,b;
	double xx,yy;

	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo_ellip !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,8,fp) !=NULL){
		if(strcmp(cbff,"##Ellip")==0){
			fscanf(fp,"%lf %lf %lf %lf\n",xc,xc+1,&a,&b);
			printf("%lf %lf %lf %lf\n",xc[0],xc[1],a,b);

			for(i=0;i<Ndiv[0];i++){
				xcod[0]=Xa[0]+dx[0]*(i+0.5);	
			for(j=0;j<Ndiv[1];j++){
				xcod[1]=Xa[1]+dx[1]*(j+0.5);	
				xx=(xcod[0]-xc[0])/a;
				yy=(xcod[1]-xc[1])/b;
				if(xx*xx+yy*yy <1.0) kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
};
//  ----------- DOMAIN PERFORATION  --------------
void Dom2D::perfo(char *fname){
	int i,j;
	double xcod[2];
	bool io;
	Circ cdat;

	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Circ")==0){
			fscanf(fp,"%lf %lf %lf\n",cdat.xc,cdat.xc+1, &cdat.radi);

			for(i=0;i<Ndiv[0];i++){
				xcod[0]=Xa[0]+dx[0]*(i+0.5);	
			for(j=0;j<Ndiv[1];j++){
				xcod[1]=Xa[1]+dx[1]*(j+0.5);	
				io=cdat.isin(xcod);	
				if(io==true){
				 kcell[i][j]=1;
				}
			}
			}
		}
	};
	fclose(fp);
};
void Dom2D::polygon(char *fname){
	FILE *fp;
	char cbff[128];
	int i,j,k,np;
	double *x,*y;
	double x1[2],x2[2],tx[2],nx[1],xf[2];
	double xmin,ymin,xmax,ymax,r12,h;
	double tmp;
	int imin,imax,jmin,jmax;


	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::polygon !");
		puts(fname);
		puts("--> process terminated");
		exit(-1);
	}

	i=0;
	while(fgets(cbff,7,fp)!=NULL){
		if(strcmp(cbff,"##Poly")==0){
			printf("Flag ##Poly found !!\n");
			fscanf(fp,"%d\n",&np);	// number vertices
			x=(double *)malloc(sizeof(double)*np);
			y=(double *)malloc(sizeof(double)*np);

			for(j=0;j<np;j++){
				fscanf(fp,"%lf %lf\n",x+j,y+j);
				if(j==0){
				       	xmin=x[0]; ymin=y[0];
				       	xmax=x[0]; ymax=y[0];
				}else{
					if(x[j]< xmin) xmin=x[j];
					if(x[j]> xmax) xmax=x[j];
					if(y[j]< ymin) ymin=y[j];
					if(y[j]> ymax) ymax=y[j];
				};
			}
			printf("xmin=%lf, xmax=%lf\n",xmin,xmax);
			printf("ymin=%lf, ymax=%lf\n",ymin,ymax);

			imin=floor((xmin-Xa[0])/dx[0]);
			imax=floor((xmax-Xa[0])/dx[0]);
			jmin=floor((ymin-Xa[1])/dx[1]);
			jmax=floor((ymax-Xa[1])/dx[1]);

			if(imin <0) imin=0;
			if(jmin <0) jmin=0;
			if(imax>=Ndiv[0]) imax=Ndiv[0]-1;
			if(jmax>=Ndiv[1]) jmax=Ndiv[1]-1;
			printf("imin,jmin=%d %d\n",imin,jmin);
			printf("imax,jmax=%d %d\n",imax,jmax);

			for(i=imin; i<=imax; i++){
				xf[0]=Xa[0]+dx[0]*(i+0.5);
			for(j=jmin; j<=jmax; j++){
				xf[1]=Xa[1]+dx[1]*(j+0.5);
				int iin=0;
				for(k=0;k<np;k++){
					x1[0]=x[k%np];
					x2[0]=x[(k+1)%np];
					x1[1]=y[k%np];
					x2[1]=y[(k+1)%np];

					tx[0]=x2[0]-x1[0];
					tx[1]=x2[1]-x1[1];
					r12=sqrt(tx[0]*tx[0]+tx[1]*tx[1]);

					tx[0]/=r12;
					tx[1]/=r12;

					nx[0]=-tx[1];
					nx[1]= tx[0];

					h=(xf[0]-x1[0])*nx[0]+(xf[1]-x1[1])*nx[1];
					if(h>0.e0) iin++;
				}

				if(iin==np) kcell[i][j]=1;
			}
			}
		}
	};

};
void Dom2D::angled_slit(char *fname){
	double xc[2],e1[2],e2[2],xf[2];
	double a,b,alph,x1,x2;
	double xs[4],ys[4];
	double xmin,xmax,ymin,ymax;
	double pi=4.0*atan(1.0);
	int i,j,i1,i2;
	int imin,imax,jmin,jmax;
	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::slit !");
		puts(fname);
		exit(-1);
	}

	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Slit")==0){
			fscanf(fp,"%lf %lf\n",xc,xc+1);
			fscanf(fp,"%lf %lf %lf\n",&a,&b,&alph);
			alph=alph/180.*pi;

			e1[0]=cos(alph);
			e1[1]=sin(alph);

			e2[0]=-e1[1];
			e2[1]=e1[0];

			xs[0]=xc[0]+a*e1[0]+b*e2[0];
			ys[0]=xc[1]+a*e1[1]+b*e2[1];

			xs[1]=xc[0]-a*e1[0]+b*e2[0];
			ys[1]=xc[1]-a*e1[1]+b*e2[1];

			xs[2]=xc[0]-a*e1[0]-b*e2[0];
			ys[2]=xc[1]-a*e1[1]-b*e2[1];

			xs[3]=xc[0]+a*e1[0]-b*e2[0];
			ys[3]=xc[1]+a*e1[1]-b*e2[1];

			xmin=xs[0]; xmax=xs[0]; 
			ymin=ys[0]; ymax=ys[0]; 
			for(i=1;i<4;i++){
				if(xs[i]<xmin) xmin=xs[i];
				if(xs[i]>xmax) xmax=xs[i];
				if(ys[i]<ymin) ymin=ys[i];
				if(ys[i]>ymax) ymax=ys[i];
			}

			
			imin=floor((xmin-Xa[0])/dx[0]);
			imax=floor((xmax-Xa[0])/dx[0]);
			jmin=floor((ymin-Xa[1])/dx[1]);
			jmax=floor((ymax-Xa[1])/dx[1]);

			if(imin <0) imin=0;
			if(jmin <0) jmin=0;
			if(imax>=Ndiv[0]) imax=Ndiv[0]-1;
			if(jmax>=Ndiv[1]) jmax=Ndiv[1]-1;

			for(i=imin; i<=imax; i++){
				xf[0]=Xa[0]+dx[0]*i-xc[0];
			for(j=jmin; j<=jmax; j++){
				xf[1]=Xa[1]+dx[1]*j-xc[1];


				x1=xf[0]*e1[0]+xf[1]*e1[1];
				x2=xf[0]*e2[0]+xf[1]*e2[1];
				i1=0; i2=0;
				if(fabs(x1)<a) i1=1;
				if(fabs(x2)<b) i2=1;
				if(i1+i2==2) kcell[i][j]=1;
			}
			}
/*
			for(i=0;i<2;i++){
				ixs[i]=floor((xs[i]-Xa[i])/dx[i]);
				ixe[i]=floor((xe[i]-Xa[i])/dx[i]);
				if( ixs[i] <0) ixs[i]=0;
				if( ixe[i] >= Ndiv[i]) ixe[i]=Ndiv[i]-1;
			}

			for(i=ixs[0]; i<=ixe[0]; i++){
			for(j=ixs[1]; j<=ixe[1]; j++){
				 kcell[i][j]=1;
			}
			}
*/
		}
	};
	

}
void Dom2D::slit(char *fname){
	int i,j;
	double xs[2],xe[2];
	int ixs[2],ixe[2];

	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::slit !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Rect")==0){
			fscanf(fp,"%lf %lf\n",xs,xs+1);
			fscanf(fp,"%lf %lf\n",xe,xe+1);

			for(i=0;i<2;i++){
				ixs[i]=floor((xs[i]-Xa[i])/dx[i]);
				ixe[i]=floor((xe[i]-Xa[i])/dx[i]);
				if( ixs[i] <0) ixs[i]=0;
				if( ixe[i] >= Ndiv[i]) ixe[i]=Ndiv[i]-1;
			}

			for(i=ixs[0]; i<=ixe[0]; i++){
			for(j=ixs[1]; j<=ixe[1]; j++){
				 kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
}
//-------------EXPORT GEOMETRY DATA  ------------------
void Dom2D :: out_kcell(){

	FILE *fp=fopen("kcell.dat","w");
	int i,j;
	double Ya[2],Yb[2];


	fprintf(fp,"# xa[2], xb[2] (computational domian)\n");
	fprintf(fp," %lf %lf\n %lf %lf\n",xa[0],xa[1],xb[0],xb[1]);

	fprintf(fp,"# Xa[2], Xb[2] (physical domain) \n");
	fprintf(fp," %lf %lf\n %lf %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);

	fprintf(fp,"# Ya[2], Yb[2] (imaging area) \n");
	for(i=0;i<2;i++){
		Ya[i]=Xa[i]+dx[i]*(iYa[i]+0.5);
		Yb[i]=Xa[i]+dx[i]*(iYb[i]+0.5);
	}
	fprintf(fp,"%lf %lf\n %lf %lf\n",Ya[0],Ya[1],Yb[0],Yb[1]);

	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# kcell[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp, "%d\n",kcell[i][j]);
	}
	}

	fflush(fp);
};
void Dom2D:: set_cplim(){
	int i,j;
	cmin=cp[0][0];
	cmax=cmin;
	double cc;
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		cc=cp[i][j];
		if(cmax < cc) cmax=cc;
		if(cmin > cc) cmin=cc;
	}
	}
	printf("(cmax,cmin)=%lf, %lf\n",cmax,cmin);
};
void Dom2D :: out_cp(){

	FILE *fp=fopen("cp.dat","w");
	int i,j;


	fprintf(fp,"# Xa[2], Xb[2] (physical domain) \n");
	fprintf(fp," %lf %lf\n %lf %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);
	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# cp[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp, "%lf\n",cp[i][j]);
	}
	}
	fflush(fp);
	fclose(fp);
};

//-------------COURANT NUMBER  ------------------
void Dom2D :: CFL(){

	cL=cmax;
	cfl0=cL*dt/dx[0]*sqrt(2.0);
	printf("CFL=%lf\n",cfl0);
	printf("cmax=%lf\n",cL);
	printf("dt=%f dx=(%f, %f)\n",dt,dx[0],dx[1]);
}

//------------- GRID NUMBER  ------------------

//	SET NUMBER OF GRIDS OF TYPE itype to Nx[2]
//      (identical with the function Fld2D:: gridNum) 
void Dom2D ::gridNum(int ityp){
	switch(ityp){
	case 0:	// v1-grid
		Nx[0]=Ndiv[0]+1;
		Nx[1]=Ndiv[1];
		break;
	case 1:	// v2-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1]+1;
		break;
	case 2:	// s11
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1];
		break;
	case 3:	// s12-grid
		Nx[0]=Ndiv[0]+1;
		Nx[1]=Ndiv[1]+1;
		break;
	case 4: //s22-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1];
		break;
	};
	Ng=Nx[0]*Nx[1];
};

//------------- Circ CLASS MEMBER FUNCTION  ------------------
bool Circ :: isin(double *x){
	double dist;

	dist=(x[0]-xc[0])*(x[0]-xc[0]);
	dist+=(x[1]-xc[1])*(x[1]-xc[1]);
	dist=sqrt(dist);

	if(dist <= radi){
		return true;
	}{
		return false;
	}
};

//------------- Rayleigh Wave Phase Velocity --------------
double Rs(double sL, double sT, double s){

	double s2=s*s;
	double sL2=sL*sL, sT2=sT*sT;
	double R1,RL,RT;

	
	R1=(2.0*s2-sT2);
	R1*=R1;
	RL=sqrt(s2-sL2);
	RT=sqrt(s2-sT2);

	return(R1-4.0*s2*RL*RT);
};

double Rwave(double cL, double cT){
	double s0,cR;
	double sR,sL=1.0/cL, sT=1.0/cT; // Slownesses
	double tol=1.e-06;
	int it=0,it_max=100;
	double alph=1.5;
	double sa=sT, sb=alph*sT;
	double err,err0=fabs(Rs(sL,sT,sa));

	do{
		it++;
		s0=0.5*(sa+sb);
		err=Rs(sL,sT,s0);

		if( err < 0.0){
			sb=s0;
		}else{
			sa=s0;
		}
		if(it>= it_max) break;
	}while(fabs(err)/err0 > tol);

	sR=s0;
	cR=1.0/sR;
	printf("iteration=%d, relative error=%lf\n",it,fabs(err)/err0);

	return(cR);
};


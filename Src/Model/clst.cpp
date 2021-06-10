#include<stdio.h>
#include<stdlib.h>
#include<math.h>

class Pin{
	public:
		double x;
		double y;
		double xg,yg;
		double dist(int ix, int iy);
		void set(double a, double b);
		int npix;
		void init();
	private:
	protected:
};
void Pin::init(){
	npix=0;
	xg=0.0;
	yg=0.0;
	x=0.0;
	y=0.0;
};

void Pin::set(double a, double b){
	x=a;
	y=b;
};
double Pin::dist(int ix, int iy){
	double rx=(ix-x);
	double ry=(iy-y);
	return(sqrt(rx*rx+ry*ry));
};

class IMG{
	public:
		int **M;
		int Ndiv[2];
		int npix;
		int count();
		int nB,nQ,nK,nN;
		int *B,*Q,*K,*N;
		void load(char *fn);
		int ipix(int k);
		int jpix(int k);
	private:
	protected:
};
int IMG::ipix(int k){
	return(int(k/Ndiv[0]));
};
int IMG::jpix(int k){
	return(int(k%Ndiv[1]));
};
void IMG::load(char *fname){
	char cbff[128];
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("File %s cannot be found !\n",fname);
		printf("--> abort process.");
		exit(-1);
	};
	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",Ndiv,Ndiv+1);
	printf("Ndiv=%d %d\n",Ndiv[0],Ndiv[1]);
	fgets(cbff,128,fp);
	npix=Ndiv[0]*Ndiv[1];

	int *ptmp=(int *)malloc(sizeof(int)*npix);
	int i,j;
	for(i=0;i<npix;i++) ptmp[i]=0;
	M=(int **)malloc(sizeof(int *)*Ndiv[0]);
	for(i=0;i<Ndiv[0];i++) M[i]=ptmp+i*Ndiv[1];

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fscanf(fp,"%d\n",M[i]+j);
	}
	}
	fclose(fp);
};

int IMG::count(){
	int i,mnrl,ix,iy;
	nB=0; nQ=0; nK=0; nN=0;
	for(i=0;i<npix;i++){
		ix=i/Ndiv[0];
		iy=i%Ndiv[1];
		mnrl=M[ix][iy];
		if(mnrl==0) nB++;
		if(mnrl==1) nQ++;
		if(mnrl==2) nK++;
		if(mnrl==3) nN++;
	};
	puts("--------------------");
	printf("nB=%d\n",nB);
	printf("nQ=%d\n",nQ);
	printf("nK=%d\n",nK);
	printf("nN=%d\n",nN);
	printf("total=%d\n",nB+nQ+nK+nN);
	puts("--------------------");

	B=(int *)malloc(sizeof(int)*nB);
	Q=(int *)malloc(sizeof(int)*nQ);
	K=(int *)malloc(sizeof(int)*nK);
	N=(int *)malloc(sizeof(int)*nN);
	nB=0; nQ=0; nK=0; nN=0;
	for(i=0;i<npix;i++){
		ix=i/Ndiv[0];
		iy=i%Ndiv[1];
		mnrl=M[ix][iy];
		if(mnrl==0) B[nB++]=i;
		if(mnrl==1) Q[nQ++]=i;
		if(mnrl==2) K[nK++]=i;
		if(mnrl==3) N[nN++]=i;
	};
	return(npix);
};
int main(){
	char fname[128]="minmap.out";

	IMG im3;
	im3.load(fname);
	im3.count();


	int npx=10, npy=10;
	int npk=npx*npy;
	Pin *pks;
	pks=(Pin *)malloc(sizeof(Pin)*npk);

	int dnpx,dnpy;

	dnpx=int(im3.Ndiv[0]/npx);
	dnpy=int(im3.Ndiv[1]/npy);

	int i,j,k;
	int ipx,ipy;
	k=0;
	for(i=0;i<npx;i++){
		ipx=int((i+0.5)*dnpx);
	for(j=0;j<npy;j++){
		ipy=int((j+0.5)*dnpy);
		//printf("%d %d\n",ipx,ipy);
		pks[k].init();
		pks[k++].set(ipx,ipy);
	}
	}	

	double rmin,r;
	int imin;
	k=im3.Q[503];
	ipx=im3.ipix(k);
	ipy=im3.jpix(k);
	printf("ipx,ipy=%d %d\n",ipx,ipy);
	for(i=0;i<npk;i++){
		r=pks[i].dist((float)ipx,(float)ipy);
		//printf("r=%lf\n",r);
		if(i==0){
			rmin=r;
			imin=0;
		}
		if(rmin>r){
		       rmin=r;
		       imin=i;
		}
	}
	pks[imin].xg+=ipx;
	pks[imin].yg+=ipy;
	pks[imin].npix++;
	printf("imin=%d, rmin=%lf\n",imin,rmin);
	printf("xg,yg=%lf %lf\n",pks[imin].xg,pks[imin].yg);
	printf("x,y=%lf %lf\n",pks[imin].x,pks[imin].y);
	printf("nclst=%d\n",pks[imin].npix);

	return(0);
};

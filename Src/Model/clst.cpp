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
		void clear_counter();
	private:
	protected:
};
void Pin::clear_counter(){
	npix=0;
	xg=0.0;
	yg=0.0;
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
//-----------------------------------------
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
		Pin *pks;
		int npk;
		void init_Pins(int n, int m);
		void cluster(int mtyp);
		void write_Pins(int i,int mtyp);
		int find_Pin(int ix, int iy, int mtyp);
		int **kcell_Q,**kcell_K,**kcell_N,**kcell_B;
		void deploy(int mtyp);
		void out_kcell(int mtyp);
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


	ptmp=(int *)malloc(sizeof(int)*npix);
	for(i=0;i<npix;i++) ptmp[i]=-1;
	kcell_B=(int **)malloc(sizeof(int *)*Ndiv[0]);
	for(i=0;i<Ndiv[0];i++) kcell_B[i]=ptmp+i*Ndiv[1];

	ptmp=(int *)malloc(sizeof(int)*npix);
	for(i=0;i<npix;i++) ptmp[i]=-1;
	kcell_Q=(int **)malloc(sizeof(int *)*Ndiv[0]);
	for(i=0;i<Ndiv[0];i++) kcell_Q[i]=ptmp+i*Ndiv[1];

	ptmp=(int *)malloc(sizeof(int)*npix);
	for(i=0;i<npix;i++) ptmp[i]=-1;
	kcell_K=(int **)malloc(sizeof(int *)*Ndiv[0]);
	for(i=0;i<Ndiv[0];i++) kcell_K[i]=ptmp+i*Ndiv[1];

	ptmp=(int *)malloc(sizeof(int)*npix);
	for(i=0;i<npix;i++) ptmp[i]=-1;
	kcell_N=(int **)malloc(sizeof(int *)*Ndiv[0]);
	for(i=0;i<Ndiv[0];i++) kcell_N[i]=ptmp+i*Ndiv[1];

};
int IMG::find_Pin(int ix, int iy,int mtyp){
	int j,imin;
	double rmin,r;
	rmin=pks[0].dist((float)ix,(float)iy);
	imin=0;
	int n0=(mtyp-1)*npk;
	for(j=1;j<npk;j++){
		r=pks[j+n0].dist((float)ix,(float)iy);
		if(rmin>r){
			rmin=r;
			imin=j+n0;
		}
	}
	return(imin);	// return global address
};
void IMG::cluster(int mtyp){

	int i,j,k,imin;
	int ipx,ipy;
	double r,rmin;
	int *Ms,nM;

	if(mtyp==0){ nM=nB; Ms=B;}
	if(mtyp==1){ nM=nQ; Ms=Q;}
	if(mtyp==2){ nM=nK; Ms=K;}
	if(mtyp==3){ nM=nN; Ms=N;}

	int n0=(mtyp-1)*npk;
	for(i=0;i<npk;i++) pks[n0+i].clear_counter();

	for(i=0;i<nM;i++){
		k=Ms[i];
		ipx=IMG::ipix(k);
		ipy=IMG::jpix(k);

		imin=IMG::find_Pin(ipx,ipy,mtyp);
		pks[imin].xg+=ipx;
		pks[imin].yg+=ipy;
		pks[imin].npix++;

	};

	int id;
	for(i=0;i<npk;i++){
		id=i+n0;
		if(pks[id].npix>0){
			pks[id].xg/=pks[id].npix;
			pks[id].yg/=pks[id].npix;
			pks[id].x=pks[id].xg;
			pks[id].y=pks[id].yg;
		}
	};

};
void IMG::init_Pins(int n, int m){
	npk=n*m;
	pks=(Pin *)malloc(sizeof(Pin)*npk*3);

	int dnpx,dnpy;
	dnpx=int(Ndiv[0]/n);
	dnpy=int(Ndiv[1]/m);

	int npx=n;
	int npy=m;

	int i,j,k=0;
	int ipx,ipy;
	for(i=0;i<npx;i++){
		ipx=int((i+0.5)*dnpx);
	for(j=0;j<npy;j++){
		ipy=int((j+0.5)*dnpy);
		//printf("%d %d\n",ipx,ipy);
		pks[k].init();
		pks[k+npk].init();
		pks[k+2*npk].init();

		pks[k].set(ipx,ipy);
		pks[k+npk].set(ipx,ipy);
		pks[k+2*npk].set(ipx,ipy);
		k++;
	}
	}	
};

void IMG::deploy(int mtyp){
	int i,j,k;
	int ipx,ipy,ipin;

	int *Ms,nM,**kcell;
	if(mtyp==0){
		nM=nB; Ms=B; kcell=kcell_B;
	}
	if(mtyp==1){
		nM=nQ; Ms=Q; kcell=kcell_Q;
	}
	if(mtyp==2){
		nM=nK; Ms=K; kcell=kcell_K;
	}
	if(mtyp==3){
		nM=nN; Ms=N; kcell=kcell_N;
	}

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		kcell[i][j]=-1;
	}
	}
	int n0=(mtyp-1)*npk;

	for(i=0;i<nM;i++){
		k=Ms[i];
		ipx=IMG::ipix(k);
		ipy=IMG::jpix(k);
		ipin=IMG::find_Pin(ipx,ipy,mtyp);
		//kcell[ipx][ipy]=ipin;	// global Pin No.
		kcell[ipx][ipy]=ipin-n0;// local Pin No.
	}
};
void IMG::write_Pins(int i,int mtyp){
	char fname[128];

	if(mtyp==0); sprintf(fname,"p%dB.out",i);
	if(mtyp==1); sprintf(fname,"p%dQ.out",i);
	if(mtyp==2); sprintf(fname,"p%dK.out",i);
	if(mtyp==3); sprintf(fname,"p%dN.out",i);
	FILE *fp=fopen(fname,"w");

	int n0=(mtyp-1)*npk;
	for(int j=0;j<npk;j++){
		fprintf(fp,"%lf, %lf\n",pks[j+n0].x,pks[j+n0].y);
	}

	fclose(fp);
};
void IMG::out_kcell(int mtyp){
	char fname[128];
	if(mtyp==0) sprintf(fname,"%s","grainB.out");
	if(mtyp==1) sprintf(fname,"%s","grainQ.out");
	if(mtyp==2) sprintf(fname,"%s","grainK.out");
	if(mtyp==3) sprintf(fname,"%s","grainN.out");

	int **kcell;
	if(mtyp==0) kcell=kcell_B;
	if(mtyp==1) kcell=kcell_Q;
	if(mtyp==2) kcell=kcell_K;
	if(mtyp==3) kcell=kcell_N;


	FILE *fp=fopen(fname,"w");
	fprintf(fp,"#Nx, Ny\n");	
	fprintf(fp,"%d, %d\n",Ndiv[0],Ndiv[1]);	
	fprintf(fp,"# grain No.\n");	
	int i,j;
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp,"%d\n",kcell[i][j]);
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
//-------------------------------------------------------
int main(){
	char fname[128]="minmap.out";

	IMG im3;
	im3.load(fname);
	im3.count();


	int i,itr=3;
	int npx=10, npy=10;

	for(int mtyp=1;mtyp<4;mtyp++){
		im3.init_Pins(npx,npy);
		im3.write_Pins(0,mtyp);
		for(i=0;i<itr;i++){
			im3.cluster(mtyp);
			im3.write_Pins(i+1,mtyp);
		}
		im3.deploy(mtyp);
		im3.out_kcell(mtyp);
	};
	return(0);
};

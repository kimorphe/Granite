#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"fdm2d.h"

#define DEBUG 0
/*
class cp_list{
	public:
		int load(char fn[128]);
		char fname[128];
		double df;
		int Nf,nfile;
		int nfreq,ndat;
		int nf1,nf2;
		double *cp_bank,*freq;
		void print_cp_bank();
	private:
		void read_header();
		void read_all();
		int retain(double f1, double f2);
	protected:
};
*/
int cp_list::load(char fn[128]){
	strcpy(fname,fn);
	int n=cp_list::retain(0.5,1.5);
	return(n);
};

void cp_list::print_cp_bank(){
	for(int i=0;i<ndat;i++){
		printf("%d, %lf, %lf\n",i, freq[i],cp_bank[i]);
	};
};
void cp_list::read_header(){

	FILE *fp=fopen(fname,"r");
	double row[4];
	char cbff[128];

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",&nfile, &Nf);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&df);
	fclose(fp);
};
int cp_list::retain(double f1, double f2){

	char cbff[128];
	FILE *fp=fopen(fname,"r");
	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",&nfile, &Nf);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&df);
	printf("df=%lf\n",df);
	printf("nfile=%d, Nf=%d\n",nfile, Nf);
	double *cp=(double *)malloc(sizeof(double)*Nf);
	double *frq=(double *)malloc(sizeof(double)*Nf);
	double row[4];


	nf1=int(f1/df);
	nf2=ceil(f2/df);
	nfreq=(nf2-nf1+1);
	ndat=nfreq*nfile;
	cp_bank=(double *)malloc(sizeof(double)*ndat);
	freq=(double *)malloc(sizeof(double)*ndat);


	int i,j,k=0;
	for(i=0;i<nfile;i++){
		fgets(cbff,128,fp);
		for(j=0;j<Nf;j++){
			fscanf(fp,"%lf, %lf, %lf, %lf\n",row,row+1,row+2,row+3);
			frq[j]=row[0];
			cp[j]=row[3];
		}
		for(j=nf1;j<=nf2;j++){
		       	cp_bank[k]=cp[j];
		       	freq[k++]=frq[j];
		};
	}

	printf("File=%s\n",fname);
	printf("nf1=%d, nf2=%d, nfreq=%d, nfile=%d\n",nf1,nf2,nfreq,nfile);
	printf("ndat=%d\n",ndat);
	printf("-------------------------\n");
	fclose(fp);
	return(ndat);

};

double cp_list::stats(){
	int i;
	double cp1=0.0;
	double cp2=0.0;
	double dat;
	for(i=0;i<ndat;i++){
		dat=cp_bank[i];
		cp1+=dat;
		cp2+=(dat*dat);
	};
	cp_mean=cp1/ndat;
	cp_var=cp2/ndat-cp_mean*cp_mean;;
	printf("%s\n",fname);
	printf("<cp>=%lf, ",cp_mean);
	printf("var{cp}=%lf, ",cp_var);
	printf("std{cp}=%lf\n",sqrt(cp_var));
	return(cp2);
};
void cp_list::read_all(){
	int i,j;

	char cbff[128];
	FILE *fp=fopen(fname,"r");
	double row[4];
	int Nf,nfile;

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",&nfile, &Nf);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&df);
	printf("df=%lf\n",df);
	printf("nfile=%d, Nf=%d\n",nfile, Nf);
	double *cp=(double *)malloc(sizeof(double)*Nf);
	double *freq=(double *)malloc(sizeof(double)*Nf);

	for(i=0;i<nfile;i++){
		fgets(cbff,128,fp);
	for(j=0;j<Nf;j++){
		fscanf(fp,"%lf, %lf, %lf, %lf\n",row,row+1,row+2,row+3);
		freq[j]=row[0];
		cp[j]=row[3];
	}
	}

	fclose(fp);
};

#if DEBUG==1
int main(){
	cp_list cQt, cK, cNa;
	char fn1[128]="../../PPT_Figs/phi_Qt.out";
	char fn2[128]="../../PPT_Figs/phi_K.out";
	char fn3[128]="../../PPT_Figs/phi_Na.out";

	cQt.load(fn1);
	cK.load(fn2);
	cNa.load(fn3);
	//cQt.print_cp_bank();
	double c2Q, c2K, c2N;
	c2Q=cQt.stats();
	c2K=cK.stats();
	c2N=cNa.stats();

	double c_mean,c_var;
	c_mean =cQt.cp_mean*cQt.ndat;
	c_mean+= cK.cp_mean*cK.ndat;
	c_mean+=cNa.cp_mean*cNa.ndat;

	double c2;
	c2=c2Q+c2K+c2N;
	int ndat;
	ndat=cQt.ndat+cK.ndat+cNa.ndat;
	c_mean/=ndat;
	c_var=c2/ndat-c_mean*c_mean;
	printf("c_mean=%lf, var{c}=%lf, std{c}=,%lf\n",c_mean,c_var,sqrt(c_var));


	return(0);
};
#endif

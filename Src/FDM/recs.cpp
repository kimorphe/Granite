#include<stdio.h>
#include<stdlib.h>
#include"fdm2d.h"

void rec_array::set_fd_grid(double xa[2],double del, int ng[2], double dtau, int nt ){
	Xa[0]=xa[0];
	Xa[1]=xa[1];
	dh=del;
	Ndiv[0]=ng[0]; 
	Ndiv[1]=ng[1];
	Nt=nt;
	dt=dtau;
};
void rec_array::set_array( int ixy, double u1, double u2, double vcod, int np){
	idir=ixy;
	int i,iad;
	double uu,du=0;
	if(np >1) du=(u2-u1)/(np-1);
	irecs=(int *)malloc(sizeof(int)*np);
	npnt=0;
	for(i=0;i<np;i++){
		uu=u1+i*du;
		if(idir==0){
			iad=int((uu-Xa[0])/dh);
			if(iad>=Ndiv[0]) continue;
		}

		if(idir==1){
			iad=int((uu-Xa[1])/dh);
			if(iad>=Ndiv[1]) continue;
		};
		if(iad<0) continue;
		irecs[npnt]=iad;
//		printf("irecs=%d\n",irecs[npnt]);
		npnt++;
	};

	if(idir==0){
		jrec=int((vcod-Xa[1])/dh);
		if(jrec>=Ndiv[1]) npnt=0;
	}
	if(idir==1){
	       	jrec=int((vcod-Xa[0])/dh);
		if(jrec>=Ndiv[0]) npnt=0;
	}
	if(jrec<0) npnt=0;
	printf("jrec=%d npnt=%d\n",jrec,npnt);
	printf("Ndiv==%d %d\n",Ndiv[0],Ndiv[1]);

	v1=(double *)malloc(sizeof(double)*npnt*Nt);
	v2=(double *)malloc(sizeof(double)*npnt*Nt);
	s=(double *)malloc(sizeof(double)*npnt*Nt);

	cntr=0;
};
void rec_array::write2(char fn[128]){
	FILE *fp=fopen(fn,"w");
	int i,j,jd;

	double xx,yy;
	fprintf(fp,"# idir(0:x-aligned, 1:y-alinged\n");
	fprintf(fp,"%d\n",idir);
	fprintf(fp,"# npnt\n");
	fprintf(fp,"%d\n",npnt);
	fprintf(fp,"# dt, Nt\n");
	fprintf(fp,"%lf, %d\n",dt,Nt);
	for(i=0;i<npnt;i++){
		if(idir==0){
			xx=Xa[0]+dh*(irecs[i]+0.5);
			yy=Xa[1]+dh*(jrec+0.5);
		};
		if(idir==1){
			xx=Xa[1]+dh*(jrec+0.5);
			yy=Xa[1]+dh*(irecs[i]+0.5);
		};
		fprintf(fp,"#xf, yf\n");
		fprintf(fp,"%lf, %lf\n",xx,yy);
		jd=i;
		for(j=0;j<Nt;j++){
			fprintf(fp," %lf, %lf, %lf\n",v1[jd],v2[jd],s[jd]);
			jd+=npnt;
		};
	};
	fclose(fp);
};


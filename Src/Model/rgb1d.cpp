#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main(){
	FILE *fin=fopen("rgb1d.dat","r");
	int ndat;

	fscanf(fin,"%d\n",&ndat);
	printf("ndat=%d\n",ndat);
	int *rx=(int *)malloc(sizeof(int)*ndat);
	int *gx=(int *)malloc(sizeof(int)*ndat);
	int *bx=(int *)malloc(sizeof(int)*ndat);

	int i;
	for(i=0;i<ndat;i++){
		fscanf(fin,"%d, %d, %d\n",rx+i,gx+i,bx+i);
		printf("%d  %d, %d, %d\n",i,rx[i],gx[i],bx[i]);
	}
	fclose(fin);

	FILE *fout=fopen("rgb1d.out","w");
	double y,dy,dr,dg,db;
	for(i=0;i<ndat;i++){
		y=(rx[i]+gx[i]+bx[i])/3.0;
		dr=rx[i]-y;
		dg=gx[i]-y;
		db=bx[i]-y;
		fprintf(fout,"%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",y,dr,dg,db,dr/y,dg/y,db/y);
	};
	fclose(fout);

	return(0);
};

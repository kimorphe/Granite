#include<stdio.h>
#include<math.h>
#include<stdlib.h>


int main(){
	FILE *fp=fopen("rgb.dat","r");

	int *R,*G,*B;
	int *pt;


	int i,ndat;
	fscanf(fp,"%d\n",&ndat);

	R=(int *)malloc(sizeof(int)*ndat);
	G=(int *)malloc(sizeof(int)*ndat);
	B=(int *)malloc(sizeof(int)*ndat);

	int r,g,b;
	for(i=0;i<ndat;i++){
		fscanf(fp,"%d, %d, %d\n",R+i,G+i,B+i);
	};
	printf("%d, %d, %d\n",R[ndat-1],G[ndat-1],B[ndat-1]);

	fclose(fp);
	return(0);
};

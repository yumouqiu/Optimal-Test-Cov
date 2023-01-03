#include <stdio.h>
#include <R.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Rmath.h>

void P(int *m, int *n, double *s)
{
	int i;

	for(i=0; i<*m; i++)
{
	*s=(*s)*(*n-i);
			
}
				
}

double PC1(int m, int n)
{
	int i;
	double s=1;
	for(i=0; i<m; i++)
{
	s=(s)*(n-i);
			
}
 return s;				
}

double dot(double x[], double y[], int l){
	double result=0;
	for(int i=0; i<l; i++){
		result+=x[i]*y[i];
}
	return result;
}

void adT(double *x, double *r, double *c, int *a, double *A, double *B,
	 double *mf22, double *mf1, double *mf2, double *mf3, double *res){
		
	double s1=0;
	int cval= *c;
	double temp1[cval],temp2[cval];
	int i=0, j=0;
	for(i = 0; i<*r; i++){

		double s2=0;
		for(int y=0,z=i; y<*c; y++,z+=*r){
		temp1[y]=x[z];		
		}
	
	for(j = 0; j<*r; j++){

	for(int y=0,z=j; y<*c; y++,z+=*r){
		temp2[y]=x[z];		
		}

s2+=*mf1 * dot(temp1,temp2,*c) * dot(temp1,temp2,*c) +
    *mf2 * *mf22 * dot(temp1,temp2,*c) - 
    *mf3 * dot(temp1,temp2,*c) * dot(temp2,A,*c) +
    *mf3 * dot(temp1,temp2,*c) * dot(temp2,temp2,*c);

}

s1+=s2 - *mf1 * dot(temp1,temp1,*c) * dot(temp1,temp1,*c) -
    *mf2 * *mf22 * dot(temp1,temp1,*r) +
    *mf3 * dot(temp1,temp1,*c) * dot(temp1,A,*c) -
    *mf3 * dot(temp1,temp1,*c) * dot(temp1,temp1,*c);

}
*res=s1;
}

int pullhelp(int i, int j, int r, int c){
	
	int p;
	
	if(i==0 && j>0){
		p=r*j;
}
	else if(j==0 && i>0){
		p=i;
}
	else {
		p=i*j+r*j-(j-1)*i;	
}
	return p;
}

double pull(double x[], int i, int j, int r, int c){
	
	double p = x[pullhelp(i, j, r, c)];
	return p;
}

double adb1(double x[],int r,int c,int i,int j,int q){
	double mfb1 = 1/PC1(2,r) + 2/PC1(3,r) + 2/PC1(4,r);
	double mfb2 = 2/PC1(3,r) + 3/PC1(4,r);
	double mfb3 = 1/PC1(4,r);
	double sband1 = 0;
	
	for(int n=0;n<c-q;n++){
		double sumnq=0,sumn=0;

	for(int y=0,z=n*r; y<r; y++,z++){
		sumn+=x[z];
}
	
	for(int a=0,b=n*r+r*q; a<r; a++,b++){
		sumnq+=x[b];
		
}

sband1+=mfb1 * pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c)*pull(x,i,n+q,r,c) + 
	mfb2 * (pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c) - 
	        pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c)*sumnq) +
	mfb3 * (pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumn*sumnq + 
		pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c) - 
		pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumnq - 
		pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)*sumn);

}
	return sband1;
}

double adb2(double x[], int r, int c, int i, int j, int k){

	double sband2 = 0;
	for(int m=0;m<k+1;m++){

	sband2+=2*adb1(x,r,c,i,j,m);
}
	sband2=sband2-adb1(x,r,c,i,j,0);
	return sband2;
}


void adb(double *x, int *r, int *c, int *k, double *res){

	double sband3 = 0;

	for(int i=0; i<*r; i++){
	
	double sband4 = 0;

	for(int j=0; j<*r; j++){

	sband4+=adb2(x, *r, *c, i, j, *k);
}
	sband3+=sband4-adb2(x, *r, *c, i, i, *k);
}

	*res = sband3;

}

void weight(int *khv, int *kv, double *w){

	int k = *kv, kh = *khv;
	double d=0;
	for(int i=0;i<k+1;i++){
	d=i;
	if(kh>d){
		w[i]=1;
		}
	if(kh<=d && kh>(d/2)){
		w[i]=2-(d / kh);
		}
	if(kh<=(d/2)){
		w[i]=0;
		}
}

}

void bandpen1call(double *x, int *rv, int *cv, int *iv, int *jv, int *qv, double *res){

int r = *rv, c = *cv, i = *iv, j= *jv, q = *qv;
double sband1, sum1, penalty1;
double mfb1 = (1 / PC1(2,r)) + (2 / PC1(3,r)) + (2 / PC1(4,r));
double mfb2 = (2 / PC1(3,r)) + (3 / PC1(4,r));
double mfb3 =  1 / PC1(4,r);
double c1 = (1 / PC1(2,r)) + (2 / PC1(3,r)) + (1 / PC1(4,r));
double c2 = (1 / PC1(3,r)) + (1 / PC1(4,r));
double c3 = (1 / PC1(3,r)) + (2 / PC1(4,r));
double c4 =  1 / PC1(4,r);

	for(int n=0;n<c-q;n++){
		double sumnq=0,sumn=0;
		double tempn[r], tempnq[r];
	for(int y=0,z=n*r; y<r; y++,z++){
		tempn[y]=x[z];
		sumn+=x[z];
}
	
	for(int a=0,b=n*r+r*q; a<r; a++,b++){
		tempnq[a]=x[b];		
		sumnq+=x[b];
		
}

sum1 = dot(tempn,tempnq,r);
sband1+= mfb1 * pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c)*pull(x,i,n+q,r,c) + 
	 mfb2 * (pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c) - 
	        pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c)*sumnq) +
	 mfb3 * (pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumn*sumnq + 
		pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c) - 
		pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumnq - 
		pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)*sumn);

penalty1+= c1*(pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c))-
	   c2*(pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)*sumn +
	       pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumnq) +
	   c3*(pull(x,i,n,r,c)*pull(x,j,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)+
	       pull(x,i,n+q,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c)*pull(x,j,n,r,c)) +
	   c4*(pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumn*sumnq +
	       pull(x,i,n,r,c)*pull(x,i,n+q,r,c)*pull(x,j,n,r,c)*pull(x,j,n+q,r,c)-
	       pull(x,i,n,r,c)*pull(x,j,n,r,c)*pull(x,j,n+q,r,c)*sumnq-
               pull(x,i,n,r,c)*pull(x,i,n+q,r,c)*pull(x,j,n+q,r,c)*sumn-
	       pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sum1);

}

double store[2] = {sband1, penalty1};

for(int a=0; a<2; a++){
res[a] = store[a];
}

}


double bandpen1(double x[], int r, int c, int i, int j, int q){
double sum1, penalty1;
double c1 = (1 / PC1(2,r)) + (2 / PC1(3,r)) + (1 / PC1(4,r));
double c2 = (1 / PC1(3,r)) + (1 / PC1(4,r));
double c3 = (1 / PC1(3,r)) + (2 / PC1(4,r));
double c4 =  1 / PC1(4,r);

	for(int n=0;n<c-q;n++){
		double sumnq=0,sumn=0;
		double tempn[r], tempnq[r];
	for(int y=0,z=n*r; y<r; y++,z++){
		tempn[y]=x[z];
		sumn+=x[z];
}
	
	for(int a=0,b=n*r+r*q; a<r; a++,b++){
		tempnq[a]=x[b];		
		sumnq+=x[b];
		
}


sum1 = dot(tempn,tempnq,r);
penalty1+= c1*(pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c))-
	   c2*(pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)*sumn +
	       pull(x,i,n,r,c)*pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumnq) +
	   c3*(pull(x,i,n,r,c)*pull(x,j,n,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n+q,r,c)+
	       pull(x,i,n+q,r,c)*pull(x,j,n+q,r,c)*pull(x,j,n,r,c)*pull(x,j,n,r,c)) +
	   c4*(pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sumn*sumnq +
	       pull(x,i,n,r,c)*pull(x,i,n+q,r,c)*pull(x,j,n,r,c)*pull(x,j,n+q,r,c)-
	       pull(x,i,n,r,c)*pull(x,j,n,r,c)*pull(x,j,n+q,r,c)*sumnq-
               pull(x,i,n,r,c)*pull(x,i,n+q,r,c)*pull(x,j,n+q,r,c)*sumn-
	       pull(x,i,n,r,c)*pull(x,j,n+q,r,c)*sum1);

}

return penalty1;
}

void sbandpen2call(double *x, int *rv, int *cv, int *qv, double *res){
	int r = *rv, c = *cv, q = *qv;
	double sband2=0;
	for(int i=0; i<r; i++){

	double sband3=0;

		for(int j=0; j<r; j++){
		
		sband3+=adb1(x,r,c,i,j,q);
	}
	
	sband2+=sband3-adb1(x,r,c,i,i,q);
}

	*res=sband2;

}

double sbandpen2(double x[], int r, int c, int q){
	double sband2=0;
	for(int i=0; i<r; i++){

	double sband3=0;

		for(int j=0; j<r; j++){
		
		sband3+=adb1(x,r,c,i,j,q);
	}
	
	sband2+=sband3-adb1(x,r,c,i,i,q);
}

	return sband2;

}


void pbandpen2call(double *x, int *rv, int *cv, int *qv, double *res){
	int r = *rv, c = *cv, q = *qv;
	double pband2=0;
	for(int i=0; i<r; i++){

	double pband3=0;

		for(int j=0; j<r; j++){
		
		pband3+=bandpen1(x,r,c,i,j,q);
					}
	
	pband2+=pband3-bandpen1(x,r,c,i,i,q);
}

	*res=pband2;

}

double pbandpen2(double x[], int r, int c, int q){
	double pband2=0;
	for(int i=0; i<r; i++){

	double pband3=0;

		for(int j=0; j<r; j++){
		
		pband3+=bandpen1(x,r,c,i,j,q);
					}
	
	pband2+=pband3-bandpen1(x,r,c,i,i,q);
}

	return pband2;

}

void bandpen(double *x, int *rv, int *cv, int *kv, double *store){

	int r = *rv, c = *cv, k = *kv;
	
	for(int z=1;z<k+1;z++){
	int a=2*z;
	int b=2*z+1;
	store[a] = 2*sbandpen2(x,r,c,z);
	store[b] = 2*pbandpen2(x,r,c,z);	
	
}

}


void bandloss(int *rv, int *cv, int *kv, double *adtv,
 double *tempres, double *res){

	int k = *kv;
	double r = *rv, c =*cv, adt = *adtv;
	double sumi1,sumi2;

	for(int i=0;i<k+1;i++){

	sumi1=0;
	sumi2=0;

	for(int j=0;j<i+1;j++){
	sumi1+=tempres[j];
	}	
		
	for(int l=k+1;l<k+1+i+1;l++){
	sumi2+=tempres[l];
	}

	res[i] = (1 / c) * (adt-sumi1) + (1 / (r*c)) * sumi2;

}

}


double tapehelp(double tempres[],double temp1[], double temp2[],
		int i, int k, int boo){
	double w[k+1],wgttape[k+1],wgtpen[k+1];
	double ret;
	
	for(int z=0;z<k+1;z++){
		double d=z;
		if(i>z){
			w[z]=1;
			}
		if(i<=d && i>(d / 2)){
			w[z]=2-(d / i);
			}
		if(i<=(d / 2)){
			w[z]=0;
			}
		}
		
		for(int y = 0; y<k+1;y++){
			wgttape[y] = 1-(1-w[y])*(1-w[y]);
			wgtpen[y] = w[y]*w[y];	
			}
	if(boo==1){
	ret = dot(temp1,wgttape,k+1);	
	}
	else{
	ret = dot(temp2,wgtpen,k+1);
	}

	return ret;
}


void tapeloss(int *rv, int *cv, int *kv, double *adtv,
 double *tempres, double *res){

	int k = *kv;
	double kf=*kv, f = floor(kf / 2);
	double r = *rv, c =*cv, adt = *adtv;
	double dt1,dt2;
	double temp1[k+1],temp2[k+2];
	
	for(int j=0,l=k+1;j<k+1;j++,l++){
	temp1[j]=tempres[j];
	temp2[j]=tempres[l];		
	}	


	res[0] = (1 / c) * (adt-tempres[0]) + (1 / (r*c)) * tempres[k+1];
	
	for(int i=1;i<f+1;i++){	
	dt1 = tapehelp(tempres,temp1,temp2,i,k,1);
	dt2 = tapehelp(tempres,temp1,temp2,i,k,2);
	
	res[i] = (1 / c) * ((adt-dt1) + (dt2 / r));

	}
	
}




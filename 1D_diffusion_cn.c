/*
Solve the 1D diffusion equation using the implicit Crank-Nicolson method
*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define PI 3.14159
#define E 2.71828
#define C 100
#define DIFFUSIVITY 0.02
#define LENGTH 2
#define TIME_SOL 10
#define M 100
#define TIME_STEPS 1000

//Function at t=0
void f_init(double** x, double** tmp_old){
	int i;
	for(i=0;i<=M;i++){
	  (*tmp_old)[i] = sin(PI * (*x)[i]);
	}
}

//Analytical solution
void f_analytic(double** x, double t, double** tmp_analytic){
	int i;
	for(i=0;i<=M;i++){
	  (*tmp_analytic)[i] = pow(E,-1*(PI*PI*t)) * sin(PI * (*x)[i]);
	}
}

//Crank Nicolson
void f_cn(double** tmp_old, double** tmp_new, double z, void (*tri)(double**,double**,double**,double**,double**,double**)){
	double *a, *b, *c, *d, *cn, *dn;
	int i,j,k;
	
	a = (double*)malloc((M+1)*sizeof(double));
	b = (double*)malloc((M+1)*sizeof(double));
	c = (double*)malloc((M+1)*sizeof(double));
	d = (double*)malloc((M+1)*sizeof(double));
	cn = (double*)malloc((M+1)*sizeof(double));
	dn = (double*)malloc((M+1)*sizeof(double));
	
	//update boundary points
	(*tmp_new)[0] = (*tmp_old)[0];
	(*tmp_new)[M] = (*tmp_old)[M];
	
	//populate coefficient vectors for tridiagonal system
	for(i=0;i<=M;i++){
		if(i==0 || i==M){
			a[i] = 0;
			b[i] = 0;
			c[i] = 0;
			d[i] = 0;
			cn[i] = 0;
			dn[i] = 0;
		}
		else if(i==1){
			a[i] = 0;
			b[i] =  2 + 2 * z;
			c[i] = -1 * z;			
			d[i] = z*(*tmp_old)[i-1] + (2-2*z)*(*tmp_old)[i] + z*(*tmp_old)[i+1] + z*(*tmp_new)[i-1];
		}
		else if(2<=i<=M-2){
			a[i] = -1 * z;
			b[i] = 2 + 2 * z;
			c[i] = -1 * z;
			d[i] = z*(*tmp_old)[i-1] + (2-2*z)*(*tmp_old)[i] + z*(*tmp_old)[i+1];										
		}
		if(i==M-1){
			a[i] = -1 * z;
			b[i] = 2 + 2 * z;
			c[i] = 0;
			d[i] = z*(*tmp_old)[i-1] + (2-2*z)*(*tmp_old)[i] + z*(*tmp_old)[i+1] + z*(*tmp_new)[i+1]; 						
		}
	}

	//solve the tridiagonal system
	tri(&a, &b, &c, &d, &cn, &dn);
	
	//update tmp_new
	for(j=M-1;j>=1;j--){
		if(j==M-1){
			(*tmp_new)[j] = dn[j];
		}
		else{
			(*tmp_new)[j] = dn[j] - cn[j]*(*tmp_new)[j+1];
		}
	}
}

//Tridiagonal
void tri(double** a, double** b, double** c, double** d, double** cn, double** dn){
	int i;
		
	for(i=1;i<=M-1;i++){
		if(i==1){
			(*cn)[i] = (*c)[i]/(*b)[i];
			(*dn)[i] = (*d)[i]/(*b)[i];
		}
		if(1<i<M-1){
			(*cn)[i] = (*c)[i]/((*b)[i]-(*a)[i]*(*cn)[i-1]);
			(*dn)[i] = ((*d)[i]-(*a)[i]*(*dn)[i-1])/((*b)[i]-(*a)[i]*(*cn)[i-1]);
		}
		if(i==M-1){
			(*cn)[i] = 0;
			(*dn)[i] = ((*d)[i]-(*a)[i]*(*dn)[i-1])/((*b)[i]-(*a)[i]*(*cn)[i-1]);
		}
	}
}

//update tmp_old
void tmp_update(double** tmp_old, double** tmp_new){
	int i;
	for(i=0;i<=M;i++){
	  (*tmp_old)[i] = (*tmp_new)[i];
	}
}
//calculate error
double error(double** tmp_cn, double** tmp_analytic){
	double err=0.0,sqsum=0.0;
	int i;
	for(i=0;i<=M;i++){
		sqsum += ((*tmp_cn)[i]-(*tmp_analytic)[i]) * ((*tmp_cn)[i]-(*tmp_analytic)[i]);
	}
	err = (1.0/(M+1))*sqrt(sqsum);
	return err;
}


//main
int main(){
  int i,j,k=0;
	double *x, *tmp_old, *tmp_new, *tmp_analytic;
	double z, dt, dx, t=0.0;
	double err;
	
	FILE *ft20 = fopen("ft20.dat","w+");
	FILE *ft22 = fopen("ft22.dat","w+");
	FILE *ft24 = fopen("ft24.dat","w+");
	FILE *ft26 = fopen("ft26.dat","w+");
	FILE *ft28 = fopen("ft28.dat","w+");
	FILE *ft210 = fopen("ft210.dat","w+");
	
	dt = (0.05)/(TIME_STEPS);
	dx = (1.0)/(M);

	//memory allocation
	x = (double*)malloc((M+1) * sizeof(double));
	tmp_old = (double*)malloc((M+1) * sizeof(double));
	tmp_new = (double*)malloc((M+1) * sizeof(double));
	tmp_analytic = (double*)malloc((M+1) * sizeof(double));

	//calculate diffusivity for dimensionless equation	
	z = (dt)/(dx*dx);
	
	//initialization
	//x t=0
	for(i=0;i<=M;i++){
		x[i] = dx * i;
	}
	//tmp at all points t=0
	f_init(&x, &tmp_old);


	//time marching CN
	for(k=0;k<=TIME_STEPS;k++){
		t = (k * dt);
		if(t!=0){
			f_cn(&tmp_old, &tmp_new, z, tri);
			tmp_update(&tmp_old, &tmp_new);
		}
		
		//write solution at even time intervals
		if(k==0){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
			  fprintf(ft20,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j]);
			}
		}
		if(k==(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft22,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==2*(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft24,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==3*(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft26,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==4*(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft28,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==(TIME_STEPS)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft210,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		
	}
	
	fclose(ft20);
	fclose(ft22);
	fclose(ft24);
	fclose(ft26);
	fclose(ft28);
	fclose(ft210);
	
	return 0;
}


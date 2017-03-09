/*
Solving the 1D diffusion equation using the Fowrard Time Central Space explicit finite difference scheme
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

//FTCS
void f_ftcs(double z, double** tmp_old, double** tmp_new){
	int i;
	//grid points other than boundary
	for(i=1;i<=M-1;i++){
	  (*tmp_new)[i] = (*tmp_old)[i] + z * ((*tmp_old)[i + 1] - 2 * (*tmp_old)[i] + (*tmp_old)[i - 1]);	
	}
	//boundary points
	(*tmp_new)[0] = (*tmp_old)[0];
	(*tmp_new)[M] = (*tmp_old)[M];
}

//update tmp_old
void tmp_update(double** tmp_old, double** tmp_new){
	int i;
	for(i=0;i<=M;i++){
	  (*tmp_old)[i] = (*tmp_new)[i];
	}
}

//calculate error
double error(double** tmp_ftcs, double** tmp_analytic){
	double err=0.0,avg=0.0;
	int i;
	for(i=0;i<=M;i++){
		avg += ((*tmp_ftcs)[i]- (*tmp_analytic)[i]) * ((*tmp_ftcs)[i]- (*tmp_analytic)[i]);
	}
	err = (1.0/(M+1))*sqrt(avg);
	return err;
}

//main
int main(){
	int i,j,k=0;
	double *x, *tmp_old, *tmp_new, *tmp_analytic;
	double z, dt, dx, t=0.0;
	double err;
	
	FILE *ft0 = fopen("ft0.dat","w+");
	FILE *ft2= fopen("ft2.dat","w+");
	FILE *ft4 = fopen("ft4.dat","w+");
	FILE *ft6 = fopen("ft6.dat","w+");
	FILE *ft8 = fopen("ft8.dat","w+");
	FILE *ft10 = fopen("ft10.dat","w+");
	
	dt = (0.05)/(TIME_STEPS);
	printf("dt = %lf",dt);
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


	//time marching FTCS
	for(k=0;k<=TIME_STEPS;k++){
		//printf("tinit = %lf",t);
		t = (k*dt);
		if(k!=0){
			f_ftcs(z, &tmp_old, &tmp_new);
			tmp_update(&tmp_old, &tmp_new);
			
		}
		
		//write solution at even time intervals
		if(k==0){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
			  fprintf(ft0,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			//printf("t = %lf",t);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft2,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==2*(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			//printf("t = %lf",t);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft4,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==3*(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			//printf("t = %lf",t);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft6,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==4*(TIME_STEPS/5)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			//printf("t = %lf",t);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft8,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		if(k==(TIME_STEPS)){
			//calculate analytic solution
			f_analytic(&x, t, &tmp_analytic);
			//printf("t = %lf",t);
			err = error(&tmp_old,&tmp_analytic);
			for(j=0;j<=M;j++){
				fprintf(ft10,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", x[j],tmp_old[j],tmp_analytic[j],err);
			}
		}
		
	}
	
	fclose(ft0);
	fclose(ft2);
	fclose(ft4);
	fclose(ft6);
	fclose(ft8);
	fclose(ft10);
	
	return 0;
}


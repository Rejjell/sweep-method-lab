#define _CRT_SECURE_NO_WARNINGS


#include "stdafx.h"
#include "stdio.h";
#include "conio.h";
#include "stdlib.h";
#include <time.h>
#include <omp.h>
#include <iostream>
#include <math.h>

float ** A;

double *a, *b, *c, *f, *x, *func;

int n = 50;
double left = -10;
double right = 10;
double step;
double halfstep;

/*void generateData(double ** &A,int n)
{
	A=new double*[n];
	copyA=new double*[n];
	for (int i = 0; i < n; i++)
	{
		A[i]=new double[n];
		copyA[i]=new double[n];
		for (int j = 0; j  < n; j ++)
		{
			A[i][j]=0;
			copyA[i][j]=0;
		}
	}

	x = new double[n];
	f = new double[n];
	copyf = new double[n];

	srand (time(NULL));

	for (int i = 0; i < n; i++)
	{
		A[i][i]=1000+rand()%1000;
		if (i>0) A[i][i-1]=1 + rand()%100;
		if (i<n-1) A[i][i+1]=1 + rand()%100;
		
		f[i]=rand()%100;
		
		copyA[i][i]=A[i][i];
		 copyA[i][i-1]= A[i][i-1];
		 copyA[i][i+1]= A[i][i+1];
		 copyf[i]=f[i];

		 x[i]=0;
	}

}*/

void deleteData(double **A,int n)
{
	for (int i = 0; i < n+1; i++)
		delete A[i];

	delete A;
	delete f;
	delete c;
}

void printMatrix(double **A, double *f, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("%0*.*f ", 4, 2,A[i][j]);

		printf("    %0*.*f", 4, 2, f[i]);
		printf("\n");
	}
}

double calcFunction(double value)
{
	return sin(value);
}

double calcDer(double value)
{
	return cos(value);
}

double calcSecondDer(double value)
{
	return -sin(value);
}

double calcSpline(double value)
{
	int k;
	double x0;

	if (value < left + halfstep) 
	{
		k = 0;
		x0 = x[1] - halfstep;
	}
	else if (value >= right - halfstep) 
	{
		k = n;
		x0 = x[n];
	}
	else 
	{
		k = (value - (left + halfstep))/step + 1;
		x0 = x[k] + halfstep;
	}

	return a[k] + b[k]*(value - x0) + c[k]*pow(value - x0,2);
}

double calcSplineDer(double value)
{
	int k;
	double x0;

	if (value < left + halfstep) 
	{
		k = 0;
		x0 = x[1] - halfstep;
	}
	else if (value >= right - halfstep) 
	{
		k = n;
		x0 = x[n];
	}
	else 
	{
		k = (value - (left + halfstep))/step + 1;
		x0 = x[k] + halfstep;
	}

	return b[k] + 2*c[k]*(value - x0);
}

double deltaf(int k)
{
	return func[k+1] - func[k];
}

double delta2f(int k)
{
	return deltaf(k+1) - deltaf(k);
}

void makeMatrix(float ** &A, int n)
{
	A = new float*[n+1];

	for (int i = 0; i < n+1; i++)
	{
		A[i] = new float[n+1];
		for (int j = 0; j  < n+1; j ++)
			A[i][j]=0;
	}

	A[0][0] = 2;
	A[n][n] = 2;

	for (int i=1; i<n; i++)
	{
		A[i][i] = 6;
		if (i>1) A[i][i-1] = 1;
		if (i<n-1) A[i][i+1] = 1;
	}

	f = new double[n+1];

	f[0] = calcSecondDer(left);
	f[n] = calcSecondDer(right);

	for (int i=1; i<n; i++)
	{
		f[i] = delta2f(i-1)/(halfstep*halfstep);
		if (i == 1) f[i] -= calcSecondDer(left)/2;
		if (i == n-1) f[i] -= calcSecondDer(right)/2;
	}

}

void serialSweepMethod(float** &A, double*&f, int size, double* &x )
{
	int n = size;

	double *alpha = new double[n];
	double * beta = new double[n];

	alpha[1]=-A[0][1]/A[0][0];
	beta[1]=f[0]/A[0][0];

	for (int i = 1; i < n-1; i++)
	{
		alpha[i+1]= -A[i][i+1]/(A[i][i-1]*alpha[i]+A[i][i]);	
		beta[i+1] = (f[i]-A[i][i-1]*beta[i])/(A[i][i-1]*alpha[i]+A[i][i]);
		
	}

	x[n-1]=(f[n-1] - A[n-1][n-2]*beta[n-1])/(A[n-1][n-2]*alpha[n-1] + A[n-1][n-1]);
	for (int i = n-2; i >=0; i--)
	{
		x[i]=alpha[i+1]* x[i+1] + beta[i+1];
	}
	 
	delete alpha; 
	delete beta;

}

int main(int argc, char **argv)
{
	n = 50000; //Размерность сетки
	step = (right - left)/n;
	halfstep = step/2;
	//generateData(A,n);;	

	x = new double[n+2];
	func = new double[n+2];

	for (int i=0; i<n+1; i++)
	{
		x[i] = left + i*step;
		func[i] = calcFunction(x[i]);
	}

	makeMatrix(A, n);

	a = new double[n+1];
	b = new double[n+1];
	c = new double[n+1];

	//printMatrix(A, f, n+1);

	serialSweepMethod(A, f, n+1, c);

	for (int i=0; i<n; i++)
	{
		b[i] = (halfstep/2)*(c[i] - c[i+1]) + deltaf(i)/step;
		a[i] = halfstep*(b[i] - halfstep*c[i]) + func[i];
	}

	b[n] = b[n-1] + calcSecondDer(right)*halfstep;
	a[n] = func[n];


	/*for (int i = 0; i < n+1; i++)
	{
		printf("%*.*f ", 4, 10, a[i]);
		printf("%*.*f ", 4, 10, b[i]);
		printf("%*.*f ", 4, 10, c[i]);
		printf("\n");
	}*/

	/*FILE *F; 
	F = fopen("E:\\Run.txt","w"); 
	
	for (int i=0; i<n+1; i++)
	{
		fprintf(F,"(%f;%f)", x[i], calcSpline(x[i]));
		printf("%*.*f ", 4, 10, calcSpline(x[i]));
	}

	fclose(F);*/

	double inaccuracy = 0; //Погрешность
	double inaccuracyDer = 0; //Погрешность производной
	double littleStep = step/4;

	for (double x = left; x <=right; x += littleStep)
	{
		double tmp = abs(calcFunction(x) - calcSpline(x));
		if (tmp > inaccuracy) inaccuracy = tmp;

		tmp = abs(calcDer(x) - calcSplineDer(x));
		if (tmp > inaccuracyDer) inaccuracyDer = tmp;
	}

	printf("%*.*f \n", 4, 10, inaccuracy);
	printf("%*.*f ", 4, 10, inaccuracyDer);


	/*double s=0;
	for (int i = 0; i < n+1; i++)
	{
		s=0;
		for (int j = 0; j < n+1; j++)
		{
			s+=A[i][j]*c[j];
		}
		f[i]-=s;
	}*/
		
	/*printf("\n Nevyazka \n");

	for (int i = 0; i < n+1; i++)
	{
		printf("%*.*f ", 4, 5, f[i]);
	}*/

	//deleteData(A,n);
	system("pause");
	return 0;
}


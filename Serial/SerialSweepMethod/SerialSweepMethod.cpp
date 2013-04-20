#define _CRT_SECURE_NO_WARNINGS


#include "stdafx.h"
#include "stdio.h";
#include "conio.h";
#include "stdlib.h";
#include <time.h>
#include <omp.h>
#include <iostream>
#include <math.h>

double *a, *b, *c, *f, *x, *func;

int n = 50;
double left = -10;
double right = 10;
double step;
double halfstep;
float mainDiag;
float sideDiag;

void deleteData(int n)
{
	delete[] f;
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] x;
	delete[] func;
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

void makeMatrix(int n)
{
	mainDiag = 6;
	sideDiag = 1;

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

void serialSweepMethod(double*&f, int size, double* &x )
{
	int n = size;

	double *alpha = new double[n];
	double * beta = new double[n];

	alpha[1]=-sideDiag/mainDiag;
	beta[1]=f[0]/mainDiag;

	for (int i = 1; i < n-1; i++)
	{
		alpha[i+1]= -sideDiag/(sideDiag*alpha[i] + mainDiag);	
		beta[i+1] = (f[i] - sideDiag*beta[i])/(sideDiag*alpha[i] + mainDiag);
		
	}

	x[n-1]=(f[n-1] - sideDiag*beta[n-1])/(sideDiag*alpha[n-1] + mainDiag);

	for (int i = n-2; i >=0; i--)
	{
		x[i]=alpha[i+1]*x[i+1] + beta[i+1];
	}
	 
	delete[] alpha; 
	delete[] beta;

}

int main(int argc, char **argv)
{
	n = 1000000; //Размерность сетки
	step = (right - left)/n;
	halfstep = step/2;
	

	x = new double[n+2];
	func = new double[n+2];

	for (int i=0; i<n+1; i++)
	{
		x[i] = left + i*step;
		func[i] = calcFunction(x[i]);
	}

	makeMatrix(n);

	a = new double[n+1];
	b = new double[n+1];
	c = new double[n+1];

	//printMatrix(A, f, n+1);

	serialSweepMethod(f, n+1, c);

	for (int i=0; i<n; i++)
	{
		b[i] = (halfstep/2)*(c[i] - c[i+1]) + deltaf(i)/step;
		a[i] = halfstep*(b[i] - halfstep*c[i]) + func[i];
	}

	b[n] = b[n-1] + calcSecondDer(right)*halfstep;
	a[n] = func[n];

	float inaccuracy = 0; //Погрешность
	float inaccuracyDer = 0; //Погрешность производной
	float littleStep = step/4;

	for (double x = left; x <=right; x += littleStep)
	{
		double tmp = abs(calcFunction(x) - calcSpline(x));
		if (tmp > inaccuracy) inaccuracy = tmp;

		tmp = abs(calcDer(x) - calcSplineDer(x));
		if (tmp > inaccuracyDer) inaccuracyDer = tmp;
	}

	printf("%*.*f \n", 4, 10, inaccuracy);
	printf("%*.*f ", 4, 10, inaccuracyDer);


	deleteData(n);
	system("pause");
	return 0;
}


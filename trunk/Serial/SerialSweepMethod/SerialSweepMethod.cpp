// SerialSweepMethod.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h";
#include "conio.h";
#include "stdlib.h";
#include <time.h>
#include <omp.h>
#include <iostream>


double ** A, **copyA ;

double *x,*f , *copyf;

void generateData(double ** &A,int n)
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

}

void deleteData(double **A,int n)
{
	for (int i = 0; i < n; i++)
	{
		delete A[i];
		delete copyA[i];
		
	}
	delete A;
	delete copyA;
	delete f;
	delete copyf;
	delete x;
}

void printMatrix(double **A,int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%0*.*f ", 4, 2,A[i][j]);
		}
		printf("\n");
	}
}

void sub(double** &A,int size ,double coef,int j,int k)//вычитание j строки, умноженной на coef из k
{
			for (int i = 0; i < size; i++)
			{
				A[k][i]-=coef * A[j][i];
			}
			f[k]-=coef * f[j];
}

void deleteUnderDiag(double** &A,int size ,int k)//удаление поддиагонального эелемента k строки
{
	sub(A,size , A[k+1][k] /A[k][k] ,k,k+1);
}

void deleteUpperDiag(double** &A,int size ,int k)//удаление наддиагонального эелемента k строки
{
	sub(A,size , A[k-1][k] /A[k][k] ,k,k-1);
}



void serialSweepMethod(double** &A, double*&f ,int size , double* &x )
{
	int n=size;

	double *alpha = new double[n];
	double * beta = new double[n];

	alpha[1]=-A[0][1]/A[0][0];
	beta[1]=f[0]/A[0][0];

	for (int i = 1; i < n-1; i++)
	{
		alpha[i+1]= -A[i][i+1]/(A[i][i-1]*alpha[i]+A[i][i]);	
		beta[i+1] = (f[i]-A[i][i-1]*beta[i])/(A[i][i-1]*alpha[i]+A[i][i]);
		
	}

	x[n-1]=(f[n-1]-A[n-1][n-2]*beta[n]) /(A[n-1][n-2]*alpha[n-1]+A[n-1][n-1]);
	for (int i = n-2; i >=0; i--)
	{
		x[i]=alpha[i+1]* x[i+1] + beta[i+1];
	}
	 
	delete alpha; 
	delete beta;

}

int main(int argc, char **argv)
{
	int n=100;

	generateData(A,n);;	
	

	serialSweepMethod(A, f ,n,x);


		double s=-0;
		for (int i = 0; i < n; i++)
		{
			s=0;
			for (int j = 0; j < n; j++)
			{
				s+=copyA[i][j]*x[j];
			}
			copyf[i]-=s;
		}


		
		printf("\n Nevyazka \n");


		for (int i = 0; i < n; i++)
		{
			printf("%0*.*f ", 4, 2,copyf[i]);
		}
		


	deleteData(A,n);
	system("pause");
	return 0;
}


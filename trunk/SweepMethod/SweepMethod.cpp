// SweepMethod.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h";
#include "conio.h";
#include "stdlib.h";
#include <time.h>

double ** A;

double * alpha ,*beta,*x,*f;

void generateData(double ** &A,int n)
{
	A=new double*[n];
	for (int i = 0; i < n; i++)
	{
		A[i]=new double[n];
		for (int j = 0; j  < n; j ++)
		{
			A[i][j]=0;
		}
	}

	alpha = new double[n];
	beta = new double[n];
	x = new double[n];
	f = new double[n];


	srand (time(NULL));

	for (int i = 0; i < n; i++)
	{
		/*A[i][i]=1000+rand()%1000;
		if (i>0) A[i][i-1]=1 + rand()%100;
		if (i<n-1) A[i][i+1]=1 + rand()%100;*/
		A[i][i]=10;
		if (i>0) A[i][i-1]=1;
		if (i<n-1) A[i][i+1]=2;

		f[i]=10;

		alpha[i]=0;
		beta[i]=0;
	}

}

void deleteData(double **A,int n)
{
	for (int i = 0; i < n; i++)
	{
		delete A[i];
	}
	delete A;
	delete alpha;
	delete beta;
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
	for (int l = 0; l < size; l++)
	{
		if (l!=j)
		{
			double c= A[l][j]/A[j][j];
			for (int i = 0; i < size; i++)
			{
				A[l][i]-=c * A[j][i];
			}
			f[k]-=f[j];
		}
	}
}

void deleteUnderDiag(double** &A,int size ,int k)//удаление поддиагонального эелемента k строки
{
	sub(A,size , A[k+1][k] /A[k][k] ,k,k+1);
}

void deleteUpperDiag(double** &A,int size ,int k)//удаление наддиагонального эелемента k строки
{
	sub(A,size , A[k-1][k] /A[k][k] ,k,k-1);
}


void forwardStep(double** &A,int size )
{
	
}

void reverseStep(double** &A,int size )
{
	
}


int _tmain(int argc, _TCHAR* argv[])
{
	int n=3;
	generateData(A,n);
	printMatrix(A,n);

	alpha[1]=-A[0][1]/A[0][0];
	beta[1]=f[0]/A[0][0];

	for (int i = 1; i < n-1; i++)
	{
		alpha[i+1]= -A[i][i+1]/(A[i][i-1]*alpha[i]+A[i][i]);	
		beta[i+1] = (f[i]-A[i][i-1]*beta[i])/(A[i][i-1]*alpha[i]+A[i][i]);
		
	}
	printf("\n");
	for (int i = 0; i < n; i++)
	{
		printf("%0*.*f ", 4, 2,alpha[i]);
	}
	printf("\n");

	for (int i = 0; i < n; i++)
	{
		printf("%0*.*f ", 4, 2,beta[i]);
	}
	printf("\n");

	x[n-1]=(f[n-1]-A[n-1][n-2]*beta[n]) /(A[n-1][n-2]*alpha[n-1]+A[n-1][n-1]);
	for (int i = n-2; i >=0; i--)
	{
		x[i]=alpha[i+1]* x[i+1] + beta[i+1];
	}

	for (int i = 0; i < n; i++)
	{
		printf("%0*.*f ", 4, 2,x[i]);
	}
	printf("\n");

	deleteData(A,n);
	system("pause");
	return 0;
}


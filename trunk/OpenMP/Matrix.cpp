#include "Matrix.h"
#include <time.h>
#include "stdio.h";
#include "conio.h";
#include "stdlib.h";

Matrix::Matrix(int n)
{	A=new double*[n];
	for (int i = 0; i < n; i++)
	{
		A[i]=new double[n];
		for (int j = 0; j  < n; j ++)
		{
			A[i][j]=0;
		}
	}

	//alpha = new double[n];
	//beta = new double[n];
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


	}
}


Matrix::~Matrix(void)
{
}

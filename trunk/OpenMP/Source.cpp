#include "stdio.h";
#include "conio.h";
#include "stdlib.h";
#include <time.h>
#include <omp.h>
#include <iostream>
#include "RMatrix.h"

double *a, *b, *c, *f, *copyf, *x, *func;
double **A, **copyA;

int shift;

int n;
double left = -10;
double right = 10;
double step;
double halfstep;

RMatrix* matrix; 

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

	//alpha = new double[n];s
	//beta = new double[n];
	x = new double[n];
	f = new double[n];
	copyf = new double[n];

	srand (time(NULL));

	for (int i = 0; i < n; i++)
	{
		if(i>=shift) A[i][i]=1000+rand()%1000;
		if (i>shift) A[i][i-1]=1 + rand()%100;
		if ((i<n-1)&&(i>=shift)) A[i][i+1]=1 + rand()%100;
		
		if(i>=shift) matrix->Set(i,i,A[i][i]);   
		if (i>shift) matrix->Set(i,i-1,A[i][i-1]);
		if ((i<n-1)&&(i>=shift)) matrix->Set(i,i+1,A[i][i+1]);

		/*A[i][i]=10;
		if (i>0) A[i][i-1]=1;
		if (i<n-1) A[i][i+1]=2;
		f[i]=10;*/
		f[i]=rand()%100;
		
		copyA[i][i]=A[i][i];
		 copyA[i][i-1]= A[i][i-1];
		 copyA[i][i+1]= A[i][i+1];
		 copyf[i]=f[i];

		 x[i]=0;

		/*copyA[i][i]=10;
		if (i>0) copyA[i][i-1]=1;
		if (i<n-1) copyA[i][i+1]=2;
		copyf[i]=10;*/

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
	matrix = new RMatrix(n);

	matrix->Set(0, 0, 2);
	matrix->Set(n-1, n-1, 2);

	for (int i=1; i<n-1; i++)
	{
		matrix->Set(i, i, 6);

		if (i>1)
			matrix->Set(i, i-1, 1);

		if (i<n-2)
			matrix->Set(i, i+1, 1);
	}

	f = new double[n];

	f[0] = calcSecondDer(left);
	f[n-1] = calcSecondDer(right);

	for (int i=1; i<n-1; i++)
	{
		f[i] = delta2f(i-1)/(halfstep*halfstep);
		if (i == 1) f[i] -= calcSecondDer(left)/2;
		if (i == n-2) f[i] -= calcSecondDer(right)/2;
	}

}

/*
void sub(double** &A,int size ,double coef,int j,int k)//вычитание j строки, умноженной на coef из k
{
			for (int i = 0; i < size; i++)
			{
				A[k][i]-=coef * A[j][i];

				matrix->Set(k,i, matrix->Get(k,i) - coef*matrix->Get(j,i));


			}
			f[k]-=coef * f[j];
}*/

void subMatrix(int size,double coef,int j,int k)
{
			for (int i = 0; i < size; i++)
			{
				matrix->Set(k,i, matrix->Get(k,i) - coef*matrix->Get(j,i));
			}
			f[k]-=coef * f[j];
}


void deleteUnderDiag(int size ,int k)//удаление поддиагонального эелемента k строки
{
	subMatrix(size,   matrix->Get(k+1,k) /matrix->Get(k,k) ,k,k+1);

}

void deleteUpperDiag(int size ,int k)//удаление наддиагонального эелемента k строки
{
	subMatrix(size,   matrix->Get(k-1,k) /matrix->Get(k,k) ,k,k-1);
}


void serialSweepMethod(double** &A, double*&f,int size , double* &x )
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

	x[n-1]=(f[n-1]-A[n-1][n-2]*beta[n-1]) /(A[n-1][n-2]*alpha[n-1]+A[n-1][n-1]);
	for (int i = n-2; i >=0; i--)
	{
		x[i]=alpha[i+1]* x[i+1] + beta[i+1];
	}
	 
	delete alpha; 
	delete beta;

}

//int main(int argc, char **argv)
void sweepMethod(int size)
{
	int n = size;
	int p=2; // количество потоков
	int m=n/p;

	shift = n%p;

	shift =2;

	generateData(A,n);
	
	printMatrix(A,n);
	printf("\n");

	serialSweepMethod(A, f ,n,x);

	omp_set_num_threads(p);

    #pragma omp parallel 
    {
		int k=omp_get_thread_num();
		int startIndex=k*m;
		int	endIndex=(k+1)*m-1;

		if (k==0) startIndex = shift;

		for (int i = startIndex; i <= endIndex ; i++)
		{
			if (i!=n-1) deleteUnderDiag(n,i);
			//#pragma omp barrier
		}

		for (int i = endIndex-1; i >= startIndex-1 ; i--)
		{
			if (i>shift) deleteUpperDiag(n,i);
			//#pragma omp barrier
		}
    }

	printMatrix(A,n);

		double** B = new double *[p]; 
		double * fB=new double[p];
		double * xB=new double[p];

		for (int i = 0; i < p; i++)
		{
			B[i]=new double[p];
		}

		for (int i = 0; i < p; i++)
		{
			for (int j = 0; j  < p; j ++)
				B[i][j]=matrix->Get((i+1)*m-1,(j+1)*m-1);//A[(i+1)*m-1][(j+1)*m-1];

			fB[i]=f[(i+1)*m-1];
		}

		serialSweepMethod(B,fB,p,xB);

		for (int i = 0; i < p; i++)
		{
			x[(i+1)*m-1] = xB[i];
		}

		printMatrix(B,p);
		
    #pragma omp parallel 
    {
		int k=omp_get_thread_num();
		int startIndex=k*m;
		int	endIndex=(k+1)*m-1;

		double x1=xB[ (k*m-1)/p];
		double x2 =xB [ ((k+1)*m-1)/p];


		int index1= k*m-1;
		int index2 = (k+1)*m-1;

		//#pragma omp barrier
		for (int i = startIndex; i < endIndex ; i++)
		{
			x[i]= -matrix->Get(i,index2)*x2+f[i];
			if (index1>=0) x[i]-= matrix->Get(i,index2)*x1;
			if(matrix->Get(i,i)!=0) x[i]/=(double) matrix->Get(i,i);
		}
	}

		for (int i = 0; i < p; i++)
		{
			delete B[i];
		}

		delete B;
		delete xB;

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

}


int main(int argc, char **argv)
{
	n = 4;

	step = (right - left)/n;
	halfstep = step/2;

	x = new double[n+2];
	func = new double[n+2];

	for (int i=0; i<n+1; i++)
	{
		x[i] = left + i*step;
		func[i] = calcFunction(x[i]);
	}

	int mSize = n + 1;

	//int n = size;
	//int p = 2; // количество потоков
	//int m = mSize/p;

	//shift = n%p;

	makeMatrix(mSize);

	a = new double[n+1];
	b = new double[n+1];
	c = new double[n+1];

	for (int i=0;i<mSize;i++)
	{
		for (int j=0;j<6;j++)
			printf("%4.2f  ", matrix->Get(i, j));

		printf("%4.2f\n", f[i]);
	}

	sweepMethod(mSize);
}
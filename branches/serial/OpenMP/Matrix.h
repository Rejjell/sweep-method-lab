#pragma once
class Matrix
{
public:

	
	double ** A;
	double *x,*f;

	Matrix(int n);
	~Matrix(void);
};


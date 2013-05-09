#include "RMatrix.h"
#include <time.h>
#include "stdio.h";
#include "conio.h";
#include "stdlib.h";


RMatrix::RMatrix(int n)
{
	lines = new List[n];
	size = n;
}

RMatrix::~RMatrix()
{
	delete[] lines;
}

Element* RMatrix::Find(int i, int j)
{
	List &currentLine = lines[i];
	ListIterator itr = currentLine.begin();

	if (!currentLine.empty())
	{
		for (ListIterator itr = currentLine.begin(); itr != currentLine.end(); itr++)
			if (itr->GetColumn() == j)
				return &*itr;
	}

	return nullptr;
}

double RMatrix::Get(int i, int j)
{
	Element* findResult = Find(i, j);

	if (findResult != nullptr)
		return findResult->GetValue();

	else
		return 0;
}

void RMatrix::Set(int i, int j, double value)
{
	Element* findResult = Find(i, j);
		
	if (findResult != nullptr)
		findResult->SetValue(value);
	else
	{
		Element newElement(j, value);
		lines[i].push_back(newElement);
	}
}

void RMatrix::SetZero(int i, int j)
{
	Element* findResult = Find(i, j);
		
	if (findResult != nullptr)
		lines[i].remove(*findResult);
}



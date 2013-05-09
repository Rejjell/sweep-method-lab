#include <list>
#include "Element.h"

typedef std::list<Element>::iterator ListIterator;
typedef std::list<Element> List;

class RMatrix
{
public:
	RMatrix(int n);
	~RMatrix(void);

	double Get(int i, int j);
	void Set(int i, int j, double value);
	void SetZero(int i, int j);

private:
	List* lines;
	int size;

	Element* RMatrix::Find(int i, int j);

};


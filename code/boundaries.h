#pragma once

/* class Boundary
{
protected:
	unsigned int size = 0;
	unsigned int *indexesOfBoundaries = (unsigned int *)malloc(sizeof(unsigned int) * size);

public:
	Boundary(unsigned int *indexes, unsigned int s = 0);
	~Boundary();
	bool isOnBoundary(unsigned int index);
	virtual void streamBoundary(float *f_new, float *f_old, unsigned int index_new, unsigned int index_old, int dir);
}; */

class BounceBack /* : public Boundary */
{
private:
	unsigned int NX = 0;
	unsigned int NY = 0;
	bool *map;

public:
	BounceBack(bool *boolMap, unsigned int x, unsigned int y)
		/* : Boundary(indexes, s){} */;
	~BounceBack();

	bool isOnBoundary(unsigned int index);
};

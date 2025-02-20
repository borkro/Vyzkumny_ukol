#include "boundaries.h"

/* Boundary::Boundary(unsigned int *indexes, unsigned int s = 0)
{
	size = s;
	indexesOfBoundaries = indexes;
}

Boundary::~Boundary()
{
	free(indexesOfBoundaries);
}

bool Boundary::isOnBoundary(unsigned int index)
{
	for (unsigned int j = 0; j < size; j++)
	{
		if (indexesOfBoundaries[j] == index)
		{
			return true;
		}
	}
	return false;
} */

BounceBack::BounceBack(bool *boolMap, unsigned int x, unsigned int y)
	: map(boolMap), NX(x), NY(y)
{
}

BounceBack::~BounceBack()
{
	NX = NY = 0;
	free(map);
}

bool BounceBack::isOnBoundary(unsigned int index)
{
	return map[index];
}

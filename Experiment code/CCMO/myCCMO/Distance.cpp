/* shift-based similarity degree calculation routine */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

// Euclidean distance for minimization problems
double Euclidean_Distance(individual *ind1,individual *ind2)
{
	int i;
	double d = 0.0;
	double diff;    //Auxiliar var

	for (i=0; i<nobj; i++)
	{
		diff = ind1->obj[i] - ind2->obj[i];
		d += pow(diff, 2.0);
	}

	return (sqrt(d));
}
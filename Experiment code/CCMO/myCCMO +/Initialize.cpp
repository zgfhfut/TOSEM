/* Data initializtion routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Function to initialize a population randomly */
void initialize_pop (population *pop, int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        initialize_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to initialize an individual randomly */
void initialize_ind (individual *ind)
{
	int j;
	double sum = 0,low_T=0;
	for (j = 0; j<nreal; j++)
	{
		
		ind->xreal[j] = rndreal(min_realvar[j], max_realvar[j]);
		sum += ind->xreal[j];
		low_T += min_realvar[j];
		if (ind->xreal[j]<0)
			exit(-1);
	}
	if (sum > T_MAX_RESOURCE)
	{
	   double d = randomperc();
		for (j = 0; j < nreal; j++)
		{
         ind->xreal[j] = min_realvar[j] + (ind->xreal[j] - min_realvar[j]) * (T_MAX_RESOURCE - low_T) / (sum - low_T) * d;
		}
	}
	return;
} 
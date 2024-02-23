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

	for (j = 0; j<nreal; j++)
	{
		ind->xreal[j] = rndreal(min_realvar[j], max_realvar[j]);

		if (ind->xreal[j]<0)
			exit(-1);
	}

	return;
} 
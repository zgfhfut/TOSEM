/* environmental selection routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Routine to perform environmental selection */
void environmental_selection (population *mixed_set, population *archive_set)
{
	int i, j;
	int dom_size, nondom_size;
	int *dom_set, *nondom_set;
	double min;
	int mark, temp;

	dom_set = (int*)malloc((popsize+archsize)*sizeof(int));
	nondom_set = (int*)malloc((popsize+archsize)*sizeof(int));

	dom_size = 0;
	nondom_size = 0;
	for(i=0; i<popsize+archsize; i++)
	{
		if(mixed_set->ind[i].fitness < 1.0)
		{
			nondom_set[nondom_size] = i;
			nondom_size += 1;
		}
		else
		{
			dom_set[dom_size] = i;
			dom_size += 1;
		}
	}

	// put some dominated individuals into the archive set when the number of the nondominated individuals is not larger than archive size
	if (nondom_size <= archsize)
	{
		for (i=0; i<nondom_size; i++)
			copy_ind(&(mixed_set->ind[nondom_set[i]]), &(archive_set->ind[i]));

		for (i=0; i<archsize-nondom_size; i++)
		{
			min = mixed_set->ind[dom_set[i]].fitness;
			mark = i;
			for (j=i+1; j<dom_size; j++)
			{
				if (min > mixed_set->ind[dom_set[j]].fitness)
				{
					mark = j;
					min = mixed_set->ind[dom_set[j]].fitness;
				}
			}
			temp = dom_set[i];
			dom_set[i] = dom_set[mark];
			dom_set[mark] = temp;

			copy_ind(&(mixed_set->ind[dom_set[i]]), &(archive_set->ind[i+nondom_size]));
		}
	}
	else// truncate the non-dominated individuals when their number is larger than archive size
	{
		truncation (mixed_set, nondom_set, nondom_size, archive_set);
	}

	free(dom_set);
	free(nondom_set);

	return;
}


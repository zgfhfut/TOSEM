/* Routine for mergeing two populations */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Routine to merge two populations into one */
void merge(population *pop1, population *pop2, population *pop3, int size1, int size2)
{
	int i, k;
	for (i = 0; i<size1; i++)
	{
		copy_ind(&(pop1->ind[i]), &(pop3->ind[i]));
	}
	for (i = 0, k = size1; i<size2; i++, k++)
	{
		copy_ind(&(pop2->ind[i]), &(pop3->ind[k]));
	}
	return;
}

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void copy_ind (individual *ind1, individual *ind2)
{
	int i;

	ind2->constr_violation = ind1->constr_violation;


	for (i = 0; i<nreal; i++)
		ind2->xreal[i] = ind1->xreal[i];

	for (i = 0; i<nobj; i++)
		ind2->obj[i] = ind1->obj[i];

	if (ncon != 0)
	{
		for (i = 0; i<ncon; i++)
			ind2->constr[i] = ind1->constr[i];
	}

	ind2->rank = ind1->rank;

	ind2->crowd_dist = ind1->crowd_dist;
   
    return;
}

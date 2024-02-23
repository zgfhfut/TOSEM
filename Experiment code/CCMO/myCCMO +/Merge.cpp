/* Routine for mergeing two populations */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Routine to merge two populations into one */
void merge(population *pop1, population *pop2, population *pop3, population *pop4, int size1, int size2, int size3)
{
    int i, k;
    for (i=0; i<size1; i++)
    {
        copy_ind (&(pop1->ind[i]), &(pop4->ind[i]));
    }
    for (i=0, k=size1; i<size2; i++, k++)
    {
        copy_ind (&(pop2->ind[i]), &(pop4->ind[k]));
    }
	for (i = 0, k = size1+size2; i<size3; i++, k++)
	{
		copy_ind(&(pop3->ind[i]), &(pop4->ind[k]));
	}

    return;
}

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void copy_ind (individual *ind1, individual *ind2)
{
	int i;

	ind2->constr_violation = ind1->constr_violation;

	ind2->fitness = ind1->fitness;

	for (i = 0; i<nreal; i++)
		ind2->xreal[i] = ind1->xreal[i];

	for (i = 0; i<nobj; i++)
		ind2->obj[i] = ind1->obj[i];

	if (ncon != 0)
	{
		for (i = 0; i<ncon; i++)
			ind2->constr[i] = ind1->constr[i];
	}

   
    return;
}

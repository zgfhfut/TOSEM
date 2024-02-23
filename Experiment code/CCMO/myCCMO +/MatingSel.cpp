/* Mating selection routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"


/* Routine for tournament selection and crossover, it creates new population from the archive set by performing tournament selection and the crossover */
void mating_selection(population *old_pop, population *new_pop)
{
	int *a1, *a2;
	int temp;
	int i;
	int rand;
	individual *parent1, *parent2;
	a1 = (int *)malloc(popsize*sizeof(int));
	a2 = (int *)malloc(popsize*sizeof(int));
	for (i = 0; i < popsize; i++)
	{
		a1[i] = a2[i] = i;
	}
	for (i = 0; i < popsize; i++)
	{
		rand = rnd(i, popsize - 1);
		temp = a1[rand];
		a1[rand] = a1[i];
		a1[i] = temp;
		rand = rnd(i, popsize - 1);
		temp = a2[rand];
		a2[rand] = a2[i];
		a2[i] = temp;
	}
	for (i = 0; i < popsize/2; i += 4)
	{
		parent1 = tournament(&old_pop->ind[a1[i]], &old_pop->ind[a1[i + 1]]);
		parent2 = tournament(&old_pop->ind[a1[i + 2]], &old_pop->ind[a1[i + 3]]);
		crossover(parent1, parent2, &new_pop->ind[i], &new_pop->ind[i + 1]);
		parent1 = tournament(&old_pop->ind[a2[i]], &old_pop->ind[a2[i + 1]]);
		parent2 = tournament(&old_pop->ind[a2[i + 2]], &old_pop->ind[a2[i + 3]]);
		crossover(parent1, parent2, &new_pop->ind[i + 2], &new_pop->ind[i + 3]);
	}
	free(a1);
	free(a2);

	return;
}

/* Routine for binary tournament */
individual* tournament(individual *ind1, individual *ind2)
{
	if (ind1->fitness<ind2->fitness)
	{
		return (ind1);
	}
	if (ind1->fitness>ind2->fitness)
	{
		return (ind2);
	}
	
	if ((randomperc()) <= 0.5)
	{
		return(ind1);
	}
	else
	{
		return(ind2);
	}
}

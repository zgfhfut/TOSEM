/* fitness assignment routine */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Function to assign fitness to a population */
void assign_fitness(population *pop,int size)
{
	int i, j;
	int *strength, *rawfitness;
	double **distance;
	double *D;
	int Kmin;

	strength = (int*)malloc(size*sizeof(int));
	rawfitness = (int*)malloc(size*sizeof(int));
	D = (double*)malloc(size*sizeof(double));
	distance = (double**)malloc(size*sizeof(double*));
	for(i=0; i<size; i++)
	{
		distance[i] = (double*)malloc(size*sizeof(double));
	}

	for(i=0; i<size; i++)
	{
		strength[i] = 0;
		rawfitness[i] = 0;
	}

	/*-------calculate the strength value------start---------*/
	for(i=0; i<size ;i++)
	{
		for(j=0; j<size; j++)
		{
			if (check_dominance(&(pop->ind[i]), &(pop->ind[j])) == 1)
				strength[i] += 1;
		}
	}
	/*-------calculate the strength value------end---------*/

	/*-------calculate the raw fitness------start---------*/
	for(i=0; i<size ;i++)
	{
		for(j=0; j<size; j++)
		{
			if (check_dominance(&(pop->ind[i]), &(pop->ind[j])) == -1)
				rawfitness[i] += strength[j];
		}
	}
	/*-------calculate the raw fitness------end---------*/

	/*-------calculate the distance between individuals in the population------start---------*/
	for(i=0; i<size; i++)
	{		
		for(j=0; j<size; j++)
		{
			if (i==j)
				distance[i][j] = INF;
			else
				distance[i][j] = Euclidean_Distance(&(pop->ind[i]), &(pop->ind[j]));
		}
	}
	/*-------calculate the distance between individuals in the population------end---------*/

	/*-------calculate the density D------start---------*/
	Kmin = (int)floor(sqrt((double)size));		// determine k value
	for(i=0; i<size; i++)
	{
		D[i] = 1.0 / (2.0 + findKmin(distance[i], Kmin, size));
	}
	/*-------calculate the density D------end---------*/

	for(i=0; i<size; i++)
		pop->ind[i].fitness = (double)rawfitness[i] + D[i];		// calculate the fitness of individuals in the population

	free(strength);
	free(rawfitness);
	free(D);
	for(i=0; i<size; i++)
		free(distance[i]);
	free(distance);

	return;
}

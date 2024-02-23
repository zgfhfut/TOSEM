/* Routine for evaluating population members  */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop(population *pop, int size)
{
	int i;
	for (i = 0; i<size; i++)
	{
		evaluate_ind(&(pop->ind[i]));
		currenteval++;
	}
	return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind(individual *ind)
{
	int j;

	test_problem(ind->xreal, ind->obj, ind->constr);
	
	ind->constr_violation = 0.0;
	for (j = 0; j<ncon; j++)
	{
		if (ind->constr[j]<0.0)
		{
			ind->constr_violation += ind->constr[j];
		}
	}

	ind->rank = 0;
	ind->crowd_dist = 0;

	return;
}

double GetFeasible_Ratio(population *pop, int size)
{
	double result = 0.0;
	for (int i = 0; i < size; i++)
	{
		if (pop->ind[i].constr_violation== 0.0) {
			result += 1.0;
		}
	}

	return result / size; 



}
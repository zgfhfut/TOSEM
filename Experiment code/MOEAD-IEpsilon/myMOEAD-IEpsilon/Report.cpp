/* Routines for putting population data into files */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"
# include <cassert>

/* Output of individuals in the objective space */
void report_objective(population *pop, int size, FILE *fpt)
{
	int i, j;
	for (i = 0; i<size; i++)
	{
		if (pop->ind[i].constr_violation==0.0)
		{
			for (j = 0; j < nobj; j++)
			{
				
				fprintf(fpt, "%6.6f\t", pop->ind[i].obj[j]);
			}
			fprintf(fpt, "\n");
		}
		
	}
	//fprintf(fpt, "\n");
	return;
}

/* Output of individuals in the decision space */
void report_variable (population *pop, int size, FILE *fpt)
{
	int i, j;
	for (i = 0; i<size; i++)
	{
		for (j = 0; j<nreal; j++)
		{
			fprintf(fpt, "%6.6f\t", pop->ind[i].xreal[j]);
		}
		fprintf(fpt, "\n");
	}
	fprintf(fpt, "\n");

    return;
}

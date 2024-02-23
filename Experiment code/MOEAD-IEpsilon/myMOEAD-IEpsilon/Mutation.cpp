/* Mutation routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Routine for real polynomial mutation of an individual */
void real_mutate_ind(individual *ind)
{
	int j;
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;

	for (j = 0; j<nreal; j++)
	{
		rnd = randomperc();
		if (rnd <= pmut_real)
		{
			y = ind->xreal[j];
			yl = min_realvar[j];
			yu = max_realvar[j];
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);
			rnd = randomperc();
			mut_pow = 1.0 / (eta_m + 1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd + (1.0 - 2.0*rnd)*(pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0*(1.0 - rnd) + 2.0*(rnd - 0.5)*(pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			y = y + deltaq * (yu - yl);
			if (y<yl)
			{
				y = yl;
			}
			if (y>yu)
			{
				y = yu;
			}
			ind->xreal[j] = y;
		}
	}

	return;
}

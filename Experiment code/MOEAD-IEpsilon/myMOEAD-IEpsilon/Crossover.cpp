/* Crossover routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

void DECrossover(individual * child, individual *parents) {


	int jrand;

	jrand = rnd(0, nreal - 1);

	// STEP 4. Checking the DE variant
	for (int j = 0; j < nreal; j++) {
		if (randomperc() < pcross_real || j == jrand) {
			double value;
			value = parents[2].xreal[j] + F_ * (parents[0].xreal[j] - parents[1].xreal[j]);

			if (value < min_realvar[j]) {
				value = min_realvar[j];
			}
			if (value > max_realvar[j]) {
				value = max_realvar[j];
			}

			child->xreal[j] = value;
		}
		else {
			double value;
			value = parents[3].xreal[j];
			child->xreal[j] = value;
		} // else

		  //	printf_s("child->xreal[j]=%6.6f\n", child->xreal[j]);

	} // for


} // execute
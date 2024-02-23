/* Crossover routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

void DECrossover(individual * child, individual *parents) {

	
	int jrand;
	double rand;
	double sum=0, low_T=0;
	jrand = rnd(0, nreal - 1);
	double lim[nreal];
	int count[nreal];
	// STEP 4. Checking the DE variant
	for (int j = 0; j < nreal; j++) {
		lim[j] = 0;
		count[j] = 0;
		if (randomperc() < pcross_real || j == jrand) {
			count[j] = 1;
			double value;
			value = parents[2].xreal[j] + F_ * (parents[0].xreal[j] - parents[1].xreal[j]);

			if (value < min_realvar[j]||value>max_realvar[j]) {
				value = repair_s(min_realvar[j], max_realvar[j], parents[3].xreal[j], value);
			}
		 	child->xreal[j] = value;
			if (child->xreal[j] > parents[3].xreal[j]) lim[j] = parents[3].xreal[j];
			else if (child->xreal[j] < parents[3].xreal[j]) lim[j] = min_realvar[j];
	}	
			else {
				double value;
				value = parents[3].xreal[j];
				child->xreal[j] = value;
				
			} 		
		} 

	repair_all(child, lim, count);
} // execute
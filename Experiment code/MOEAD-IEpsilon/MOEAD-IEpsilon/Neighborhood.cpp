
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/**
* initNeighborhood
*/
void initNeighborhood() {
	double * x = (double *)malloc(popsize*sizeof(double));
	int * idx = (int *)malloc(popsize*sizeof(int));

	for (int i = 0; i < popsize; i++) {
		// calculate the distances based on weight vectors
		for (int j = 0; j < popsize; j++) {
			x[j] = distVector(lambda_[i], lambda_[j],nobj);
			idx[j] = j;
			// cout << "x[" << j << "]: " << x[j] << ". idx[" << j << "]: " <<
			//    idx[j] << endl ;
		} // for

		// find 'niche' nearest neighboring subproblems
		minFastSort(x, idx, popsize, T_);
	//	neighborhood_[i] = new int[T_];
		for (int k = 0; k < T_; k++) {
			neighborhood_[i][k] = idx[k];
			//cout << "neg[ << i << "," << k << "]: " << neighborhood_[i][k] << endl;
		}
	} // for

	free(x);
	free(idx);

} // initNeighborhood


double distVector(double * vector1, double * vector2, int dim) {
	//int dim = vector1.size();
	double sum = 0;
	for (int n = 0; n < dim; n++) {
		sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
	}
	return sqrt(sum);
} // distVector


void minFastSort(double * x, int * idx, int n, int m) {

	for (int i = 0; i < m; i++) {
		for (int j = i + 1; j < n; j++) {
			if (x[i] > x[j]) {
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				int id = idx[i];
				idx[i] = idx[j];
				idx[j] = id;
			} // if
		}
	} // for
} // minFastSort


void randomPermutation(int * perm, int size) {
	int * index = (int *)malloc(size*sizeof(int));
	bool * flag = (bool *)malloc(size*sizeof(bool));

	for (int n = 0; n < size; n++) {
		index[n] = n;
		flag[n] = true;
	}

	int num = 0;
	while (num < size) {
		int start = rnd(0, size - 1);	
		while (true) {
			if (flag[start]) {
				perm[num] = index[start];
				flag[start] = false;
				num++;
				break;
			}
			if (start == (size - 1)) {
				start = 0;
			}
			else {
				start++;
			}
		}
	} // while

	free(index);
	free(flag);

} // randomPermutation
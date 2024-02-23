#include "stdafx.h"

# include "Global.h"
# include "Random.h"

#include <string>
#include <sstream>
#include <iostream> 
#include <fstream>

using namespace std;

void initUniformWeight() {
	if ((nobj == 2) && (popsize <= 300)) {
		for (int n = 0; n < popsize; n++) {
		//	lambda_[n] = (double *)malloc(nobj*sizeof(double));
			double a = 1.0 * n / (popsize - 1);
			lambda_[n][0] = a;
			lambda_[n][1] = 1 - a;
		} // for
	} // if
	else {
		ostringstream os;
		os << "data/Weight/" << "W" << nobj << "D_"
			<< popsize << ".dat";
		string dataFileName;
		dataFileName = os.str();

		// Open the file
		std::ifstream in(dataFileName.c_str());
		if (!in) {
			cout << "initUniformWeight: failed when reading from file: : " <<
				dataFileName << endl;
			//     exit(-1);
		} // if

		//int numberOfObjectives = 0;
		int i = 0;
		int j = 0;
		string aux;
		while (getline(in, aux)) {
			istringstream iss(aux);
			j = 0;
			// TODO: Check if number of tokens per line is equal to number of
			//       objectives
			//lambda_[i] = new double[nobj];
			while (iss) {
				string token;
				iss >> token;
				if (token.compare("") != 0) {
					double value = atof(token.c_str());
					lambda_[i][j] = value;
					//cout << "lambda[" << i << "," << j << "] = " << value << endl;
					j++;
				} // if
			} // while
			i++;
		} // while
		in.close();
	} // else
} // initUniformWeight

void initIdealPoint(population *pop) {
	for (int i = 0; i < nobj; i++) {
		z_[i] = 1.0e+30;		
	} // for

	for (int i = 0; i < popsize; i++) {
		updateReference(&pop->ind[i]);
	} // for
} // initIdealPoint

void updateReference(individual *ind) {
	for (int n = 0; n <nobj; n++) {
		if (ind->obj[n] < z_[n]) {
			z_[n] = ind->obj[n];			
		}
	}
} // updateReference
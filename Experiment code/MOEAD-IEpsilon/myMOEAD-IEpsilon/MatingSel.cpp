/* Mating selection routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include <vector>

# include "Global.h"
# include "Random.h"

using namespace std;
/**
* matingSelection
*/
void matingSelection(vector<int> &list, int cid, int size, int type) {

	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss;
	int r;
	int p;

	//ss = neighborhood_[cid].length;
	ss = T_;
	while ((int)list.size() < size) {
		if (type == 1) {
			r = rnd(0, ss - 1);
			p = neighborhood_[cid][r];
			//p = population[cid].table[r];
		}
		else {
			p = rnd(0, popsize - 1);
		}
		bool flag = true;
		for (int i = 0; i < (int)list.size(); i++) {
			if (list[i] == p) // p is in the list
			{
				flag = false;
				break;
			}
		}
		
		if (flag) {
			list.push_back(p);
		}
	}
} // matingSelection
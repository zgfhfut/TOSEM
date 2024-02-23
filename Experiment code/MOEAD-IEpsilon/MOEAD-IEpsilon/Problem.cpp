/* Test problem definitions */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include <string.h>

# include "Global.h"
# include "Random.h"

#define maxvar nreal
double Mis_tim[nreal];
double Y[maxvar];
double tmid[maxvar];
double c_min[maxvar];
double t_min[maxvar];
double C;
double dd = 0.0000001;

extern double u[maxvar];
extern double reli_thresh;
void Cal_res_lim(double r) {
	
	double sum_u = 0;/*所有组件预期访问次数之和*/
	//计算总预期访问次数
	for (int i = 0; i < nreal; i++)
	{
		sum_u += u[i];
	}

	//计算单个组件任务时间
	for (int i = 0; i < nreal; i++)
	{
		Mis_tim[i] = MISSION_TIME * u[i] / sum_u;
	}

}

void Cal_c()
{


	double C_min = 0.0;
	double max = DBL_MAX;
	double min = 0.0;
	double c[maxvar];
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;


	while (fabs(sum - log(Max_R)) > dd && (max - min) > dd)
	{

		sum = 0.0;
		for (int i = 0; i < maxvar; ++i)
		{

			sum += para_module[i][4] * Mis_tim[i] / (para_module[i][2] - para_module[i][3] - mid * para_module[i][1] * Mis_tim[i]);

		}

		double aa = log(Max_R);
		if (sum < log(Max_R))
		{
			min = mid;
			mid = (min + max) / 2;
		}
		if (sum > log(Max_R))
		{
			max = mid;
			mid = (min + max) / 2;
		}

	}

	for (int j = 0; j < maxvar; ++j)
	{
		Y[j] = para_module[j][4] * Mis_tim[j] / (para_module[j][2] - para_module[j][3] - mid * para_module[j][1] * Mis_tim[j]);
		tmid[j] = -1 / para_module[j][1] * log(-Y[j] / (para_module[j][0] * para_module[j][1] * Mis_tim[j]));
		c[j] = (para_module[j][2] - para_module[j][3]) * para_module[j][0] * (1 - exp(-para_module[j][1] * tmid[j])) + para_module[j][3] * para_module[j][0] + para_module[j][4] * tmid[j];
		C_min += c[j];

	}
	C = C_min * 1.20;

}

double cal_reliab(int ind_module, double m_res) {

	double rel = 0, rel1 = 0;
	rel = para_module[ind_module][0] * para_module[ind_module][1] * exp(-para_module[ind_module][1] * m_res);
	rel1 = exp(-rel * Mis_tim[ind_module]);


	return rel1;
}

double cal_cost(int ind_module, double m_res) {

	double cost = 0;

	cost = (para_module[ind_module][2] - para_module[ind_module][3]) * para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res)) + para_module[ind_module][3] * para_module[ind_module][0] + para_module[ind_module][4] * m_res;
	return cost;
}

void test_problem(double* xreal, double* obj, double* constr)
{

	int j, k, p;
	/********************************************************************/
	j = 0;
	double sum1 = 1.0;
	double sum2=1.0;

	for (k = 0; k < nreal; k++)
	{
		if (xreal[k] < 0.0)
		{
			printf_s("error xreal[j] < 0.0,%d---%6.6f\n", k, xreal[k]);
		}
		sum2 = sum2 * cal_reliab(k, xreal[k]);
	}
	sum1 = sum1 * sum2;/*系统可靠性*/
	obj[0] = 1.0 - sum1;/*here we use 1.0-sum1 to minimize obj[0]*///reliability
	/********************************************************************/
	double sum3 = 0;/*系统初始成本为50.0*/
	double sum4 = 0;

	for (j = 0; j < nreal; j++)
	{
		sum3 = sum3 + cal_cost(j, xreal[j]);
		sum4 = sum4 + xreal[j];
	}
	obj[1] = sum3 ;
	obj[2] = sum4;
	/********************************************************************/

	//violate the constraint
	/********************************************************************/
	constr[0] = sum1 - reli_thresh;//reliability constraint
	constr[1] = (T_MAX_RESOURCE - sum4) / T_MAX_RESOURCE; // time constraint
	constr[2] = (C - sum3) / C;
	/********************************************************************/

	return;
}

/**
* @param individual
* @param id
* @param type
*/
void updateProblem(population *pop, individual * indiv, int id, int type) {
	// indiv: child solution
	// id: the id of current subproblem
	// type: update solutions in neighborhood (1) or whole population (otherwise)
	int size;
	int time;

	time = 0;

	if (type == 1) {
		//size = neighborhood_[id].length;
		size = T_;
	}
	else {
		//size = population_.size();
		size = popsize;
	}
	int * perm = (int *)malloc(size*sizeof(int));

	randomPermutation(perm, size);

	for (int i = 0; i < size; i++) {

		int k;

		if (type == 1) {
			k = neighborhood_[id][perm[i]];
		}
		else {
			k = perm[i];      // calculate the values of objective function regarding
			// the current subproblem
		}

		double f1, f2, con1, con2;

		f1 = fitnessFunction(&pop->ind[k], lambda_[k]);
		f2 = fitnessFunction(indiv, lambda_[k]);

		con1 = abs(pop->ind[k].constr_violation);
		con2 = abs(indiv->constr_violation);

		//printf_s("con1 = %6.6f;; con2 = %6.6f;; epsilon_k_ = %6.6f\n", con1, con2, epsilon_k_);
		// use epsilon constraint method

		if (con1 <= epsilon_k_ && con2 <= epsilon_k_) {
			if (f2 < f1) {
				copy_ind(indiv, &pop->ind[k]);
				time++;
			}
		}
		else if (con2 == con1) {
			if (f2 < f1) {
				copy_ind(indiv, &pop->ind[k]);
				time++;
			}
		}
		else if (con2 < con1) {
			copy_ind(indiv, &pop->ind[k]);
			time++;
		}

		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= nr_) {
			//delete indiv;
			free(perm);
			return;
		}
	}

	//delete indiv;
	free(perm);

} // updateProblem

/**
* fitnessFunction
*/
double fitnessFunction(individual * indiv, double * lambda) {
	double fitness;
	fitness = 0.0;

	
	double maxFun = -1.0e+30;

	for (int n = 0; n < nobj; n++) {
		double diff = fabs(indiv->obj[n] - z_[n]);

		double feval;
		if (lambda[n] == 0) {
			feval = 0.0001 * diff;
		}
		else {
			feval = diff * lambda[n];
		}
		if (feval > maxFun) {
			maxFun = feval;
		}
	} // for

	fitness = maxFun;
	
	return fitness;
} // fitnessEvaluation
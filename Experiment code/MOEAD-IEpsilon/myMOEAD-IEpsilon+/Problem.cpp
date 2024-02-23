/* Test problem definitions */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include <string.h>
#include<iostream>
# include "Global.h"
# include "Random.h"
using namespace std;
#define maxvar nreal
#define nvar nreal

double C;
double reli_thresh;
double Mis_tim[nreal];/*构件分配到的任务时间*/
double lim_r[nreal][2];/*Limits of variable in array*/
double lim_t[nreal][2];
double M;
double T_low_sys;
double zz, zz2, zz1;
int q1, q2;/*系统测试资源总量下限时单个构件测试资源大于0的个数*/
int Z_0[nreal];/*记录*/
double Z2[nreal];/*系统测试资源总量下限时单个构件测试资源*/
double e[nreal];
double Zmax[nreal];
double reciprocal_sum = 0.0;/*故障检测率倒数之和*/
double Y[nreal];
double tmid[nreal];//测试资源中间值；
double c_min[nreal];//每个模块成本的极小值点
double t_min[nreal];//每个模块成本的极小值点对应的测试资源

void Cal_x();
void Cal_res_lim(double r);/*计算构件测试资源上下限*/
void Cal_res_lim3(int e, double r);/*计算构件测试资源的中间参数T_middle,为了求最终的构件测试资源上下限*/
double Cal_res_lim4(double left, double right, int module, double r);/*求构件测试资源上限或下限*/
double fun4(double x, int module, double r);
double fun3(double x);

double cal_z(int  g, double f);//除去第g个组件的z值
double cal_t1(int i);//每个模块对应c最小值所需的测试资源
double cal_t2(int i);//每个模块可靠性达到MAX_R所需要的测试资源
double cal_sumcost(int  g, double f);//求除去g模块的n-1个模块所花费的成本最小值
double cal_lowtime(int g);//成本和可靠性约束下每个模块测试资源下限
double cal_hightime(int g);//成本和可靠性约束下每个模块测试资源上限
void Cal_res_t();

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
void Cal_x()
{
	double sum_u = 0.0;/*所有组件预期访问次数之和*/
	//计算总预期访问次数
	for (int i = 0; i < maxvar; i++)
	{
		sum_u += u[i];
	}
	//计算单个组件任务时间
	for (int i = 0; i < maxvar; i++)
	{
		Mis_tim[i] = MISSION_TIME * u[i] / sum_u;
	}
}


double cal_z(int  g, double f)//除去第g个组件的z值
{
	double f1 = cal_reliab(g, f);
	double max = DBL_MAX;
	double min = 0.0;
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;
	double f2 = log(Max_R) - log(f1);
	while (fabs(sum - f2) > EPSINON && (max - min) > EPSINON)
	{

		sum = 0.0;
		for (int i = 0; i < maxvar; ++i)
		{
			if (i == g) continue;
			sum += para_module[i][4] * Mis_tim[i] / (para_module[i][2] - para_module[i][3] - mid * para_module[i][1] * Mis_tim[i]);

		}


		if (sum < f2)
		{
			min = mid;
			mid = (min + max) / 2;
		}
		if (sum > f2)
		{
			max = mid;
			mid = (min + max) / 2;
		}

	}
	return mid;
}

double cal_t1(int i)//每个模块对应c最小值的测试资源
{
	double t;
	t = -1 / para_module[i][1] * log(para_module[i][4] / (para_module[i][0] * para_module[i][1] * (para_module[i][3] - para_module[i][2])));

	return t;
}

double cal_t2(int i)//每个模块可靠性达到MAX_R所需要的测试资源
{
	double a = -log(Max_R) / (para_module[i][0] * para_module[i][1] * Mis_tim[i]);
	double t = -1 / para_module[i][1] * log(a);
	return t;
}
double cal_sumR(int  g, double f)//求n-1个模块总可靠性
{
	double c;
	double sum = 0.0, sum1 = 1.0;
	double mid = cal_z(g, f);
	for (int k = 0; k < maxvar; ++k)
	{
		if (k == g) continue;
		Y[k] = para_module[k][4] * Mis_tim[k] / (para_module[k][2] - para_module[k][3] - mid * para_module[k][1] * Mis_tim[k]);
		sum1 = sum1 * exp(Y[k]);

	}

	return sum1;
}

double cal_sumcost(int  g, double f)//求n-1个模块所花费的成本最小值
{
	double c;
	double sum = 0.0, sum1 = 1.0;
	double mid = cal_z(g, f);
	for (int k = 0; k < maxvar; ++k)
	{
		if (k == g) continue;
		Y[k] = para_module[k][4] * Mis_tim[k] / (para_module[k][2] - para_module[k][3] - mid * para_module[k][1] * Mis_tim[k]);
		sum1 = sum1 * exp(Y[k]);
		c = para_module[k][2] * para_module[k][0] + (para_module[k][2] - para_module[k][3]) * Y[k] / (para_module[k][1] * Mis_tim[k])
			+ para_module[k][4] * (-1 / para_module[k][1] * log(-Y[k] / (para_module[k][0] * para_module[k][1] * Mis_tim[k])));
		sum += c;
	}

	return sum;
}

double cal_lowtime(int g)//成本和可靠性约束下每个模块测试资源下限
{
	double t = cal_t1(g);
	double t1 = cal_t2(g);
	double min, max, mid;
	if (t < t1)//当模块成本极小值对应的测试资源小于模块可靠性达到MAX_R对应的测试资源时,模块的测试资源下限在(t1,tmid[g])
	{
		min = t1;
		max = tmid[g];
		mid = (min + max) / 2;
		while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > EPSINON && max - min > EPSINON)
		{
			if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
			if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

		}
	}

	if (t > t1)//当模块成本极小值对应的测试资源大于模块可靠性达到MAX_R对应的测试资源时
	{
		
			min = t1;
			max = t;
			mid = (min + max) / 2;
			while (fabs(C - cal_cost(g, mid) + cal_sumcost(g, mid)) > EPSINON && max - min > EPSINON)
			{
				if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
				if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

			}
		
		if (cal_cost(g, t) + cal_sumcost(g, t) > C)//判断t对应的成本加其他模块的总最低成本是否大于C，如果大于则下限范围在(t,T)
		{
			min = t;
			max = tmid[g];
			mid = (min + max) / 2;
			while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > EPSINON && max - min > EPSINON)
			{
				if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
				if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

			}

		}
	}

	return  mid;

}
double cal_hightime(int g)//成本和可靠性约束下每个模块测试资源上限
{
	double t = cal_t1(g);
	double min, max, mid;
	if (t <= tmid[g])
	{
		min = tmid[g];
		max = T_MAX_RESOURCE;
		mid = (min + max) / 2;
		while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > EPSION && max - min > EPSION)
		{
			if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  min = mid; mid = (max + min) / 2;
			if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  max = mid; mid = (max + min) / 2;

		}
	}
	;
	if (t > tmid[g])
	{
		if (cal_cost(g, t) + cal_sumcost(g, t) < C)
		{
			min = t;
			max = T_MAX_RESOURCE;
			mid = (min + max) / 2;
			while (fabs(C - cal_cost(g, mid) + cal_sumcost(g, mid)) > EPSION && max - min > EPSION)
			{

				if (cal_cost(g, mid) + cal_sumcost(g, mid) < C) min = mid; mid = (max + min) / 2;
				if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  max = mid; mid = (max + min) / 2;

			}

		}
		if (cal_cost(g, t) + cal_sumcost(g, t) > C)
		{
			min = tmid[g];
			max = t;
			mid = (min + max) / 2;
			while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > EPSION && max - min > EPSION)
			{
				if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  min = mid; mid = (max + min) / 2;
				if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  max = mid; mid = (max + min) / 2;

			}

		}
	}

	return  mid;

}
void Cal_res_t()
{
	Cal_x();
	
	double C_min = 0.0;
	double max = DBL_MAX;
	double min = 0.0;
	double c[maxvar];
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;


	while (fabs(sum - log(Max_R)) > EPSINON && (max - min) > EPSINON)
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

	for (int k = 0; k < maxvar; ++k)
	{

		lim_r[k][0] = cal_lowtime(k);//lim_t[k][0] 是成本和可靠性约束的下限
		if (lim_r[k][0] < 0) lim_r[k][0] = 0;
		lim_r[k][1] = cal_hightime(k);//lim_t[k][1] 是成本和可靠性约束的上限
		
	}
	for (int i = 0; i < maxvar; ++i)
	{
		cout << "(" << lim_r[i][0] << "," << lim_r[i][1] << ")" << endl;
	}
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
	obj[1] = sum3;
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

double repair_s(double low, double high, double l, double r)
{
	double rand = randomperc();
	double res;
	if (r < low)
		res = rand * (l - low) + low;
	if (r > high)
		res = rand * (high - l) + l;
	return res;
}

void  repair_all(individual* ind, double lim[], int count[])
{

	double sum = 0, sum1 = 0;
	double low_T = 0;
	for (int j = 0; j < nreal; j++)
	{
		if (count[j] == 1)
		{
			low_T += lim[j];
			sum1 += ind->xreal[j];
		}
		sum += ind->xreal[j];
	}

	double	T = T_MAX_RESOURCE - (sum - sum1);

	if (sum > T_MAX_RESOURCE)
	{


		for (int k = 0; k < nvar; k++)
		{
			double rnd = randomperc();
			if (count[k] == 1)
				ind->xreal[k] = lim[k] + (ind->xreal[k] - lim[k]) * rnd * ((T - low_T) / (sum1 - low_T));
		}

	}
	sum = 0.0;
	for (int j = 0; j < nreal; j++)
	{
		sum += ind->xreal[j];
	}
	if (sum > T_MAX_RESOURCE + EPSION) {
	
		printf("error：修复越界！！！");
		system("pause");
	}
	
}
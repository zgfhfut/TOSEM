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
double Mis_tim[nreal];/*�������䵽������ʱ��*/
double lim_r[nreal][2];/*Limits of variable in array*/
double lim_t[nreal][2];
double M;
double T_low_sys;
double zz, zz2, zz1;
int q1, q2;/*ϵͳ������Դ��������ʱ��������������Դ����0�ĸ���*/
int Z_0[nreal];/*��¼*/
double Z2[nreal];/*ϵͳ������Դ��������ʱ��������������Դ*/
double e[nreal];
double Zmax[nreal];
double reciprocal_sum = 0.0;/*���ϼ���ʵ���֮��*/
double Y[nreal];
double tmid[nreal];//������Դ�м�ֵ��
double c_min[nreal];//ÿ��ģ��ɱ��ļ�Сֵ��
double t_min[nreal];//ÿ��ģ��ɱ��ļ�Сֵ���Ӧ�Ĳ�����Դ

void Cal_x();
void Cal_res_lim(double r);/*���㹹��������Դ������*/
void Cal_res_lim3(int e, double r);/*���㹹��������Դ���м����T_middle,Ϊ�������յĹ���������Դ������*/
double Cal_res_lim4(double left, double right, int module, double r);/*�󹹼�������Դ���޻�����*/
double fun4(double x, int module, double r);
double fun3(double x);

double cal_z(int  g, double f);//��ȥ��g�������zֵ
double cal_t1(int i);//ÿ��ģ���Ӧc��Сֵ����Ĳ�����Դ
double cal_t2(int i);//ÿ��ģ��ɿ��ԴﵽMAX_R����Ҫ�Ĳ�����Դ
double cal_sumcost(int  g, double f);//���ȥgģ���n-1��ģ�������ѵĳɱ���Сֵ
double cal_lowtime(int g);//�ɱ��Ϳɿ���Լ����ÿ��ģ�������Դ����
double cal_hightime(int g);//�ɱ��Ϳɿ���Լ����ÿ��ģ�������Դ����
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
	double sum_u = 0.0;/*�������Ԥ�ڷ��ʴ���֮��*/
	//������Ԥ�ڷ��ʴ���
	for (int i = 0; i < maxvar; i++)
	{
		sum_u += u[i];
	}
	//���㵥���������ʱ��
	for (int i = 0; i < maxvar; i++)
	{
		Mis_tim[i] = MISSION_TIME * u[i] / sum_u;
	}
}


double cal_z(int  g, double f)//��ȥ��g�������zֵ
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

double cal_t1(int i)//ÿ��ģ���Ӧc��Сֵ�Ĳ�����Դ
{
	double t;
	t = -1 / para_module[i][1] * log(para_module[i][4] / (para_module[i][0] * para_module[i][1] * (para_module[i][3] - para_module[i][2])));

	return t;
}

double cal_t2(int i)//ÿ��ģ��ɿ��ԴﵽMAX_R����Ҫ�Ĳ�����Դ
{
	double a = -log(Max_R) / (para_module[i][0] * para_module[i][1] * Mis_tim[i]);
	double t = -1 / para_module[i][1] * log(a);
	return t;
}
double cal_sumR(int  g, double f)//��n-1��ģ���ܿɿ���
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

double cal_sumcost(int  g, double f)//��n-1��ģ�������ѵĳɱ���Сֵ
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

double cal_lowtime(int g)//�ɱ��Ϳɿ���Լ����ÿ��ģ�������Դ����
{
	double t = cal_t1(g);
	double t1 = cal_t2(g);
	double min, max, mid;
	if (t < t1)//��ģ��ɱ���Сֵ��Ӧ�Ĳ�����ԴС��ģ��ɿ��ԴﵽMAX_R��Ӧ�Ĳ�����Դʱ,ģ��Ĳ�����Դ������(t1,tmid[g])
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

	if (t > t1)//��ģ��ɱ���Сֵ��Ӧ�Ĳ�����Դ����ģ��ɿ��ԴﵽMAX_R��Ӧ�Ĳ�����Դʱ
	{
		
			min = t1;
			max = t;
			mid = (min + max) / 2;
			while (fabs(C - cal_cost(g, mid) + cal_sumcost(g, mid)) > EPSINON && max - min > EPSINON)
			{
				if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
				if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

			}
		
		if (cal_cost(g, t) + cal_sumcost(g, t) > C)//�ж�t��Ӧ�ĳɱ�������ģ�������ͳɱ��Ƿ����C��������������޷�Χ��(t,T)
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
double cal_hightime(int g)//�ɱ��Ϳɿ���Լ����ÿ��ģ�������Դ����
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

		lim_r[k][0] = cal_lowtime(k);//lim_t[k][0] �ǳɱ��Ϳɿ���Լ��������
		if (lim_r[k][0] < 0) lim_r[k][0] = 0;
		lim_r[k][1] = cal_hightime(k);//lim_t[k][1] �ǳɱ��Ϳɿ���Լ��������
		
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
	sum1 = sum1 * sum2;/*ϵͳ�ɿ���*/
	obj[0] = 1.0 - sum1;/*here we use 1.0-sum1 to minimize obj[0]*///reliability
	/********************************************************************/
	double sum3 = 0;/*ϵͳ��ʼ�ɱ�Ϊ50.0*/
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
	
		printf("error���޸�Խ�磡����");
		system("pause");
	}
	
}
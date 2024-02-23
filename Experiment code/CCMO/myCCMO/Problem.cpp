/* Test problem definitions */
#include "stdafx.h"
#include<iostream>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include <string.h>
;
# include "Global.h"
# include "Random.h"
using namespace std;
#define maxvar nreal
#define nvar nreal
double C;
double dd = 0.0000001;

extern double u[maxvar];
double Mis_tim[maxvar];/*�������䵽������ʱ��*/
double lim_r[maxvar][2];/*Limits of variable in array*/
double lim_t[maxvar][2];
double T_low_sys;
double zz, zz2, zz1;
int q1, q2;/*ϵͳ������Դ��������ʱ��������������Դ����0�ĸ���*/
int Z_0[maxvar];/*��¼*/
double Z2[maxvar];/*ϵͳ������Դ��������ʱ��������������Դ*/
double e[maxvar];
double Zmax[maxvar];
double reciprocal_sum = 0.0;/*���ϼ���ʵ���֮��*/

double Y[maxvar];
double tmid[maxvar];//������Դ�м�ֵ��
double c_min[maxvar];//ÿ��ģ��ɱ��ļ�Сֵ��
double t_min[maxvar];//ÿ��ģ��ɱ��ļ�Сֵ���Ӧ�Ĳ�����Դ
void Cal_x();


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
	/*cost = para_module[ind_module][2] * para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res)) +
		para_module[ind_module][3] * (para_module[ind_module][0] - para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res))) +
		para_module[ind_module][4] * pow(m_res, para_module[ind_module][5]);*/
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
	while (fabs(sum - f2) > dd && (max - min) > dd)
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
		while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > dd && max - min > dd)
		{
			if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
			if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

		}
	}

	if (t > t1)//��ģ��ɱ���Сֵ��Ӧ�Ĳ�����Դ����ģ��ɿ��ԴﵽMAX_R��Ӧ�Ĳ�����Դʱ
	{
		if (cal_cost(g, t) + cal_sumcost(g, t) < C)//�ж�t��Ӧ�ĳɱ�������ģ�������ͳɱ��Ƿ�С��C�����С�������޷�Χ�ڣ�0��t)
		{
			min = t1;
			max = t;
			mid = (min + max) / 2;
			while (fabs(C - cal_cost(g, mid) + cal_sumcost(g, mid)) > dd && max - min > dd)
			{
				if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
				if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

			}
		}
		if (cal_cost(g, t) + cal_sumcost(g, t) > C)//�ж�t��Ӧ�ĳɱ�������ģ�������ͳɱ��Ƿ����C��������������޷�Χ��(t,T)
		{
			min = t;
			max = tmid[g];
			mid = (min + max) / 2;
			while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > dd && max - min > dd)
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
	double sum2 = 1.0;

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
	constr[0] = sum1 - Max_R;//reliability constraint
	constr[1] = (T_MAX_RESOURCE - sum4) / T_MAX_RESOURCE; // time constraint
	constr[2] = (C - sum3) / C;
	/********************************************************************/

	return;
}

void test2_problem(double* xreal, double* obj, double* constr)
{

	int j, k, p;
	/********************************************************************/
	j = 0;
	double sum1 = 1.0;
	double sum2 = 1.0;

	for (k = 0; k < nreal; k++)
	{
		sum2 = sum2 * cal_reliab(k, xreal[k]);

	}
	sum1 = sum1 * sum2;
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

	constr[0] = 0;//reliability constraint
	constr[1] = 0; // time constraint
	constr[2] = 0;
	return;
}


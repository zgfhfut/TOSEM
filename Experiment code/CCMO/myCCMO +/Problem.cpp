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
double Mis_tim[maxvar];/*构件分配到的任务时间*/
double lim_r[maxvar][2];/*Limits of variable in array*/
double lim_t[maxvar][2];
double T_low_sys;
double zz, zz2, zz1;
int q1, q2;/*系统测试资源总量下限时单个构件测试资源大于0的个数*/
int Z_0[maxvar];/*记录*/
double Z2[maxvar];/*系统测试资源总量下限时单个构件测试资源*/
double e[maxvar];
double Zmax[maxvar];
double reciprocal_sum = 0.0;/*故障检测率倒数之和*/

double Y[maxvar];
double tmid[maxvar];//测试资源中间值；
double c_min[maxvar];//每个模块成本的极小值点
double t_min[maxvar];//每个模块成本的极小值点对应的测试资源
void Cal_x();
void Cal_res_lim(double r);
double fun4(double x, int module, double r);
double Cal_res_lim4(double left, double right, int module, double r);
double fun3(double x);
void Cal_res_lim3(int e, double r);

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
	/*cost = para_module[ind_module][2] * para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res)) +
		para_module[ind_module][3] * (para_module[ind_module][0] - para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res))) +
		para_module[ind_module][4] * pow(m_res, para_module[ind_module][5]);*/
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

void Cal_res_lim(double r)
{
	double Tmin[maxvar];

	//double Tmin_sum = 0;/*组件测试资源下限之和*/
	double sum_u = 0.0;/*所有组件预期访问次数之和*/
	double left1[maxvar], right1[maxvar];
	reciprocal_sum = 0.0;
	for (int i = 0; i < maxvar; i++) {
		if (para_module[i][0] > EPSION)reciprocal_sum += 1.0 / para_module[i][1];
	}

	//计算总预期访问次数
	Cal_x();

	T_low_sys = 0.0; zz2 = 0.0; q1 = 0;
	for (int i = 0; i < maxvar; i++) {
		Z_0[i] = 1;
	}
	//计算系统测试资源下限
	Cal_res_lim3(nvar, r);
	for (int i = 0; i < maxvar; i++)
	{//因Cal_res_lim2()计算中已假设z2[i]计算中的负数为0，此时对其赋值0
		if (Z2[i] > 0.0)T_low_sys = T_low_sys + Z2[i];/*系统测试资源下限*/
	}
	if (T_low_sys > T_MAX_RESOURCE) {
		printf("\n%f测试资源总量不够", T_low_sys);
		//system("pause");
	}
	//计算单个组件达到可靠性约束时的测试资源
	for (int i = 0; i < maxvar; i++)
	{
		//double rrr = exp(-para_module[i][0] * 0.04 * Mis_tim[i]*exp(-0.04*400));
		if (para_module[i][0] < EPSION)Tmin[i] = 0.0;
		else Tmin[i] = -1 / para_module[i][1] * log(-log(r) / (para_module[i][1] * para_module[i][0] * Mis_tim[i]));

		/*if (Tmin[i] < 0.0)
		{
			Tmin[i] = 0.0;
		}*/
		/*Tmax[i] = 0.0; zz2 = 0.0; q1 = 0;
		for (int j = 0; j < maxvar; j++) {
		if (i == j)
		Z_0[j] = 0;
		else
		Z_0[j] = 1;
		}
		Cal_res_lim3(i);
		for (int j = 0; j < maxvar; j++)
		{
		if (Z2[j] < 0.0)
		Z2[j] = 0.0;
		Tmax[i] = Tmax[i] + Z2[j];
		}*/
	}
	//M = -log(MAX_R);

	//zz1 = 0;	q1 = 0;
	//for (int i = 0; i < maxvar; i++) {
	//	Z_0[i] = 1;
	//}
	//double www = 0.0;/*用来查看得到最大可靠性的分配之和*/
	//Cal_res_lim2();/*计算组件测试资源上限*/
	//for (int i = 0; i < maxvar; i++)
	//{
	//	if (Zmax[i] < 0.0)
	//		Zmax[i] = 0.0;
	//	www += Zmax[i];
	//}
	//求组件下限
	double time0 = 0.0;
	for (int i = 0; i < maxvar; i++) {
		zz1 = 0;	q1 = 1;
		left1[i] = Tmin[i] + EPSION;
		right1[i] = Z2[i];
		for (int j = 0; j < maxvar; j++) {
			if (i == j)
				Z_0[j] = 0;
			else
				Z_0[j] = 1;
		}
		if (para_module[i][0] < EPSION || Z2[i] < 0.0)lim_r[i][0] = 0.0;
		else lim_r[i][0] = Cal_res_lim4(left1[i], right1[i], i, r);
		if (lim_r[i][0] < 0.0)lim_r[i][0] = 0.0;
		time0 += lim_r[i][0];
	}
	/*for (int i = 0; i < maxvar; i++) {
		re_parameter(lim_r[i][0], time0,i);
	}*/
	//求组件上限
	for (int i = 0; i < maxvar; i++) {
		zz1 = 0;	q1 = 1;
		if (Z2[i] < 0.0)Z2[i] = 0.0;
		left1[i] = Z2[i];
		right1[i] = T_MAX_RESOURCE + EPSION;
		for (int j = 0; j < maxvar; j++) {
			if (i == j)
				Z_0[j] = 0;
			else
				Z_0[j] = 1;
		}
		if (para_module[i][0] < EPSION)lim_r[i][1] = 0.0;
		else lim_r[i][1] = Cal_res_lim4(left1[i], right1[i], i, r);
	}

	double sum3 = 0.0;
	for (int i = 0; i < maxvar; i++) {
		sum3 += lim_r[i][0];
	}
	if (sum3 > T_MAX_RESOURCE + EPSION)printf("ERROR:下限之和越界！！！");

	double sum = 0.0;
	for (int i = 0; i < maxvar; i++) {
		sum += lim_r[i][0];
		if (lim_r[i][1] > T_MAX_RESOURCE + EPSION || lim_r[i][0] < 0.0 - EPSION) {
			printf("ERROR:Cal_res_lim4中出错！！！");
		}
	}

}

double fun4(double x, int module, double r)
{
	int i;
	double aa = T_MAX_RESOURCE - x;
	//小于0的测试资源设为0，即后续不考虑测试资源小于0的组件
	for (i = 0; i < maxvar; i++) {
		if (i == module || para_module[i][0] < EPSION)aa -= 0.0;

		if (Z_0[i] == 1) {
			aa -= 1.0 / para_module[i][1] * log(-(para_module[i][1] * para_module[i][1] * para_module[i][0] * Mis_tim[i] * (reciprocal_sum - 1.0 / para_module[module][1]))
				/ (log(r) + para_module[module][1] * para_module[module][0] * Mis_tim[module] * exp(-para_module[module][1] * x)));
		}
	}
	//double www = para_module[module][1] * para_module[module][0] * Mis_tim[module] * exp(-para_module[module][1] * x);
	return aa;
}
double Cal_res_lim4(double left, double right, int module, double r) {
	double dd = 1.E-10;/*精确度*/
	q2 = q1; q1 = 0;
	double result_value;
	double test_R[maxvar];
	double _left = left, _right = right;
	//double 类型的输入输出格式是lf，判断条件是_left,_right处的函数值异号
	/*if ((fun4(_left, module, r, nop) * fun4(_right, module, r, nop) <= 0.0) && (_left < _right))
	{*/
	while (_right - _left >= EPSION)
		//判断解在哪个区间是应该是区间端点的函数值异号，而不是中点处函数值的大小
	{
		if (fun4(_left, module, r) * fun4((_left + _right) / 2.0, module, r) > 0)
			//if判断后执行的不是一条语句就要用大括号括起来
		{
			_left = (_left + _right) / 2.0;
			continue;
		}
		if (fun4((_left + _right) / 2.0, module, r) * fun4(_right, module, r) > 0)
		{
			_right = (_left + _right) / 2.0;
			continue;
		}
		if (fabs(_left - _right) <= EPSION) {
			break;
		}
	}
	//当某一个端点的值足够，退出循环

	/*if (fun4(_left, module, r, nop) * fun4(_right, module, r, nop) > 0.0)
	{
		printf("\n不在二分法区间内");
	}*/
	/*printf("%lf\t%lf", _left, _right);*/
	if (_left < _right) {
		result_value = _right;
	}
	else
		result_value = _left;
	/*}*/
	//else printf("不在区间内！！！！！！！！！");
	//else result_value= _right;
	//zz = log(MAX_R);
	double etete = T_low_sys;
	for (int i = 0; i < maxvar; i++)
	{
		/*for (int j = 0; j < maxvar; j++)
		{*/
		if (i == module)test_R[i] = 0.0;
		else {
			test_R[i] = 1.0 / para_module[i][1] * log(-(para_module[i][1] * para_module[i][1] * para_module[i][0] * Mis_tim[i] * (reciprocal_sum - 1.0 / para_module[module][1]))
				/ (log(r) + para_module[module][1] * para_module[module][0] * Mis_tim[module] * exp(-para_module[module][1] * result_value)));
		}

		//如果得到的组件测试资源小于0，则q1加1，随后继续递归
		if (test_R[i] <= 0) {
			Z_0[i] = 0;
			q1 += 1;

		}
		else
			Z_0[i] = 1;

	}

	zz1 += 1;/*记录循环次数*/
	//zz2 = zz;

	//如果组件测试资源下限值小于0的个数增加，则继续递归
	if (q1 > q2) {
		return  Cal_res_lim4(left, right, module, r);
	}

	return result_value;
}
double fun3(double x)
{
	int i;
	double aa = T_MAX_RESOURCE;
	for (i = 0; i < maxvar; i++) {
		if (Z_0[i] == 1)
		{
			aa += 1 / para_module[i][1] * log(x / (para_module[i][1] * para_module[i][1] * para_module[i][0] * Mis_tim[i]));
		}
	}
	return aa;
}

void Cal_res_lim3(int e, double r) {

	q2 = q1; q1 = 0;

	double ff = 0.0;
	zz = log(r);
	for (int i = 0; i < nvar; i++) {
		if (Z_0[i] == 1)
			ff -= 1.0 / para_module[i][1];
	}
	ff = ff * 1.0 / (zz - zz2);
	zz2 = 0.0;
	for (int i = 0; i < nvar; i++)
	{
		if (i == e)
			Z2[i] = 0.0;
		else
			Z2[i] = 1.0 / para_module[i][1] * log(para_module[i][1] * para_module[i][1] * para_module[i][0] * Mis_tim[i] * ff);

		if (Z2[i] < 0.0) {
			Z_0[i] = 0;
			q1 += 1;
			zz2 += -Mis_tim[i] * para_module[i][1] * para_module[i][0];
		}
	}
	//如果构件测试资源下限值小于0的个数增加，则返回
	if (q1 > q2) {
		return  Cal_res_lim3(e, r);
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
		while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > dd && max - min > dd)
		{
			if (cal_cost(g, mid) + cal_sumcost(g, mid) < C)  max = mid; mid = (max + min) / 2;
			if (cal_cost(g, mid) + cal_sumcost(g, mid) > C)  min = mid; mid = (max + min) / 2;

		}
	}

	if (t > t1)//当模块成本极小值对应的测试资源大于模块可靠性达到MAX_R对应的测试资源时
	{
		if (cal_cost(g, t) + cal_sumcost(g, t) < C)//判断t对应的成本加其他模块的总最低成本是否小于C，如果小于则下限范围在（0，t)
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
		if (cal_cost(g, t) + cal_sumcost(g, t) > C)//判断t对应的成本加其他模块的总最低成本是否大于C，如果大于则下限范围在(t,T)
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
	//Cal_res_lim(Max_R);
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
	double sum2 = 1.0;

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
	double sum3 = 0;/*系统初始成本为50.0*/
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

double repair_s(double low,double high,double l,double r )
{  double rand = randomperc();
  double res;
  if (r < low)
	  res = rand * (l- low) + low;
  if (r> high)
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
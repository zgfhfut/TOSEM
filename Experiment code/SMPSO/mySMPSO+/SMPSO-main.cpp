




// SMPSO-main.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "gobel.h"
#include "Problem.h"
#include "Algorithm.h"
#include "Solution.h"
#include "Operator.h"
#include "ProblemFactory.h"
#include "SMPSO.h"
#include "PolynomialMutation.h"
#include <iostream>
#include <time.h>

using namespace std;

/**
* This main executes the SMPSO algorithm described in:
* A.J. Nebro, J.J. Durillo, J. Garcia-Nieto, C.A. Coello Coello, F. Luna and E. Alba
* "SMPSO: A New PSO-based Metaheuristic for Multi-objective Optimization".
* IEEE Symposium on Computational Intelligence in Multicriteria Decision-Making
* (MCDM 2009), pp: 66-73. March 2009
*/

double para_module[nreal][5];
double u[nreal];
double Mis_tim[nreal];
double C;
double Y[nreal];
double tmid[nreal];//测试资源中间值；
double lim_r[nreal][2];

void  Cal_res_lim() {

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

double cal_z(int  g, double f)
{
	double f1 = cal_reliab(g, f);
	double max = DBL_MAX;
	double min = 0.0;
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;
	double f2 = log(Max_R) - log(f1);
	while (fabs(sum - f2) > EPSION && (max - min) > EPSION)
	{

		sum = 0.0;
		for (int i = 0; i < nreal; ++i)
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
	for (int k = 0; k < nreal; ++k)
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
	for (int k = 0; k < nreal; ++k)
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
		while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > EPSION && max - min > EPSION)
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
			while (fabs(C - cal_cost(g, mid) + cal_sumcost(g, mid)) > EPSION && max - min > EPSION)
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
			while (fabs(cal_cost(g, mid) + cal_sumcost(g, mid) - C) > EPSION && max - min > EPSION)
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
void Cal_C()
{
	double C_min = 0.0;
	double max = DBL_MAX;
	double min = 0.0;
	double c[nreal];
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;


	while (fabs(sum - log(Max_R)) > EPSION && (max - min) > EPSION)
	{

		sum = 0.0;
		for (int i = 0; i < nreal; ++i)
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
	for (int j = 0; j < nreal; ++j)
	{
		Y[j] = para_module[j][4] * Mis_tim[j] / (para_module[j][2] - para_module[j][3] - mid * para_module[j][1] * Mis_tim[j]);
		tmid[j] = -1 / para_module[j][1] * log(-Y[j] / (para_module[j][0] * para_module[j][1] * Mis_tim[j]));
		c[j] = (para_module[j][2] - para_module[j][3]) * para_module[j][0] * (1 - exp(-para_module[j][1] * tmid[j])) + para_module[j][3] * para_module[j][0] + para_module[j][4] * tmid[j];
		C_min += c[j];

	}
	C = C_min * 1.20;
}
int _tmain(int argc, _TCHAR* argv[])
{
	clock_t t_ini, t_fin;

	
	Problem   * problem;   // The problem to solve
	Algorithm * algorithm; // The algorithm to use
	Operator  * mutation;  // "Turbulence" operator

	map<string, void *> parameters; // Operator parameters

	//TODO: QualityIndicator * indicators; // Object to get quality indicators
	for (int o = 1; o < 31; o++)
	{
		FILE* end_ptr2; FILE* end_ptr1;
		char file3[500]; char file2[500];

		int i;
		sprintf_s(file3, "D:\\论文代码数据\\Parameters\\single-input and single-output\\SMPSO\\组件参数\\50\\%d.txt",o);

		fopen_s(&end_ptr2, file3, "r+");
		//在file3中写上每次迭代产生的不同fitness文件的路径

		for (i = 0; i < nreal; i++)
		{

			fscanf_s(end_ptr2, "%lf %lf %lf %lf %lf", &para_module[i][0], &para_module[i][1], &para_module[i][2], &para_module[i][3], &para_module[i][4]);// 循环读

		}

		fclose(end_ptr2); //关闭文件
		sprintf_s(file2, "D:\\论文代码数据\\Parameters\\single-input and single-output\\SMPSO\\预期访问次数\\50\\%d.txt",o);

		fopen_s(&end_ptr1, file2, "r+");
		for (i = 0; i < nreal; i++)
			fscanf_s(end_ptr1, "%lf", &u[i]); // 循环读
		fclose(end_ptr1); //关闭文件
		for (int int_exp = 1; int_exp <= 30; int_exp++)
		{
			printf_s("Run experiment %d!\n", int_exp);
			Cal_res_lim();
			Cal_C();
			for (int i = 0; i < nreal; ++i)
			{
				lim_r[i][0] = cal_lowtime(i);
				lim_r[i][1] = cal_hightime(i);
			}
			problem = ProblemFactory::getProblem(const_cast<char*>("Tanaka"));


			algorithm = new SMPSO(problem);

			// Algorithm parameters
			int swarmSizeValue = 300;//主种群规模
			int archiveSizeValue = 300;//扩展种群规模
			int maxIterationsValue = 200;//最大迭代次数
			algorithm->setInputParameter("swarmSize", &swarmSizeValue);
			algorithm->setInputParameter("archiveSize", &archiveSizeValue);
			algorithm->setInputParameter("maxIterations", &maxIterationsValue);

			// Mutation operator
			double probabilityParameter = 1.0 / (problem->getNumberOfVariables());//变异概率
			double distributionIndexParameter = 20.0;
			parameters["probability"] = &probabilityParameter;
			parameters["distributionIndex"] = &distributionIndexParameter;
			mutation = new PolynomialMutation(parameters);

			// Add the operators to the algorithm
			algorithm->addOperator("mutation", mutation);

			// Add the indicator object to the algorithm
			//algorithm->setInputParameter("indicators", indicators) ;

			// Execute the Algorithm
			t_ini = clock();
			SolutionSet* population = algorithm->execute();
			t_fin = clock();
			double secs = (double)(t_fin - t_ini);
			secs = secs / CLOCKS_PER_SEC;

			// Result messages
			cout << "Total execution time: " << secs << "s" << endl;
			//cout << "Variables values have been written to file VAR" << endl;
			//population->printVariablesToFile("VAR");
			population->printObjectivesToFile("D:\\论文代码数据\\experimental data\\single-input and single-output\\mySMPSO+\\50\\Data of 30 instances\\test" + to_string(o) + "\\fitness" + to_string(int_exp) + ".txt"); 

			delete mutation;
			delete population;
			delete algorithm;
		}
		printf("\n Routine successfully exited \a\n");
	}

	return 0;
}


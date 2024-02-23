




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
	for (int o =1; o < 31; o++)
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
				lim_r[i][0] =0;
				lim_r[i][1] =T_MAX_RESOURCE;
		;
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
			population->printObjectivesToFile("D:\\论文代码数据\\experimental data\\single-input and single-output\\SMPSO\\50\\Data of 30 instances\\test" + to_string(o) + "\\fitness" + to_string(int_exp) + ".txt"); 

			delete mutation;
			delete population;
			delete algorithm;
		}
		printf("\n Routine successfully exited \a\n");
	}

	return 0;
}


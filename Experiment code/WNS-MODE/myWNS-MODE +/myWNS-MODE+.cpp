// WNS-MODE.cpp : 定义控制台应用程序的入口点。
//


#include "stdafx.h"
#include<iostream>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#include "string.h"
/********************************************************************/
#include "random.h"       /*Random Number Generator*/
/********************************************************************/
using namespace std;

/********************************************************************/
#define maxpop   300  /*Max population */

//#define MISSION_TIME  100.0
#define MISSION_TIME  200.0  //complex system
//#define MISSION_TIME  200.0  //large system

/*Max value of testing time*/
//#define T_MAX_RESOURCE  8000.0 
#define T_MAX_RESOURCE  150000.0  //complex system
//#define T_MAX_RESOURCE  70000.0  //large system

/*Max no. of variables*/

//#define maxvar    10
#define maxvar    100 //complex system
//#define maxvar    50 ////large system

#define MAX_R  0.9//the reliability contraint

//#define Cost_4  250000.0//the potential cost caused by software failure for large system

#define EPSION  0.000001
#define maxfun    3  /*Max no. of functions */
#define maxcons   3  /*Max no. of Constraints*/
/********************************************************************/
/********************************************************************/
#define MAX_EXPERIMENT  30  /*Experimental times*/
#define MAX_EVALUATION  90000/*Max no. of evaluation*/

/*****************************************a_i,b_i,c_1^i,c_2^i,c_3^i,\gama_i***************************/
double para_module[maxvar][5];/*parameters for component*/
double  C;
double u[maxvar] ; //Expected number of visits

double Mis_tim[maxvar]; //Expected number of visits
/********************************************************************/


double lim_r[maxvar][2];/*Limits of variable in array*/
double lim_t[maxvar][2];
double reciprocal_sum = 0.0;
double Y[maxvar];
double tmid[maxvar];//测试资源中间值；
double c_min[maxvar];//每个模块成本的极小值点
double t_min[maxvar];//每个模块成本的极小值点对应的测试资源

double T_low_sys;
double zz, zz2, zz1;
int q1, q2;/*系统测试资源总量下限时单个组件测试资源大于0的个数*/
int Z_0[maxvar];/*记录*/
double Z2[maxvar];/*系统测试资源总量下限时单个组件测试资源*/


# define INF 1.0e14
# define EPS 1.0e-14
# define PI 3.14159265358979

/********************************************************************/

#define BLOCK_INTERVAL 90  //
const double weight[maxfun] = { 0.1, 0.4, 0.5 };
/********************************************************************/

/********************************************************************/

/********************************************************************/
int gener;       /*No of generations*/
int nvar;          /*No of variables*/
int ncons;         /*No of Constraints*/

int eva_num = 0;
/********************************************************************/

/********************************************************************/
double seed;      /*Random Seed*/
double pcross;        /*Cross-over Probability*/
double pmut_r;          /*Mutation Probability*/

double di;       /*Distribution Index for the Cross-over*/
double dim;          /*Distribution Index for the Mutation*/


int nfunc;       /*No of functions*/
static int popsize;/*Population Size*/
/********************************************************************/

/********************************************************************/
typedef struct       /*individual properties*/
{
	double xreal[maxvar]; /*list of real variables*/
	double fitness[maxfun];/*Fitness values */
	double constr[maxcons];     /*Constraints values*/
	double error;              /* overall constraint violation for the individual*/

}individual;        /*Structure defining individual*/

typedef struct
{
	individual ind[maxpop]; /*Different Individuals*/
	individual* ind_ptr;
}population;             /*Population Structure*/


/********************************************************************/

/********************************************************************/

/********************************************************************/
typedef struct
{
	int maxrank;   /*Max rank of the global population*/
	int rank[2 * maxpop],            /*rank of different individuals*/
		flag[2 * maxpop];            /*Setting the flag */

	double fitness[2 * maxpop][maxfun], /*Fitness function values for the different individuals*/
		xreal[2 * maxpop][maxvar],      /*value of the decoded variables for different individuals */
		error[2 * maxpop],               /*Error Values of the individuals*/
		constr[2 * maxpop][maxcons];

	double norm_fitness[2 * maxpop][maxfun];
	double weight_fitness[2 * maxpop];

}globpop;/*Population structure for the pool having both the old as well as new population*/


globpop arch_globalpop, * arch_global_pop_ptr, arch_globaltemp;
/********************************************************************/

/********************************************************************/
void input(); /*Input Parameters from user*/

void Cal_x();
void realinit(population* pop_ptr);/*initializes the population*/
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
double repair_s(double low, double high, double l, double r);
void  repair_all(individual* ind, double lim[], int count[]);

double cal_reliab(int ind_module, double m_res);//calculate the reliability of each module
double cal_cost(int ind_module, double m_res);//calculate the cost of each module
void func(population* pop_ptr);/*evaluate the value of the function & errors*/

void realcross(population* new_pop_ptr, population* mate_pop_ptr);/*simulated binary crossover for Real Coded GA*/
void real_mutate(population* old_pop_ptr, population* new_pop_ptr);/*formulate the mutation routine*/


void addInArchive(population* pop_ptr, population* arch2pop_ptr);

void DivideBlock(population* pop_ptr1, population* pop_ptr2, population* pop_ptr3);
/********************************************************************/
void sleep(clock_t wait)
{
	clock_t goal;
	goal = wait + clock();
	while (goal > clock());
}

int main(int argc, _TCHAR* argv[])
{
	int pp = 0;
	for (pp = 1; pp < 31; pp++)
	{
		printf_s("test Instance %d!\n", pp);

		FILE* end_ptr1;
		char file1[500];
		sprintf_s(file1, "C:\\Users\\Administrator\\Desktop\\论文代码数据\\Parameters\\multi-input and multi-output\\WNS-MOEAD\\组件参数\\100\\%d.txt", pp);
		fopen_s(&end_ptr1, file1, "r+");
		for (int i = 0; i < maxvar; i++)
			fscanf_s(end_ptr1, "%lf %lf %lf %lf %lf", &para_module[i][0], &para_module[i][1], &para_module[i][2], &para_module[i][3], &para_module[i][4]); // 循环读
		fclose(end_ptr1); //关闭文件


		FILE* end_ptr2;
		char file2[500];
		sprintf_s(file2, "C:\\Users\\Administrator\\Desktop\\论文代码数据\\Parameters\\multi-input and multi-output\\WNS-MOEAD\\预期访问次数\\100\\%d.txt", pp);
		fopen_s(&end_ptr2, file2, "r+");
		for (int i = 0; i < maxvar; i++)
			fscanf_s(end_ptr2, "%lf", &u[i]); // 循环读
		fclose(end_ptr2); //关闭文件

		int int_exp = 0;
	
		input();
		Cal_res_t();
		for (int_exp = 0; int_exp < MAX_EXPERIMENT; int_exp++)
		{
			
			printf_s("Run experiment %d!\n", int_exp + 1);

			eva_num = 0;


			FILE* end_ptr;
			char file4[500];
			
			sprintf_s(file4, "C:\\Users\\Administrator\\Desktop\\论文代码数据\\experimental data\\multi-input and multi-output\\myWNS-MOEAD+\\100\\Data of 30 instances\\test%d\\fitness%d.txt", pp, int_exp + 1);//在file1中写上每次迭代产生的不同fitness文件的路径
			
			fopen_s(&end_ptr, file4, "w+");


			int j, l, f;

			population oldpop, newpop, matepop, newpop2, * old_pop_ptr, * new_pop_ptr, * mate_pop_ptr, * new_pop_ptr2;/*Defining the population Structures*/

			population archivepop, * archive_pop_ptr;

			int random;

			do
			{
				srand((unsigned)time(NULL));
				random = rand() % 1000;

				seed = (float)random / 1000.0;

			} while (seed <= 0.0 || seed >= 1.0);


			warmup_random(seed);

			

			/***********************************************/

			arch_globalpop = arch_globaltemp;

			old_pop_ptr = NULL;
			new_pop_ptr = NULL;
			new_pop_ptr2 = NULL;

			mate_pop_ptr = NULL;
			archive_pop_ptr = NULL;

			old_pop_ptr = &(oldpop);
			new_pop_ptr = &(newpop);
			new_pop_ptr2 = &(newpop2);
			mate_pop_ptr = &(matepop);
			archive_pop_ptr = &(archivepop);

			for (j = 0; j < popsize; j++)
			{
				old_pop_ptr->ind[j].error = 0;
				new_pop_ptr->ind[j].error = 0;
				new_pop_ptr2->ind[j].error = 0;
				mate_pop_ptr->ind[j].error = 0;
				archive_pop_ptr->ind[j].error = 0;

				old_pop_ptr->ind_ptr = NULL;
				new_pop_ptr->ind_ptr = NULL;
				new_pop_ptr2->ind_ptr = NULL;
				mate_pop_ptr->ind_ptr = NULL;
				archive_pop_ptr->ind_ptr = NULL;

				for (l = 0; l < ncons; l++)
				{
					old_pop_ptr->ind[j].constr[l] = 0;
					new_pop_ptr->ind[j].constr[l] = 0;
					new_pop_ptr2->ind[j].constr[l] = 0;
					mate_pop_ptr->ind[j].constr[l] = 0;
					archive_pop_ptr->ind[j].constr[l] = 0;
				}

				for (l = 0; l < nvar; l++)
				{

					old_pop_ptr->ind[j].xreal[l] = 0;
					new_pop_ptr->ind[j].xreal[l] = 0;
					new_pop_ptr2->ind[j].xreal[l] = 0;
					mate_pop_ptr->ind[j].xreal[l] = 0;
					archive_pop_ptr->ind[j].xreal[l] = 0;
				}

				for (f = 0; f < nfunc; f++)
				{
					old_pop_ptr->ind[j].fitness[f] = 0;
					new_pop_ptr->ind[j].fitness[f] = 0;
					new_pop_ptr2->ind[j].fitness[f] = 0;
					mate_pop_ptr->ind[j].fitness[f] = 0;
					archive_pop_ptr->ind[j].fitness[f] = 0;
				}

			}
			/*******************************************************/

			old_pop_ptr = &(oldpop);
			realinit(old_pop_ptr);

			old_pop_ptr = &(oldpop);
			
			func(old_pop_ptr);/*Function Calculation*/

			archive_pop_ptr = &(archivepop);
			addInArchive(old_pop_ptr, archive_pop_ptr);

			/********************************************************************/
			/*----------------------GENERATION STARTS HERE----------------------*/
			while (eva_num < MAX_EVALUATION)
			{
				old_pop_ptr = &(oldpop);
				mate_pop_ptr = &(matepop);
				real_mutate(old_pop_ptr, mate_pop_ptr); /*Real Mutation*/
				/********************************************************************/
				mate_pop_ptr = &(matepop);
				

				old_pop_ptr = &(oldpop);
				realcross(old_pop_ptr, mate_pop_ptr);/*Real Cross-over*/
				/********************************************************************/
				mate_pop_ptr = &(matepop);
			
				func(mate_pop_ptr);/*----------FUNCTION EVALUATION-----------*/
				/************************selection********************************************/

				old_pop_ptr = &(oldpop);
				mate_pop_ptr = &(matepop);
				new_pop_ptr = &(newpop);
				DivideBlock(old_pop_ptr, mate_pop_ptr, new_pop_ptr);
				/************************selection********************************************/

				/*-------------------Update the Archive--------------*/
				new_pop_ptr = &(newpop);
				archive_pop_ptr = &(archivepop);
				new_pop_ptr2 = &(newpop2);
				DivideBlock(new_pop_ptr, archive_pop_ptr, new_pop_ptr2);
				addInArchive(new_pop_ptr2, archive_pop_ptr);

				/********************************************************************/
				/*=======Copying the new population to old population======*/
				old_pop_ptr = &(oldpop);
				new_pop_ptr = &(newpop);
				for (j = 0; j < popsize; j++)
				{

					old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
					new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);

					old_pop_ptr->ind_ptr->error = new_pop_ptr->ind_ptr->error;

					if (nvar > 0)
					{
						/*For Real Coded GA copying of the chromosomes*/
						for (l = 0; l < nvar; l++)
							old_pop_ptr->ind_ptr->xreal[l] = new_pop_ptr->ind_ptr->xreal[l];
					}

					/*Copying the fitness vector */
					for (l = 0; l < nfunc; l++)
						old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];

					for (l = 0; l < ncons; l++)
						old_pop_ptr->ind_ptr->constr[l] = new_pop_ptr->ind_ptr->constr[l];

				}   // end of j			
				/********************************************************************/



			}// end of while

			/*Printing the fitness record for last generation in a file last*/

			archive_pop_ptr = &(archivepop);
			for (f = 0; f < popsize; f++) // for printing
			{
				archive_pop_ptr->ind_ptr = &(archive_pop_ptr->ind[f]);
				if (archive_pop_ptr->ind_ptr->error <= 0.0)  // for all feasible solutions and non-dominated solutions
				{
					for (l = 0; l < nfunc; l++)
						//cout << archive_pop_ptr->ind_ptr->fitness[l] << " ";
						fprintf(end_ptr, "%5.6f\t", archive_pop_ptr->ind_ptr->fitness[l]);
						//cout << endl;
					fprintf(end_ptr, "%s", "\n");
				}
			} // end of f (printing)


			fclose(end_ptr);
		}
	}
	//fclose(run_time);

	//	printf_s("please input a char to exit the program!\n");

	//	getchar();

	return 0;
}

void input()
{
	int i;
	/*number of the variables*/
	nvar = maxvar;

	/*number of the functions*/
	nfunc = maxfun;

	/*number of the individuals in the population*/
	popsize = maxpop;

	ncons = maxcons;

	/*the crossover probability (between 0.5 and 1)*/
	pcross = 0.1;

	/*the mutation probability*/
	pmut_r = 1.4;

	/************************************************************************************/
	

}

void realinit(population* pop_ptr)
{
	int i, j;
	double d,sum;
	double low_T = 0, high_T = 0;
	for (int q = 0; q < nvar; q++)
	{
		low_T += lim_r[q][0];
		high_T += lim_r[q][1];
	}
	for (i = 0; i < popsize; i++)
	{
		sum = 0;
		for (j = 0; j < nvar; j++)
		{
			d = randomperc();
			/*limits are specified it generates the value in range of minimum and maximum value of the variable*/
			pop_ptr->ind[i].xreal[j] = d * (lim_r[j][1] - lim_r[j][0]) + lim_r[j][0];
			sum += pop_ptr->ind[i].xreal[j];
		}

		if (sum > T_MAX_RESOURCE)
		{
			double new_sum = 0;/*计算修复后的系统测试资源总量*/
			d = randomperc();
			for (j = 0; j < nvar; j++)
			{
                pop_ptr->ind[i].xreal[j] = lim_r[j][0] + ((pop_ptr->ind[i].xreal[j] - lim_r[j][0]) * (T_MAX_RESOURCE - low_T) / (sum - low_T)) * d;
				new_sum += pop_ptr->ind[i].xreal[j];
			}
			if (new_sum > T_MAX_RESOURCE)
			{
				cout << "修复失败111" << endl;
			}
		}

	}
	

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
void Cal_C()
{

	double C_min = 0.0;
	double max = DBL_MAX;
	double min = 0.0;
	double c[maxvar];
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;
	while (fabs(sum - log(MAX_R)) > EPSION && (max - min) > EPSION)
	{

		sum = 0.0;
		for (int i = 0; i < maxvar; ++i)
		{

			sum += para_module[i][4] * Mis_tim[i] / (para_module[i][2] - para_module[i][3] - mid * para_module[i][1] * Mis_tim[i]);

		}

		double aa = log(MAX_R);
		if (sum < log(MAX_R))
		{
			min = mid;
			mid = (min + max) / 2;
		}
		if (sum > log(MAX_R))
		{
			max = mid;
			mid = (min + max) / 2;
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
double cal_z(int  g, double f)//除去第g个组件的z值
{
	double f1 = cal_reliab(g, f);
	double max = DBL_MAX;
	double min = 0.0;
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;
	double f2 = log(MAX_R) - log(f1);
	while (fabs(sum - f2) > EPSION && (max - min) > EPSION)
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
	double a = -log(MAX_R) / (para_module[i][0] * para_module[i][1] * Mis_tim[i]);
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
void Cal_res_t()
{
	Cal_x();
	double C_min = 0.0;
	double max = DBL_MAX;
	double min = 0.0;
	double c[maxvar];
	double mid = (max + min) / 2;
	double sum = 0, sum1 = 0;
	

	while (fabs(sum - log(MAX_R)) > EPSION && (max - min) > EPSION)
	{

		sum = 0.0;
		for (int i = 0; i < maxvar; ++i)
		{

			sum += para_module[i][4] * Mis_tim[i] / (para_module[i][2] - para_module[i][3] - mid * para_module[i][1] * Mis_tim[i]);

		}

		double aa = log(MAX_R);
		if (sum < log(MAX_R))
		{
			min = mid;
			mid = (min + max) / 2;
		}
		if (sum > log(MAX_R))
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
	for (int j = 0; j < nvar; j++)
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
	for (int j = 0; j < nvar; j++)
	{
		sum += ind->xreal[j];
	}
	if (sum > T_MAX_RESOURCE + EPSION) {
		cout << sum << endl;
		printf("error：修复越界！！！");
		system("pause");
	}

}
void func(population* pop_ptr)
{
	double* realx_ptr, /*Pointer to the array of x values*/
		* fitn_ptr,      /*Pointer to the array of fitness function*/
		x[2 * maxvar],     /* problem variables */
		f[maxfun],     /*array of fitness values*/
		* err_ptr,      /*Pointer to the error */
		cstr[maxcons];


	double error, cc;

	int i, j, k;


	double sum1 = 0, sum2 = 0;//added by Guofu Zhang

	double sum3 = 0;
	double sum4 = 0;

	pop_ptr->ind_ptr = &(pop_ptr->ind[0]);


	for (i = 0; i < popsize; i++)
	{
		eva_num++;

		pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
		realx_ptr = &(pop_ptr->ind_ptr->xreal[0]);


		for (j = 0; j < nvar; j++)
		{ // Real-coded variables 
			x[j] = *realx_ptr++;

			if (x[j] < 0.0)
			{
				printf_s("error：x[j]=%0.4f<0.0\n", x[j]);
				system("pause");
			}
		}

		fitn_ptr = &(pop_ptr->ind_ptr->fitness[0]);
		err_ptr = &(pop_ptr->ind_ptr->error);


		/*First fitness function */
		/* calculate the reliability of each subsystem*/
		sum1 = 1.0;
		sum2 = 0.0;
		for (j = 0; j < nvar; j++)
		{
			sum1 = sum1 * (cal_reliab(j, x[j]));
			sum2 += x[j];

		}

		f[0] = 1.0 - sum1;/*here we use 1-sum2 to minimize f[0]*/

		// Second and Third Fitness Functions
		sum3 = 0;
		for (j = 0; j < nvar; j++)
		{
			sum3 += cal_cost(j, x[j]);

		}


		f[1] = sum3;
		f[2] = sum2;

		if ((f[2] - T_MAX_RESOURCE) > EPSION)
		{
			//cout << f[2] << " " << T_MAX_RESOURCE << endl;
			printf_s("error in the upper bound constraint11\n");
			system("pause");
		}


		/*=========End Your Coding Upto This Point===============*/

		/******************************************************************/
		/*              Put The Constraints Here                          */
		/******************************************************************/
		cstr[0] = sum1 - MAX_R;
		cstr[1] = (T_MAX_RESOURCE - sum2) / T_MAX_RESOURCE;
		cstr[2] = (C - sum3) / C;
		/*===========Constraints Are Coded Upto Here=============*/
		/*   DO NOT CHANGE ANYTHING BELOW  */
		for (k = 0; k < nfunc; k++)
		{
			*fitn_ptr++ = f[k];
		}

		for (k = 0; k < ncons; k++)
		{
			pop_ptr->ind_ptr->constr[k] = cstr[k];
		}
		error = 0.0;
		for (k = 0; k < ncons; k++)
		{
			cc = cstr[k];
			if (cc < 0.0)
				error = error - cc;
		}
		*err_ptr = error;


	}

}

void realcross(population* new_pop_ptr, population* mate_pop_ptr)
{
	int i, j;
	double rnd2,old;

	for (i = 0; i < popsize; i++)
	{
		double lim[maxvar];
		int count[maxvar];
		for (j = 0; j < nvar; j++)
		{
			count[j] = 0;
			lim[j] = 0
				;
			rnd2 = randomperc();
			if (rnd2 < pcross)
			{
				mate_pop_ptr->ind[i].xreal[j] = mate_pop_ptr->ind[i].xreal[j];

			}
			else {
			
			    old = mate_pop_ptr->ind[i].xreal[j];
				mate_pop_ptr->ind[i].xreal[j] = new_pop_ptr->ind[i].xreal[j];
				if (mate_pop_ptr->ind[i].xreal[j] < old)
				{
					lim[j] = lim_r[j][0];
					count[j] = 1;
				}
				if (mate_pop_ptr->ind[i].xreal[j] > old)
				{
					lim[j] = old;
					count[j] = 1;
				}
				
			
			}

		}
		repair_all(&mate_pop_ptr->ind[i], lim, count);
	}
}

void real_mutate(population* old_pop_ptr, population* new_pop_ptr)
{   
	int p1 = 0, p2 = 0, p3 = 0;
	int i;
	int j;
	double value;
	
	for (i = 0; i < popsize; i++)
	{
		double lim[maxvar];
		int  count[maxvar];
		do {
			p1 = rnd(0, popsize - 1);
			p2 = rnd(0, popsize - 1);
			p3 = rnd(0, popsize - 1);

		} while (i == p1 || i == p2 || i == p3 || p1 == p2 || p1 == p3 || p2 == p3);

		for (j = 0; j < nvar; j++)
		{
			count[j] = 0;
			lim[j] = 0;
			value = old_pop_ptr->ind[p1].xreal[j] + pmut_r * (old_pop_ptr->ind[p2].xreal[j] - old_pop_ptr->ind[p3].xreal[j]);

			if (value < lim_r[j][0]|| value > lim_r[j][1])
			{
				value = repair_s(lim_r[j][0], lim_r[j][1], old_pop_ptr->ind[p1].xreal[j],value);
			}

			new_pop_ptr->ind[i].xreal[j] = value;
			if (new_pop_ptr->ind[i].xreal[j] < old_pop_ptr->ind[p1].xreal[j])
			{
				lim[j] = lim_r[j][0];
				count[j] = 1;
			}
			if (new_pop_ptr->ind[i].xreal[j] > old_pop_ptr->ind[p1].xreal[j])
			{
				lim[j] = old_pop_ptr->ind[p1].xreal[j];
				count[j] = 1;
			}
			
			
		}
		repair_all(&new_pop_ptr->ind[i], lim, count);
	}

}




void addInArchive(population* pop_ptr, population* arch2pop_ptr)
{
	int j, l, i;

	for (j = 0; j < popsize; j++)
	{
		pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
		arch2pop_ptr->ind_ptr = &(arch2pop_ptr->ind[j]);
		if (nvar > 0)
		{
			/*For Real Coded GA copying of the chromosomes*/
			for (l = 0; l < nvar; l++)
				arch2pop_ptr->ind_ptr->xreal[l] = pop_ptr->ind_ptr->xreal[l];
		}

		arch2pop_ptr->ind_ptr->error = pop_ptr->ind_ptr->error;

		/*Copying the fitness vector */
		for (l = 0; l < nfunc; l++)
			arch2pop_ptr->ind_ptr->fitness[l] = pop_ptr->ind_ptr->fitness[l];

		for (i = 0; i < ncons; i++)
		{
			arch2pop_ptr->ind_ptr->constr[i] = pop_ptr->ind_ptr->constr[i];
		}

	}   // end of j

}

void DivideBlock(population* pop_ptr1, population* pop_ptr2, population* pop_ptr3)
{
	int i, j, k, l;

	int min_index = -1;

	int no_selected = 0;

	int sel = 0;

	int p1;

	double max_fit[maxfun], min_fit[maxfun];

	double min, Diff;

	double block_interval_value;

	for (i = 0; i < popsize; i++)//copy two population
	{
		arch_globalpop.error[i] = pop_ptr1->ind[i].error;
		arch_globalpop.error[i + popsize] = pop_ptr2->ind[i].error;

		for (l = 0; l < ncons; l++)
		{
			arch_globalpop.constr[i][l] = pop_ptr1->ind[i].constr[l];
			arch_globalpop.constr[i + popsize][l] = pop_ptr2->ind[i].constr[l];
		}

		if (nvar > 0)
		{	/*For Real Coded GA x values are copied */
			for (k = 0; k < nvar; k++)
			{
				arch_globalpop.xreal[i][k] = pop_ptr1->ind[i].xreal[k];
				arch_globalpop.xreal[i + popsize][k] = pop_ptr2->ind[i].xreal[k];
			}
		}

		/*Fitness is copied to the global pool */
		for (l = 0; l < nfunc; l++)
		{
			arch_globalpop.fitness[i][l] = pop_ptr1->ind[i].fitness[l];
			arch_globalpop.fitness[i + popsize][l] = pop_ptr2->ind[i].fitness[l];
		}

	}

	for (i = 0; i < 2 * popsize; i++)
	{
		for (l = 0; l < nfunc; l++)
		{
			arch_globalpop.norm_fitness[i][l] = 0;
		}

		arch_globalpop.weight_fitness[i] = 0;

		arch_globalpop.rank[i] = 0;

		arch_globalpop.flag[i] = 0;

	}

	for (j = 0; j < nfunc; j++)
	{
		max_fit[j] = -EPS;
		min_fit[j] = INF;

		for (i = 0; i < 2 * popsize; i++)
		{
			if (arch_globalpop.fitness[i][j] > max_fit[j])
			{
				max_fit[j] = arch_globalpop.fitness[i][j];
			}

			if (arch_globalpop.fitness[i][j] < min_fit[j])
			{
				min_fit[j] = arch_globalpop.fitness[i][j];
			}

		}

		Diff = max_fit[j] - min_fit[j];

		if (Diff < 0.0)
		{
			printf("Something wrong in DivideBlock\n");
			exit(1);
		}

		for (i = 0; i < 2 * popsize; i++)
		{
			arch_globalpop.norm_fitness[i][j] = (arch_globalpop.fitness[i][j] - min_fit[j]) / Diff;
		}
	}

	for (i = 0; i < 2 * popsize; i++)
	{
		for (j = 0; j < nfunc; j++)
		{
			arch_globalpop.weight_fitness[i] += weight[j] * arch_globalpop.norm_fitness[i][j];
		}
	}



	for (j = 0; j < nfunc; j++)
	{
		Diff = max_fit[j] - min_fit[j];

		block_interval_value = Diff / (BLOCK_INTERVAL - 1);//evenly divide the range into BLOCK_INTERVAL blocks

		for (i = 0; i < 2 * popsize; i++)
		{
			arch_globalpop.rank[i] = (int)((arch_globalpop.fitness[i][j] - min_fit[j]) / block_interval_value) + 1;
		}

		for (l = 1; l <= BLOCK_INTERVAL; l++)
		{
			min = INF;

			for (i = 0; i < 2 * popsize; i++)
			{
				if (arch_globalpop.rank[i] == l)
				{
					if (arch_globalpop.weight_fitness[i] < min)
					{
						min = arch_globalpop.weight_fitness[i];

						min_index = i;//the smallest fitness
					}
				}
			}

			arch_globalpop.flag[min_index] = 1;


		}
	}

	no_selected = 0;
	for (i = 0; i < 2 * popsize; i++)
	{
		if (arch_globalpop.flag[i] == 1)
			no_selected++;
	}

	while (no_selected > popsize)//select randomly sel individuals
	{
		do {

			p1 = rnd(0, 2 * popsize - 1);

		} while (arch_globalpop.flag[p1] == 0);

		arch_globalpop.flag[p1] = 0;

		no_selected--;
	}

	if (no_selected < popsize)
	{

		sel = popsize - no_selected;

		while (sel > 0)//select randomly sel individuals
		{
			do {

				p1 = rnd(0, 2 * popsize - 1);

			} while (arch_globalpop.flag[p1] == 1);

			arch_globalpop.flag[p1] = 1;

			sel--;
		}
	}
	j = 0;
	for (i = 0; i < 2 * popsize; i++)
	{
		if (arch_globalpop.flag[i] == 1)
		{
			if (nvar > 0)
			{	/*For Real Coded GA x values are copied */
				for (k = 0; k < nvar; k++)
				{
					pop_ptr3->ind[j].xreal[k] = arch_globalpop.xreal[i][k];

				}
			}

			/*Fitness is copied to the global pool */
			for (l = 0; l < nfunc; l++)
			{
				pop_ptr3->ind[j].fitness[l] = arch_globalpop.fitness[i][l];

			}

			pop_ptr3->ind[j].error = arch_globalpop.error[i];

			for (l = 0; l < ncons; l++)
			{
				pop_ptr3->ind[j].constr[l] = arch_globalpop.constr[i][l];
			}

			j++;

		}
	}

}


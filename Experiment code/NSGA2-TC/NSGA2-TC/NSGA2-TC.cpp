// MSD-NSGA2.cpp : 定义控制台应用程序的入口点。
//
/* This is a Multi-Stage Dynamic NSGA-II for the Software Testing Resource Allocation Problem.
**********************************************************************
*  This VS2013 program is developed by Dr. Guofu Zhang               *
*  School of Computer and Information                                *
*  Hefei University of Technology                                    *
*  zgf@hfut.edu.cn                                                   *
*The C code of nsga2 is from http://www.iitk.ac.in/kangal/codes.shtml*
**********************************************************************
*/

#include "stdafx.h"

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include<iostream>
/********************************************************************/
#include "random.h"       /*Random Number Generator*/
using namespace std;
/********************************************************************/
# define INF 1.0e14/**/
# define EPS 1.0e-14/**/
# define PI 3.14159265358979
/********************************************************************/
#define maxpop   300 /*Max population */

/*Max no. of variables*/


#define maxvar    100/*simple system*/
//#define maxvar    30/*complex system*/
//#define maxvar    50/*large system*/

#define maxfun    3  /*Max no. of functions */
#define maxcons 3 /*Max no. of Constraints*/
/********************************************************************/
#define MAX_EXPERIMENT  30  /*Experimental times*/
#define MAX_GENERATION  200 /*Max no. of generation*/
//#define MAX_GENERATION  200
//#define MAX_GENERATION  500 
								/********************************************************************/
#define MAX_R 0.9
double C;
#define T_MAX_RESOURCE  150000.0
//#define T_MAX_RESOURCE  30000.0
//#define T_MAX_RESOURCE  70000.0


#define EPSION  0.000001

//#define MISSION_TIME  200.0
#define MISSION_TIME  200.0
//#define MISSION_TIME  200.0  //large system

/**************************************************************************/
								//计算测试资源约束的变量

double sum_u = 0.0;/*所有构件预期访问次数之和*/

double time_stage; /*the value of testing time in each stage*/

				   /********************************************************************/

				   /********************************************************************/
int gener;       /*No of generations*/
int nvar;          /*No of variables*/
int ncons;         /*No of Constraints*/
				   /********************************************************************/

				   /********************************************************************/
double seed;      /*Random Seed*/
double pcross;        /*Cross-over Probability*/
double pmut_r;          /*Mutation Probability*/
double Mis_tim[maxvar];/*构件分配到的任务时间*/
double lim_r[maxvar][2];/*Limits of variable in array*/
						//double min_orel[MAX_SYSTEM];
double Y[maxvar];
double tmid[maxvar];//测试资源中间值；



double di;       /*Distribution Index for the Cross-over*/
double dim;          /*Distribution Index for the Mutation*/
double delta_fit;     /* variables required for fitness for fitness sharing */
double min_fit;
double front_ratio;

int nfunc;       /*No of functions*/
static int popsize;/*Population Size*/
				   /********************************************************************/

				   /********************************************************************/
typedef struct       /*individual properties*/
{
	int rank,              /*Rank of the individual*/
		flag;              /*Flag for ranking*/
	double xreal[maxvar]; /*list of real variables*/
	double fitness[maxfun];/*Fitness values */
	double constr[maxcons];     /*Constraints values*/
	double cub_len;            /*crowding distance of the individual*/
	double error;              /* overall constraint violation for the individual*/
}individual;        /*Structure defining individual*/

typedef struct
{
	int maxrank;            /*Maximum rank present in the population*/
	double rankrat[maxpop];  /*Rank Ratio*/
	int rankno[maxpop];     /*Individual at different ranks*/
	individual ind[maxpop]; /*Different Individuals*/
	individual* ind_ptr;
}population;             /*Population Structure*/
						 /********************************************************************/
double para_module[maxvar][6]; double u[maxvar];

/********************************************************************/
typedef struct
{
	int maxrank,   /*Max rank of the global population*/
		rankar[2 * maxpop][2 * maxpop], /*record of array of individual numbers at a particular rank */
		rankno[2 * maxpop];           /*record of no. of individuals at a particular rank*/

	int rank[2 * maxpop],            /*rank of different individuals*/
		flag[2 * maxpop];            /*Setting the flag */

	double fitness[2 * maxpop][maxfun], /*Fitness function values for the different individuals*/
		cub_len[2 * maxpop],              /*Dummy fitness*/
		xreal[2 * maxpop][maxvar],       /*value of the decoded variables for different individuals */
		error[2 * maxpop],               /*Error Values of the individuals*/
		constr[2 * maxpop][maxcons];
}globpop;/*Population structure for the pool having both the old as well as new population*/

globpop globalpop, * global_pop_ptr;
int Lastrank;
double fpara1[2 * maxpop][2];
/********************************************************************/

/********************************************************************/
void input(); /*Input Parameters from user*/

void realinit(population* pop_ptr);/*initializes the population*/
void Cal_res_lim(double r);/*计算构件测试资源上下限*/
void Cal_C();


double cal_reliab(int ind_module, double m_res);//calculate the reliability of each module
double cal_cost(int ind_module, double m_res);//calculate the cost of each module
void func(population* pop_ptr, double r);/*evaluate the value of the function & errors*/

void rankcon(population* pop_ptr);/*demarcates the different Pareto Fronts*/
int indcmp3(double* ptr1, double* ptr2);

void nselect(population* old_pop_ptr, population* pop2_ptr);/*get the different individuals selected*/

void realcross(population* new_pop_ptr, population* mate_pop_ptr);/*simulated binary crossover for Real Coded GA*/

void real_mutate(population* new_pop_ptr);/*formulate the mutation routine*/

void grankc(int gen);/*Ranking the global pool when the constraints are there*/
int indcmp1(double* ptr1, double* ptr2);/*Comparison of the variables*/
void gsort(int rnk, int sel);/*Sorting for the function values in ascending order*/
void gshare(int rnk);/*Sharing the fitness*/
void sort(int rnk);
void keepalive(population* pop1_ptr, population* pop2_ptr, population* pop3_ptr, int gen);/*keep the fronts alive (caring the end problem)*/
																						  /********************************************************************/
void sleep(clock_t wait)
{
	clock_t goal;
	goal = wait + clock();
	while (goal > clock());
}

int _tmain(int argc, _TCHAR* argv[])
{
	int o = 1;
	for (o = 1; o < 31; o++)
	{
		//读取文件数据
		FILE* end_ptr2; FILE* end_ptr1;
		char file3[500]; char file2[500];

		int i;
		sprintf_s(file3, "D:\\论文代码数据\\Parameters\\single-input and single-output\\NSGA2\\组件参数\\100\\%d.txt", o);
		fopen_s(&end_ptr2, file3, "r+");
		//在file3中写上每次迭代产生的不同fitness文件的路径

		for (i = 0; i < maxvar; i++)
			fscanf_s(end_ptr2, "%lf %lf %lf %lf %lf", &para_module[i][0], &para_module[i][1], &para_module[i][2], &para_module[i][3], &para_module[i][4]); // 循环读


		fclose(end_ptr2); //关闭文件

		sprintf_s(file2, "D:\\论文代码数据\\Parameters\\single-input and single-output\\NSGA2\\预期访问次数\\100\\%d.txt", o);
		fopen_s(&end_ptr1, file2, "r+");
		for (i = 0; i < maxvar; i++)
			fscanf_s(end_ptr1, "%lf", &u[i]); // 循环读
		fclose(end_ptr1); //关闭文件



		int random = rand() % 1000;

		seed = (float)random / 1000.0;
		srand((unsigned)time(NULL));
		while (seed <= 0.0 || seed >= 1.0)
		{
			int random = rand() % 1000;

			seed = (float)random / 1000.0;
		}

		warmup_random(seed);
		time_stage = T_MAX_RESOURCE;

		input();
		Cal_res_lim(MAX_R);
		Cal_C();
		int int_exp = 0;

		for (int_exp = 0; int_exp < MAX_EXPERIMENT; int_exp++)
		{
			printf_s("Run experiment %d!\n", int_exp + 1);
			printf("1");
			//	sleep((clock_t)3 * CLOCKS_PER_SEC);	

			FILE* end_ptr;
			char file1[500];
			sprintf_s(file1, "D:\\论文代码数据\\experimental data\\single-input and single-output\\NSGA2-TC\\100\\Data of 30 instances\\test%d\\fitness%d.txt", o, int_exp + 1);//在file1中写上每次迭代产生的不同fitness文件的路径
			fopen_s(&end_ptr, file1, "w+");



			int int_exp = 0;


			//FILE *end_ptr;
			//char file1[500];
			//sprintf_s(file1, "F:\\test1\\const_nsga2\\d%d\\fitness%d.txt",z, int_exp + 1);//在file1中写上每次迭代产生的不同fitness文件的路径
			//fopen_s(&end_ptr, file1, "w+");
			//z++;
			//	clock_t start, finish;
			//	double duration;

			int i, j, l, f, maxrank1;
			double tot; int k;

			population oldpop, newpop, matepop, * old_pop_ptr, * new_pop_ptr, * mate_pop_ptr;/*Defining the population Structures*/




			seed = (1.0 / (MAX_EXPERIMENT + 1)) * (int_exp + 1);
			warmup_random(seed);

			//	start = clock();

			old_pop_ptr = &(oldpop);
			realinit(old_pop_ptr);

			old_pop_ptr = &(oldpop);
			new_pop_ptr = &(newpop);

			for (j = 0; j < popsize; j++)
			{
				/*Initializing the Rank array having different individuals at a particular  rank to zero*/
				old_pop_ptr->rankno[j] = 0;
				new_pop_ptr->rankno[j] = 0;
			}

			old_pop_ptr = &(oldpop);
			func(old_pop_ptr, MAX_R);/*Function Calculation*/

							  /********************************************************************/
							  /*----------------------GENERATION STARTS HERE----------------------*/
			for (i = 0; i < gener; i++)
			{
				printf("1");
				old_pop_ptr = &(oldpop);
				mate_pop_ptr = &(matepop);

				nselect(old_pop_ptr, mate_pop_ptr);/*--------SELECT----------------*/
												   /*for (k = 0; k < maxpop; k++){
												   for (j = 0; j < maxvar; j++){
												   if (mate_pop_ptr->ind[k].xreal[j]>lim_r[j][1]){
												   mate_pop_ptr->ind[k].xreal[j] = lim_r[j][0] + (lim_r[j][1] - lim_r[j][0])*randomperc();
												   printf("交叉约束使值超出上界");

												   }
												   }
												   }*/
												   /********************************************************************/

				new_pop_ptr = &(newpop);
				mate_pop_ptr = &(matepop);

				realcross(new_pop_ptr, mate_pop_ptr);/*Real Cross-over*/

													 /********************************************************************/

				new_pop_ptr = &(newpop);
				real_mutate(new_pop_ptr); /*Real Mutation*/

				new_pop_ptr = &(newpop);
				func(new_pop_ptr, MAX_R);/*----------FUNCTION EVALUATION-----------*/

								  /********************************************************************/

								  /*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
				old_pop_ptr = &(oldpop);
				new_pop_ptr = &(newpop);

				mate_pop_ptr = &(matepop);

				keepalive(old_pop_ptr, new_pop_ptr, mate_pop_ptr, i + 1);/*Elitism And Sharing Implemented*/

																		 /********************************************************************/

																		 /*----------------Rank Ratio Calculation------------------------*/
				new_pop_ptr = &(matepop);
				old_pop_ptr = &(oldpop);
				if (old_pop_ptr->maxrank > new_pop_ptr->maxrank)/*Finding the greater maxrank among the two populations*/
					maxrank1 = old_pop_ptr->maxrank;
				else
					maxrank1 = new_pop_ptr->maxrank;

				for (j = 0; j < maxrank1; j++)
				{
					/*Sum of the no of individuals at any rank in old population and the new populaion*/

					tot = (old_pop_ptr->rankno[j]) + (new_pop_ptr->rankno[j]);

					/*Finding the rank ratio for new population at this rank*/

					new_pop_ptr->rankrat[j] = (new_pop_ptr->rankno[j]) / tot;

				}

				/********************************************************************/

				/*=======Copying the new population to old population======*/
				old_pop_ptr = &(oldpop);
				new_pop_ptr = &(matepop);
				for (j = 0; j < popsize; j++)
				{
					old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
					new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);
					if (nvar > 0)
					{
						/*For Real Coded GA copying of the chromosomes*/
						for (l = 0; l < nvar; l++)
							old_pop_ptr->ind_ptr->xreal[l] = new_pop_ptr->ind_ptr->xreal[l];
					}

					/*Copying the fitness vector */
					for (l = 0; l < nfunc; l++)
						old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];

					/*Copying the dummy fitness*/
					old_pop_ptr->ind_ptr->cub_len = new_pop_ptr->ind_ptr->cub_len;

					/*Copying the rank of the individuals*/
					old_pop_ptr->ind_ptr->rank = new_pop_ptr->ind_ptr->rank;

					/*Copying the error and constraints of the individual*/

					old_pop_ptr->ind_ptr->error = new_pop_ptr->ind_ptr->error;
					for (l = 0; l < ncons; l++)
					{
						old_pop_ptr->ind_ptr->constr[l] = new_pop_ptr->ind_ptr->constr[l];
					}

					/*Copying the flag of the individuals*/
					old_pop_ptr->ind_ptr->flag = new_pop_ptr->ind_ptr->flag;
				}   // end of j

				maxrank1 = new_pop_ptr->maxrank;
				/*Copying the array having the record of the individual	at different ranks */
				for (l = 0; l < popsize; l++)
				{
					old_pop_ptr->rankno[l] = new_pop_ptr->rankno[l];
				}
				/*Copying the maxrank */
				old_pop_ptr->maxrank = new_pop_ptr->maxrank;

				/********************************************************************/

				/*Printing the fitness record for last generation in a file last*/
				if (i == gener - 1)
				{  // for the last generation 
					old_pop_ptr = &(matepop);
					for (f = 0; f < popsize; f++) // for printing
					{
						old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[f]);

						if ((old_pop_ptr->ind_ptr->error <= 0.0) && (old_pop_ptr->ind_ptr->rank == 1))  // for all feasible solutions and non-dominated solutions
						{
							for (l = 0; l < nfunc; l++) {
								if (l == 0) {
									fprintf(end_ptr, "%5.6f\t", old_pop_ptr->ind_ptr->fitness[l]);
								}
								else
									fprintf(end_ptr, "%5.6f\t", old_pop_ptr->ind_ptr->fitness[l]);
							}
							fprintf(end_ptr, "%s", "\n");

						}  // feasibility check
					} // end of f (printing)
				}

			}// end of i

			 //finish = clock();
			 //duration = (double)(finish - start) / CLOCKS_PER_SEC;
			 //fprintf_s(run_time, "%f ", duration);

			fclose(end_ptr);



			//fclose(run_time);

			//printf_s("please input a char to exit the program!\n");

			//getchar();
			/*fclose(end_ptr);*/

		}
	}
	return 0;
}

void input()
{
	/*number of the variables*/
	nvar = maxvar;

	/*number of the functions*/
	nfunc = 3;

	/*number of  constraints*/
	ncons = maxcons;

	/*number of the individuals in the population*/
	popsize = maxpop;

	/*No. of generations*/
	gener = MAX_GENERATION;

	/*the crossover probability (between 0.5 and 1)*/
	pcross = 0.9;

	/*the mutation probability (between 0 and cc)*//*cc=1.0/nvar*/
	pmut_r = 1.0 / maxvar;

	/*Distribution Index for real-coded crossover between 0.5 to 100*/
	di = 100;

	/*Distribution Index for real-coded mutation between 0.5 to 500*/
	dim = 500;

}

void realinit(population* pop_ptr)
{
	int i, j;
	double d;
	double sum;
	double low_T = 0;
	for (i = 0; i < popsize; i++)
	{
		sum = 0;
		for (j = 0; j < nvar; j++)
		{
			d = randomperc();
			/*printf("%f\t", d);*/
			/*limits are specified it generates the value in range of minimum and maximum value of the variable*/
			pop_ptr->ind[i].xreal[j] = d * (lim_r[j][1] - lim_r[j][0]) + lim_r[j][0];
			sum = sum + pop_ptr->ind[i].xreal[j];

		}

		if (sum > time_stage)
		{

			for (j = 0; j < nvar; j++)
			{
				d = randomperc();

				pop_ptr->ind[i].xreal[j] = (pop_ptr->ind[i].xreal[j]) * time_stage / sum;
			}
		}

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

void Cal_res_lim(double r)
{
	//double Tmin_sum = 0;/*组件测试资源下限之和*/
	double sum_u = 0.0;/*所有组件预期访问次数之和*/
	double left1[maxvar], right1[maxvar];

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

	//计算单个组件达到可靠性约束时的测试资源

	double time0 = 0.0;
	for (int i = 0; i < maxvar; i++) {
		lim_r[i][0] = 0.0;
		lim_r[i][1] = T_MAX_RESOURCE;

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
void func(population* pop_ptr, double r)
{
	double* realx_ptr, /*Pointer to the array of x values*/
		* fitn_ptr,      /*Pointer to the array of fitness function*/
		x[2 * maxvar],     /* problem variables */
		f[maxfun],     /*array of fitness values*/
		* err_ptr,      /*Pointer to the error */
		cstr[maxcons];
	/*系统失效产生的单位成本*/
	int i, j, k;
	double error, cc;



	double sum1 = 1.0, sum2 = 0, sum3 = 0, sum4 = 0;
	//added by Guofu Zhang

	pop_ptr->ind_ptr = &(pop_ptr->ind[0]);


	for (i = 0; i < popsize; i++)
	{
		pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
		realx_ptr = &(pop_ptr->ind_ptr->xreal[0]);
		for (j = 0; j < nvar; j++)
		{ // Real-coded variables 
			x[j] = *realx_ptr++;


			if (x[j] < lim_r[j][0])
				sum3 = sum3 - 1;//Calculate the number of values which are <0.0

		}

		fitn_ptr = &(pop_ptr->ind_ptr->fitness[0]);
		err_ptr = &(pop_ptr->ind_ptr->error);

		/*First fitness function */
		/* calculate the reliability of each subsystem*/
		sum1 = 1.0;
		for (k = 0; k < maxvar; k++)
		{
			sum1 = sum1 * cal_reliab(k, x[k]);
		}

		f[0] = 1 - sum1;/*here we use 1-sum1 to minimize f[0]*/

						// Second and Third Fitness Functions
		sum2 = 0, sum4 = 0;
		for (j = 0; j < nvar; j++)
		{
			sum4 = sum4 + cal_cost(j, x[j]);
			sum2 = sum2 + x[j];
		}
		f[1] = sum4;

		f[2] = sum2;
		if (sum2 - time_stage > 0.1)
		{
			printf_s("error in the upper bound constraint\n");
			exit(1);
		}

		/******************************************************************/
		/*              Put The Constraints Here                          */
		/******************************************************************/

		cstr[0] = sum1 - r;/*if cstr[1]<0.0, we have sum3>time_stage*/
		cstr[1] = (C - f[1]) / C;/*if cstr[0]<0.0, we have at least a x[j]<0.0*/
		cstr[2] = (T_MAX_RESOURCE - f[2]) / T_MAX_RESOURCE;

		//cstr[2] = (time_stage - sum2);

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

	/*---------------------------* RANKING *------------------------------*/

	rankcon(pop_ptr);

}
void rankcon(population* pop_ptr)
{
	int i, j, k,       /*counters*/
		rnk,           /*rank*/
		val,           /*value obtained after comparing two individuals*/
		nondom,        /*no of non dominated members*/
		maxrank1,      /*Max rank of the population*/
		rankarr[maxpop], /*Array storing the individual number at a rank*/
		q;

	double* ptr1, * ptr2, * err_ptr1, * err_ptr2;

	/*------------------------------* RANKING *------------------------------*/

	/*Initializing the ranks to zero*/
	rnk = 0;

	nondom = 0;
	maxrank1 = 0;
	/*min_fit is initialize to start distributing the dummy fitness =
	popsize to the rank one individuals and keeping the record such
	that the minimum fitness of the better rank individual is always
	greater than max fitness of the relatively worse rank*/

	min_fit = popsize;

	/*Difference in the fitness of minimum dummy fitness of better rank
	and max fitness of the next ranked individuals*/

	delta_fit = 0.1 * popsize;
	/*Initializing all the flags to 2*/

	for (j = 0; j < popsize; j++)
	{
		pop_ptr->ind[j].flag = 2;
	}

	q = 0;

	for (k = 0; k < popsize; k++, q = 0)
	{
		for (j = 0; j < popsize; j++)
		{
			if (pop_ptr->ind[j].flag != 1)break;
			/*Break if all the individuals are assigned a rank*/
		}
		if (j == popsize)break;

		rnk = rnk + 1;

		for (j = 0; j < popsize; j++)
		{
			if (pop_ptr->ind[j].flag == 0) pop_ptr->ind[j].flag = 2;
			/*Set the flag of dominated individuals to 2*/
		}

		for (i = 0; i < popsize; i++)
		{
			/*Select an individual which rank to be assigned*/

			pop_ptr->ind_ptr = &(pop_ptr->ind[i]);

			if (pop_ptr->ind_ptr->flag != 1 && pop_ptr->ind_ptr->flag != 0)
			{
				ptr1 = &(pop_ptr->ind_ptr->fitness[0]);
				err_ptr1 = &(pop_ptr->ind_ptr->error);

				for (j = 0; j < popsize; j++)
				{

					/*Select the other individual which has not got a rank*/
					if (i != j)
					{
						if (pop_ptr->ind[j].flag != 1)
						{
							pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
							ptr2 = &(pop_ptr->ind_ptr->fitness[0]);
							err_ptr2 = &(pop_ptr->ind_ptr->error);

							if (*err_ptr1 < EPSION && *err_ptr2 > EPSION)
							{
								/*first ind is feasible second individual is infeasible*/
								pop_ptr->ind[j].flag = 0;
							}
							else
							{
								if (*err_ptr1 > EPSION && *err_ptr2 < EPSION)
								{
									/*first individual is infeasible and second is feasible*/
									pop_ptr->ind[i].flag = 0;
									break;
								}
								else
								{
									/*both are feasible or both are infeasible*/
									if (*err_ptr1 > * err_ptr2)
									{
										pop_ptr->ind[i].flag = 0;
										/*first individual is more infeasible*/
										break;
									}
									else
									{
										if (*err_ptr1 < *err_ptr2)
										{
											pop_ptr->ind[j].flag = 0;
											/*second individual is more
											infeasible*/
										}
										else
										{
											/*Compare the two individuals for
											fitness*/
											val = indcmp3(ptr1, ptr2);
											/*VAL = 2 for dominated individual which rank to be given*/
											/*VAL = 1 for dominating individual	which rank to be given*/
											/*VAL = 3 for non comparable individuals*/

											if (val == 2)
											{
												pop_ptr->ind[i].flag = 0;/* individual 1 is dominated */
												break;
											}

											if (val == 1)
											{
												pop_ptr->ind[j].flag = 0;/* individual 2 is dominated */
											}

											if (val == 3)
											{
												nondom++;
												/* individual 1 & 2 are	non dominated */
												if (pop_ptr->ind[j].flag != 0)
													pop_ptr->ind[j].flag = 3;
											}

										}   /*if loop ends*/
									}       /* i != j loop ends*/
								}
							}
						}
					}
				}        /*loop over j ends*/
				if (j == popsize)
				{
					/*Assign the rank and set the flag*/
					pop_ptr->ind[i].rank = rnk;
					pop_ptr->ind[i].flag = 1;
					rankarr[q] = i;
					q++;
				}
			}       /*Loop over flag check ends*/
		}           /*Loop over i ends */
		pop_ptr->rankno[rnk - 1] = q;
	}
	maxrank1 = rnk;

	/*     Find Max Rank of the population    */
	for (i = 0; i < popsize; i++)
	{
		rnk = pop_ptr->ind[i].rank;

		if (rnk > maxrank1)maxrank1 = rnk;

	}

	pop_ptr->maxrank = maxrank1;

}


int indcmp3(double* ptr1, double* ptr2)/*Routine Comparing the two individuals*/
{
	double fit1[maxfun], fit2[maxfun];
	int i, value, m, n;
	for (i = 0; i < nfunc; i++)
	{
		fit1[i] = *ptr1++;
		fit2[i] = *ptr2++;
	}
	m = 0;
	n = 0;
	while (m < nfunc && fit1[m] <= fit2[m])
	{
		if (fit1[m] == fit2[m]) n++;
		m++;
	}
	if (m == nfunc)
	{
		if (n == nfunc) value = 3;
		else value = 1;             /*value = 1 for domination*/
	}
	else
	{
		m = 0;
		n = 0;
		while (m < nfunc && fit1[m] >= fit2[m])
		{
			if (fit1[m] == fit2[m]) n++;
			m++;
		}
		if (m == nfunc)
		{
			if (n != nfunc)
				value = 2;                       /*value =  2 for dominated */
			else value = 3;
		}
		else value = 3;                   /*value = 3 for incomparable*/
	}

	return value;
}

void nselect(population* old_pop_ptr, population* pop2_ptr)
{
	int* fit_ptr1, * fit_ptr2;

	double rnd2, * f1_ptr, * f2_ptr;

	double* select_ptr_r, * s1_ptr_r, * s2_ptr_r;

	individual* j, * j1;

	int i, rnd, rnd1, k, n, j2, r;

	old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);

	pop2_ptr->ind_ptr = &(pop2_ptr->ind[0]);

	j = &(old_pop_ptr->ind[popsize - 1]);

	old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
	j2 = 0;
	r = popsize;

	for (n = 0, k = 0; n < popsize; n++, k++)
	{
		pop2_ptr->ind_ptr = &(pop2_ptr->ind[k]);
		select_ptr_r = &(pop2_ptr->ind_ptr->xreal[0]);

		rnd2 = randomperc();

		rnd2 = popsize * rnd2;

		rnd = (int)floor(rnd2);

		if (rnd == 0)
			rnd = popsize - k;

		if (rnd == popsize)
			rnd = (popsize - 2) / 2;

		/*Select first parent randomly*/
		j = &(old_pop_ptr->ind[rnd - 1]);

		rnd2 = randomperc();

		rnd2 = popsize * rnd2;

		rnd1 = (int)floor(rnd2);

		if (rnd1 == 0)
			rnd1 = popsize - n;

		if (rnd1 == popsize)
			rnd1 = (popsize - 4) / 2;

		/*Select second parent randomly*/
		j1 = &(old_pop_ptr->ind[rnd1 - 1]);

		old_pop_ptr->ind_ptr = j;

		s1_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
		fit_ptr1 = &(old_pop_ptr->ind_ptr->rank);
		f1_ptr = &(old_pop_ptr->ind_ptr->cub_len);

		old_pop_ptr->ind_ptr = j1;
		s2_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
		fit_ptr2 = &(old_pop_ptr->ind_ptr->rank);
		f2_ptr = &(old_pop_ptr->ind_ptr->cub_len);
		/*--------------------------------------------------------------------------*/

		/*------------------SELECTION PROCEDURE------------------------------------*/

		/*Comparing the fitnesses*/

		if (*fit_ptr1 > * fit_ptr2)
		{
			for (i = 0; i < nvar; i++)
				*select_ptr_r++ = *s2_ptr_r++;
		}
		else
		{
			if (*fit_ptr1 < *fit_ptr2)
			{
				for (i = 0; i < nvar; i++)
					*select_ptr_r++ = *s1_ptr_r++;
			}
			else
			{
				if (*f1_ptr < *f2_ptr)
				{
					for (i = 0; i < nvar; i++)
						*select_ptr_r++ = *s2_ptr_r++;
				}
				else
				{
					for (i = 0; i < nvar; i++)
						*select_ptr_r++ = *s1_ptr_r++;
				}
			}
		}
	}

}
void realcross(population* new_pop_ptr, population* mate_pop_ptr)
{
	int i, j, y, n;
	double rnd1, rnd3, par1, par2, chld1, chld2, betaq, beta, alpha;
	double y1, y2, yu, yl, expp;

	int k, q;
	double rnd2, low_T;
	double sum_chld1, sum_chld2, sum_cross1, sum_cross2, temp1, temp2;
	int no_cross, pos_cross[maxvar];//save the crossover points in the encoding
	double gama;//the scale factor for controlling the enlargement and shrink of gene values

	y = 0; n = 0;
	for (i = 0; i < popsize / 2; i++)
	{
		rnd1 = randomperc();
		low_T = 0;
		/*Check Whether the cross-over to be performed*/
		if (rnd1 <= pcross)
		{
			sum_chld1 = 0;
			sum_chld2 = 0;
			sum_cross1 = 0;
			sum_cross2 = 0;

			no_cross = 0;
			for (j = 0; j < nvar; j++)
				pos_cross[j] = -1;

			q = -1;

			/*Loop over no of variables*/
			for (j = 0; j < nvar; j++)
			{

				/*Selected Two Parents*/
				par1 = mate_pop_ptr->ind[y].xreal[j];
				par2 = mate_pop_ptr->ind[y + 1].xreal[j];


				yl = lim_r[j][0];
				yu = lim_r[j][1];

				rnd3 = randomperc();

				/* Check whether variable is selected or not*/
				if (rnd3 <= 0.5)
				{
					pos_cross[no_cross] = j;
					no_cross++;
					/*Variable selected*/

					if (fabs(par1 - par2) > EPS) // changed by Deb (31/10/01)
					{
						if (par2 > par1)
						{
							y2 = par2;
							y1 = par1;
						}
						else
						{
							y2 = par1;
							y1 = par2;
						}

						/*Find beta value*/
						if ((y1 - yl) > (yu - y2))
						{
							beta = 1 + (2 * (yu - y2) / (y2 - y1));
							//printf("beta = %f\n",beta);
						}
						else
						{
							beta = 1 + (2 * (y1 - yl) / (y2 - y1));
							//printf("beta = %f\n",beta);
						}

						/*Find alpha*/
						expp = di + 1.0;

						beta = 1.0 / beta;

						alpha = 2.0 - pow(beta, expp);

						if (alpha < 0.0)
						{
							printf("ERRRROR %f %d %d %f %f\n", alpha, y, n, par1, par2);
							exit(-1);
						}

						rnd2 = randomperc();

						if (rnd2 <= 1.0 / alpha)
						{
							alpha = alpha * rnd2;
							expp = 1.0 / (di + 1.0);
							betaq = pow(alpha, expp);
						}
						else
						{
							alpha = alpha * rnd2;
							alpha = 1.0 / (2.0 - alpha);
							expp = 1.0 / (di + 1.0);
							if (alpha < 0.0)
							{
								printf("ERRRORRR \n");
								exit(-1);
							}
							betaq = pow(alpha, expp);
						}

						/*Generating two children*/
						chld1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
						chld2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

					}
					else
					{

						betaq = 1.0;
						y1 = par1; y2 = par2;

						/*Generation two children*/
						chld1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
						chld2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

					}

					/*Repair encoding*/
					if (chld1 < lim_r[j][0] || chld2 > lim_r[j][1])//This means that chld1 is reduced, then it is certain that chld2 is enlarged
					{/*Select a minimum interval to generate a child to ensure the other child obeys the upper or lower constraint*/

						if ((y1 - lim_r[j][0]) <= (lim_r[j][1] - y2))
						{
							rnd2 = randomperc();
							chld1 = rnd2 * (y1 - lim_r[j][0]) + lim_r[j][0];//generate a value over the range (min_time[j],y1) to ensure chld1 obeys the lower constraint
							chld2 = y1 + y2 - chld1;
						}
						else {
							rnd2 = randomperc();
							chld2 = rnd2 * (lim_r[j][1] - y2) + y2;//generate a value over the range (y2,yu) to ensure chld2 obeys the upper constraint
							chld1 = y1 + y2 - chld2;
						}
					}

					if (chld1 > lim_r[j][1])
					{
						printf("ERRRORRR for chld1\n");
						exit(-1);
						//	rnd2 = randomperc();
						//	chld1 = rnd2 * (yu - y1) + y1; 
						//	chld2 = y1 + y2 - chld1;

					}
					sum_cross2 = sum_cross2 + chld2;//the sum of variables in the crossover points
					sum_cross1 = sum_cross1 + chld1;//the sum of variables in the crossover points

				}
				else
				{

					/*Copying the children to parents*/
					chld1 = par1;
					chld2 = par2;
				}
				new_pop_ptr->ind[n].xreal[j] = chld1;
				new_pop_ptr->ind[n + 1].xreal[j] = chld2;

				sum_chld2 = sum_chld2 + chld2;//the sum of all variables
				sum_chld1 = sum_chld1 + chld1;
			}

			/*Repair encoding*/
			if (sum_chld1 > time_stage)/*chld1 violates the upper constraint,then it is certain than sum_chld2<yu2 based on the average property*/
			{

				if (sum_cross1 - time_stage + sum_chld2 <= 0.0)
					gama = 0;
				else {
					if (time_stage - sum_chld1 + sum_cross1 == 0)
						gama = 1.0;
					else
						gama = (sum_cross1 - time_stage + sum_chld2) / (time_stage - sum_chld1 + sum_cross1);
					if (gama > 1.0)
						gama = 1.0;	//if (time_stage - sum_chld1 + sum_cross1)==0 							
				}
				for (k = 0; k < no_cross; k++)//Only repair the mutation genes to maintain and inherit the parent's information
				{
					q = pos_cross[k];

					temp1 = new_pop_ptr->ind[n].xreal[q];//chld1
					temp2 = new_pop_ptr->ind[n + 1].xreal[q];//chld2
					new_pop_ptr->ind[n].xreal[q] = (temp1 - lim_r[q][0]) * ((1.0 - gama) + gama) * ((time_stage - sum_chld1 + sum_cross1 - low_T) / (sum_cross1 - low_T)) + lim_r[q][0];
					new_pop_ptr->ind[n + 1].xreal[q] = temp1 + temp2 - new_pop_ptr->ind[n].xreal[q];
				}


			}

			if (sum_chld2 > time_stage)/*chld2 violates the upper constraint,then it is certain than sum_chld1<yu1 based on the average property*/
			{

				if (sum_cross2 - time_stage + sum_chld1 <= 0.0)
					gama = 0;
				else {
					if (time_stage - sum_chld2 + sum_cross2 == 0)
						gama = 1.0;
					else
						gama = (sum_cross2 - time_stage + sum_chld1) / (time_stage - sum_chld2 + sum_cross2);

					if (gama > 1.0)
						gama = 1.0;// if (yu - sum_chld2 + sum_cross2)==0							
				}

				for (k = 0; k < no_cross; k++)//Only repair the mutation genes to maintain and inherit the parent's information
				{
					q = pos_cross[k];

					temp1 = new_pop_ptr->ind[n + 1].xreal[q];//chld2
					temp2 = new_pop_ptr->ind[n].xreal[q];//chld1
					new_pop_ptr->ind[n + 1].xreal[q] = (temp1 - lim_r[q][0]) * ((1.0 - gama) + gama) * ((time_stage - sum_chld2 + sum_cross2 - low_T) / (sum_cross2 - low_T)) + lim_r[q][0];
					new_pop_ptr->ind[n].xreal[q] = temp1 + temp2 - new_pop_ptr->ind[n + 1].xreal[q];
				}

			}

		}
		else
		{
			for (j = 0; j < nvar; j++)
			{
				par1 = mate_pop_ptr->ind[y].xreal[j];
				par2 = mate_pop_ptr->ind[y + 1].xreal[j];
				chld1 = par1;
				chld2 = par2;
				new_pop_ptr->ind[n].xreal[j] = chld1;
				new_pop_ptr->ind[n + 1].xreal[j] = chld2;
			}
		}
		n = n + 2; y = y + 2;
	}
}

void real_mutate(population* new_pop_ptr)
{
	int i, j;

	int k, q;

	double rnd1, rnd3, delta, indi, deltaq;

	double y, yl, yu, val, xy, low_T;

	double rnd2, chld, sum, sum_muta;
	int no_muta, pos_muta[maxvar];//save the mutation points in the encoding
	for (j = 0; j < popsize; j++)
	{
		low_T = 0;
		sum = 0;
		sum_muta = 0;

		no_muta = 0;
		for (i = 0; i < nvar; i++)
			pos_muta[i] = -1;

		q = -1;

		for (i = 0; i < nvar; i++)
		{
			rnd1 = randomperc();

			/*For each variable find whether to do mutation or not*/
			if (rnd1 <= pmut_r)
			{
				pos_muta[no_muta] = i;//save the positions of mutation genes
				no_muta++;

				y = new_pop_ptr->ind[j].xreal[i];
				yl = lim_r[i][0];
				yu = lim_r[i][1];


				/*Calculate delta*/

				if ((y - yl) < (yu - y))
					delta = (y - yl) / (yu - yl);
				else
					delta = (yu - y) / (yu - yl);

				rnd3 = randomperc();

				indi = 1.0 / (dim + 1.0);

				if (rnd3 <= 0.5)
				{
					xy = 1.0 - delta;
					val = 2 * rnd3 + (1 - 2 * rnd3) * (pow(xy, (dim + 1)));
					deltaq = pow(val, indi) - 1.0;
				}
				else
				{
					xy = 1.0 - delta;
					val = 2.0 * (1.0 - rnd3) + 2.0 * (rnd3 - 0.5) * (pow(xy, (dim + 1)));
					deltaq = 1.0 - (pow(val, indi));
				}

				/*Change the value for the parent */
				chld = y + deltaq * (yu - yl);

				/*Repair encoding*/
				if (chld < lim_r[i][0])
				{
					rnd2 = randomperc();
					chld = rnd2 * (y - lim_r[i][0]) + lim_r[i][0];

				}

				if (chld > lim_r[i][1])
				{
					rnd2 = randomperc();
					chld = rnd2 * (lim_r[i][1] - y) + y;
				}

				new_pop_ptr->ind[j].xreal[i] = chld;

				sum_muta = sum_muta + new_pop_ptr->ind[j].xreal[i];//the sum of variables in the mutation points
			}//end of mutation
			sum = sum + new_pop_ptr->ind[j].xreal[i];//the sum of all variables
		}//end of i

		/*Repair encoding*/

		if (sum > time_stage)
		{

			for (k = 0; k < no_muta; k++)//Only repair the mutation genes to maintain and inherit the parent's information
			{
				q = pos_muta[k];
				new_pop_ptr->ind[j].xreal[q] = new_pop_ptr->ind[j].xreal[q] * ((time_stage - sum + sum_muta) / sum_muta);
			}

		}

	}

}


void keepalive(population* pop1_ptr, population* pop2_ptr, population* pop3_ptr, int gen)
{
	int i, j, jj, k, m, l, rec;

	int st, pool, poolf, sel;

	double* gene3_ptr, * gene4_ptr;

	/*Forming the global mating pool*/

	for (i = 0; i < popsize; i++)
	{
		if (nvar > 0)
		{	/*For Real Coded GA x values are copied */
			for (k = 0; k < nvar; k++)
			{
				globalpop.xreal[i][k] = pop1_ptr->ind[i].xreal[k];
				globalpop.xreal[i + popsize][k] = pop2_ptr->ind[i].xreal[k];
			}
		}

		/*Fitness is copied to the global pool */
		for (l = 0; l < nfunc; l++)
		{
			globalpop.fitness[i][l] = pop1_ptr->ind[i].fitness[l];
			globalpop.fitness[i + popsize][l] = pop2_ptr->ind[i].fitness[l];
		}

		/*Initial;ising the dummyfitness to zero */
		globalpop.cub_len[i] = 0;
		globalpop.cub_len[i + popsize] = 0;
		globalpop.error[i] = pop1_ptr->ind[i].error;
		globalpop.error[i + popsize] = pop2_ptr->ind[i].error;
		for (jj = 0; jj < ncons; jj++)
		{
			globalpop.constr[i][jj] = pop1_ptr->ind[i].constr[jj];
			globalpop.constr[i + popsize][jj] = pop2_ptr->ind[i].constr[jj];
		}
	}


	global_pop_ptr = &(globalpop);

	/*Finding the global ranks */
	grankc(gen);

	m = globalpop.maxrank;

	/* Sharing the fitness to get the dummy fitness */
	for (i = 0; i < m; i++)
	{
		gshare(i + 1);
	}

	poolf = popsize;
	pool = 0;


	/*Initializing the flags of population to zero */
	for (i = 0; i < 2 * popsize; i++)
	{
		globalpop.flag[i] = 0;
	}
	// decide which all solutions belong to the pop3 
	rec = 0;
	st = 0;
	for (i = 0; i < m; i++)
	{
		/*    Elitism Applied Here     */
		st = pool;
		pool += globalpop.rankno[i];

		if (pool <= popsize)
		{
			for (k = 0; k < 2 * popsize; k++)
			{
				if (globalpop.rank[k] == i + 1)
					globalpop.flag[k] = 1;
			}
			pop3_ptr->rankno[i] = globalpop.rankno[i];
		}
		else
		{
			sel = popsize - st;
			Lastrank = i + 1;
			pop3_ptr->rankno[i] = sel;
			gsort(i + 1, sel);
			break;
		}
	}

	k = 0;
	for (i = 0, k = 0; i < 2 * popsize && k < popsize; i++)
	{
		if (nvar > 0)
		{
			if (globalpop.flag[i] == 1)
			{
				gene3_ptr = &(globalpop.xreal[i][0]);
				pop3_ptr->ind_ptr = &(pop3_ptr->ind[k]);
				gene4_ptr = &(pop3_ptr->ind_ptr->xreal[0]);

				for (j = 0; j < nvar; j++)
				{
					*gene4_ptr++ = *gene3_ptr++;
				}
			}
		}
		if (globalpop.flag[i] == 1)
		{
			for (j = 0; j < nfunc; j++)
				pop3_ptr->ind[k].fitness[j] = globalpop.fitness[i][j];
			pop3_ptr->ind[k].cub_len = globalpop.cub_len[i];
			if (ncons != 0)
				pop3_ptr->ind[k].error = globalpop.error[i];
			for (jj = 0; jj < ncons; jj++)
				pop3_ptr->ind[k].constr[jj] = globalpop.constr[i][jj];
			pop3_ptr->ind[k].rank = globalpop.rank[i];
			k++;  // increment the pop3 counter
		}
	}

	pop3_ptr->maxrank = Lastrank;

}

void grankc(int gen)
{
	int i, j, k, rnk, val, nondom, popsize1, gflg[2 * maxpop], q;
	double* ptr1, * ptr2;
	double* err_ptr1, * err_ptr2;

	/*----------------------------* RANKING *---------------------------------*/
	rnk = 0;
	nondom = 0;
	popsize1 = 2 * popsize;
	min_fit = popsize1;
	delta_fit = 0.1 * popsize1;
	for (i = 0; i < popsize1; i++)
	{
		gflg[i] = 2;
	}
	for (k = 0; k < popsize1; k++)
	{
		q = 0;
		for (j = 0; j < popsize1; j++)
		{
			if (gflg[j] != 1) break;
		}
		if (j == popsize1) break;
		rnk = rnk + 1;
		for (j = 0; j < popsize1; j++)
		{
			if (gflg[j] == 0) gflg[j] = 2;
		}
		for (i = 0; i < popsize1; i++)
		{
			if (gflg[i] != 1 && gflg[i] != 0)
			{
				ptr1 = &(global_pop_ptr->fitness[i][0]);
				err_ptr1 = &(global_pop_ptr->error[i]);
				for (j = 0; j < popsize1; j++)
				{
					if (i != j)
					{
						if (gflg[j] != 1)
						{
							ptr2 = &(global_pop_ptr->fitness[j][0]);
							err_ptr2 = &(global_pop_ptr->error[j]);

							if (*err_ptr1 < 1.0e-6 && *err_ptr2 > 1.0e-6)
							{/* first feasible second individual is infeasible*/
								gflg[j] = 0;
							}
							else
							{
								if (*err_ptr1 > 1.0e-6 && *err_ptr2 < 1.0e-6)
								{/*first individual is infeasible and second is feasible*/
									gflg[i] = 0;
									break;
								}
								else
								{/*both feasible or both infeasible*/
									if (*err_ptr1 > * err_ptr2)
									{
										gflg[i] = 0;
										/*first individual is more infeasible*/
										break;
									}
									else
									{
										if (*err_ptr1 < *err_ptr2)
											gflg[j] = 0;
										/*second individual is more infeasible*/

										else
										{
											val = indcmp1(ptr1, ptr2);
											if (val == 2)
											{
												gflg[i] = 0;
												/* individual 1 is dominated */
												break;
											}
											if (val == 1)
											{
												gflg[j] = 0;
												/* individual 2 is dominated */
											}
											if (val == 3)
											{
												nondom++;/* individual 1 & 2 are non dominated */
												if (gflg[j] != 0) gflg[j] = 3;
											}
										}
									}
								}
							}
						}
					}
				}
				if (j == popsize1)
				{
					global_pop_ptr->rank[i] = rnk;
					gflg[i] = 1;
					global_pop_ptr->rankar[rnk - 1][q] = i;
					q++;
				}
			}
		}
		global_pop_ptr->rankno[rnk - 1] = q;
	}
	global_pop_ptr->maxrank = rnk;

}

int indcmp1(double* ptr1, double* ptr2)
{
	double fit1[maxfun], fit2[maxfun];
	int i, value, m, n;
	for (i = 0; i < nfunc; i++)
	{
		fit1[i] = *ptr1++;
		fit2[i] = *ptr2++;
	}
	m = 0; n = 0;
	while (m < nfunc && fit1[m] <= fit2[m])
	{
		if ((fit2[m] - fit1[m]) < 1e-7) n++;
		m++;
	}
	if (m == nfunc)
	{
		if (n == nfunc) value = 3;
		else value = 1;                    /*value = 1 for dominating*/
	}
	else
	{
		m = 0; n = 0;
		while (m < nfunc && fit1[m] >= fit2[m])
		{
			if ((fit1[m] - fit2[m]) < 1e-7) n++;
			m++;
		}
		if (m == nfunc)
		{
			if (n != nfunc)
				value = 2;                       /*value =  2 for dominated */
			else value = 3;
		}
		else value = 3;                   /*value = 3 for incomparable*/
	}
	return value;
}

void gsort(int rnk, int sel)/* sort the dummy fitness arrays */
{
	int i, j, a, q;
	double array[2 * maxpop][2], temp, temp1;

	q = globalpop.rankno[rnk - 1];

	for (i = 0; i < q; i++)
	{
		array[i][0] = globalpop.rankar[rnk - 1][i];
		a = globalpop.rankar[rnk - 1][i];
		array[i][1] = globalpop.cub_len[a];
	}
	for (i = 0; i < q; i++)
	{
		for (j = i + 1; j < q; j++)
		{
			if (array[i][1] < array[j][1])
			{
				temp = array[i][1];
				temp1 = array[i][0];
				array[i][1] = array[j][1];
				array[i][0] = array[j][0];

				array[j][1] = temp;
				array[j][0] = temp1;
			}
		}
	}

	for (i = 0; i < sel; i++)
	{
		a = (int)array[i][0];
		globalpop.flag[a] = 1;
	}

}

void gshare(int rnk)
{
	double length[2 * maxpop][2], max;
	int i, j, m1, a;
	double min, Diff;  // Added 18.08.2003

	m1 = globalpop.rankno[rnk - 1];

	for (j = 0; j < nfunc; j++)
	{
		for (i = 0; i < m1; i++)
		{
			fpara1[i][0] = 0;
			fpara1[i][1] = 0;
		}

		for (i = 0; i < m1; i++)
		{
			a = globalpop.rankar[rnk - 1][i];
			fpara1[i][0] = (double)a;
			fpara1[i][1] = globalpop.fitness[a][j];
		}

		sort(m1); /*Sort the arrays in ascending order of the fitness*/

		max = fpara1[m1 - 1][1];
		min = fpara1[0][1];  // Added 18.08.2003
		Diff = max - min;      // Added 18.08.2003 and 5 subsequent lines
		if (Diff < 0.0)
		{
			printf("Something wrong in keepaliven.h\n");
			exit(1);
		}
		for (i = 0; i < m1; i++)
		{
			if (i == 0 || i == (m1 - 1))
			{
				length[i][0] = fpara1[i][0];
				length[i][1] = 100 * max;
			}
			else
			{
				length[i][0] = fpara1[i][0];
				length[i][1] = fabs(fpara1[i + 1][1] - fpara1[i - 1][1]) / Diff; // crowding distances are normalized 18.08.2003
			}
		}
		for (i = 0; i < m1; i++)
		{
			a = (int)length[i][0];
			globalpop.cub_len[a] += length[i][1];
		}
	}

}

void sort(int m1)
{
	double temp, temp1;
	int i1, k1;
	for (k1 = 0; k1 < m1 - 1; k1++)
	{
		for (i1 = k1 + 1; i1 < m1; i1++)
		{
			if (fpara1[k1][1] > fpara1[i1][1])
			{
				temp = fpara1[k1][1];
				temp1 = fpara1[k1][0];
				fpara1[k1][1] = fpara1[i1][1];
				fpara1[k1][0] = fpara1[i1][0];
				fpara1[i1][1] = temp;
				fpara1[i1][0] = temp1;
			}
		}
	}

}




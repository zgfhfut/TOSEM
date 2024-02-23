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


/********************************************************************/
#define maxpop   300  /*Max population */

#define MISSION_TIME  200.0
//#define MISSION_TIME  150.0  //complex system
//#define MISSION_TIME  200.0  //large system

/*Max value of testing time*/
#define T_MAX_RESOURCE  150000.0
//#define T_MAX_RESOURCE  30000.0  //complex system
//#define T_MAX_RESOURCE  70000.0  //large system

/*Max no. of variables*/
//#define maxvar    10
//#define maxvar    30 //complex system
#define maxvar    100 ////large system

#define MAX_R  0.9//the reliability contraint

//#define Cost_4  250000.0//the potential cost caused by software failure for large system

#define EPSION  0.000001
#define maxfun    3  /*Max no. of functions */
#define maxcons   3  /*Max no. of Constraints*/
/********************************************************************/
/********************************************************************/
#define MAX_EXPERIMENT  30  /*Experimental times*/
#define MAX_EVALUATION 90000/*Max no. of evaluation*/

/*****************************************a_i,b_i,c_1^i,c_2^i,c_3^i,\gama_i***************************/
double para_component[maxvar][5] = { 0 };/*parameters for component*/
double  C;
double mu[maxvar] = { 0 }; //Expected number of visits

double theta[maxvar] = { 0 }; //Expected number of visits
/********************************************************************/
double lim_r[maxvar][2] = { 0 };/*Limits of variable in array*/
double Y[maxvar];
double tmid[maxvar];//测试资源中间值；
const  double EPSINON = 0.0000001;    //minimal precision

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
void Cal_C();
void realinit(population* pop_ptr);/*initializes the population*/

double cal_reliab(int ind_module, double m_res);//calculate the reliability of each module
double cal_cost(int ind_module, double m_res);//calculate the cost of each module
void func(population* pop_ptr);/*evaluate the value of the function & errors*/

void realcross(population* new_pop_ptr, population* mate_pop_ptr);/*simulated binary crossover for Real Coded GA*/
void real_mutate(population* old_pop_ptr, population* new_pop_ptr);/*formulate the mutation routine*/

void checking(population* pop_ptr);

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
			fscanf_s(end_ptr1, "%lf %lf %lf %lf %lf", &para_component[i][0], &para_component[i][1], &para_component[i][2], &para_component[i][3], &para_component[i][4]); // 循环读
		fclose(end_ptr1); //关闭文件


		FILE* end_ptr2;
		char file2[500];
		sprintf_s(file2, "C:\\Users\\Administrator\\Desktop\\论文代码数据\\Parameters\\multi-input and multi-output\\WNS-MOEAD\\预期访问次数\\100\\%d.txt", pp);
		fopen_s(&end_ptr2, file2, "r+");
		for (int i = 0; i < maxvar; i++)
			fscanf_s(end_ptr2, "%lf", &mu[i]); // 循环读
		fclose(end_ptr2); //关闭文件

		int int_exp = 0;
		input();
		Cal_C();
		for (int_exp = 0; int_exp < MAX_EXPERIMENT; int_exp++)
		{
			printf_s("Run experiment %d!\n", int_exp + 1);

			eva_num = 0;


			FILE* end_ptr;
			char file4[500];

			sprintf_s(file4, "C:\\Users\\Administrator\\Desktop\\论文代码数据\\experimental data\\multi-input and multi-output\\WNS-MOEAD\\100\\Data of 30 instances\\test%d\\fitness%d.txt", pp, int_exp + 1);//在file1中写上每次迭代产生的不同fitness文件的路径
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
			checking(old_pop_ptr);
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
				checking(mate_pop_ptr);

				old_pop_ptr = &(oldpop);
				realcross(old_pop_ptr, mate_pop_ptr);/*Real Cross-over*/
				/********************************************************************/
				mate_pop_ptr = &(matepop);
				checking(mate_pop_ptr);
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
						fprintf(end_ptr, "%5.6f\t", archive_pop_ptr->ind_ptr->fitness[l]);

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
	double sum_u = 0.0;/*所有组件预期访问次数之和*/

	//计算总预期访问次数
	for (i = 0; i < nvar; i++)
	{
		sum_u += mu[i];
	}

	//计算单个组件任务时间
	for (i = 0; i < nvar; i++)
	{
		theta[i] = MISSION_TIME * mu[i] / sum_u;
	}

	for (i = 0; i < nvar; i++)
	{
		/*Specify the limits of the variables*/
		lim_r[i][0] = 0;
		lim_r[i][1] = T_MAX_RESOURCE;
	}


}

void realinit(population* pop_ptr)
{
	int i, j;
	double d;

	for (i = 0; i < popsize; i++)
	{
		for (j = 0; j < nvar; j++)
		{
			d = randomperc();
			/*limits are specified it generates the value in range of minimum and maximum value of the variable*/
			pop_ptr->ind[i].xreal[j] = d * (lim_r[j][1] - lim_r[j][0]) + lim_r[j][0];

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

			sum += para_component[i][4] * theta[i] / (para_component[i][2] - para_component[i][3] - mid * para_component[i][1] * theta[i]);

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
		Y[j] = para_component[j][4] * theta[j] / (para_component[j][2] - para_component[j][3] - mid * para_component[j][1] * theta[j]);
		tmid[j] = -1 / para_component[j][1] * log(-Y[j] / (para_component[j][0] * para_component[j][1] * theta[j]));
		c[j] = (para_component[j][2] - para_component[j][3]) * para_component[j][0] * (1 - exp(-para_component[j][1] * tmid[j])) + para_component[j][3] * para_component[j][0] + para_component[j][4] * tmid[j];
		C_min += c[j];

	}
	C = C_min * 1.20;

}

double cal_reliab(int ind_module, double m_res) {

	double rel = 0, rel1 = 0;
	rel = para_component[ind_module][0] * para_component[ind_module][1] * exp(-para_component[ind_module][1] * m_res);
	rel1 = exp(-rel * theta[ind_module]);
	return rel1;
}

double cal_cost(int ind_module, double m_res) {

	double cost = 0;
	/*cost = para_module[ind_module][2] * para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res)) +
		para_module[ind_module][3] * (para_module[ind_module][0] - para_module[ind_module][0] * (1 - exp(-para_module[ind_module][1] * m_res))) +
		para_module[ind_module][4] * pow(m_res, para_module[ind_module][5]);*/
	cost = (para_component[ind_module][2] - para_component[ind_module][3]) * para_component[ind_module][0] * (1 - exp(-para_component[ind_module][1] * m_res)) + para_component[ind_module][3] * para_component[ind_module][0] + para_component[ind_module][4] * m_res;
	return cost;
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

		if ((f[2] - T_MAX_RESOURCE) > EPSINON)
		{
			printf_s("error in the upper bound constraint11\n");
			system("pause");
		}


		/*=========End Your Coding Upto This Point===============*/

		/******************************************************************/
		/*              Put The Constraints Here                          */
		/******************************************************************/
		cstr[0] = sum1 - MAX_R;
		cstr[1] = (T_MAX_RESOURCE - sum2)/T_MAX_RESOURCE;
		cstr[2] = (C-sum3)/C;
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
	double rnd2;

	for (i = 0; i < popsize; i++)
	{
		for (j = 0; j < nvar; j++)
		{
			rnd2 = randomperc();
			if (rnd2 < pcross)
			{
				mate_pop_ptr->ind[i].xreal[j] = mate_pop_ptr->ind[i].xreal[j];
			}
			else {
				mate_pop_ptr->ind[i].xreal[j] = new_pop_ptr->ind[i].xreal[j];
			}

		}

	}
}

void real_mutate(population* old_pop_ptr, population* new_pop_ptr)
{
	int p1 = 0, p2 = 0, p3 = 0;
	int i;
	int j;

	for (i = 0; i < popsize; i++)
	{
		do {
			p1 = rnd(0, popsize - 1);
			p2 = rnd(0, popsize - 1);
			p3 = rnd(0, popsize - 1);

		} while (i == p1 || i == p2 || i == p3 || p1 == p2 || p1 == p3 || p2 == p3);

		for (j = 0; j < nvar; j++)
		{
			new_pop_ptr->ind[i].xreal[j] = old_pop_ptr->ind[p1].xreal[j] + pmut_r * (old_pop_ptr->ind[p2].xreal[j] - old_pop_ptr->ind[p3].xreal[j]);

			if (new_pop_ptr->ind[i].xreal[j] < lim_r[j][0])
			{
				new_pop_ptr->ind[i].xreal[j] = lim_r[j][0];
			}

			if (new_pop_ptr->ind[i].xreal[j] > lim_r[j][1])
			{
				new_pop_ptr->ind[i].xreal[j] = lim_r[j][1];
			}

		}

	}

}


void checking(population* pop_ptr)
{
	double sum;
	int i;
	int j;
	double rnd3;

	for (i = 0; i < popsize; i++)
	{
		sum = 0;
		for (j = 0; j < nvar; j++)
		{
			if (pop_ptr->ind[i].xreal[j] < 0.0)
			{
				printf("Something wrong in xreal[j]\n");
				exit(1);

			}
			sum += pop_ptr->ind[i].xreal[j];
		}

		if (sum > T_MAX_RESOURCE)
		{

			for (j = 0; j < nvar; j++)
			{
				rnd3 = randomperc();
				pop_ptr->ind[i].xreal[j] = pop_ptr->ind[i].xreal[j] * rnd3 * (T_MAX_RESOURCE / sum);
			}

		}
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


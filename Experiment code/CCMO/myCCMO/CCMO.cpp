#include "stdafx.h"
#include<iostream>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include<iostream>
# include "global.h" 
# include "Random.h"
using namespace std;
#define MAX_R 0.9
#define MAX_NUM_EXP 30
double para_module[nreal][5];
double u[nreal];
double reli_thresh;
extern double lim_r[nreal][2];
int _tmain(int argc, _TCHAR* argv[])
{
	for (int o =1; o < 31; o++)
	{
		
	FILE* end_ptr2; FILE* end_ptr1;
	char file3[500]; char file2[500];

	int i;
	sprintf_s(file3, "C:\\Users\\shao\\Desktop\\论文代码数据\\Parameters\\single-input and single-output\\CCMO\\组件参数\\100\\%d.txt", o);
	
	fopen_s(&end_ptr2, file3, "r+");
	//在file3中写上每次迭代产生的不同fitness文件的路径

	for (i = 0; i < nreal; i++)
	{
	
		fscanf_s(end_ptr2, "%lf %lf %lf %lf %lf", &para_module[i][0], &para_module[i][1], &para_module[i][2], &para_module[i][3], &para_module[i][4]);// 循环读

	}

	fclose(end_ptr2); //关闭文件
	sprintf_s(file2, "C:\\Users\\shao\\Desktop\\论文代码数据\\Parameters\\single-input and single-output\\CCMO\\预期访问次数\\100\\%d.txt", o);

	fopen_s(&end_ptr1, file2, "r+");
	for (i = 0; i < nreal; i++)
		fscanf_s(end_ptr1, "%lf", &u[i]); // 循环读
	fclose(end_ptr1); //关闭文件

	
	
	    Cal_res_t();
		for (int int_exp = 1; int_exp <= MAX_NUM_EXP; int_exp++)
		{
			printf_s("Run experiment %d!\n", int_exp);

			FILE* end_ptr;
			char file1[500];

		//	sprintf_s(file1, "E:\\new_CCMO\\fitness\\100\\2\\R=0.95\\fitness%d.txt", int_exp);
			sprintf_s(file1, "C:\\Users\\shao\\Desktop\\论文代码数据\\experimental data\\single-input and single-output\\myCCMO\\100\\Data of 30 instances\\test%d\\fitness%d.txt", o, int_exp);

			fopen_s(&end_ptr, file1, "w+");

			int i;

			//	int t = 0;

			int random;


			population* pop_set;

			population* child1_pop;

			population* archive_set;

			population* child2_pop;

			population* mixed1_set;

			population* mixed2_set;

			srand((unsigned)time(NULL));
			random = rand() % 1000;

			seed = (float)random / 1000.0;
			while (seed <= 0.0 || seed >= 1.0)
			{
				srand((unsigned)time(NULL));
				random = rand() % 1000;

				seed = (float)random / 1000.0;
				
			}


			/*-------test problem setting------start---------*/
			// for the moersa instance
			nobj = 3;						// number of objectives

			ncon = 3;			            // number of constraints

			// set the range of the variables
			// set the range of the variables
			min_realvar = (double*)malloc(nreal * sizeof(double));
			max_realvar = (double*)malloc(nreal * sizeof(double));
			for (i = 0; i < nreal; i++)
			{
				min_realvar[i] = lim_r[i][0];
				max_realvar[i] = lim_r[i][1];
			}
			/*-------test problem setting------end---------*/
			//大规模变量，1000个变量
			/*-------parameter setting in the algorithm------start---------*/
			popsize = 304;					// population size
			archsize = 304;					// archive size

			neval = 90000;

			pcross_real = 0.9;				// crossover probability
			pmut_real = 1.0 / nreal;			// mutation probability
			eta_c = 20.0;					// distribution index of simulated binary crossover (SBX)
			eta_m = 20.0;					// distribution index of polynomial mutation
			/*-------parameter setting in the algorithm------end---------*/

			currenteval = 0;

			pop_set = (population*)malloc(sizeof(population));
			child1_pop = (population*)malloc(sizeof(population));

			archive_set = (population*)malloc(sizeof(population));
			child2_pop = (population*)malloc(sizeof(population));

			mixed1_set = (population*)malloc(sizeof(population));
			mixed2_set = (population*)malloc(sizeof(population));

			allocate_memory_pop(pop_set, popsize);

			allocate_memory_pop(child1_pop, popsize / 2);

			allocate_memory_pop(archive_set, archsize);
			allocate_memory_pop(child2_pop, archsize / 2);


			allocate_memory_pop(mixed1_set, popsize + archsize);
			allocate_memory_pop(mixed2_set, popsize + archsize);

			randomize();

			initialize_pop(pop_set, popsize);
			evaluate_pop(pop_set, popsize);
			assign_fitness(pop_set, popsize);

			initialize_pop(archive_set, archsize);
			evaluate2_pop(archive_set, archsize);
			assign_fitness(archive_set, archsize);


			do
			{
				mating_selection(pop_set, child1_pop);
				mutation_pop(child1_pop);

				mating_selection(archive_set, child2_pop);
				mutation_pop(child2_pop);

				merge(pop_set, child1_pop, child2_pop, mixed1_set, popsize, popsize / 2, popsize / 2);
				merge(archive_set, child1_pop, child2_pop, mixed2_set, archsize, archsize / 2, archsize / 2);

				evaluate_pop(mixed1_set, popsize * 2);
				assign_fitness(mixed1_set, popsize * 2);
				environmental_selection(mixed1_set, pop_set);	// environmental selection

				evaluate2_pop(mixed2_set, archsize * 2);
				assign_fitness(mixed2_set, archsize * 2);
				environmental_selection(mixed2_set, archive_set);	// environmental selection


			} while (currenteval < neval);



			report_objective(pop_set, popsize, end_ptr);
			fflush(stdout);
			fflush(end_ptr);
			fclose(end_ptr);


			free(min_realvar);
			free(max_realvar);

			deallocate_memory_pop(pop_set, popsize);
			deallocate_memory_pop(archive_set, archsize);
			deallocate_memory_pop(child1_pop, popsize / 2);
			deallocate_memory_pop(child2_pop, archsize / 2);
			deallocate_memory_pop(mixed1_set, popsize + archsize);
			deallocate_memory_pop(mixed2_set, popsize + archsize);

			free(pop_set);
			free(archive_set);
			free(child1_pop);
			free(child2_pop);
			free(mixed1_set);
			free(mixed2_set);

		}

		printf("\n Routine successfully exited \a\n");

		//system("pause");
	}
	return 0;
}


// MOEAD-IEpsilon.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

#include <iostream>
#include <algorithm>

#include <vector>

using namespace std;

# include "global.h" 
# include "Random.h"

#define MAX_R 0.9

#define MAX_NUM_EXP 30
double para_module[nreal][5];
double u[nreal];

int _tmain(int argc, _TCHAR* argv[])
{
	for (int o =1; o <31; o++)
	{
		reli_thresh = MAX_R;

	FILE* end_ptr2; FILE* end_ptr1;
	char file3[500]; char file2[500];

	int i;

	sprintf_s(file3, "C:\\Users\\shao\\Desktop\\论文代码数据\\Parameters\\multi-input and multi-output\\MOEAD-IEpsilon\\组件参数\\50\\%d.txt", o);

	fopen_s(&end_ptr2, file3, "r+");
	//在file3中写上每次迭代产生的不同fitness文件的路径

	for (i = 0; i < nreal; i++)
		fscanf_s(end_ptr2, "%lf %lf %lf %lf %lf", &para_module[i][0], &para_module[i][1], &para_module[i][2], &para_module[i][3], &para_module[i][4]); // 循环读


	fclose(end_ptr2); //关闭文件

	sprintf_s(file2, "C:\\Users\\shao\\Desktop\\论文代码数据\\Parameters\\multi-input and multi-output\\MOEAD-IEpsilon\\预期访问次数\\50\\%d.txt", o);
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

			//sprintf_s(file1, "D:\\MOEAD-IEpsilon\\fitness\\100\\3\\R=0.95\\fitness%d.txt", int_exp);
			sprintf_s(file1, "C:\\Users\\shao\\Desktop\\论文代码数据\\experimental data\\multi-input and multi-output\\myMOEAD-IEpsilon+\\50\\Data of 30 instances\\test%d\\fitness%d.txt", o, int_exp);

			fopen_s(&end_ptr, file1, "w+");

			int i;


			int random;


			population* pop_set;
			population* archive_set;
			population* mixed_set;

			srand((unsigned)time(NULL));
			random = rand() % 1000;

			seed = (float)random / 1000.0;
			while (seed <= 0.0 || seed >= 1.0)
			{
				int random = rand() % 1000;

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
			//
			/*-------parameter setting in the algorithm------start---------*/
			popsize = 300;					// population size
			archsize = 300;					// archive size

			niter = 50; // maximal number of iterations

			//neval = 15000;
			neval =90000;
			currenteval = 0;

			pcross_real = 0.9;				// crossover probability
			pmut_real = 1.0 / nreal;			// mutation probability		
			eta_m = 20.0;					// distribution index of polynomial mutation

			phi_max_ = -1e+30;
			T_ = 30;
			delta_ = 0.9;
			nr_ = 2;

			F_ = 0.5;

			neighborhood_ = (int**)malloc(popsize * sizeof(int*));
			for (i = 0; i < popsize; i++)
			{
				neighborhood_[i] = (int*)malloc(T_ * sizeof(int));
			}

			z_ = (double*)malloc(nobj * sizeof(double));

			lambda_ = (double**)malloc(popsize * sizeof(double*));
			for (i = 0; i < popsize; i++)
			{
				lambda_[i] = (double*)malloc(nobj * sizeof(double));
			}

			pop_set = (population*)malloc(sizeof(population));
			archive_set = (population*)malloc(sizeof(population));
			mixed_set = (population*)malloc(sizeof(population));

			allocate_memory_pop(pop_set, popsize);
			allocate_memory_pop(archive_set, archsize);
			allocate_memory_pop(mixed_set, popsize + archsize);

			randomize();

			initUniformWeight();

			initNeighborhood();

			initialize_pop(pop_set, popsize);
			evaluate_pop(pop_set, popsize);

			initialize_pop(archive_set, archsize);
			evaluate_pop(archive_set, archsize);

			// Initialize the epsilon_zero_
			double* constraints = (double*)malloc(popsize * sizeof(double));
			for (int i = 0; i < popsize; i++) {
				constraints[i] = pop_set->ind[i].constr_violation;
				//	  cout << "constraints[i]=" << constraints[i] << endl;
			}

			sort(constraints, constraints + popsize); // each constraints is less or equal than zero;

			double epsilon_zero_ = abs(constraints[(int)ceil(0.05 * popsize)]);


			if (phi_max_ < abs(constraints[0])) {
				phi_max_ = abs(constraints[0]);
			}
			int tc_ = (int)(0.8 * niter / popsize);
			double r_k_ = GetFeasible_Ratio(pop_set, popsize);
			double tao_ = 0.05;

			//cout << "r_k_=" << r_k_ << endl;

			// STEP 1.3. Initialize z_
			initIdealPoint(pop_set);

			int gen = 0;
			epsilon_k_ = epsilon_zero_;

			do
			{
				printf("1");
				// update the epsilon level
				if (gen >= tc_) {
					epsilon_k_ = 0;
				}
				else {
					if (r_k_ < 0.95) {
						epsilon_k_ = (1 - tao_) * epsilon_k_;
					}
					else {
						epsilon_k_ = phi_max_ * (1 + tao_);
					}
				}

				int* permutation = (int*)malloc(popsize * sizeof(int));
				randomPermutation(permutation, popsize);
				for (int i = 0; i < popsize; i++) {
					int n = permutation[i];
					int type;
					double rnd = randomperc();

					// STEP 2.1. Mating selection based on probability
					if (rnd < delta_) // if (rnd < realb)
					{
						type = 1;   // neighborhood
					}
					else {
						type = 2;   // whole population
					}
					vector<int> p;
					matingSelection(p, n, 2, type);

					// STEP 2.2. Reproduction
					individual* child = (individual*)malloc(sizeof(individual));
					allocate_memory_ind(child);

					individual* parents = (individual*)malloc(4 * sizeof(individual));

					parents[0] = pop_set->ind[p[0]];
					parents[1] = pop_set->ind[p[1]];
					parents[2] = pop_set->ind[n];
					parents[3] = pop_set->ind[n];
					// Apply DE crossover			
					DECrossover(child, parents);

					real_mutate_ind(child);

					evaluate_ind(child);
					currenteval++;

					//update phi_max_
					if (phi_max_ < abs(child->constr_violation)) {
						phi_max_ = abs(child->constr_violation);
					}

					// STEP 2.4. Update z_
					updateReference(child);

					// STEP 2.5. Update of solutions
					updateProblem(pop_set, child, n, type);

					free(parents);

					deallocate_memory_ind(child);
					free(child);
				}

				r_k_ = GetFeasible_Ratio(pop_set, popsize);
				//	cout << "r_k_=" << r_k_ << endl;

				merge(pop_set, archive_set, mixed_set, popsize, archsize);

				fill_nondominated_sort(mixed_set, archive_set);

				//cout << "gen=" << gen << endl;

				gen++;

				free(permutation);

			} while (currenteval < neval);

			report_objective(archive_set, archsize, end_ptr);
			fflush(stdout);
			fflush(end_ptr);

			fclose(end_ptr);


			free(min_realvar);
			free(max_realvar);

			deallocate_memory_pop(pop_set, popsize);
			deallocate_memory_pop(archive_set, archsize);
			deallocate_memory_pop(mixed_set, popsize + archsize);

			free(pop_set);
			free(archive_set);
			free(mixed_set);


			free(z_);
			for (i = 0; i < popsize; i++)
			{
				free(neighborhood_[i]);
			}
			free(neighborhood_);

			for (i = 0; i < popsize; i++)
			{
				free(lambda_[i]);
			}
			free(lambda_);

			free(constraints);


		}

		printf("\n Routine successfully exited \a\n");
	}
	//system("pause");
	return 0;
}


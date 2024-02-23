/* archive truncation routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

// truncate the archive set based on the k-th nearest neighbor method
void truncation (population *pop, int *front, int front_size, population *archive)
{
	int i, j;
	double **distSort;
	int **distIndex, *record;
	double min;
	int mark, count, current_size;

	distSort = (double**)malloc(front_size*sizeof(double*));
	distIndex = (int**)malloc(front_size*sizeof(int*));
	record = (int*)malloc(front_size*sizeof(int));
	for (i=0; i<front_size; i++)
	{
		distSort[i] = (double*)malloc(front_size*sizeof(double));
		distIndex[i] = (int*)malloc(front_size*sizeof(int));
	}
	for (i=0; i<front_size; i++)
	{
		for (j=0; j<front_size; j++)
		{
			if (i==j)
				distSort[i][j] = INF;
			else
				distSort[i][j] = Euclidean_Distance(&(pop->ind[front[i]]), &(pop->ind[front[j]]));
				
			distIndex[i][j] = j;
		}
	}

	for (i=0; i<front_size; i++)
	{
		q_sort_distance_index (distSort[i], distIndex[i], 0, front_size-1);
	}

	for (i=0; i<front_size; i++)
		record[i] = 1; 

	current_size = front_size;
	while (current_size > archsize)
	{
		min = INF;
		for (i=0; i<front_size; i++)
		{
			if (record[i] == 1)
			{
				if (min > distSort[i][0])
				{
					min = distSort[i][0];
					mark = i;
				}
				else
				{
					if (min == distSort[i][0])
					{
						if (check_repeat(&pop->ind[front[mark]], &pop->ind[front[i]]) == 1)
						{
							count = 0;
							while (distSort[mark][count] == distSort[i][count] && count < current_size-1)
								count++;
							if (distSort[mark][count] > distSort[i][count])
								mark = i;
						}
					}
				}
			}
		}
		record[mark] = 0;
		for (i=0; i<front_size; i++)
		{
			if (record[i] == 1)
			{
				for (j=0; j<current_size-1; j++)
				{
					if (distIndex[i][j] == mark)
					{
						while (j != current_size - 2)
						{
							distSort[i][j] = distSort[i][j+1];
							distIndex[i][j] = distIndex[i][j+1];
							j++;
						}
						distSort[i][current_size-2] = INF;
						distIndex[i][current_size-2] = mark;
						break;
					}
				}
			}
		}
		current_size -= 1;
	}

	j = 0;
	for (i=0; i<front_size; i++)
	{
		if (record[i] == 1)
		{
			copy_ind(&pop->ind[front[i]], &archive->ind[j]);
			j++;
		}
	}

	free(record);
	for (i=0; i<front_size; i++)
	{
		free(distSort[i]);
		free(distIndex[i]);
	}
	free(distSort);
	free(distIndex);

	return;
}

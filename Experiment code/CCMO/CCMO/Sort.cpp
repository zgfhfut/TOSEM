/* Routines for finding the distance between the considered individual and its k-th nearest neighbor */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

double findKmin (double *a, int Kmin, int size)
{
	q_sort_distance(a, 0, size-1);
	return a[Kmin-1];
}

/* Randomized quick sort routine to sort an array */
void q_sort_distance (double *a, int left, int right)
{
	int index;
    double temp;
    int i, j;
    double pivot;

    if (left < right)
    {
        index = rnd (left, right);
		temp = a[right];
        a[right] = a[index];
        a[index] = temp;
        pivot = a[right];
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (a[j] <= pivot)
            {
                i += 1;
                temp = a[j];
                a[j] = a[i];
                a[i] = temp;
            }
        }
        index = i+1;
        temp = a[index];
        a[index] = a[right];
        a[right] = temp;
        q_sort_distance (a, left, index-1);
        q_sort_distance (a, index+1, right);
    }
    return;
}


void q_sort_distance_index (double *a, int *b, int left, int right)
{
	int index;
    double temp_a;
	int temp_b;
    int i, j;
    double pivot;

    if (left < right)
    {
        index = rnd (left, right);
		temp_a = a[right];
        a[right] = a[index];
        a[index] = temp_a;
		temp_b = b[right];
		b[right] = b[index];
		b[index] = temp_b;
        pivot = a[right];
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (a[j] <= pivot)
            {
                i += 1;
                temp_a = a[j];
                a[j] = a[i];
                a[i] = temp_a;
				temp_b = b[j];
				b[j] = b[i];
				b[i] = temp_b;
            }
        }
        index = i+1;
        temp_a = a[index];
        a[index] = a[right];
        a[right] = temp_a;
		temp_b = b[index];
		b[index] = b[right];
		b[right] = temp_b;
        q_sort_distance_index (a, b, left, index-1);
        q_sort_distance_index (a, b, index+1, right);
    }
    return;
}
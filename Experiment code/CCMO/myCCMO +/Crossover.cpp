/* Crossover routines */
#include "stdafx.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "Global.h"
# include "Random.h"

/* Function to cross two individuals */
void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
	SBX_cross (parent1, parent2, child1, child2);
    return;
}

/* Routine for real variable SBX crossover */
void SBX_cross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{	
	int i;
	double rand;
	double sum_chld1=0, sum_chld2=0,low_T=0;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double lim1[nreal];
	double lim2[nreal];
	int count[nreal];

	if (randomperc() <= pcross_real)
	{
		
		for (int i = 0; i<nreal; i++)
		{
			lim1[i] = 0;
			lim2[i] = 0;
			count[i] = 0;
			if (randomperc() <= 0.5)
			{
				
		
				if (parent1->xreal[i] < parent2->xreal[i])
				{
					y1 = parent1->xreal[i];
					y2 = parent2->xreal[i];
				}
				else
				{
					y1 = parent2->xreal[i];
					y2 = parent1->xreal[i];
				}
				if (fabs(parent1->xreal[i] - parent2->xreal[i]) > EPS)
				{
					count[i] = 1;
					yl = min_realvar[i];
					yu = max_realvar[i];
					rand = randomperc();
					beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand * alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
					}
					c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
					beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand * alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
					}
					c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

					 if (randomperc() <= 0.5)
					{
						child1->xreal[i] = c2;
					    child2->xreal[i] = c1;
						
					}
					else
					{
						child1->xreal[i] = c1;
						child2->xreal[i] = c2;
					
					}
				
				 if (child1->xreal[i]<yl)
					{  
						child1->xreal[i] =repair_s(yl,yu, parent1->xreal[i],child1->xreal[i]);
						
					}
					if (child1->xreal[i]>yu)
					{
						child1->xreal[i] = repair_s(yl, yu, parent1->xreal[i], child1->xreal[i]);
						
					}
					if (child2->xreal[i]<yl)
					{
						child2->xreal[i] =  repair_s(yl, yu, parent2->xreal[i], child2->xreal[i]);
						
					}
					if (child2->xreal[i]>yu)
					{
						child2->xreal[i] = repair_s(yl, yu, parent2->xreal[i], child2->xreal[i]);
				
					}
				}
				else
				{
					child1->xreal[i] = parent1->xreal[i];
				
					child2->xreal[i] = parent2->xreal[i];
					
				}
				if (child1->xreal[i] < parent1->xreal[i]) lim1[i] = min_realvar[i];
				if (child1->xreal[i] > parent1->xreal[i]) lim1[i] = parent1->xreal[i];
			}
			else
			{
				child1->xreal[i] = parent1->xreal[i];
				child2->xreal[i] = parent2->xreal[i];
				
			}
			
		}
 		repair_all(child1,lim1,count);
		repair_all(child2,lim2,count);
	}
	else
	{
		for (i = 0; i<nreal; i++)
		{
			child1->xreal[i] = parent1->xreal[i];
			child2->xreal[i] = parent2->xreal[i];
		}
	}
	
    return;
}

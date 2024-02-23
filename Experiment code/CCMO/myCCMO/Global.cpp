#include "stdafx.h"

# include "Global.h"

const  double EPSINON = 0.000001;    //minimal precision

//i
extern double para_module[nreal][5]; 


int nobj;
int ncon;
int popsize;
int archsize;
double pcross_real;
double pmut_real;
double eta_c;
double eta_m;

int currenteval;
int neval;

int niter;

double *min_realvar;
double *max_realvar;

double seed;
double oldrand[55];
int jrand;

int minvalue(int a, int b)
{
	if (a>=b)
		return b;
	else
		return a;

}

int maxintvalue(int a, int b)
{
	if (a >= b)
		return a;
	else
		return b;

}

double maxvalue(double a, double b)
{
	if (a >= b)
		return a;
	else
		return b;

}



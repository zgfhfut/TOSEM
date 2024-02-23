#pragma once

#define MISSION_TIME  200.0  //workable life
#define nreal   50
#define T_MAX_RESOURCE 70000.0  //simple system
#define EPSION  0.000001
#define Max_R 0.9
extern double C;

extern double para_module[nreal][5];
extern double u[nreal];
extern double Mis_tim[nreal];
extern double lim_r[nreal][2];/*Limits of variable in array*/
extern double lim_t[nreal][2];
extern double M;
extern double T_low_sys;
extern double zz, zz2, zz1;
extern int q1, q2;/*系统测试资源总量下限时单个构件测试资源大于0的个数*/
extern int Z_0[nreal];/*记录*/
extern double Z2[nreal];/*系统测试资源总量下限时单个构件测试资源*/
extern double e[nreal];
extern double Zmax[nreal];
extern double reciprocal_sum;/*故障检测率倒数之和*/
extern double Y[nreal];
extern double tmid[nreal];//测试资源中间值；
extern double c_min[nreal];//每个模块成本的极小值点
extern double t_min[nreal];//每个模块成本的极小值点对应的测试资源
extern double lim_r[nreal][2];
double cal_lowtime(int g);
double cal_hightime(int g);
double cal_reliab(int ind_module, double m_res);
double cal_cost(int ind_module, double m_res);
double cal_z(int  g, double f);//除去第g个组件的z值
double cal_t1(int i);//每个模块对应c最小值所需的测试资源
double cal_t2(int i);//每个模块可靠性达到MAX_R所需要的测试资源
double cal_sumcost(int  g, double f);//求除去g模块的n-1个模块所花费的成本最小值
double cal_sumR(int  g, double f);



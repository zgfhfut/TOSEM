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
extern int q1, q2;/*ϵͳ������Դ��������ʱ��������������Դ����0�ĸ���*/
extern int Z_0[nreal];/*��¼*/
extern double Z2[nreal];/*ϵͳ������Դ��������ʱ��������������Դ*/
extern double e[nreal];
extern double Zmax[nreal];
extern double reciprocal_sum;/*���ϼ���ʵ���֮��*/
extern double Y[nreal];
extern double tmid[nreal];//������Դ�м�ֵ��
extern double c_min[nreal];//ÿ��ģ��ɱ��ļ�Сֵ��
extern double t_min[nreal];//ÿ��ģ��ɱ��ļ�Сֵ���Ӧ�Ĳ�����Դ
extern double lim_r[nreal][2];
double cal_lowtime(int g);
double cal_hightime(int g);
double cal_reliab(int ind_module, double m_res);
double cal_cost(int ind_module, double m_res);
double cal_z(int  g, double f);//��ȥ��g�������zֵ
double cal_t1(int i);//ÿ��ģ���Ӧc��Сֵ����Ĳ�����Դ
double cal_t2(int i);//ÿ��ģ��ɿ��ԴﵽMAX_R����Ҫ�Ĳ�����Դ
double cal_sumcost(int  g, double f);//���ȥgģ���n-1��ģ�������ѵĳɱ���Сֵ
double cal_sumR(int  g, double f);



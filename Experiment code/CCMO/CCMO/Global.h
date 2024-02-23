/* This file contains the variable and function declarations */

# define INF 1.0e14
# define EPS 1.0e-14
# define E  2.71828182845905
# define PI 3.14159265358979

extern const double EPSINON;

#define nreal   100
//#define nreal   30//complex system 
//#define nreal   50//complex system 

//#define reli_thresh   0.95//complex system 

//#define T_MAX_RESOURCE  150000.0  //complex system

//#define T_MAX_RESOURCE  230000.0  //complex system

//#define T_MAX_RESOURCE  230000.0  //complex system

#define T_MAX_RESOURCE  150000.0  //complex system
//#define T_MAX_RESOURCE  30000.0  //complex system
//#define T_MAX_RESOURCE  70000.0  //large system
#define Max_R 0.9
//#define MISSION_TIME  150.0  //workable life
//#define MISSION_TIME  150.0  //workable life
#define MISSION_TIME  200.0  //workable life


/*Max no. of subsystems*/
//#define MAX_SYSTEM   11//complex system 

//#define MAX_SYSTEM   16

//#define MAX_SYSTEM   10

//extern int no_subsys[MAX_SYSTEM];//complex system 
extern int currenteval;
extern double para_module[nreal][5];



/*-------global variables------start---------*/
typedef struct
{
	double *xreal;              // individual's variable value 
    double *obj;				// individual's objective value 

	double *constr;				// individual's constraint
	double constr_violation;	// individual's constraint violation degree

	double fitness;				// individual's fitness

} individual;

typedef struct
{
    individual *ind;
} population;

typedef struct lists
{
	int index;
	struct lists *parent;
	struct lists *child;
} list;

typedef struct ind_lists
{
	double *obj;
	struct ind_lists *parent;
	struct ind_lists *child;
} ind_list;

extern int ncon;
extern int nobj;
extern int popsize;
extern int archsize;
extern double pcross_real;
extern double pmut_real;
extern double eta_c;
extern double eta_m;


extern int niter;

extern int currenteval;
extern int neval;

extern double *min_realvar;
extern double *max_realvar;
/*-------global variables------end---------*/
extern double u[nreal];
void Cal_c();

/*-------global function declarations------start---------*/
void allocate_memory_pop (population *pop, int size);
void allocate_memory_ind (individual *ind);
void deallocate_memory_pop (population *pop, int size);
void deallocate_memory_ind (individual *ind);

void evaluate_pop(population *pop, int size);
void evaluate2_pop(population *pop, int size);

void evaluate_ind(individual *ind);
void evaluate2_ind(individual *ind);

void test_problem(double *xreal, double *obj, double *constr);
void test2_problem(double *xreal, double *obj, double *constr);

void Cal_res_lim(double r);/*用来计算组件任务时间*/
void initialize_pop (population *pop, int size);
void initialize_ind (individual *ind);

void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2);
void SBX_cross (individual *parent1, individual *parent2, individual *child1, individual *child2);

void mating_selection(population *old_pop, population *new_pop);
individual* tournament(individual *ind1, individual *ind2);
int check_dominance(individual *a, individual *b);

void mutation_pop (population *pop);
void real_mutate_ind (individual *ind);

void merge(population *pop1, population *pop2, population *pop3, population *pop4, int size1, int size2, int size3);
void copy_ind(individual *ind1, individual *ind2);

void assign_fitness(population *pop, int size);

double Euclidean_Distance(individual *ind1, individual *ind2);
double findKmin(double *a, int Kmin, int size);
void q_sort_distance(double *a, int left, int right);
void q_sort_distance_index(double *a, int *b, int left, int right);

void truncation(population *pop, int *front, int front_size, population *archive);
int check_repeat(individual *a, individual *b);
void environmental_selection(population *mixed_set, population *archive_set);


void report_objective(population *pop, int size, FILE *fpt);
void report_variable(population *pop, int size, FILE *fpt);
/*-------global function declarations------end---------*/
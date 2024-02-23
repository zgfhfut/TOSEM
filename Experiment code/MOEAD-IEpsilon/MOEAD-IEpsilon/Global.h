/* This file contains the variable and function declarations */
#include <vector>
using namespace std;

# define INF 1.0e14
# define EPS 1.0e-14
# define E  2.71828182845905
# define PI 3.14159265358979

extern const double EPSINON;
extern double C;

//#define nreal    10
//#define nreal   50//complex system 
#define nreal   50//large system 

//#define reli_thresh   0.95//complex system 

//#define T_MAX_RESOURCE  8000.0  //simple system
#define T_MAX_RESOURCE  70000.0  //complex system
//#define T_MAX_RESOURCE  70000.0//large system

//#define MISSION_TIME  100.0  //workable life
#define MISSION_TIME  150.0  //workable life
//#define MISSION_TIME  200.0  //workable life
#define Max_R 0.9
extern double para_module[nreal][5];

/*-------global variables------start---------*/
typedef struct
{
	double *xreal;              // individual's variable value 
    double *obj;				// individual's objective value 

	double *constr;				// individual's constraint
	double constr_violation;	// individual's constraint violation degree	

	int rank;
	double crowd_dist;

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
extern double eta_m;


extern int niter;

extern int currenteval;
extern int neval;

extern double *min_realvar;
extern double *max_realvar;


extern double * z_;

extern double ** lambda_;

extern int T_;

extern double F_;

extern int ** neighborhood_;

extern double delta_;

extern int nr_;

extern double phi_max_;
extern double epsilon_k_;


extern double u[nreal];


void Cal_c();
/*-------global variables------end---------*/

void insert(list *node, int x);
list* del(list *node);
void insert_ind(ind_list *node, double *x);
ind_list* del_ind(ind_list *node);


/*-------global function declarations------start---------*/
void allocate_memory_pop (population *pop, int size);
void allocate_memory_ind (individual *ind);
void deallocate_memory_pop (population *pop, int size);
void deallocate_memory_ind (individual *ind);

void evaluate_pop(population *pop, int size);
void evaluate_ind(individual *ind);
void test_problem(double *xreal, double *obj, double *constr);


void initialize_pop (population *pop, int size);
void initialize_ind (individual *ind);

void DECrossover(individual *child, individual *parents);

void Cal_res_lim(double r);/*���������������ʱ��*/

void matingSelection(vector<int> &list, int cid, int size, int type);

int check_dominance(individual *a, individual *b);

void real_mutate_ind (individual *ind);

void merge(population *pop1, population *pop2, population *pop3, int size1, int size2);
void copy_ind(individual *ind1, individual *ind2);


void report_objective(population *pop, int size, FILE *fpt);
void report_variable(population *pop, int size, FILE *fpt);

void initUniformWeight();
void initNeighborhood();
double distVector(double * vector1, double * vector2, int dim);
void minFastSort(double * x, int * idx, int n, int m);
double GetFeasible_Ratio(population *pop, int size);

void initIdealPoint(population *pop);
void updateReference(individual *ind);
void randomPermutation(int * perm, int size);

void updateProblem(population *pop, individual * indiv, int id, int type);

double fitnessFunction(individual * indiv, double * lambda);

void fill_nondominated_sort(population *mixed_pop, population *new_pop);
void crowding_fill(population *mixed_pop, population *new_pop, int count, int front_size, list *elite);

void assign_crowding_distance_list(population *pop, list *lst, int front_size);

void assign_crowding_distance_indices(population *pop, int c1, int c2);

void assign_crowding_distance(population *pop, int *dist, int **obj_array, int front_size);

void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size);
void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right);

void quicksort_dist(population *pop, int *dist, int front_size);

void q_sort_dist(population *pop, int *dist, int left, int right);


/*-------global function declarations------end---------*/
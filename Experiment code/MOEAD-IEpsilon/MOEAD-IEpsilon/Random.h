/* Declaration for random number related variables and routines */

/* Variable declarations for the random number generator */
extern double seed;
extern double oldrand[55];
extern int jrand;

/* Function declarations for the random number generator */
void randomize();
void warmup_random (double seed);
void advance_random ();
double randomperc();
int rnd (int low, int high);
double rndreal (double low, double high);

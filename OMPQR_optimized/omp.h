#include <math.h>
#include <stdlib.h>

#define FABS(X) (((X)>0)?(X):(-(X)))

#define M 256
#define L_MAX 128


void omp(float D[L_MAX][M], float X[M], int l, int k, float epsilon,
			float V[L_MAX], int supp[L_MAX], int &supp_len);



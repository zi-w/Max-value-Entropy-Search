/*  Copyright (c) 2014, J.M. Hernandez-Lobato, M.W. Hoffman, Z. Ghahramani
 This function is from the code for the paper
 Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
 Predictive Entropy Search for Efficient Global Optimization of Black-box
 Functions, In NIPS, 2014.
 https://bitbucket.org/jmh233/codepesnips2014
 */
#include "mex.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

/* Function that computes the inverse of a matrix using its Cholesky decomposition */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	#define L_matlab prhs[0]
	#define ret_matlab plhs[0]

	int i, j, n, m;
	double *L, *ret;
	gsl_matrix *result;

	n = mxGetN(L_matlab);
	m = mxGetM(L_matlab);
	
	/* We ask for memory to store the matrix in the gnu library format */

	result = gsl_matrix_alloc(n, m);

	/* We copy the matrix from the matlab format to the gnu library format */

	L = mxGetPr(L_matlab);
	for (i = 0 ; i < n ; i++)
		for (j = 0 ; j < m ; j++)
			gsl_matrix_set(result, i, j, L[ i + n * j ]);

	gsl_linalg_cholesky_decomp(result);

	/* We obtain the inverse of the matrix */

	gsl_linalg_cholesky_invert(result);

	/* We ask for memory to return the solution */

	ret_matlab = mxCreateDoubleMatrix(n, m, mxREAL);

	/* We copy the solution from the gnu library representation to the matlab representation */

	ret = mxGetPr(ret_matlab);
	for (i = 0 ; i < n ; i++)
		for (j = 0 ; j < m ; j++)
			ret[ i + j * n ] = gsl_matrix_get(result, i, j);

	/* We are done */

	gsl_matrix_free(result);

	return;
}

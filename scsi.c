#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#define H_COLUMNS 3																// number of columns of H matrix
#define DBL_EPSILON 0.00001														// to compare two double
#define IS_EQUAL(x, y) 			(((fabs(x-y)) < (DBL_EPSILON)) ? (1) : (0))		// to compare two double

void ThomasAlgorithm (int n, double *xi, double *yi, double **m);
void CubicSplineInterpolation (double step, double *xi, double *yi, double *m, double **x, double **fx, int interval);


int main (int argc, char *argv[]) {

    int n;								// number of input elements

    struct timeval begin, total_time, time_without_writing;


    FILE *source = NULL;
	char *filename_input = "./input.txt";		            // name of input file (default: "./input.txt")
	char *filename_output = "./sequential_output.txt";		// name of output file (default: "./sequential_output.txt")
	
	double step = 0.1;				// step size (default: 0.1)

    double *xi = NULL;				// xi-coordinates array
	double *yi = NULL;				// yi-coordinates array

    double *m = NULL;				// H × m = r - unknowns array

	int sigma;						// value of total number generated by Cubic Spline Interpolation method
	double *x = NULL;				// array of x-coordinate of Cubic Spline Interpolation
	double *fx = NULL;				// array of y-coordinate of Cubic Spline Interpolation

    /* Command-line options */
	for (int i=0; i<argc; i++) {
		if (strcmp(argv[i], "-fi") == 0) {
			filename_input = argv[i+1];
		}
		if (strcmp(argv[i], "-s") == 0) {
			step = atof(argv[i+1]);
		}
		if (strcmp(argv[i], "-fo") == 0) {
			filename_output = argv[i+1];
		}
	}

    gettimeofday(&begin, 0);

    if(step <= 0) {
		fprintf(stderr, "ERROR: Step-size can't be negative or equal to zero!\n");
		return -1;
	}

	// Open the input file
	if((source = fopen(filename_input, "rt")) == NULL) {
		fprintf(stderr, "Error with input fopen!\n");
		return -1;
	}

	if(fscanf(source, "%d", &n) != 1) {		// Read the number of elements n
		fprintf(stderr, "ERROR: Expected one number as input for n value!\n");
		return -1;
	}

    if(n < 2) {
		fprintf(stderr, "ERROR: n can't be less than 2!\n");
		return -1;
	}

    xi = (double*) malloc(n * sizeof(double));
	if(xi == NULL) { fprintf(stderr, "Malloc error (xi variable)!\n"); return -1; }

	yi = (double*) malloc(n * sizeof(double));
	if(yi == NULL) { fprintf(stderr, "Malloc error (yi variable)!\n"); return -1; }

    for(int i = 0; i < n; i++) {			// Read xi-coordinates
		if(fscanf(source, "%lf", &xi[i]) != 1) {
			fprintf(stderr, "ERROR: Expected one number as input for xi array!\n");
			return -1;
		}
	}
	for(int i = 0; i < n; i++) {			// Read yi-coordinates
		if(fscanf(source, "%lf", &yi[i]) != 1) {
			fprintf(stderr, "ERROR: Expected one number as input for yi array!\n");
			return -1;
		}
	}
	fclose(source);
	source = NULL;

	if(step >= 1){							// Calculate total number generated by Cubic Spline Interpolation method
		sigma = ((xi[n-1]-xi[0]) * step) + 1;
	} else {
		sigma = ((xi[n-1]-xi[0]) / step) + 1;
	}

    /* Calculate m */
	ThomasAlgorithm(n, xi, yi, &m);

	/* Cubic Spline Interpolation method */
	CubicSplineInterpolation (step, xi, yi, m, &x, &fx, sigma);

    free(xi); xi = NULL;
	free(yi); yi = NULL;

    gettimeofday(&time_without_writing, 0);

    /* Write the ouput data */
	// Open the output file
	if((source = fopen(filename_output, "wt")) == NULL) {
		fprintf(stderr, "Error with output fopen!\n");
		return -1;
	}
	fprintf(source,"%d\n", sigma);							// write value of n
	for(int i = 0; i < sigma; i++) {
		fprintf(source,"%lf\t%lf\n", x[i], fx[i]);		// write io_xi and io_yi
	}
	fclose(source);

	free(x); x = NULL;
	free(fx); fx = NULL;

    gettimeofday(&total_time, 0);

    printf("\nWALL CLOCK TIME:\n");
	printf("\tTotal time:\t%f\n", (total_time.tv_sec - begin.tv_sec) + ((total_time.tv_usec - begin.tv_usec)*1e-6));
	printf("\tComputation:\t%f\n\n", (time_without_writing.tv_sec - begin.tv_sec) + ((time_without_writing.tv_usec - begin.tv_usec)*1e-6));

    return 0;
}


void ThomasAlgorithm (
	int n,							/* IN - Total elements */
	double *xi,						/* IN - Array of distributed xi */
	double *yi,						/* IN - Array of distributed yi */
	double **m						/* OUT - H × m = r - unknowns array */
	)
{
	double **H = NULL;
	double *HStorage = NULL;            // H × m = r - coefficients matrix
	double *r = NULL;					// H × m = r - results array

	// coefficients Thomas Algorithm:
	double *alpha = NULL;
	double *beta = NULL;
	double *gamma = NULL;


	// Allocate memory for partial_H matrix
	HStorage = (double *) calloc(n * H_COLUMNS, sizeof(double));
	if(HStorage == NULL) { fprintf(stderr, "Calloc error (HStorage variable)!\n"); return; }
	
	H = (double **) malloc (n * sizeof(double *));
	if(H == NULL) { fprintf(stderr, "Malloc error (H variable)!\n"); return; }
	for(int i = 0; i < n; i++) {
		H[i] = &HStorage[i*H_COLUMNS];
	}

    // Allocate memory for partial_r vector
	r = (double*) calloc(n, sizeof(double));
	if(r == NULL) { fprintf(stderr, "Calloc error (r variable)!\n"); return; }


	// Calculate matrix H and vector r
	H[0][0] = 0;
	H[0][1] = 1;
	H[0][2] = 0;
    H[n-1][0] = 0;
	H[n-1][1] = 1;
	H[n-1][2] = 0;
	r[0] = 0;
	r[n-1] = 0;
	for(int i = 1; i < n - 1; i++) {
		H[i][0] = xi[i] - xi[i-1];
		H[i][1] = 2 * ((xi[i] - xi[i-1]) + (xi[i+1] - xi[i]));
		H[i][2] = xi[i+1] - xi[i];
		r[i] = 6 * (((yi[i+1]-yi[i])/(xi[i+1]-xi[i]))-((yi[i]-yi[i-1])/(xi[i]-xi[i-1])));
	}

	/* Allocate memory for coefficients vectors  */
	alpha = (double*) calloc(n, sizeof(double));
	if(alpha == NULL) { fprintf(stderr, "Calloc error (alpha variable)!\n"); return; }

	beta  = (double*) calloc(n, sizeof(double));
	if(beta == NULL) { fprintf(stderr, "Calloc error (beta variable)!\n"); return; }

	gamma = (double*) calloc(n, sizeof(double));
	if(gamma == NULL) { fprintf(stderr, "Calloc error (gamma variable)!\n"); return; }

	

	/* Calculate alpha, beta and gamma */
	alpha[0] = 1;
	beta[0] = 0;
	gamma[0] = 0;
	for(int i = 1; i < n; i++) {
		alpha[i] = H[i][1] - H[i][0] * beta[i-1];
		beta[i] = H[i][2] / alpha[i];
		gamma[i] = (r[i] - H[i][0] * gamma[i-1])/alpha[i];
	}

	free(H); free(HStorage); H = NULL; HStorage = NULL;
	free(alpha); alpha = NULL;
	free(r); r = NULL;


	/* Allocate memory and calculate total m */
	*m = (double*) calloc(n, sizeof(double));
	if(*m == NULL) { fprintf(stderr, "Calloc error (m variable)!\n"); return; }

	(*m)[n-1] = gamma[n-1];
	for(int i = n-2; i >= 0; i--) {
		(*m)[i] = gamma[i] - beta[i]*(*m)[i+1];
	}

	free(beta); beta = NULL;
	free(gamma); gamma = NULL;
}

void CubicSplineInterpolation (
	double step,						/* IN - Step to generate numbers */
	double *xi,							/* IN - Array of distributed xi */
	double *yi,							/* IN - Array of distributed yi */
	double *m,							/* IN - Unknowns array */
	double **x,							/* OUT - Coordinate x of Cubic Spline Interpolation */
	double **fx,						/* OUT - Coordinate y of Cubic Spline Interpolation */
	int sigma							/* OUT - Intervals generated by processes */
	)
{

	// Cubic Spline Interpolation coefficients:
	double b = 0;
	double d = 0;
	
	// The following three assignments hold when p is equal to 1:
	int k = 0;								// Index to scroll xi and yi

	/* Allocate memory and calculate x and fx */
	*x = (double*) calloc(sigma, sizeof(double));
	if(*x == NULL) { fprintf(stderr, "Calloc error (x variable)!\n"); return; }

	*fx = (double*) calloc(sigma, sizeof(double));
	if(*fx == NULL) { fprintf(stderr, "Calloc error (gx variable)!\n"); return; }
	
	(*x)[0] = xi[0];
	(*fx)[0] = yi[0];
	b = (((yi[k+1] - yi[k]) / (xi[k+1] - xi[k])) - (((xi[k+1] - xi[k])/2) * m[0]) - (((xi[k+1] - xi[k]) / 6) * (m[1] - m[0])));
	d = (m[1]-m[0])/(6*(xi[k+1] - xi[k]));

	for(int i = 1; i < sigma; i++) {
		(*x)[i] = (*x)[i-1] + step;
		if(IS_EQUAL((*x)[i], xi[k+1])) {
			k++;
			b = (((yi[k+1] - yi[k]) / (xi[k+1] - xi[k])) - (((xi[k+1] - xi[k])/2) * m[k]) - (((xi[k+1] - xi[k]) / 6) * (m[k+1] - m[k])));
			d = (m[k+1]-m[k])/(6*(xi[k+1] - xi[k]));
			(*fx)[i] = yi[k];
			continue;
		}
		(*fx)[i] = (yi[k] + b*((*x)[i] - xi[k]) + (m[k]/2) * pow((*x)[i] - xi[k], 2) + d * pow((*x)[i] - xi[k], 3));
		//         - a[k] -                       --c[k]--
	}
}

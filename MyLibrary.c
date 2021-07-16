#include "MyLibrary.h"


void create_counts_and_displs_with_replications (
	int rank,     					/* IN - Process rank */
	int p,        					/* IN - Number of processes */
	int n,       					/* IN - Total number of elements */
	int replication,				/* IN - To replicate or not data between the processes */
	int **count,  					/* OUT - Array of counts */
	int **disp   					/* OUT - Array of displacements */
	)
{
	int buff_size = n/p;
	int more = buff_size + n - buff_size * p;
	int extra = 0;      // If n is not a multiple of p
	int extra1 = 0;     // To duplicate some elements between process
	int extra2 = 0;     // To duplicate some elements between process
	if(more > n/p){
		more -= n/p;
		extra = 1;
	}
	if(replication == 1) {
		extra1 = 1;
		extra2 = 2;
	}
	
	*count = (int*) malloc(p * sizeof(int));
	if(*count == NULL) { fprintf(stderr, "Malloc error (count variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
	*disp = (int*) malloc(p * sizeof(int));
	if(*disp == NULL) { fprintf(stderr, "Malloc error (disp variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	(*count)[0] = n/p + extra + extra1;
	(*disp)[0] = 0;
	more -= extra;
	if(p>1) {
		for (int i = 1; i < p-1; i++) {
			if(more <= 0) {
				extra = 0;
			}
			(*count)[i] = n/p + extra + extra2;
			(*disp)[i] = (*disp)[i-1] + (*count)[i-1] - extra2;
			more -= extra;
		}
		(*disp)[p-1] = (*disp)[p-2] + (*count)[p-2] - extra2;
		(*count)[p-1] = n/p + extra + extra1;
	}
}





void parallelThomasAlgorithm (
	int rank,						/* IN - Process rank */
	int p,							/* IN - Number of processes */
	int n,							/* IN - Total elements */
	int total_partial_size,			/* IN - Total number of elements for partition */
	int effective_partial_size,		/* IN - Effective number of elements for partition */
	double *xi,						/* IN - Array of distributed xi */
	double *yi,						/* IN - Array of distributed yi */
	double **m						/* OUT - H × m = r - unknowns array */
	)
{

	double **partial_H = NULL;
	double *partial_HStorage = NULL;    // H × m = r - partial coefficients matrix
	double **H = NULL;
	double *HStorage = NULL;            // H × m = r - coefficients matrix

	double *r = NULL;					// H × m = r - results array
	double *partial_r = NULL;			// H × m = r - partial results array
	double *total_m = NULL;             // H × m = r - temporary unknowns array

	// coefficients Thomas Algorithm:
	double *alpha = NULL;
	double *beta = NULL;
	double *gamma = NULL;

	int *count;
	int *displ;
	int *displ_H;                        // displacement to matrix H

	int m_size;                          // dimension of partial m (vector that will go out the function)


	// Allocate memory for partial_H matrix
	partial_HStorage = (double *) calloc(effective_partial_size * H_COLUMNS, sizeof(double));
	if(partial_HStorage == NULL) { fprintf(stderr, "Calloc error (partial_HStorage variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
	
	partial_H = (double **) malloc (effective_partial_size * sizeof(double *));
	if(partial_H == NULL) { fprintf(stderr, "Malloc error (partial_H variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
	for(int i = 0; i < effective_partial_size; i++) {
		partial_H[i] = &partial_HStorage[i*H_COLUMNS];
	}

    // Allocate memory for partial_r vector
	partial_r = (double*) calloc(effective_partial_size, sizeof(double));
	if(partial_r == NULL) { fprintf(stderr, "Calloc error (partial_r variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }



	{   // Calculate matrix partial_H and vector partial_r
		int start  = 0;
		int bound  = effective_partial_size;
		if(rank == 0) {					// Process rank 0 (first process)
			partial_H[0][0] = 0;
			partial_H[0][1] = 1;
			partial_H[0][2] = 0;
			partial_r[0] = 0;
			start = 1;
		} else if(rank == p - 1 ) {	    // Process rank p - 1 (last process)
			partial_H[effective_partial_size-1][0] = 0;
			partial_H[effective_partial_size-1][1] = 1;
			partial_H[effective_partial_size-1][2] = 0;
			partial_r[effective_partial_size-1] = 0;
			bound = effective_partial_size - 1;
		}
		for(int i = start, j = 1; i < bound || j < effective_partial_size; i++, j++) {
			partial_H[i][0] = xi[j] - xi[j-1];
			partial_H[i][1] = 2 * ((xi[j] - xi[j-1]) + (xi[j+1] - xi[j]));
			partial_H[i][2] = xi[j+1] - xi[j];
			partial_r[i] = 6 * (((yi[j+1]-yi[j])/(xi[j+1]-xi[j]))-((yi[j]-yi[j-1])/(xi[j]-xi[j-1])));
		}
	}


	/* Calculate count and displ, allocate memory for vector r and exchange arrays partial_r to obtain all elements */
	create_counts_and_displs_with_replications (rank, p, n, 0, &count, &displ);							// without replications
	r = (double *) calloc(n, sizeof(double));															// Allocate memory for vector r (with all the elements)
	if(r == NULL) { fprintf(stderr, "Calloc error (r variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
    
	MPI_Allgatherv(partial_r, count[rank], MPI_DOUBLE, r, count, displ, MPI_DOUBLE, MPI_COMM_WORLD);	// Echange vectors partial_r



	/* Calculate count and displ, allocate memory for matrix H and exchange arrays partial_HStorage to obtain all elements */
	displ_H = (int *) malloc(p * sizeof(int));
	for(int i = 0; i < p; i++) {
		count[i] = count[i] * H_COLUMNS;
		displ_H[i] = displ[i] * H_COLUMNS;
	}

	HStorage = (double *) calloc(n * H_COLUMNS, sizeof(double));			// Allocate memory for matrix H (with all the elements)
	if(partial_HStorage == NULL) { fprintf(stderr, "Calloc error (HStorage variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
    
	H = (double **) malloc (n * sizeof(double *));
	if(partial_H == NULL) { fprintf(stderr, "Malloc error (H variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
	for(int i = 0; i < n; i++) {
		H[i] = &HStorage[i*H_COLUMNS];
	}
	
	MPI_Allgatherv(partial_HStorage, count[rank], MPI_DOUBLE, HStorage, count, displ_H, MPI_DOUBLE, MPI_COMM_WORLD);		// Exchange matrix partial_H



	/* Allocate memory for coefficients vectors  */
	alpha = (double*) calloc(n, sizeof(double));
	if(alpha == NULL) { fprintf(stderr, "Calloc error (alpha variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	beta  = (double*) calloc(n, sizeof(double));
	if(beta == NULL) { fprintf(stderr, "Calloc error (beta variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	gamma = (double*) calloc(n, sizeof(double));
	if(gamma == NULL) { fprintf(stderr, "Calloc error (gamma variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	

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
	free(partial_H); free(partial_HStorage); partial_H = NULL; partial_HStorage = NULL;
	free(partial_r); partial_r = NULL;
	free(alpha); alpha = NULL;
	free(r); r = NULL;
	free(displ_H); displ_H = NULL;



	/* Allocate memory and calculate total m */
	total_m = (double*) calloc(n, sizeof(double));
	if(total_m == NULL) { fprintf(stderr, "Calloc error (total_m variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	total_m[n-1] = gamma[n-1];
	for(int i = n-2; i >= 0; i--) {
		total_m[i] = gamma[i] - beta[i]*total_m[i+1];
	}


	MPI_Barrier (MPI_COMM_WORLD);

	/* Calculate size of partial m and allocate memory */
	if(rank == 0) {
		m_size = total_partial_size;
	} else {
		m_size = total_partial_size - 1;	
	}

	*m = (double*) calloc(m_size, sizeof(double));
	if(*m == NULL) { fprintf(stderr, "Calloc error (*m variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	// Assign the values of total_m to m (this reduce the space in memory because each process will need only a few elements)
	for(int i = 0; i < m_size; i++) {
		(*m)[i] = total_m[displ[rank] + i];
	}


	free(beta); beta = NULL;
	free(gamma); gamma = NULL;
	free(count); count = NULL;
	free(displ); displ = NULL;
}





void parallelCubicSplineInterpolation (
	int rank,							/* IN - Process rank */
	int p,								/* IN - Number of processes */
	double step,						/* IN - Step to generate numbers */
	int effective_partial_size,			/* IN - Effective number of elements for partition */
	double *xi,							/* IN - Array of distributed xi */
	double *yi,							/* IN - Array of distributed yi */
	double *m,							/* IN - Unknowns array */
	double **x,							/* OUT - Coordinate x of Cubic Spline Interpolation */
	double **fx,						/* OUT - Coordinate y of Cubic Spline Interpolation */
	int *sigma_i						/* OUT - Intervals generated by processes */
	)
{

	// Cubic Spline Interpolation coefficients:
	double b;
	double d;
	
	// The following three assignments hold when p is equal to 1:
	int k = 0;								// Index to scroll xi and yi
	double first_x = xi[0];
	double first_y = yi[0];
	
	if(p > 1) {
		if(rank == 0) {
			k = 0;
			first_x = xi[0];
			first_y = yi[0];
			if(step >= 1){
				*sigma_i = (xi[effective_partial_size] - xi[0]) * step;
			} else {
				*sigma_i = (xi[effective_partial_size] - xi[0]) / step;
			}
		}
		if(rank != p - 1 && rank != 0) {
			k = 1;
			first_x = xi[1];
			first_y = yi[1];
			if(step >= 1){
				*sigma_i = (xi[effective_partial_size+1] - xi[1]) * step;
			} else {
				*sigma_i = (xi[effective_partial_size+1] - xi[1]) / step;
			}
		}
		if(rank == p - 1) {
			k = 1;
			first_x = xi[1];
			first_y = yi[1];
			if(step >= 1){
				*sigma_i = ((xi[effective_partial_size] - xi[1]) * step) + 1;
			} else {
				*sigma_i = ((xi[effective_partial_size] - xi[1]) / step) + 1;
			}
		}
	}
	

	/* Allocate memory and calculate x and fx */
	*x = (double*) calloc(*sigma_i, sizeof(double));
	if(*x == NULL) { fprintf(stderr, "Calloc error (x variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }

	*fx = (double*) calloc(*sigma_i, sizeof(double));
	if(*fx == NULL) { fprintf(stderr, "Calloc error (gx variable)!\n"); MPI_Abort(MPI_COMM_WORLD, -1); return; }
	
	(*x)[0] = first_x;
	(*fx)[0] = first_y;
	b = (((yi[k+1] - yi[k]) / (xi[k+1] - xi[k])) - (((xi[k+1] - xi[k])/2) * m[0]) - (((xi[k+1] - xi[k]) / 6) * (m[1] - m[0])));
	d = (m[1]-m[0])/(6*(xi[k+1] - xi[k]));

	for(int i = 1, j = 0; i < *sigma_i; i++) {
		(*x)[i] = (*x)[i-1] + step;
		if(IS_EQUAL((*x)[i], xi[k+1])) {
			k++;
			j++;
			b = (((yi[k+1] - yi[k]) / (xi[k+1] - xi[k])) - (((xi[k+1] - xi[k])/2) * m[j]) - (((xi[k+1] - xi[k]) / 6) * (m[j+1] - m[j])));
			d = (m[j+1]-m[j])/(6*(xi[k+1] - xi[k]));
			(*fx)[i] = yi[k];
			continue;
		}
		(*fx)[i] = (yi[k] + b*((*x)[i] - xi[k]) + (m[j]/2) * pow((*x)[i] - xi[k], 2) + d * pow((*x)[i] - xi[k], 3));
		//         - a[k] -                       --c[k]--
	}
}

/* PX390 Assignment 5: This program solves a system for reaction-diffusion equation for two molecules A and B in the late-time, steady state solution, with the assumption that the diffusion rate and source rate are smooth. The programm reads in a sereis of input parameters from 'input.txt' and reads in the function K(x) and S(x) from 'coeffecients.txt' and outputs the disaplacement x, A(x), B(x) to the file 'output.txt'. PDE's are solved through using a banded system of linear equation, discretised using finite difference methods, and then solved through using LAPACKE's dgbsv.f. These routines are adapted from the lectures notes of Prof. McMillan's PX390 module. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

void read_input (double *L, int *N, double *vel, double *tau);

/* Routines for using banded matrices from Prof. McMillan's PX390: L.10 */

struct band_mat{
    long ncol;              /* Number of columns in band matrix            */
    long nbrows;            /* Number of rows (bands in original matrix)   */
    long nbands_up;         /* Number of bands above diagonal           */
    long nbands_low;        /* Number of bands below diagonal           */
    double *array;          /* Storage for the matrix in banded format  */
    long nbrows_inv;        /* Number of rows of inverse matrix   */
    double *array_inv;      /* Store the matrix decomposition used to calculate inv_matrix*/
    int *ipiv;              /* Additional inverse information         */
};
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory, and set the parameters.  */
int init_band_mat (band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
    bmat->nbrows     = nbands_lower + nbands_upper + 1;
    bmat->ncol       = n_columns;
    bmat->nbands_up  = nbands_upper;
    bmat->nbands_low = nbands_lower;
    bmat->array      = (double *) malloc (sizeof (double)*bmat->nbrows*bmat->ncol);
    bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
    bmat->array_inv  = (double *) malloc (sizeof (double)*(bmat->nbrows + bmat->nbands_low)*bmat->ncol);
    bmat->ipiv       = (int *) malloc (sizeof (int)*bmat->ncol);

    /* Initialise array to zero, prevents intialising errors when running throug valgrind */
    long i;
    for (i = 0; i < bmat->nbrows*bmat->ncol; i++) {
        bmat->array[i] = 0.0;
    }
    for (i = 0; i < (bmat->nbrows + bmat->nbands_low)*bmat->ncol; i++) {
	bmat->array_inv[i] = 0.0;
    }
    for (i = 0; i < bmat->ncol; i++) {
	bmat->ipiv[i] = 0.0;
    }
    if (bmat->array == NULL || bmat->array_inv == NULL || bmat->ipiv == NULL) {
        printf ("Allocation of bmat->array, bmat->array_inv, or bmat->ipiv failed.\n");
        return 1;
    }
    return 0;
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full (unbanded)matrix.           */
double *getp (band_mat *bmat, long row, long column) {
    int bandno = bmat->nbands_up + row - column;
    if (row < 0 || column < 0 || row >= bmat->ncol || column >= bmat->ncol ) {
        printf ("Indexes out of bounds in getp: %ld %ld %ld \n", row, column, bmat->ncol);
        exit (1);
    }
    return &bmat->array[bmat->nbrows*column + bandno];
}

/* Return the value of a location in the band matrix, using
   the row and column indexes of the full (unbanded) matrix.           */
double getv (band_mat *bmat, long row, long column) {
    return *getp (bmat, row, column);
}

double setv (band_mat *bmat, long row, long column, double val) {
    *getp(bmat, row, column) = val;
    return val;
}

/* Solve the equation Ax = b for a matrix a stored in banded format, where x and b real arrays */
int solve_Ax_eq_b (band_mat *bmat, double *x, double *b) {
    /* Copy bmat array into the temporary store */
    int i, bandno;
    for (i = 0; i < bmat->ncol; i++) {
        for (bandno = 0; bandno < bmat->nbrows; bandno++) {
            bmat->array_inv[bmat->nbrows_inv*i + (bandno + bmat->nbands_low)] = bmat->array[bmat->nbrows*i + bandno];
        }
    x[i] = b[i];
    }
    long nrhs = 1;
    long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
    int info = LAPACKE_dgbsv (LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
    return info;
}

int printmat (band_mat *bmat) {
    long i, j;
    for (i = 0; i < bmat->ncol; i++) {
        for (j = 0; j < bmat->nbrows; j++) {
            printf ("%ld %ld %g \n", i, j, bmat->array[bmat->nbrows*i + j]);
        }
    }
    return 0;
}
/* free_struct frees all the internal dependancies in the bmat structure that are assigned on the malloc 
each time the band_mat data type is used*/
void free_struct (band_mat *bmat) {
	free (bmat->array);
	free (bmat->array_inv);
	free (bmat->ipiv);
}

#define DIAGS 0 // Diagnostic routines: set to 1 to print banded matrices for A & B (Prevents double comparison)
int main(void) {
    /* INPUT PARAMETERS: Read in from 'input.txt' */
    double L;                 /* Right Boundary of X Domain */
    int N;                    /* Number of Grid Points */
    double vel;               /* Advection Velocity */
    double tau;               /* Decay Rate */
    read_input (&L, &N, &vel, &tau);
    /* COEFFICIENTS: Read in from 'coefficients.txt' */
    double *K, *S;
    K = malloc (sizeof (double)*N);
    S = malloc (sizeof (double)*N);
    for (int i = 0; i < N; i++) {
        K[i] = 0;
        S[i] = 0;
    }
    if ( K == NULL || S == NULL) {
        printf ("Memory allocation failure of K or S. \n");
        return 1;
    }
    FILE *co_infile;
    co_infile = fopen ("coefficients.txt", "r");
    if (co_infile == NULL) {
        printf ("Error opening 'coefficients.txt'.\n");
        return 1;
    }
    for (int i= 0; i < N && (fscanf (co_infile, "%lf %lf", &K[i], &S[i]) == 2); i++);
    fclose (co_infile);

    long ncols = N-1;
    double delta_x = L/(N-1);

    /* Intialisation of each of the banded matrices for A & B respectively */
    band_mat bmat_A;
    band_mat bmat_B;
    /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1 */
    long nbands_low = 1;
    long nbands_up  = 1;

    /* Intisiation of each of the banded matrices for A & B*/
    init_band_mat (&bmat_A, nbands_low, nbands_up, ncols);
    init_band_mat (&bmat_B, nbands_low, nbands_up, ncols);
     
    /* Arrays for the linear system of arrays in the form a*x_{particle}=b for Particles A & B Respectively */
    double *x_A,*x_B, *b_A, *b_B;
    x_A = malloc (sizeof (double)*ncols);
    x_B = malloc (sizeof (double)*ncols);
    b_A = malloc (sizeof (double)*ncols);
    b_B = malloc (sizeof (double)*ncols);
    for (int j = 0; j < ncols; j++) {
        x_A[j] = 0;
        x_B[j] = 0;
        b_A[j] = 0;
        b_B[j] = 0;
    }
    if (x_A == NULL || x_B == NULL || b_A == NULL || b_B == NULL) {
        printf ("Memory allocation failure of x_A, x_B, b_A or b_B. \n");
        return 1;
    }
     
    /* Setting of bmat_A according to FDE analysis of PDE */
    long i;
    for (i = 0; i < ncols; i++) {
        if (i > 1) {
            setv (&bmat_A, i, i-1, (K[i] + vel*delta_x));
        }
        setv (&bmat_A, i, i, (-K[i+1] - K[i] - vel*delta_x - tau*delta_x*delta_x)); //i,j coeff
        }
    for (i = 1; i < ncols; i++) {
        if (i < ncols - 1) {
            setv (&bmat_A, i, i+1, (K[i+1]));
        }
        /* Source term for linear system of equation according to particle A  */
        b_A[i] = -S[i]*delta_x*delta_x;
    }
    /* A(X = 0) BOUNDARY CONDITIONS */
    setv(&bmat_A, 0,0, -(vel*vel*delta_x*delta_x)/K[0] - K[1] - vel*delta_x - tau*delta_x*delta_x);
    setv(&bmat_A, 0, 1, K[1]);
    
    solve_Ax_eq_b (&bmat_A, x_A, b_A);

    /* Setting of bmat_b according to FDE analysis of PDE */
    for (i = 0; i < ncols; i++) {
        if (i > 0) {
            setv (&bmat_B, i, i-1, (K[i]+delta_x*vel));
        }
        setv (&bmat_B, i, i, (-K[i+1] - K[i]-vel*delta_x));
    }
    for (i = 1; i < ncols; i++) {
        if (i < ncols - 1) {
            setv (&bmat_B, i, i+1, (K[i+1]));
        }
        /* Source term for linear system of equation according to particle B  */
        b_B[i] = -tau*x_A[i]*delta_x*delta_x;
    }
    /* B(X = 0) BOUNDARY CONDITIONS */
    setv(&bmat_B, 0, 0, -(vel*vel*delta_x*delta_x)/K[0] - K[1] - vel*delta_x);
    setv(&bmat_B, 0, 1, K[1]);

    solve_Ax_eq_b (&bmat_B, x_B, b_B);
    
    /*  Print matrix for debugging, to turn on set DIAGS = 1 */
    if (DIAGS) {
        printmat (&bmat_A);
        printmat (&bmat_B);
    }
    
    /* Writing of x, A(x), B(x) to the file 'output.txt' */
    FILE *outfile;
    outfile = fopen ("output.txt", "w");
    if (outfile == NULL) {
        printf ("Error opening 'output.txt'.\n");
        return 1;
    }
    for (i = 0; i < ncols; i++) {
        fprintf (outfile, "%g %g %g \n", i*delta_x, x_A[i], x_B[i]);
    }
    /* Fixing of boundary conditions at x = L, where A = B = 0 */
    fprintf (outfile, "%g %g %g \n", L, 0.0, 0.0);
    fclose (outfile);
    
    /* Freeing of all variables used */
    free (x_A);
    free (x_B);
    free (b_A);
    free (b_B);
    free (K);
    free (S);
    free_struct (&bmat_A);
    free_struct (&bmat_B);
    return 0;
}

void read_input (double *L, int *N, double *vel, double *tau) {
   FILE *infile;
   if (!(infile = fopen ("input.txt","r"))) { 
       printf ("Error opening 'input.txt'\n");
       exit (1);
   }
   if (4!=fscanf (infile,"%lf %d %lf %lf ", L, N, vel, tau)) {
       printf ("Error reading parameters from 'input.txt'\n");
       exit (1);
   }
   fclose (infile);
}

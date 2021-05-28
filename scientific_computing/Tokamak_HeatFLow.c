/* PX390 Assignment 5: This program models the heat transfer around the toroidal shell of a tokamak in two separate modes: one in the steady-state (i.e. at the late time, after the temperature distribution has relaxed) and in the time-evolution mode for the initial condition of T=0.The program reads in a series of input parameters from 'input.txt' and reads in a series of position-variant functions from 'coefficients.txt' and outputs the temperature profile for all grid points in one single column to the file 'output.txt'. For the case of the steady state, the PDE is  solved through using a banded system of linear equations, discretised using finite difference methods, and then solved through using LAPACKE's dgbsv.f. While for the case of the time-evolution mode, an implicit solver is employed in combination with a finite difference method. These routines are adapted from the lectures notes of Prof. McMillan's PX390 module. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
void read_input (long *N_theta, long *N_zeta, int *E, double *time_final, double *I_min);

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

    /* Initialise array to zero, prevents initialising errors when running through valgrind */
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
   the row and column indexes of the full (unbanded) matrix.           */
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

void setv (band_mat *bmat, long row, long column, double val) {
    *getp(bmat, row, column) = val;
}

/* Solve the equation A x = b for a matrix a stored in banded format, where x and b real arrays */
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
/* free_struct frees all the internal dependencies in the bmat structure that are assigned on the malloc each time the band_mat data type is used*/
void free_struct (band_mat *bmat) {
	free (bmat->array);
	free (bmat->array_inv);
	free (bmat->ipiv);
}
/*Check that a grid point has valid coordinates 
within the specified domain of the problem */
int is_valid (long t, long z, long T, long Z) {
  return (t >= 0) && (t < T) && (z >= 0) && (z < Z);
}

/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long t:  The theta grid point index
   long z:  The zeta grid point index
   long N_theta:  The number of theta points
   long N_zeta:  The number of Xi  points
*/
long indx (long t, long z, long N_zeta) {
    // L = z + t*N_zeta; As defined in specification
    return t*N_zeta +z;
    }
/* reindx folds all of the indices so that the 2D matrix 
can be banded into a 1D matrix.*/
long reindx (long L , long N_theta, long N_zeta) {
    long N = N_theta*N_zeta;
    if (L < 0.5*N) {
        return 2*L;
    }
    else {
    return 2*(N-L)-1;
    }
}
/* Return the 2D point corresponding to a particular 1D grid index */
void gridp (long indx, long Z, long *t, long *z) {
  *t = indx%Z;
  *z = indx - (*t)*Z;
}

#define DIAGS  0  // Set to 1 to turn on diagnostic conditions
int main (void) {
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* INPUT PARAMETERS: Read in from 'input.txt' */
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    long N_theta;           // Number of theta grid points
    long N_zeta;            // Number of zeta  grid points
    int E;                  // if E = 0, solve for steady state, else time evolution mode
    double time_final;      // Final time for time evolution mode
    double I_min;           // Minimum number of time steps to take, if in time evolution mode
    read_input (&N_theta, &N_zeta, &E, &time_final, &I_min);

    /* Initialising Banded Matrix*/
    band_mat bmat;
    /* ngridpoints:  Total size of problem is number of grid points on 2D plane (theta, zeta) */
    long  ngridpoints = N_theta*N_zeta;
    /* The function fmax deals with cases for N_zeta != N_theta, and thus ensure the banded
    matrix has the correct number of bands. While the factor of 3 accounts for the maximum 
    separation in the case of an odd number of bands */
    long nbands_low = 3*fmax(N_zeta, N_theta)-1; 
    long nbands_up  = nbands_low;
    init_band_mat (&bmat, nbands_low, nbands_up,  ngridpoints);

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* COEFFICIENTS: Read in from 'coefficients.txt'*/
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* This fill contains ngridpoints (N_theta*N_zeta) lines that are all indexed according to the indexing 
    function of L = t*N_zeta +z. The first line in the coefficients file corresponds to L=0, and subsequent 
    lines are increasing values of L*/
    long double *Q_11, *Q_22;
    double *H, *S, *R;
    Q_11 = malloc (sizeof (long double)* ngridpoints);
    Q_22 = malloc (sizeof (long double)* ngridpoints);
    H    = malloc (sizeof (double)* ngridpoints);
    S    = malloc (sizeof (double)* ngridpoints);
    R    = malloc (sizeof (double)* ngridpoints);
    /*
    Q_11 & Q_22: represent combination of geometric quantities; these capture the
    shape of the shell multiplied by the thermal conductivity.
    H: mass per unit area divided by the Jacobian of  the coordinate system.
    R: represent the connection of the shell to external heat sinks.
    S: normalised energy source due to fusion plasma heating at walls.
    */
    long i = 0;
    for (i = 0; i <  ngridpoints; i++) {
        Q_11[i] = 0;
        Q_22[i] = 0;
        H[i] = 0;
        S[i] = 0;
        R[i] = 0;
    }
    if (Q_11 == NULL || Q_22 == NULL || H == NULL || S == NULL || R == NULL) {
        printf ("Memory allocation failure of Q_11, Q_22, H, S or R. \n");
        return 1;
    }
    FILE *co_infile;
    co_infile = fopen ("coefficients.txt", "r");
    if (co_infile == NULL) {
        printf ("Error opening 'coefficients.txt'.\n");
        return 1;
    }
    for (i = 0; i <  ngridpoints && (fscanf (co_infile, "%Lf %Lf %lf %lf %lf", &Q_11[i], &Q_22[i], &H[i], &S[i], &R[i]) == 5); i++);
    fclose (co_infile);

    

    double *T   = malloc (sizeof (double)* ngridpoints);
    double *b   = malloc (sizeof (double)* ngridpoints);
    double *T_n = malloc (sizeof (double)* ngridpoints);
    /* Initialisation of arrays, before they are read in from 'coefficients.txt' */
    for (i = 0; i <  ngridpoints; i++) {
        T[i]   = 0;
        b[i]   = 0;
        T_n[i] = 0;
    }
    if (T == NULL || T_n == NULL || b == NULL) {
        printf ("Memory allocation failure of T, T_n or b.\n");
    }
    /* theta and xi both are bounded in the domain of 0 to 2*PI.
    Thus the interval of steps between each of these respective points 
    in this domain divided by the number of grid points */
    /* Defining pi exactly as 4*atan(1) */
    double d_theta = 2*(4*atan(1))/N_theta;
    double d_zeta    = 2*(4*atan(1))/N_zeta;
    /* Time step dt is set from the condition that our system must finish at t = time_final 
    and I_min is the minimum number of steps to get from 0 to time_final. Using I_min + 1 and not I_min so that the 
    timestep is just smaller than the absolute minimum required This is only used in the time-stepping mode. */
    double dt = time_final/(I_min+1);
    /* t & x both run from 0 to N_theta or N_zeta respectively: they define the 
    grind points for theta and zeta respectively*/
    
    long t = 0, z = 0;
    if (E == 0) {
        for (t = 0; t < N_theta; t++) {
            for (z = 0; z < N_zeta; z++) {
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                // LATE TIME STEADY STATE SOLUTION
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                /* Enforces periodicity of indices */
                long tp = (t+1)%N_theta;
                long tm = (t+N_theta-1)%N_theta;
                long zp = (z+1)%N_zeta;
                long zm = (z+N_zeta-1)%N_zeta;
                long index = indx (t, z, N_zeta);
                long unknown_reindx = reindx (index,  N_theta, N_zeta);

                /* Coefficients produced through using finite difference method with a backwards scheme  */
                double C_A = H[index] * Q_11[index] / (d_theta*d_theta);
                double C_B = H[index] * Q_11[indx (tp, z, N_zeta)] / (d_theta*d_theta);
                double C_C = H[index] * Q_22[index] / (d_zeta*d_zeta);
                double C_D = H[index] * Q_22[indx (t, zp, N_zeta)] / (d_zeta*d_zeta);
                double C_E = -H[index] * (Q_11[index] + Q_11[indx (tp, z, N_zeta)]) / (d_theta*d_theta) - H[index] * (Q_22[index] + Q_22[indx (t, zp, N_zeta)]) / (d_zeta*d_zeta) - R[index];

                /* Setting of the coefficients of the banded matrix before solving using solve_Ax_eq_b function */
                setv (&bmat, unknown_reindx, reindx (indx (tm, z, N_zeta), N_theta, N_zeta), C_A);
                setv (&bmat, unknown_reindx, reindx (indx (tp, z, N_zeta), N_theta, N_zeta), C_B);
                setv (&bmat, unknown_reindx, reindx (indx (t, zm, N_zeta), N_theta, N_zeta), C_C);
                setv (&bmat, unknown_reindx, reindx (indx (t, zp, N_zeta), N_theta, N_zeta), C_D);
                setv (&bmat, unknown_reindx, reindx (index, N_theta, N_zeta), C_E);

                /* Source term (normalised energy source due to fusion)*/
                b[unknown_reindx] = -S[index];
            }
        }
        solve_Ax_eq_b(&bmat, T, b);
    }

    else {
        // %%%%%%%%%%%%%%%%%%%%%%%
        // TIME EVOLUTION SOLUTION
        // %%%%%%%%%%%%%%%%%%%%%%%
        double c_time = 0;
        /* Setting of Initial Condition: T(t=0) = 0. */
        for (i = 0; i < ngridpoints; i++) {
            T[i] = 0.0;
        }
        while (c_time < time_final) {
            // Looping over all points in 2d space
            for (t = 0; t < N_theta; t++) {
                for (z = 0; z < N_zeta; z++) {
                    /* Enforces periodicity of indices */
                    long tp = (t+1)%N_theta;
                    long tm = (t+N_theta-1)%N_theta;
                    long zp = (z+1)%N_zeta;
                    long zm = (z+N_zeta-1)%N_zeta;
                    long index = indx (t, z, N_zeta);                        
                    long unknown_reindx = reindx (index,  N_theta, N_zeta)
                    ;
                    /* Coefficients produced through using finite difference 
                    method with a backwards scheme  */
                    double C_A = H[index] * Q_11[index] / (d_theta*d_theta);
                    double C_B = H[index] * Q_11[indx (tp, z, N_zeta)] / (d_theta*d_theta);
                    double C_C = H[index] * Q_22[index] / (d_zeta*d_zeta);
                    double C_D = H[index] * Q_22[indx (t, zp, N_zeta)] / (d_zeta*d_zeta);
                    double C_E = H[index] * (Q_11[index]  + Q_11[indx(tp, z, N_zeta)]) /
                    (d_theta*d_theta) + H[index] * (Q_22[index] + Q_22[indx(t, zp, N_zeta)]) / (d_zeta*d_zeta) + R[index];

                    /* Setting of the coefficients of the banded matrix before 
                    solving using solve_Ax_eq_b function */
                    setv (&bmat, unknown_reindx, reindx (indx (tm, z, N_zeta), N_theta, N_zeta), -dt*C_A);
                    setv (&bmat, unknown_reindx, reindx (indx (tp, z, N_zeta), N_theta, N_zeta), -dt*C_B);
                    setv (&bmat, unknown_reindx, reindx (indx (t, zm, N_zeta), N_theta, N_zeta), -dt*C_C);       
                    setv (&bmat, unknown_reindx, reindx (indx (t, zp, N_zeta), N_theta, N_zeta), -dt*C_D);
                    setv (&bmat, unknown_reindx, reindx (index, N_theta, N_zeta), 1+dt*C_E);
                    /* Source term now equal to b = T+dt*s*/
                    b[unknown_reindx] = T[reindx (index, N_theta, N_zeta)] + dt*S[index];
                }
            }
            solve_Ax_eq_b (&bmat, T_n, b);
            /* Pointer Swap using temporary variable tmp_var */
            double *tmp_var;
            tmp_var = T_n;
            T_n = T;
            T = tmp_var;
            
            /* Incrementing of time up from c_time = 0 up to the final value of c_time = time_final */
            c_time += dt;
        }
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* Writing of T(theta, xi) to the file 'output.txt' */
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE *outfile;
    outfile = fopen ("output.txt", "w");
    if (outfile == NULL) {
        printf ("Error opening 'output.txt'.\n");
        return 1;
    }
    for (t = 0; t < N_theta; t++) {
        for (z = 0; z < N_zeta; z++) {
            long final_index = reindx (indx(t, z, N_zeta), N_theta, N_zeta);
            /* Must re-index T so that is it in the indexing defined in specification
            i.e. the function L as opposed to the re-indexing folding function used */
            fprintf (outfile, "%g\n", T[final_index]);
        }
    }
    fclose (outfile);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* DIAGNOSTIC TOOLS: Set DIAGS = 1, to turn on */
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (DIAGS) {
        printf ("N_theta, N_zeta, E, time_final, I_min\n");
        printf ("%li %li %d %lf %lf \n", N_theta, N_zeta, E, time_final, I_min);
        printf ("d_zeta = %g; d_theta = %g; dt = %g\n", d_zeta, d_theta, dt);
        printf ("Q_11[i], Q_22[i], H[i], S[i], R[i]\n");
        for (i = 0; i < ngridpoints; i++) {
            printf ("%Lf %Lf %lf %lf %lf\n", Q_11[i], Q_22[i], H[i], S[i], R[i]);
        }
        printf ("Theta, Zeta, T\n");
        for (t=0; t < N_theta; t++) {
            for (z = 0; z < N_zeta; z++) {
                long final_index = reindx (indx(t, z, N_zeta), N_theta, N_zeta);
                printf("%li %li %g \n", t, z, T[final_index]);
            }
        }
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* Freeing of all variables assigned to the malloc */
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    free (Q_11);
    free (Q_22);
    free (H);
    free (S);
    free (R);
    free (T);
    free (T_n);
    free (b);
    free_struct (&bmat);
    return 0;
}
void read_input (long *N_theta, long *N_zeta, int *E, double *time_final, double *I_min) {
   FILE *infile;
   if (!(infile = fopen ("input.txt","r"))) { 
       printf ("Error opening 'input.txt'\n");
       exit (1);
   }
   if (5!=fscanf (infile,"%li %li %d %lf %lf", N_theta, N_zeta, E, time_final, I_min)) {
       printf ("Error reading parameters from 'input.txt'\n");
       exit (1);
   }
   fclose (infile);
}

/* PX390 Assignment 4: This program solves the non linear diffusion equation as described in the 
specification (d_dx(T) + d_dx(gamma) = 0). This is subject to both initial and boundary conditions.
The program reads in a series of parameters from 'input.txt' and outputs the time, space and temperature 
profile solution to the PDE desribed above. This is particularly useful when looking at the heat transport 
in Tokamaks.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void read_input (double *L, int *N, double *sim_time, double *timestep_d, 
                double *T_0, double *amp, double *omega, double *P);

int main (void) {
    
    /* INPUT PARAMETERS: These are all read in from 'input.txt' */

    double L;               /* Right Boundary of X Domain */
    int N;                  /* Number of Grid Points */
    double sim_time;        /* Length of Time to run Simulation */
    double timestep_d;      /* Time-step for diagnostic output */
    double T_0;             /* Initial Condition Parameter */
    double amp;             /* Amplitude of Oscillation ( Boundary Parameter ) */
    double omega;           /* Period of Oscillation ( Boundary Parameter ) */
    double P;               /* Non Linearity Exponent ( This can be a Float ) */

    read_input (&L, &N, &sim_time, &timestep_d, &T_0, &amp, &omega, &P);
    /* Printed outputs that check input variables are sensible. */
    if (T_0 < 0) {
        printf ("The value of T_0 is negative. Consider changing input parameters.\n");
    }
    if (omega < 0) {
        printf ("The value frequency of oscillation is negative. Consider changing input parameters.\n");
    }

    double dx = L/ (N-1);
    /* Time Step set through Van Neumann Analysis; 
    dt < dx^2 /(2*P* (abs (B))^ (P-1); B ~ 14 */ 
    int B = 14;
    double dt = (dx*dx)/ (3*P* pow (abs (B), (P-1)));       
    
    /* MALLOC ALLOCATION */
    
    /* T is the temperature of the system as a function of space and time,
    while T_n is this same function but forward in time by one time grid space */
    double *T, *T_n;
    T      = malloc (sizeof (double)*N);
    T_n    = malloc (sizeof (double)*N);
    if (T == NULL || T_n == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    /*  INTIALISATION */

    /* Intialising of arrays assigned on malloc to prevent
    intialisation error when running through valrgind */
    for (int j = 0; j < N; j++) {
       T[j] = 0.0;
       T_n[j] = 0.0;
    }
    
    int i;      /* Space Index */
    double x;   /* Displacement 'x': x = i*dx */

    for (i = 0; i < N; i++) {
        x = i*dx;
        /* Initial condition at t = 0 for all x */
        T[i] = T_0* (1- (x/L));
    }
    double c_time = 0.0;     /* Current time (in the loop) */

    /* Output file containing the time, x coordinate, and T(x,t): spaced on interval of t_d.
    This is the result of solving the PDE as specified.*/
    FILE *outfile;
    outfile = fopen ("output.txt", "w");
    if (outfile == NULL) {
        printf ("Error opening file.\n");
        return 1;
    }

    for (i = 0; i < N; i++) {
        x = i*dx;
        fprintf (outfile, "%g %g %g \n", c_time, x, T[i]);
    }
    double next_output_time = timestep_d;

    /* INTERATING TO SOLVE PDE */
    
    while (c_time < sim_time) {
        double dt0 = dt;
        int output = 0;
        /* Prevents overshooting of time-step */
        if (c_time + dt0 > next_output_time) {
            dt0 = next_output_time - c_time;
            output = 1;
        }

        // LOOPING OVER POINTS
        for (i = 0; i < N; i++) {
            x = i*dx;
            int im = (i - 1);
            int ip = (i + 1);
            if (i == 0) {
                /* Boundary condition set at x =L for all t */
                T[i] = amp* (1 - cos (omega*c_time)) + T_0;
            }
            
            else if (i == N-1) { // dT/dx | (x=L) = -T_0/L
                /* Boundary condition set at x = L for all t */
                T[i] = T[im] - (T_0*dx)/L;
            } 
            
            else {
                /* The PDE here is now begin solved in across the indices where boundary 
                conditions do not apply. This has been achieved through expanding out 
                derivative using chain rule and then discretising the equation using 
                a central difference scheme. */
                double dT_dx = fabs ((T[ip] - T[im]) / (2*dx));
                double dT2_dx2 = (T[im] - 2*T[i] + T[ip]) / (dx*dx);

                T_n[i] = T[i] + (dt0*P*pow (dT_dx, (P-1))*dT2_dx2);
            }
        }
        
        /* REALLOCATING ARRAYS VIA POINTER SWAP */

        double *tmp;
        tmp = T_n;
        T_n = T;
        T = tmp;

        c_time += dt0;
        if (output) {
            for (i = 0; i < N; i++) {
                x = i*dx;
                fprintf (outfile, "%g %g %g \n", c_time, x, T[i]);
            }
            next_output_time += timestep_d;
        }
    }
    free (T);
    free (T_n);
    fclose (outfile);
    return 0;
}

void read_input (double *L, int *N, double *sim_time, double *timestep_d, 
                double *T_0, double *amp, double* omega, double *P) {
   FILE *infile;
   if (!(infile=fopen ("input.txt","r"))) { 
       printf ("Error opening file\n");
       exit(1);
   }
   if (8!=fscanf (infile,"%lf %d %lf %lf %lf %lf %lf %lf", L, N, sim_time, timestep_d, T_0, amp, omega, P)) {
       printf ("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}


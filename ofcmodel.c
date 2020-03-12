/*======================================================/
/              PX425 2019 Assignment 5                  /
/            Self organised criticality                 /
/                      OF                               /
/                  EARTHQUAKES!                         /
/======================================================*/
/*=====================================================*/

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <time.h>
#include "mt19937ar.h"
#include "makePNG.h"

int read_args(int argc, char ** argv, double * alpha, int * nalpha, int * Nsteps, int * Nx, int * Ny);

int isprint(int opt);

int main(int argc, char ** argv) {

  /* Timing Variables */
  double start_time, end_time;

  /* Number of processor and current rank */
  int p, my_rank;

  /* Critical threshold for an earthquake */
  const double threshold = 1.0;

  /* 'conservative' parameter (also 'redistribution' parameter) */
  double alpha, alpha_of;

  /* Whether we are still resolving an active earthquake */
  int quake;

  /* Size of grid in x and y directions */
  int Nx, Ny;

  /* Number of steps to run for */
  int Nsteps;

  /* Number of alpha values to use (counting from alpha0 to 0.25) */
  int nalpha;

  /* Interval at which to print images - reduce to print more often */
  const int Nimage = 10000000;

  /* Interval at which to print counters */
  const int Ncount = 10000;

  /* loop counters */
  int ix, iy, istep, isnap = 0;
  int ichain, maxchain;

  /* local variables */
  double xi;
  double maxforce;

  /* filename for output images and histograms */
  char filename[25];

  /* Seed random number generator */
  unsigned long seed = 20340293;
  init_genrand(seed);

  /* Output file pointer */
  FILE *outfile;

  /* Default sizes of grid, conservative parameter, step count and alpha count */
  Nx = 15;
  Ny = 15;
  double alpha0 = 0.25;
  Nsteps = 100000;
  nalpha = 1;

  /* Read command line arguments */
  int ra = read_args (argc, argv, &alpha0, &nalpha, &Nsteps, &Nx, &Ny);

  /* Echo command line arguments to stdout */
  if (!ra) {
    printf("# Command line arguments successfully processed\n");
  }
  if (ra) {
    printf("Error Processing command line arguments\n");
    return EXIT_FAILURE;
  }
  /* Allocate memory for a grid this size */
  double *force_mem = (double *) calloc(Ny * Nx, sizeof(double));
  double *oldforce_mem = (double *) calloc(Ny * Nx, sizeof(double));
  double **force = (double **) malloc(Ny * sizeof(double * ));
  double **oldforce = (double **) malloc(Ny * sizeof(double * ));
  for (iy = 0; iy < Ny; iy++) {
    force[iy] = & force_mem[iy * Nx];
    oldforce[iy] = & oldforce_mem[iy * Nx];
  }

  /* Generate the random initial grid */
  for (iy = 0; iy < Ny; iy++) {
    for (ix = 0; ix < Nx; ix++) {
      xi = genrand();
      force[iy][ix] = xi;
    }
  }

  /* Allocate memory for a histogram of the earthquake sizes */
  int histsize = 100000;
  int *hist = (int *) calloc(histsize, sizeof(double));
  int *hist_allranks = (int *) calloc(histsize, sizeof(double));
  int *histval = (int *) calloc(histsize, sizeof(double));
  /* Last 20% of the histogram is "overflow" points on a logarithmic scale up to 10^8 */
  int ih = histsize - 1;
  histval[ih] = pow(10, 8);
  double histp = 0.8;
  double histb = log(histp * histval[ih] / (double)(histsize)) / (double)(histsize) / (1.0 - histp);
  double histA = histval[ih] / exp(histb * (double)(histsize));
  for (ih = 0; ih < histsize; ih++) {
    hist[ih] = 0;
    hist_allranks[ih] = 0;
    if (ih < histsize * histp) {
      /* Linearly-spaced bit of the histogram */
      histval[ih] = ih;
    } else {
      /* Exponential bit */
      histval[ih] = (int)(histA * exp(histb * (ih + 1)));
    }
  }
  maxchain = 0;

  int ialpha;
  double dalpha = (0.25 - alpha0);
  if (fabs(dalpha) < DBL_EPSILON) {
    nalpha = 1;
    dalpha = 0.0;
  } else {
    dalpha = dalpha / ((double)(nalpha - 1));
  }

  /* -------------- */
  /* Initialise MPI */
  /* -------------- */
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  /* ----------------------------/
  / Start loop over alpha values /
  /-----------------------------*/
  for (ialpha = 0; ialpha < nalpha; ialpha++) {
    alpha = alpha0 + ((double) ialpha) * dalpha;

    /*-------------/
    / Start Timer  /
    /-------------*/
    start_time = MPI_Wtime();
    /*----------------/
    / BEGIN MAIN LOOP /
    /----------------*/
    istep = 0;
    /* Current index in istep reached */
    int istep_current = 0; 
    /* Current index in chain reached */
    int ichain_current = 0;
    /* Iterations reaming before either output */
    int istep_out = 0; 
    /* Iterations remaining before Ncount ouput*/
    int istep_Ncount = 0; 
    /* Nsteps completed this round */
    int istep_Nstep = 0; 
    while (istep_current <= Nsteps) {
      istep_out = Nimage - istep_current % Nimage;
      istep_Ncount = Ncount - istep_current % Ncount;
      istep_Nstep = 1 + Nsteps - istep_current;
      if (istep_out > istep_Ncount) {
        istep_out = istep_Ncount;
      }
      if (istep_out > istep_Nstep) {
        istep_out = istep_Nstep;
      }
      if (my_rank < istep_out) {
       /* Local istep inside each rank */
        istep = istep_current + my_rank;
        /*-------------------------------------------/
        /  Check for sites over threshold
        /-------------------------------------------*/
        quake = 0;
        maxforce = 0.0;
        for (iy = 0; iy < Ny; iy++) {
          for (ix = 0; ix < Nx; ix++) {
            /* Find maximum value of force */
            if (force[iy][ix] > maxforce) {maxforce = force[iy][ix];}
            /* Is any site over the threshold? */
            if (force[iy][ix] >= threshold) {quake = 1;}
          } /* ix */
        } /* iy */

        if (!quake) {
          for (iy = 0; iy < Ny; iy++) {
            for (ix = 0; ix < Nx; ix++) {
              /* Advance time to when highest site becomes critical */
              force[iy][ix] += (1.0 - maxforce);
            } /* ix */
          } /* iy */
          quake = 1;
        }

        /*------------------------/
        / Evolution of Earthquake
        /------------------------*/
        ichain = 0;
        while (quake) {
          quake = 0;
          for (iy = 0; iy < Ny; iy++) {
            for (ix = 0; ix < Nx; ix++) {
              oldforce[iy][ix] = 0.0;
              /* Is this site above the threshold ? */
              if (force[iy][ix] >= threshold) {
                ichain++; /* increment counter for earthquake size */
                oldforce[iy][ix] = force[iy][ix]; /* make note of old force */
                force[iy][ix] = 0.0;
              } /* if(force[iy][ix] >= threshold) */
            } /* iy */
          } /* ix */

          for (iy = 0; iy < Ny; iy++) {
            for (ix = 0; ix < Nx; ix++) {
              if (oldforce[iy][ix] >= threshold) {
                alpha_of = alpha * oldforce[iy][ix];
                if (ix != 0)      {force[iy][ix-1] += alpha_of;}
                if (ix != Nx - 1) {force[iy][ix+1] += alpha_of;}
                if (iy != 0)      {force[iy-1][ix] += alpha_of;}
                if (iy != Ny - 1) {force[iy+1][ix] += alpha_of;}
              } /* if (oldforce) */
            } /* iy */
          } /* ix */

          for (iy = 0; iy < Ny; iy++) {
            for (ix = 0; ix < Nx; ix++) {
              if (force[iy][ix] >= threshold) {quake = 1;}
            } /* ix */
          } /* iy */

        }

        /* --------------------------------------------------------/
        / Find the histogram bin into which to deposit this result /
        /---------------------------------------------------------*/
        ih = ichain;
        if (ih >= histsize * histp) {
          for (ih = histsize * histp - 1; ih < histsize; ih++) {
            if (histval[ih] > ichain) {
              break;
            }
          }
        }
        hist[ih]++; /* increment histogram bin */

        /*---------------------------------/
        /  Print images every Nimage steps /
        /---------------------------------*/
        if (istep % Nimage == 0) {
          sprintf(filename, "snapshot%08d.png", isnap);
          writePNG(filename, force, Ny, Nx);
          isnap++;
        }
        /*-------------------------------------/
        /  Print counters every Ncount steps   /
        /-------------------------------------*/
        MPI_Allreduce (&ichain, &ichain_current, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (ichain_current > maxchain) {
          maxchain = ichain_current;
        };
        if (istep % Ncount == 0) {
          printf("# Step : %8d %8d\n", istep, maxchain);
          maxchain = 0;
        }
      } /* if(my_rank<istep_out) */
      else {
        MPI_Allreduce (&ichain, &ichain_current, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      } 
      if (istep_out > p) {
        istep_current += p;
      } else {
        istep_current += istep_out;
      }

    } /* istep */

    MPI_Reduce (hist, hist_allranks, histsize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
      /*------------------*/
      /* Open output file */
      /*------------------*/
      sprintf(filename, "hist_%5.3lf.dat", alpha);
      outfile = fopen(filename, "w");
      if (outfile == NULL) {
        printf("Error opening output file!\n");
        exit(EXIT_FAILURE);
      }
      /*--------------------------
      / Print histogram to file  /
      /-------------------------*/
      maxchain = 0;
      for (ih = 0; ih < histsize; ih++) {
        if (hist_allranks[ih] > 0) {maxchain = ih + 1;}
      }
      for (ih = 0; ih < maxchain; ih++) {
        fprintf(outfile, "%10d %10d\n", histval[ih], hist_allranks[ih]);
      }
      fclose(outfile);

      /*------------/
      / Stop Timer  /
      /------------*/
      end_time = MPI_Wtime();
      //printf("Total time elapsed on alpha %f: %lf seconds.\n", alpha, end_time - start_time);
    } /* if (my_rank == 0) */
  }
  if (my_rank==0) {
    printf("Total time elapsed on whole program: %lf seconds.\n", end_time - start_time);
  }
  /* Shutdown MPI */
  MPI_Finalize();

  exit(EXIT_SUCCESS);

  /* Release memory */
  for (iy = 0; iy < Ny; iy++) {
    free(force[iy]);
  }
  free(force);
  free(force_mem);
  free(oldforce);
  free(oldforce_mem);
  free(hist);
  free(histval);
  free(hist_allranks);

}

/* Parse Command line arguments */
int read_args(int argc, char ** argv, double * alpha, int * nalpha, int * Nsteps, int * Nx, int * Ny) {
  int index;
  int c;
  opterr = 0;

  /* Process all flags found in argc */
  while ((c = getopt(argc, argv, "a:n:N:X:Y:")) != -1)
    switch (c) {
    case 'a':
      *
      alpha = atof(optarg);
      if (( * alpha <= 0.0) || ( * alpha > 0.25)) {
        opterr = 1;
        fprintf(stderr, "Alpha argument could not be read: %s\n", optarg);
      }
      break;
    case 'n':
      *
      nalpha = atoi(optarg);
      if (( * nalpha <= 0)) {
        opterr = 1;
        fprintf(stderr, "nalpha argument could not be read: %s\n", optarg);
      }
      break;
    case 'N':
      *
      Nsteps = atoi(optarg);
      if ( * Nsteps <= 0) {
        opterr = 1;
        fprintf(stderr, "Nsteps argument could not be read: %s\n", optarg);
      }
      break;
    case 'X':
      *
      Nx = atoi(optarg);
      if ( * Nx <= 0) {
        opterr = 1;
        fprintf(stderr, "Nx argument could not be read: %s\n", optarg);
      }
      break;
    case 'Y':
      *
      Ny = atoi(optarg);
      if ( * Ny <= 0) {
        opterr = 1;
        fprintf(stderr, "Nx argument could not be read: %s\n", optarg);
      }
      break;
    case '?':
      if ((optopt == 'X') || (optopt == 'Y') || (optopt == 'N') || (optopt == 'a') || (optopt == 'n'))
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint(optopt))
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
      return 1;
    default:
      abort();
    }

  /* List unrecognised arguments */
  for (index = optind; index < argc; index++)
    printf("Non-option argument %s\n", argv[index]);

  return opterr;
}

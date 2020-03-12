/*==========================================================//
//  Evolves a function u in 2D via finite differences       //
//  of the Fisher-Kolmogorov equation                       //
//                                                          //
//  d u      d^2 u    d^2 u                                 //
//  ---  - D ----- -D -----  = (alpha/k) u (1 - u^q)        //
//  d t      d x^2    d y^2                                 //
//                                                          //
//  Fisher-Kolmogorov version by N. Hine - October 2019     //
//  Based on code created by D.Quigley - October 2015       //
//==========================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "makePNG.h"     /* For visualisation */
#include "mt19937ar.h"   /* Random number generator */

/* Function prototypes for memory management routines */
void allocate2d(double ***a,int num_rows,int num_cols);
void free2d(double ***a,int num_rows);

int main () {

  /* Function phi on new and current grid */
  double **u_new, **u;

  /* Approximate Laplacian */
  double Laplu;

  /* Number of grid points */
  register int Nx = 512;
  register int Ny = 1024;

  /* Loop counters */
  register int ix,iy, istep;

  /* Filename to which the grid is drawn */
  register int  isnap=0;
  char filename[25];

  /*--------------*/
  /* Initial time */
  /*--------------*/
  clock_t t1 = clock();

  /*------------------------------------*/
  /* Initialise random number generator */
  /*------------------------------------*/
  unsigned long seed = 120492783972;
  init_genrand(seed);

  /*--------------------------*/
  /* Set grid spacing         */
  /*--------------------------*/
  double dx = 1;
  double dy = 1;
  double dx2 = 1 / (dx*dx);
  double dy2 = 1 / (dy*dy);

  /*---------------------------------*/
  /* Set timestep, D, alpha, k and q */
  /*---------------------------------*/
  double dt = 1;
  double D  = 0.2;
  double alpha = 2;
  double k = 4;
  double alphak = alpha/k;
  double q = 1;

  /*--------------------------*/
  /* Number of steps to run   */
  /*--------------------------*/
  register int nstep = 2000;

  /*--------------------------------------*/
  /* Allocate memory for a bunch of stuff */
  /*--------------------------------------*/
  allocate2d(&u,Nx,Ny);
  allocate2d(&u_new,Nx,Ny);

  /*--------------------------------------*/
  /* Initialise with random numbers       */
  /*--------------------------------------*/
  for(ix=0;ix<Nx;ix++) {
    for(iy=0;iy<Ny;iy++) {
      u[ix][iy] = genrand();
    }
  }
      
  /*------------------------------------*/
  /* Write an image of the initial grid */
  /*------------------------------------*/
  int stepincr = 10; 
  sprintf(filename,"snapshot%08d.png",isnap); 
  writePNG(filename,u,Nx,Ny); 
  isnap++;

  /* setup time */
  clock_t t2 = clock();
  printf("Setup time                    : %15.6f seconds\n",(double)(t2-t1)/(double)CLOCKS_PER_SEC);
  t1 = t2;

  /*===============================*/
  /* BEGIN SECTION TO BE OPTIMISED */
  /*===============================*/
 
  /*--------------------------------------*/
  /* Loop over number of output timesteps */
  /*--------------------------------------*/
  for (istep=nstep;istep--;) {
    
    /*------------------------*/
    /* Loops over grid points */
    /*------------------------*/
    //both Nx and Ny are divisable by 8
    for(ix=Nx;ix--;) {
      for(iy=Ny;iy--;) {
            Laplu = 0; /* Initialise approximation to Laplacian */
            if (ix==Nx-1) {
            Laplu += (-2*u[ix][iy] + u[0][iy] + u[ix-1][iy])*dx2;
            } else if (ix==0){
            Laplu +=( -2*u[ix][iy] + u[ix+1][iy] + u[Nx-1][iy])*dx2;
            } else {
            Laplu += (-2*u[ix][iy] + u[ix+1][iy]+ u[ix-1][iy])*dx2;
            }
            if (iy==Ny-1) {
            Laplu += (-2*u[ix][iy-1] + u[ix][0] + u[ix][iy-1])*dy2;
            } else if (iy==0) {
             Laplu += (-2*u[ix][iy] + u[ix][iy+1] + u[ix][Ny-1])*dy2;
            } else {
            Laplu += (-2*u[ix][iy] + u[ix][iy+1] + u[ix][iy-1])*dy2; 
        }
        /* compute new value at this grid point */	
        u_new[ix][iy] = u[ix][iy] + dt*(alphak*(u[ix][iy] - pow(u[ix][iy],q+1)) + D*Laplu);
      }
    }
   
    //Pointer swap to transfer u_new into u
    double **tmp;
    tmp = u_new;
    u_new = u;
    u = tmp;

    /*-------------------------------*/
    /* Snapshots of grid to png file */
    /*-------------------------------*/
    if (istep==isnap)  { 
        sprintf(filename,"snapshot%08d.png",isnap); 
        writePNG(filename,u,Nx,Ny); 
        if (stepincr==1){ isnap += stepincr; } else { isnap *= stepincr;}
    }  
  }

  /*===============================*/
  /* END SECTION TO BE OPTIMISED */
  /*===============================*/
   
  /* calculation time */
  t2 = clock();
  printf("Time taken for %8d steps : %15.6f seconds\n",nstep,(double)(t2-t1)/(double)CLOCKS_PER_SEC);
    
  /*-----------------------------------*/
  /* Write png image of the final grid */
  /*-----------------------------------*/
  sprintf(filename,"snapshot%08d.png",istep); 
  writePNG(filename,u,Nx,Ny); 

  /*--------------------------------------------*/
  /* Write final time-evolved solution to file. */
  /*--------------------------------------------*/
  FILE *fp = fopen("final_grid.dat","w");
  if (fp==NULL) printf("Error opening final_grid.dat for output\n");

  for(ix=0;ix<Nx-1;ix++) {
    for(iy=0;iy<Ny-1;iy++) {
      /* x and y at the current grid points */
      double x = dx*(double)ix;
      double y = dy*(double)iy;
      fprintf(fp,"%8.4f %8.4f %8.4e\n",x,y,u[ix][iy]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  /* Release memory */
  free2d(&u,Nx);
  free2d(&u_new,Nx);
  
  return 0;
  
}
  


/*===========================================*/
/* Auxilliary routines for memory management */ 
/*===========================================*/
void allocate2d(double ***a,int Nx,int Ny) {

  double **b_loc; 

  b_loc = (double **)calloc(Nx,sizeof(double *));
  if (b_loc==NULL) printf("malloc error in allocate2d\n"); 

  int iy;
  for (iy=0;iy<Nx;iy++) {
    
    b_loc[iy] = (double *)calloc(Ny,sizeof(double));
    if (b_loc[iy]==NULL) printf("malloc error for row %d of %d in allocate2d\n",iy,Nx);

  }

  *a = b_loc;

}

void free2d(double ***a,int Nx) {

  int iy;

  double **b_loc = *a;

  /* Release memory */
  for (iy=0;iy<Nx;iy++) { 
    free(b_loc[iy]);
  }
  free(b_loc);
  *a = b_loc;

}









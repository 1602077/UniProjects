/*=========================================================/
/  Grid functions, PX425 2018 assignment 4. Contains all   /
/  functions which define the local (to each processor)    /
/  grid of spins within an Ising model.                    /
/=========================================================*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"

/* Integer constants defining coordinates and directions */
/* We take the left/right direction to be -x/+x */
/* We take the down/up directions to be   -y/+y */
const int x=0, y=1;
const int left=0, right=1, down=2, up=3;

/* Integers defining the local (to this rank) domain */
/* GLOBAL coordinates                                */
int grid_domain_size;
int grid_domain_start[2], grid_domain_end[2];

/* Integer arrays defining the 4 subgrids of this domain */
/* LOCAL coordinates                                     */
int grid_subdom_size;
int grid_subid[2];
int grid_subdom_start[4][2];
int grid_subdom_end[4][2];
  
/* Local grid of spins */
int **grid_spin;

/* Local grid of field values */
double **grid_field;

/* Global grid of spins */
int **global_grid_spin;

/* Four halo arrays containing copies of spins which  */
/* are neighbours to my boundary spins                */
int **grid_halo;

/* Array of the above structs to avoid repeat lookups */
neighbour_ptr_t **grid_neighbours;


void grid_initialise_local(int Ngrid,int p,int my_rank,int *my_rank_coords) {
  /*==================================================================/
  / Function to identify the section of the global grid on which      /
  / the current MPI task will operate, create four subdomains within  /
  / that section, and allocate memory accordingly.                    /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/

  /* Loop counters */
  int isub,iy,idir,ih;

  /* Domain in each dimension */
  grid_domain_size     = Ngrid/(int)sqrt((double)p+0.5);                          
  grid_domain_start[x] = my_rank_coords[x]*grid_domain_size;
  grid_domain_start[y] = my_rank_coords[y]*grid_domain_size;
 
  grid_domain_end[x]   = grid_domain_start[x] + grid_domain_size - 1;
  grid_domain_end[y]   = grid_domain_start[y] + grid_domain_size - 1;

  if (my_rank == 0) {
    printf("Size of each local processor grid : %5d x %5d\n",
	   grid_domain_size,grid_domain_size);
  }
  
  /* Set up the four subdomains on this rank */
  grid_subdom_size = grid_domain_size/2;
  for (isub = 0 ; isub < 4 ; isub++) {
   
    /* X and Y coordinates of the current sub-domain */
    grid_subid[x] = (isub)/2;
    grid_subid[y] = (isub+1)%2;
       
    /* Start and end LOCAL coordinates of the subdomain */
    grid_subdom_start[isub][x] =  grid_subid[x]*grid_subdom_size;   
    grid_subdom_start[isub][y] =  grid_subid[y]*grid_subdom_size;   
    
    grid_subdom_end[isub][x] = grid_subdom_start[isub][x] + grid_subdom_size - 1; 
    grid_subdom_end[isub][y] = grid_subdom_start[isub][y] + grid_subdom_size - 1; 
      
  }
    
  /* Allocate the array of spins local to this MPI task */
  grid_spin = (int **)malloc(grid_domain_size*sizeof(int *));
  if (grid_spin==NULL) { printf("Error allocating grid_spin pointers\n"); exit(EXIT_FAILURE); }
  for (iy = 0 ; iy < grid_domain_size ; iy++) {
    grid_spin[iy] = (int *)malloc(grid_domain_size*sizeof(int));
    if (grid_spin[iy]==NULL) { printf("Error allocating grid_spin array\n"); exit(EXIT_FAILURE); }
  }

  /* Allocate the array of random field values local to this MPI task */
  grid_field = (double **)malloc(grid_domain_size*sizeof(double *));
  if (grid_field==NULL) { printf("Error allocating grid_field pointers\n"); exit(EXIT_FAILURE); }
  for (iy = 0 ; iy < grid_domain_size ; iy++) {
    grid_field[iy] = (double *)malloc(grid_domain_size*sizeof(double));
    if (grid_field[iy]==NULL) { printf("Error allocating grid_field array\n"); exit(EXIT_FAILURE); }
  }

   
  /* Allocate the four halo arrays */
  grid_halo = (int **)malloc(4*sizeof(int *));
  if (grid_halo==NULL) { printf("Error allocating grid_halo pointers\n"); exit(EXIT_FAILURE); }
  for (idir = 0 ; idir < 4 ; idir++) {
    grid_halo[idir]=(int*)malloc(grid_domain_size*sizeof(int));
    if (grid_halo[idir]==NULL) { printf("Error allocating grid_halo array\n"); exit(EXIT_FAILURE); }
  }

  /* Fill halos with +1 to ensure no uninitialised data is encountered */
  for (idir = 0 ; idir < 4 ; idir++ ) {
    for (ih = 0 ; ih < grid_domain_size ; ih ++ ) {
      grid_halo[idir][ih] = +1;
    }
  }

  /* Allocate neighbour pointers for each local grid point */
  grid_neighbours = (neighbour_ptr_t **)malloc(grid_domain_size*sizeof(neighbour_ptr_t *));
  if (grid_neighbours==NULL) { printf("Error allocating grid_neighbour pointers\n"); exit(EXIT_FAILURE); }
  for (iy = 0 ; iy < grid_domain_size ; iy++) {
    grid_neighbours[iy] = (neighbour_ptr_t *)malloc(grid_domain_size*sizeof(neighbour_ptr_t));
    if (grid_neighbours[iy]==NULL) { printf("Error allocating grid_neighbour array\n"); exit(EXIT_FAILURE); }
  }
  
  /* Populate the above */
  grid_set_neighbour_spins();

  
}

void grid_initialise_global(int Ngrid,int my_rank, int p) {
  /*==================================================================/
  / Allocates memory to store the global grid. Normally called only   /
  / on MPI rank 0 before printing the final image.                    /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/
  int iy,ix;   /* loop counters */

  /* Only rank zero needs to allocate memory for the global grid  */
  /* and there's no point if P=1 as we can reuse existing storage */
  if (my_rank == 0 && p > 1 ) {

    /* Allocate the array of spins global to entire grid */
    global_grid_spin = (int **)malloc(Ngrid*sizeof(int *));
    if (global_grid_spin==NULL) { printf("Error allocating global_grid_spin pointers\n"); exit(EXIT_FAILURE); }
    for (iy = 0 ; iy < Ngrid ; iy++) {
      global_grid_spin[iy] = (int *)malloc(Ngrid*sizeof(int));
      if (global_grid_spin[iy]==NULL) { printf("Error allocating global_grid_spin array\n"); exit(EXIT_FAILURE); }
      /* Set elements so can safely generate an image even if missing contributions from other tasks */
      for (ix = 0 ; ix < Ngrid ; ix++ ) {
	global_grid_spin[iy][ix] = 0;
      }
    }
  }

}

void grid_shutdown_local() {
  /*==================================================================/
  / Release all memory allocated in grid_initialise_local.            /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/
  int iy;   /* loop counter */

  /* Release memory used in storing the grid of spins and fields */
  for (iy = 0 ; iy < grid_domain_size ; iy++) {
    free(grid_spin[iy]);
    free(grid_field[iy]);
  }
  free(grid_spin);
  free(grid_field);

  /* Release memory used in storing halo of data from neighbours */
  for (iy = 0 ; iy < 4 ; iy++) {
    free(grid_halo[iy]);
  }
  free(grid_halo);

}
void grid_shutdown_global(int Ngrid,int my_rank,int p) {
  /*==================================================================/
  / Release all memory allocated in grid_initialise_global.           /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/
  int iy;   /* error flag */

  /* Only rank zero should have allocated data */
  if (my_rank == 0 && p > 1 ) {

    /* Release memory used in storing the grid of spins */
    for (iy = 0 ; iy < Ngrid ; iy++) {
      free(global_grid_spin[iy]);
    }
    free(global_grid_spin);
  }

}

void grid_set_neighbour_spins() {
  /*==================================================================/
  / Given the LOCAL (to this MPI task) coordinates of a spin ix,iy    /
  / build a list of pointers to the four neighbouring spins. This     /
  / will sometimes point to to halo rather than local data if the     /
  / neighbours 'live' on other MPI tasks.                             /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/
    
  int ix,iy; /* loop counters */

  /* loop over all grid points */
  for (iy=0;iy<grid_domain_size;iy++){
    for(ix=0;ix<grid_domain_size;ix++) {

      /* Set pointer to left hand neighbour */
      if (ix==0) {
	grid_neighbours[iy][ix].left = &grid_halo[left][iy];
      } else {
	grid_neighbours[iy][ix].left = &grid_spin[iy][ix-1];
      }
    
      /* Set pointer to right hand neighbour */
      if (ix==grid_domain_size-1) {
	grid_neighbours[iy][ix].right = &grid_halo[right][iy];
      } else {
	grid_neighbours[iy][ix].right = &grid_spin[iy][ix+1];
      }
      
      /* Set pointer to bottom neighbour */
      if (iy==0) {
	grid_neighbours[iy][ix].down = &grid_halo[down][ix];
      } else {
	grid_neighbours[iy][ix].down = &grid_spin[iy-1][ix];
      }
      
      /* Set pointer to top neighbour */
      if (iy==grid_domain_size-1) {
	grid_neighbours[iy][ix].up = &grid_halo[up][ix];
      } else {
	grid_neighbours[iy][ix].up = &grid_spin[iy+1][ix];
      }

    } /* ix */

  } /* iy */

}


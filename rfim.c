/*========================================================/
/       ANSI-C code for PX425 assignment 4 2018           /
/ Performs Monte Carlo on the random field Ising model.   /
/ To be parallelised using an MPI domain decomposition.   /
/========================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "grid.h"      /* grid variables and constants  */
#include "comms.h"     /* comms variables and constants */
#include "mt19937ar.h" /* random number generator       */
#include "makePNG.h"   /* makePNG header                */

int main(int argc, char **argv) {

  /* Ising model coupling parameter */
  double J = 1.0;

  /* Width of random field distribution */
  double zeta = 0.65;

  /* Size of grid in each direction */
  int Ngrid = 480;
  
  /* Number of MC cycles to perform */
  int Ncyc = 1000;
  
  /* Counter of accepted moves */
  int accepted_moves = 0;
  
  /* Random number */
  double xi;
  
  /* Energy before and after an MC move */
  double energy_old,energy_new;
  
  /* Local and global magnetisation */
  double local_mag,global_mag;
  
  /* 1/kT in reduced units (T = temperature, k = 1) */
  double beta = 1.0/0.6;
  
  /* Loop counters, error flag and random number seed */
  int icyc,ix,iy,isub;
  long imove;
  unsigned long seed;

  /* Old spin for MC move resets */
  int oldspin;

  /* Filename to which the grid is drawn */
  char filename[25];

  /* Number of processors along each side of grid */
  int proot;

  /*-----------------------------------/
  / Initialise MPI, get my_rank and P. /
  /-----------------------------------*/
  /* You will need to modify this function - see comms.c */
  comms_initialise(argc,argv);

  /* Square root of number of processors */
  proot = (int)sqrt((double)p+0.5);

  /* Check that the grid size is divisible by 2Poot */
  if (Ngrid%(2*proot)!=0) {
    if (my_rank==0) { 
      printf("Total number of grid points must be divisible by 2sqrt(P)!\n");
    }
    comms_finalise();
    exit(EXIT_FAILURE);      
  }

  /*---------------------------------------/
  / Allocate sections of grid to MPI tasks /
  /---------------------------------------*/
  /* You will need to modify this function - see comms.c */
  comms_processor_map();

  /*--------------------------------/
  / Create the local grid of spins  /
  /--------------------------------*/
  /* You will need to modify this function - see comms.c */
  grid_initialise_local(Ngrid,p,my_rank,my_rank_coords);

  /*------------------------------/
  / Create the global grid array  /
  /------------------------------*/
  /* You will need to modify this function - see comms.c */
  grid_initialise_global(Ngrid,my_rank,p);

  /*------------------------------------------------------------------/
  / Initialise random number generator and populate local spin array  /
  /------------------------------------------------------------------*/
  /* different seed for each MPI task */
  seed = -1 * my_rank; 
  init_genrand(seed);

  /* Populate the initial array with spin = -1 */

  for (iy = 0 ; iy < grid_domain_size ; iy++ ) {
    for (ix = 0 ; ix < grid_domain_size ; ix++ ) {

      /* Create a random field */
      xi = genrand();
      grid_field[iy][ix] = 2.0*(xi-0.5)*zeta;

      /* Set spin randomly */
      xi = genrand();
      grid_spin[iy][ix] = 1;
      if (xi<0.5) grid_spin[iy][ix] *= -1;

    }
  }

  /*----------------------------------------/
  / Begin main Metropolis Monte Carlo loop  /
  /----------------------------------------*/
  /* Initialise count of accepted moves */
  accepted_moves = 0;

  /* Loop over Ncyc cyles */
  for (icyc = 0 ; icyc <= Ncyc ; icyc++ ) { 

    /* Loop over subdomains */
    for (isub = 0 ; isub < 4 ; isub++ ) {

      /*----------------------------------------------------------/
      / Perform halo swaps to update neighbours of boundary spins /
      /----------------------------------------------------------*/
      /* You need to modify this function - see comms.c */
      comms_halo_swaps();

      /* Perform moves within the current subdomain only, such
      / that these moves cannot invalidate neighbour information held on 
      / other MPI ranks. */
      for (imove = 0 ; imove < grid_subdom_size*grid_subdom_size ; imove++ ) {

	/* Choose a random point between subdomain start and end in each direction */
	xi  = genrand();
	ix  = grid_subdom_start[isub][x] + (int)((double)grid_subdom_size*xi);
	xi  = genrand();
	iy  = grid_subdom_start[isub][y] + (int)((double)grid_subdom_size*xi);

	/* Compute the old interactions with neighbours */
	energy_old = 0.0;
	energy_old -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].left;
	energy_old -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].right;
	energy_old -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].down;
	energy_old -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].up;

	/* Compute the interaction with the field at this grid point */
	energy_old -= grid_field[iy][ix]*grid_spin[iy][ix];

	/* Change the 'spin' at this point */
	oldspin = grid_spin[iy][ix];
	grid_spin[iy][ix] = -1*oldspin;

	/* Compute the new interactions with neighbours */
	energy_new = 0.0;
	energy_new -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].left;
	energy_new -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].right;
	energy_new -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].down;
	energy_new -= J*grid_spin[iy][ix]**grid_neighbours[iy][ix].up;

	/* Compute the interaction with the field at this grid point */
	energy_new -= grid_field[iy][ix]*grid_spin[iy][ix];

	/* Spin a random number to decide if we accept this change */
	xi = genrand();
	if ( xi < exp(-beta*(energy_new-energy_old)) ) {
	  /* Move accepted */
	  accepted_moves++;
	}  else {
	  /* Restore original spin */
	  grid_spin[iy][ix] = oldspin;
	}
	
      }  /* End loop over moves within the current subdomain */
      
    }  /* End loop over subdomains */
    
    /*-------------------------------------------/
    / Compute the local and global magnetisation /
    /-------------------------------------------*/
    /* Only do this every 100 steps */
    if (icyc%100==0) {

      /* Compute the local magnetisation */
      local_mag = 0.0;      
      for (iy = 0 ; iy < grid_domain_size ; iy++ ) {
	for (ix = 0 ; ix < grid_domain_size ; ix++ ) {
	  local_mag += (double)grid_spin[iy][ix];
	}
      }
      local_mag = local_mag/(double)(grid_domain_size*grid_domain_size);

      /* You will need to modify this routine - see comms.c */
      comms_get_global_mag(local_mag,&global_mag);

      /* Print global magnetisation */
      if (my_rank==0) {
	printf("Global magnetisation at cycle %8d : %12.6f\n",icyc,global_mag);
      }

      /*------------------------------------------------------------------------/
      / Gather data from all MPI ranks and create a PNG image of the whole grid /
      /------------------------------------------------------------------------*/
      /* Collect the glocal grid onto MPI task zero          */
      /* You will need to complete this function in comms.c  */
      comms_get_global_grid();
      
      /* Call function to write the global grid on rank zero only */
      sprintf(filename,"snapshot%08d.png",icyc);
      if (my_rank==0) {

	writepng(filename,global_grid_spin,Ngrid,Ngrid); 

      }

    }

  } /* End loop over MC cycles */

  /*-----------------------/
  / Report some statistics /
  /-----------------------*/
  printf("End of simulation. Rank %5d accepted %12d moves.\n",my_rank,accepted_moves);
  
  /*----------------/
  / Release storage / 
  /----------------*/
  grid_shutdown_local();
  grid_shutdown_global(Ngrid,my_rank,p);
  
  /*---------------/
  / Shut down MPI. /          
  /---------------*/
  /* You will need to complete this function - see comms.c */
  comms_finalise();

  return(EXIT_SUCCESS);

}

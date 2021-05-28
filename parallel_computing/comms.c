/*=========================================================/
/  Comms routines for PX425 assignment 4. Contains all     /
/  routines which interact with MPI libraries. Many of     /
/  these are currently incomplete and will work only in    /
/  serial. You will need to correct this.                  /
/=========================================================*/
#include "mpi.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comms.h"

int p;               /* Number of processor       */
int my_rank;         /* Rank of current processor */
int my_rank_in_cart; /* Rank in cartesian grid */

MPI_Comm cart_comm;  /* Cartesian communicator    */

/* Coordinates of current rank in the processor grid */
int my_rank_coords[2];

/* Ranks of neighbours to the current processor (left, right, down, up) */
int my_rank_neighbours[4];

/* Time and initialisation and shutdown */
double t1,t2;

void comms_initialise(int argc, char **argv) {
  /*==================================================================/
  / Function to initialise MPI, get the communicator size p and the   /
  / rank my_rank of the current process within that communicator.     / 
  /==================================================================*/
  
  int proot;                   /* square root of p */

  /* Initial MPI and get communicator of size p and rank of current   /
  /   process within that communicator                               */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Start timer after initialiasing MPI */
  t1 = MPI_Wtime();

  /* Check that we have a square number of processors */
  proot = (int)sqrt((double)p+0.5);
  if (proot*proot!=p) {
    if (my_rank==0) {
      printf("Number of processors must be an exact square!\n");
      exit(EXIT_FAILURE); 
    }
  }

  return;

}

void comms_processor_map() {
  /*====================================================================/
  / Function to map our p processors into a 2D Cartesian grid of      /
  / dimension proot by proot where proot = sqrt(p).                   /
  /                                                                   /
  / Should populate the arrays my_rank_cooords, which contains the    /
  / location of the current MPI task within the processor grid, and   /
  / my_rank_neighbours, which contains (in the order left, right,     /
  / down and up) ranks of neighbouring MPI tasks on the grid with     /
  / which the current task will need to communicate.                  /
  /==================================================================*/
  
  /* Information for setting up a Cartesian communicator */
  int ndims = 2;
  int reorder = 0;
  int pbc [2] = {1,1};
  int dims[2];
    
  /* Local variables */
  int proot;            /* square root of p */

  /* Square root of number of processors */
  proot = (int)sqrt((double)p+0.5);
  
  /* Dimensions of Cartesian communicator */
  dims[x] = proot; dims[y] = proot;

  /* Creating Cartesian topology */
  MPI_Cart_create (MPI_COMM_WORLD, 2, dims, pbc, reorder, &cart_comm);
  /* No need to check for reordering, as reorder=0 */

  /* Assigning communicator for in the cartesian grid*/
  MPI_Comm_rank (cart_comm, &my_rank_in_cart);

  /* Find coordinates of current rank in the topology */
  MPI_Cart_coords (cart_comm, my_rank_in_cart, 2, my_rank_coords);
 
  /* Constructing halo of neighbouring data */
  MPI_Cart_shift (cart_comm, x, 1, &my_rank_neighbours[left],  &my_rank_neighbours[right]);
  MPI_Cart_shift (cart_comm, y, 1, &my_rank_neighbours[down], &my_rank_neighbours[up]);

  return;

} 

void comms_get_global_mag(double local_mag,double *global_mag) {
  /*==================================================================/
  / Function to compute the glocal magnetisation of the grid by       /
  / averaging over all values of local_mag, and storing the result    /
  / in global_mag.                                                    /
  /==================================================================*/

  /* Assigns local_mag to global_mag, valid on one processor only */
  /* i.e. when the local grid size is the same as the gloabl one */
  *global_mag = local_mag;


  /* Reduce operation sums all values of local mag and assigns them to  /
  /  global mag  on rank 0.                                            */
  MPI_Reduce(&local_mag, global_mag,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  /*In combination with MPI_Reduce operation this calculates             / 
   / this calculates the  average global magnetisation                  */
  *global_mag = *global_mag/(double)p;

} 

void comms_halo_swaps() {
  /*==================================================================/
  / Function to send boundary spins on each side of the local grid    /
  / to neighbour processors, and to receive from those processors the /
  / halo information needed to perform computations involving spins   /
  / on the boundary processors grid.                                  /
  /==================================================================*/

  /* Send and receive buffers */
  int *sendbuf, *recvbuf;
    
  /* MPI Status */
  MPI_Status status;

  int ix,iy;      /* Loop counters */

  /* If running on 1 processor copy boundary elements into opposite halo */
  if (p==1) {

    for (iy = 0 ; iy < grid_domain_size ; iy++) {
      grid_halo[right][iy] = grid_spin[iy][0]; 
      grid_halo[left][iy]  = grid_spin[iy][grid_domain_size-1];
    }

    for (ix = 0 ; ix < grid_domain_size ; ix++) {
      grid_halo[up][ix]   = grid_spin[0][ix];
      grid_halo[down][ix] = grid_spin[grid_domain_size-1][ix];
    }

    return; /* Do not do any comms */
  }

  /* Allocate buffers */
  sendbuf = (int *)malloc(grid_domain_size*sizeof(int));
  if (sendbuf==NULL) { 
    printf("Error allocating sendbuf in comms_halo_swaps\n");
    exit(EXIT_FAILURE);
  }
  recvbuf = (int *)malloc(grid_domain_size*sizeof(int));
  if (recvbuf==NULL) { 
    printf("Error allocating recvbuf in comms_halo_swaps\n");
    exit(EXIT_FAILURE);
  }

  /*Intialising status and request structures used by MPI_Isend and MPI_Irecv
   * to confirm that communications have bee successful sent and recieved   */
  MPI_Status send_status1, recv_status1;
  MPI_Request send_req1, recv_req1;


   for (iy = 0; iy < grid_domain_size; iy++) {
     /* Filling sendbuf with the left hand boundary elements of grid_spin */
     sendbuf[iy] = grid_spin[iy][0];
   }
    
  /* Send left hand boundary elements of grid_spin to my_rank_neighbours[left] 
     and receive from my_rank_neighbours[right] into the appropriate part
     of grid_halo. Remember to use the appropriate communicator. */
   MPI_Isend (sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], my_rank_neighbours[left]+888, MPI_COMM_WORLD, &send_req1);
   MPI_Irecv (recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], my_rank+888, MPI_COMM_WORLD, &recv_req1);

   MPI_Wait (&send_req1, &send_status1);
   MPI_Wait (&recv_req1, &send_status1);

   for (iy = 0; iy < grid_domain_size; iy++) {
    /* Allocating right boundary elements stored in recvbuf into grid_halo[right][iy] */
       grid_halo[right][iy] = recvbuf[iy];
   }

   MPI_Status send_status2, recv_status2;
   MPI_Request send_req2, recv_req2;

   for (iy = 0; iy < grid_domain_size; iy++) {
     /* Filling sendbuf with the right hand boundary elements of grid_spin */
     sendbuf[iy] = grid_spin[iy][grid_domain_size-1];
   }

   /* Send right hand boundary elements of grid_spin to my_rank_neighbours[right]
     and receive from my_rank_neighbours[left] into the appropriate part 
     of grid_halo. Remember to use the appropriate communicator. */
   MPI_Isend (sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], my_rank_neighbours[right]+889, MPI_COMM_WORLD, &send_req2);
   MPI_Irecv (recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], my_rank+889, MPI_COMM_WORLD, &recv_req2);


   MPI_Wait (&send_req2, &send_status2);
   MPI_Wait (&recv_req2, &send_status2);

   for (iy = 0; iy < grid_domain_size; iy++) {
    /* Allocating left boundary elements stored in recvbuf into grid_halo[left][iy] */
       grid_halo[left][iy] = recvbuf[iy];
   }

   MPI_Status send_status3, recv_status3;
   MPI_Request send_req3, recv_req3;

   for (ix = 0; ix < grid_domain_size; ix++) {
    /* Filling sendbuf with the up  boundary elements of grid_spin */
     sendbuf[ix] = grid_spin[0][ix];
   }
   /* Send top boundary elements of grid_spin to my_rank_neighbours[up]
     and receive from my_rank_neighbours[down] into the appropriate part
     of grid halo.*/
   MPI_Isend (sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], my_rank_neighbours[down]+777, MPI_COMM_WORLD, &send_req3);
   MPI_Irecv (recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], my_rank+777, MPI_COMM_WORLD, &recv_req3);

   MPI_Wait (&send_req3, &send_status3);
   MPI_Wait (&recv_req3, &send_status3);

   for (ix = 0; ix < grid_domain_size; ix++) {
        /* Allocating up boundary elements stored in recvbuf into grid_halo[up][ix] */
       grid_halo[up][ix] = recvbuf[ix];
   }
 
   MPI_Status send_status4, recv_status4;
   MPI_Request send_req4, recv_req4;

   for (ix = 0; ix < grid_domain_size; ix++) {
     /* Filling sendbuf with the down  boundary elements of grid_spin */
     sendbuf[ix] = grid_spin[grid_domain_size-1][ix];
   }
    /* Send bottom boundary elements of grid_spin to my_rank_neighbours[down]
     and receive from my_rank_neighbours[up] into the appropriate part
     of grid halo.  */
   MPI_Isend (sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], my_rank_neighbours[up]+666, MPI_COMM_WORLD, &send_req4);
   MPI_Irecv (recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], my_rank+666, MPI_COMM_WORLD, &recv_req4);

   MPI_Wait (&send_req4, &send_status4);
   MPI_Wait (&recv_req4, &send_status4);

   for (ix = 0; ix < grid_domain_size; ix++) {
       /* Allocating down boundary elements stored in recvbuf into grid_halo[down][ix] */
       grid_halo[down][ix] = recvbuf[ix];
   }

  /* Release memory */
  free(sendbuf);
  free(recvbuf);

  return;

}


void comms_get_global_grid() {
  /*==================================================================/
  / Function to collect all contributions to the global grid onto     /
  / rank zero for visualisation.                                      /
  /==================================================================*/

  /* comms buffer */
  int *combuff;

  /* MPI Status */
  MPI_Status status;

  /* Information on the remote domain */
  int remote_domain_start[2] = {0,0};

  /* Loop counters and error flag */
  int ix,iy,ixg,iyg,ip;

  /* Just point at local grid if running on one processor */
  if (p==1) {
    global_grid_spin=grid_spin;
    return;
  }

  if (my_rank==0) {

    /* Rank 0 first fills out its part of the global grid */
    for (iy = 0 ; iy < grid_domain_size ; iy++) {
      for (ix = 0 ; ix < grid_domain_size ; ix++) {
	
	/* Global indices */
	ixg = ix+grid_domain_start[x];
	iyg = iy+grid_domain_start[y];
	
	global_grid_spin[iyg][ixg] = grid_spin[iy][ix];
      }
    }

  } /* End of if: rank 0 filling global_grid_spin with its data */

  /* Allocate buffer */
  combuff = (int *)malloc(grid_domain_size*sizeof(int));
  if (combuff==NULL) { 
    printf("Error allocating combuff in comms_get_global_grid\n");
    exit(EXIT_FAILURE);
  }

  if (my_rank==0) {
    
    /* Now loops over all other ranks receiving their data */
    for (ip = 1; ip < p ; ip++) {

      /* First receive remote_domain_start from rank ip */
      MPI_Recv(&remote_domain_start, 2, MPI_INT, ip, 998, MPI_COMM_WORLD, &status);


      /* Loop over rows within a domain */
      for (iy = 0 ; iy < grid_domain_size ; iy++) {
	
	/* Receive this row of grid spin data from rank ip */ 
    MPI_Recv(combuff,grid_domain_size,MPI_INT,ip,999,MPI_COMM_WORLD, &status);

	for (ix = 0 ; ix < grid_domain_size ; ix++) {

	  /* Global indices */
	  ixg = ix+remote_domain_start[x];
	  iyg = iy+remote_domain_start[y];
	
	  /* Store in global_grid_spin */
	  global_grid_spin[iyg][ixg] = combuff[ix];


	} /* elements in row */

      } /* rows */

    } /* processors */

  } else {

    /* Sends start point (grid_domain_start) for each rank to rank 0 */
    MPI_Send (&grid_domain_start,2, MPI_INT,0,998,MPI_COMM_WORLD);

    /* Loop over rows in the domain, sending them to rank 0*/
    for (iy = 0 ; iy < grid_domain_size ; iy++) {

      /* Sending grid_spin data row by row to rank 0 */
      MPI_Send (&grid_spin[iy][0],grid_domain_size, MPI_INT,0,999,MPI_COMM_WORLD);
    }
    
  } /* End of else statement */
  
  /* Free memory */
  free(combuff);

  return;
  
}
 

void comms_finalise() {
  /*==================================================================/
  / Function to finalise MPI functionality and exit cleanly           /
  /==================================================================*/
  
  /* Measure the time t2 using MPI_Wtime() which returns a double */
  t2 = MPI_Wtime();

  if (my_rank==0 && p>=1) {
    printf("Total time elapsed since MPI initialised :  %12.6f s\n",t2-t1);
  }

  /* Shutdown MPI */
  MPI_Finalize();
  return;

}

			     
  



